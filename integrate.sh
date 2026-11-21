#!/usr/bin/env bash
#
# integrate.sh - build the lab's combined GROMACS branch
#
# WHAT THIS IS FOR
#   The lab runs GROMACS with a number of our merge requests (MRs) that are
#   still under review upstream. This script takes upstream GROMACS `main`
#   at one commit, adds those MRs on top in a fixed order, and produces
#
#     * a local branch (default `lab/main`) holding the combination, and
#     * a manifest recording exactly what went in: the base commit, and for
#       every MR its number, the head commit used, its state and title.
#
#   The last commit on the branch sets GMX_VERSION_STRING_OF_FORK, so every
#   binary built from it reports a version like
#   "2027.0-smithlab_integration_20261121_bc88adaf-dev" (upstream version,
#   lab tag, date of the run, base commit). That string can
#   never be mistaken for an upstream release. The manifest is committed
#   into the tree as lab-integration-manifest.txt and a copy, with the
#   resulting commit added, goes to manifests/.
#
#   MR heads are fetched from the upstream project itself
#   (refs/merge-requests/<iid>/head on gitlab.com/gromacs/gromacs), so the
#   script does not depend on anyone's personal fork.
#
# USAGE
#   ./integrate.sh run      fetch, combine, verify, write branch + manifest.
#                           Exit 0 = done, 2 = error, 3 = conflict.
#   ./integrate.sh check    report only: has upstream main moved, has any MR
#                           head moved, been merged or closed, since the last
#                           manifest? Needs no clone. Exit 0 = nothing moved,
#                           1 = something moved (re-run `run`), 2 = error.
#                           Meant to be called from a scheduled job.
#
#   Settings (environment variables, all optional):
#     GMX_REPO       clone to work in (default: ./gromacs next to this
#                    script; cloned on first use)
#     REFERENCE_REPO an existing GROMACS clone to borrow objects from when
#                    cloning, to save download time (the clone is still
#                    self-contained afterwards)
#     BASE_COMMIT    pin the upstream base (default: tip of upstream main)
#     OUTPUT_BRANCH  default lab/main
#     LAB_TAG        default smithlab
#
# THE MR LIST
#   mr-list.txt, next to this script, is the only place that says which MRs
#   go in and in what order. See its header. To add an MR, add a line; to
#   drop one, delete its line. MRs upstream has merged are skipped and
#   reported. A closed (not merged) MR stops the run: decide whether to
#   drop it from the list.
#
#   Stacked MRs form a series. Only the last MR of a series is merged; the
#   earlier ones must be inside it. Each one is still recorded and checked:
#   either its head is an ancestor of the result, or (when the MR branch on
#   GitLab is older than the rebased copy inside the later MR) all of its
#   commits are present as identical patches. The manifest says which.
#
#   An MR whose target is not main (e.g. release-2026) is carried forward by
#   cherry-picking its own commits, i.e. <target>..<head>, one by one, each
#   with a "(cherry picked from commit ...)" line. Their changes under
#   docs/release-notes/ are NOT carried: those notes describe the release
#   branch, the files differ between branches, and main receives them when
#   upstream merges the release branch into main. Every carried commit is
#   verified to be patch-identical to the original outside that directory.
#
# ON A CONFLICT
#   The script stops, names the MR being applied, the MR(s) already applied
#   that touched the same files, and the files. It resolves nothing and
#   leaves the output branch and the previous manifest untouched. Do not
#   resolve the conflict here. Ask the owners of the two MRs to make them
#   compatible upstream (usually by rebasing one MR onto the other), then
#   re-run. If one MR must wait, comment it out in mr-list.txt.
#
# REQUIREMENTS
#   git, curl, python3, network access to gitlab.com. No GitLab token: only
#   public data is read.

set -euo pipefail

UPSTREAM_URL=https://gitlab.com/gromacs/gromacs.git
API=https://gitlab.com/api/v4/projects/gromacs%2Fgromacs
BASE_BRANCH=main

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MR_LIST=$HERE/mr-list.txt
MANIFEST_DIR=$HERE/manifests
GMX_REPO=${GMX_REPO:-$HERE/gromacs}
OUTPUT_BRANCH=${OUTPUT_BRANCH:-lab/main}
LAB_TAG=${LAB_TAG:-smithlab}
IN_TREE_MANIFEST=lab-integration-manifest.txt
CARRY_EXCLUDE=docs/release-notes   # not carried forward from release-branch MRs

die()  { echo "ERROR: $*" >&2; exit 2; }
note() { echo "== $*"; }

# --- MR list -----------------------------------------------------------------

# Prints "series iid" per MR, in file order.
read_list() {
    [ -f "$MR_LIST" ] || die "missing $MR_LIST"
    sed -e 's/#.*//' "$MR_LIST" | awk 'NF >= 2 { print $1, $2 }'
}

# Prints the series names in order of first appearance.
series_order() { read_list | awk '!seen[$1]++ { print $1 }'; }

# --- GitLab API (public, read-only) ------------------------------------------

# Prints "state<TAB>target_branch<TAB>sha<TAB>title" for one MR.
mr_info() {
    curl -sf "$API/merge_requests/$1" | python3 -I -c '
import json, sys
m = json.load(sys.stdin)
print("\t".join([m["state"], m["target_branch"], m["sha"], " ".join(m["title"].split())]))
' || die "could not read !$1 from the GitLab API"
}

# --- check mode --------------------------------------------------------------

do_check() {
    local last=$MANIFEST_DIR/latest.txt
    [ -f "$last" ] || { echo "No previous manifest ($last): run './integrate.sh run' first."; exit 1; }
    local moved=0 old new
    old=$(awk -F= '$1 == "base_commit" { print $2 }' "$last")
    new=$(git ls-remote "$UPSTREAM_URL" "refs/heads/$BASE_BRANCH" | cut -f1)
    [ -n "$new" ] || die "cannot read upstream $BASE_BRANCH"
    if [ "$old" != "$new" ]; then
        echo "upstream $BASE_BRANCH moved: ${old:0:10} -> ${new:0:10}"; moved=1
    else
        echo "upstream $BASE_BRANCH unchanged at ${old:0:10}"
    fi
    local series iid state target sha title oldsha
    while read -r series iid; do
        IFS=$'\t' read -r state target sha title < <(mr_info "$iid")
        oldsha=$(awk -F'\t' -v i="$iid" '$1 == "mr" && $2 == i { print $3 }' "$last")
        if [ -z "$oldsha" ]; then
            echo "!$iid is new in mr-list.txt ($state)"; moved=1
        elif [ "$state" != opened ]; then
            echo "!$iid is now $state upstream"; moved=1
        elif [ "$sha" != "$oldsha" ]; then
            echo "!$iid head moved: ${oldsha:0:10} -> ${sha:0:10}  ($title)"; moved=1
        else
            echo "!$iid unchanged at ${sha:0:10}"
        fi
    done < <(read_list)
    while IFS=$'\t' read -r _ iid _; do
        read_list | awk -v i="$iid" '$2 == i { f = 1 } END { exit !f }' \
            || { echo "!$iid was dropped from mr-list.txt"; moved=1; }
    done < <(awk -F'\t' '$1 == "mr"' "$last")
    if [ $moved = 1 ]; then echo "RESULT: changed, re-run './integrate.sh run'"; exit 1; fi
    echo "RESULT: nothing moved"; exit 0
}

# --- run mode ----------------------------------------------------------------

G() { git -C "$GMX_REPO" "$@"; }

ensure_repo() {
    if [ ! -d "$GMX_REPO/.git" ]; then
        note "cloning $UPSTREAM_URL into $GMX_REPO"
        if [ -n "${REFERENCE_REPO:-}" ]; then
            git clone --reference-if-able "$REFERENCE_REPO" --dissociate "$UPSTREAM_URL" "$GMX_REPO"
        else
            git clone "$UPSTREAM_URL" "$GMX_REPO"
        fi
    fi
    [ "$(G remote get-url origin)" = "$UPSTREAM_URL" ] || die "$GMX_REPO: origin is not $UPSTREAM_URL"
    [ -z "$(G status --porcelain --untracked-files=no)" ] || die "$GMX_REPO has local changes"
    [ ! -e "$(G rev-parse --git-path MERGE_HEAD)" ] && [ ! -e "$(G rev-parse --git-path CHERRY_PICK_HEAD)" ] \
        || die "$GMX_REPO has a merge or cherry-pick in progress"
}

# File(s) a commit range touches.
files_of() { G diff --name-only "$1" "$2"; }

# The first MR of a series (in list order) whose head contains commit $1.
owner_of() {
    local c=$1; shift
    local i
    for i in "$@"; do
        if G merge-base --is-ancestor "$c" "refs/lab/mr/$i" 2>/dev/null; then echo "!$i"; return; fi
    done
    echo "(unknown MR)"
}

# Reports a conflict, puts the repository back, exits.
conflict() {
    local what=$1 series=$2 range_from=$3 range_to=$4
    local files f c owners
    files=$(G diff --name-only --diff-filter=U)
    echo
    echo "CONFLICT while applying series '$series' ($what)."
    echo "Nothing was resolved, and $OUTPUT_BRANCH and the manifests are unchanged."
    for f in $files; do
        c=$(G log -1 --format=%H "$range_from..$range_to" -- "$f")
        owners=""
        local s
        for s in "${APPLIED_SERIES[@]}"; do
            local e
            e=$(G log -1 --format=%H "$BASE..${APPLIED_TIP[$s]}" -- "$f" 2>/dev/null || true)
            if [ -n "$e" ]; then
                # shellcheck disable=SC2086
                owners+=" $(owner_of "$e" ${SERIES_MRS[$s]}) (series $s)"
            fi
        done
        echo "  $f"
        echo "      incoming: $(owner_of "${c:-$range_to}" ${SERIES_MRS[$series]})"
        echo "      already applied:${owners:- (no earlier MR touched it; conflict is with upstream $BASE_BRANCH)}"
    done
    G merge --abort 2>/dev/null || G cherry-pick --abort 2>/dev/null || true
    G switch -q --detach "$BASE"
    exit 3
}

# Patch-id of a commit, ignoring CARRY_EXCLUDE.
patch_id() { G diff "$1^" "$1" -- . ":(exclude)$CARRY_EXCLUDE" | G patch-id --stable | cut -d' ' -f1; }

# The carried copy (in BASE..HEAD) of original commit $1, if any.
carried_copy() { G log -1 --format=%H --fixed-strings --grep="(cherry picked from commit $1)" "$BASE..HEAD"; }

# Cherry-picks the own commits of MR $1 (relative to branch $2) onto HEAD,
# leaving out CARRY_EXCLUDE. Prints the number of commits carried.
carry_forward() {
    local iid=$1 from=$2 c msg carried=0
    msg=$(mktemp)
    for c in $(G rev-list --reverse "$from..refs/lab/mr/$iid"); do
        G cherry-pick --no-commit "$c" >/dev/null 2>&1 || { rm -f "$msg"; return 1; }
        if G ls-files --error-unmatch -- "$CARRY_EXCLUDE" >/dev/null 2>&1 \
           || ! G diff --cached --quiet -- "$CARRY_EXCLUDE"; then
            G restore --staged --worktree --source=HEAD -- ":(glob)$CARRY_EXCLUDE/**"
        fi
        if G diff --cached --quiet; then continue; fi   # it only touched CARRY_EXCLUDE
        { G log -1 --format=%B "$c"; echo "(cherry picked from commit $c)"; } > "$msg"
        G commit -q -F "$msg" --author="$(G log -1 --format='%an <%ae>' "$c")" \
            --date="$(G log -1 --format=%aD "$c")"
        carried=$((carried + 1))
    done
    rm -f "$msg"
    echo "$carried"
}

do_run() {
    ensure_repo
    command -v curl >/dev/null || die "curl is needed"

    # Read the list and the state of every MR.
    declare -gA SERIES_MRS STATE TARGET SHA TITLE APPLIED_TIP HOW
    declare -ga APPLIED_SERIES=()
    local series iid
    while read -r series iid; do
        SERIES_MRS[$series]="${SERIES_MRS[$series]:-} $iid"
        IFS=$'\t' read -r "STATE[$iid]" "TARGET[$iid]" "SHA[$iid]" "TITLE[$iid]" < <(mr_info "$iid")
    done < <(read_list)

    # Fetch base branch, the release branches MRs target, and every MR head.
    local refspecs=() targets
    targets=$(for iid in "${!TARGET[@]}"; do echo "${TARGET[$iid]}"; done | sort -u)
    for t in $BASE_BRANCH $targets; do refspecs+=("+refs/heads/$t:refs/remotes/origin/$t"); done
    for iid in "${!STATE[@]}"; do refspecs+=("+refs/merge-requests/$iid/head:refs/lab/mr/$iid"); done
    note "fetching ${#refspecs[@]} refs from upstream"
    G fetch -q origin "${refspecs[@]}"

    BASE=$(G rev-parse "${BASE_COMMIT:-origin/$BASE_BRANCH}^{commit}")
    note "base: $(G log -1 --format='%h %cs %s' "$BASE")"
    for iid in "${!STATE[@]}"; do
        local fetched; fetched=$(G rev-parse "refs/lab/mr/$iid")
        [ "$fetched" = "${SHA[$iid]}" ] || die "!$iid: API says head ${SHA[$iid]:0:10}, fetched ${fetched:0:10}; re-run"
        case ${STATE[$iid]} in
            opened|merged) ;;
            *) die "!$iid is ${STATE[$iid]} upstream. Drop it from mr-list.txt or find out why." ;;
        esac
    done

    G switch -q --detach "$BASE"

    # Apply each series.
    local mrs tip
    for series in $(series_order); do
        mrs=${SERIES_MRS[$series]}
        tip=""
        for iid in $mrs; do [ "${STATE[$iid]}" = opened ] && tip=$iid; done
        for iid in $mrs; do [ "${STATE[$iid]}" = merged ] && { HOW[$iid]="skipped: merged upstream"; note "!$iid is merged upstream, skipped"; }; done
        if [ -z "$tip" ]; then note "series $series: everything merged upstream"; continue; fi

        if [ "${TARGET[$tip]}" = "$BASE_BRANCH" ]; then
            note "series $series: merging tip !$tip ${SHA[$tip]:0:10}"
            G merge -q --no-ff --no-edit -m "Merge !$tip: ${TITLE[$tip]}" \
                -m "Series '$series' of the lab integration; MR head ${SHA[$tip]}." "refs/lab/mr/$tip" \
                >/dev/null 2>&1 || conflict "merge of !$tip" "$series" "$BASE" "refs/lab/mr/$tip"
            HOW[$tip]="merged (series tip)"
        else
            local from=origin/${TARGET[$tip]} n carried
            n=$(G rev-list --count "$from..refs/lab/mr/$tip")
            note "series $series: carrying !$tip forward from ${TARGET[$tip]} ($n commits)"
            carried=$(carry_forward "$tip" "$from") \
                || conflict "cherry-pick of !$tip from ${TARGET[$tip]}" "$series" "$from" "refs/lab/mr/$tip"
            HOW[$tip]="cherry-picked from ${TARGET[$tip]}: $carried of $n commits ($((n - carried)) touched only $CARRY_EXCLUDE)"
        fi
        APPLIED_SERIES+=("$series")
        APPLIED_TIP[$series]=$(G rev-parse "refs/lab/mr/$tip")
    done

    # Verify that every open MR is contained in the result.
    local missing=0 from unmatched
    for series in $(series_order); do
        for iid in ${SERIES_MRS[$series]}; do
            [ "${STATE[$iid]}" = opened ] || continue
            if G merge-base --is-ancestor "refs/lab/mr/$iid" HEAD; then
                HOW[$iid]=${HOW[$iid]:-contained (head is an ancestor)}
                continue
            fi
            if [ "${TARGET[$iid]}" != "$BASE_BRANCH" ]; then
                local c copy bad=0
                for c in $(G rev-list "origin/${TARGET[$iid]}..refs/lab/mr/$iid"); do
                    copy=$(carried_copy "$c")
                    if [ -z "$copy" ]; then
                        [ -z "$(patch_id "$c")" ] && continue   # only touched CARRY_EXCLUDE
                        echo "NOT CONTAINED: !$iid commit ${c:0:10} was not carried forward."; bad=1
                    elif [ "$(patch_id "$c")" != "$(patch_id "$copy")" ]; then
                        echo "DIFFERS: !$iid commit ${c:0:10} and its carried copy ${copy:0:10} are not the same patch."; bad=1
                    fi
                done
                [ $bad = 0 ] || missing=1
                continue
            fi
            from=$(G merge-base "origin/$BASE_BRANCH" "refs/lab/mr/$iid")
            unmatched=$(G cherry HEAD "refs/lab/mr/$iid" "$from" | grep -c '^+' || true)
            if [ "$unmatched" = 0 ]; then
                HOW[$iid]=${HOW[$iid]:-contained (identical patches; the MR head on GitLab is older than the copy merged)}
            else
                echo "NOT CONTAINED: !$iid has $unmatched commit(s) missing from the result."
                echo "  Its series tip does not include it. Check the order in mr-list.txt, or ask"
                echo "  the series owner to rebase the later MRs onto it."
                missing=1
            fi
        done
    done
    if [ $missing = 1 ]; then
        G switch -q --detach "$BASE"
        die "verification failed; $OUTPUT_BRANCH and the manifests are unchanged"
    fi

    # Version string and manifest, committed as the last commit.
    local date stamp version
    date=$(date -u +%Y%m%d)
    stamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    version="${LAB_TAG}_integration_${date}_$(G rev-parse --short=8 "$BASE")"
    sed -i "s/^set(GMX_VERSION_STRING_OF_FORK \"[^\"]*\"/set(GMX_VERSION_STRING_OF_FORK \"$version\"/" \
        "$GMX_REPO/cmake/gmxVersionInfo.cmake"
    grep -q "GMX_VERSION_STRING_OF_FORK \"$version\"" "$GMX_REPO/cmake/gmxVersionInfo.cmake" \
        || die "could not set GMX_VERSION_STRING_OF_FORK"
    {
        echo "# GROMACS lab integration manifest. Generated by integrate.sh; do not edit."
        echo "version_string_of_fork=$version"
        echo "generated=$stamp"
        echo "upstream=$UPSTREAM_URL"
        echo "base_branch=$BASE_BRANCH"
        echo "base_commit=$(G rev-parse "$BASE")"
        echo "base_subject=$(G log -1 --format=%s "$BASE")"
        echo "# mr<TAB>iid<TAB>head<TAB>state<TAB>target<TAB>series<TAB>how<TAB>title"
        while read -r series iid; do
            printf 'mr\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$iid" "${SHA[$iid]}" "${STATE[$iid]}" \
                "${TARGET[$iid]}" "$series" "${HOW[$iid]}" "${TITLE[$iid]}"
        done < <(read_list)
    } > "$GMX_REPO/$IN_TREE_MANIFEST"
    G add cmake/gmxVersionInfo.cmake "$IN_TREE_MANIFEST"
    G commit -q -m "Lab integration build $version" \
        -m "Not an upstream GROMACS release. Contents: $IN_TREE_MANIFEST."
    local result; result=$(G rev-parse HEAD)

    G branch -f "$OUTPUT_BRANCH" "$result"
    G switch -q "$OUTPUT_BRANCH"
    mkdir -p "$MANIFEST_DIR"
    local out=$MANIFEST_DIR/manifest-$date-${result:0:10}.txt
    { cat "$GMX_REPO/$IN_TREE_MANIFEST"; echo "result_branch=$OUTPUT_BRANCH"; echo "result_commit=$result"; } > "$out"
    cp "$out" "$MANIFEST_DIR/latest.txt"

    echo
    note "done: $OUTPUT_BRANCH = ${result:0:10} in $GMX_REPO"
    note "version string: $version"
    note "manifest: $out"
    awk -F'\t' '$1 == "mr" { printf "   !%s %s  %s\n", $2, substr($3, 1, 10), $7 }' "$out"
}

case ${1:-} in
    run)   do_run ;;
    check) do_check ;;
    *)     sed -n '2,/^set -euo/p' "$0" | sed -e '$d' -e 's/^# \{0,1\}//'; exit 2 ;;
esac
