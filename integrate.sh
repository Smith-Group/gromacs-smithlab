#!/usr/bin/env bash
#
# integrate.sh - build the lab's combined GROMACS branch
#
# WHAT THIS IS FOR
#   The lab runs GROMACS with a number of our merge requests (MRs) that are
#   still under review upstream. This script takes upstream GROMACS `main`
#   at one commit, merges the newest upstream release branch into it, adds
#   those MRs on top in a fixed order, and produces
#
#     * a local branch (default `lab/main`) holding the combination, and
#     * a manifest recording exactly what went in: the base commit, the
#       release-branch commit merged into it, and for every MR its number,
#       the head commit used, its state, how it went in, and its title.
#
#   The last commit on the branch sets GMX_VERSION_STRING_OF_FORK, so every
#   binary built from it reports a version like
#   "2028.0-smithlab_integration_20261123_bc88adaf-dev" (version, lab tag,
#   date of the run, upstream main commit). That string can never be
#   mistaken for an upstream release. The manifest is committed into the
#   tree as lab-integration-manifest.txt and a copy, with the resulting
#   commit added, goes to manifests/.
#
#   MR heads are fetched from the upstream project itself
#   (refs/merge-requests/<iid>/head on gitlab.com/gromacs/gromacs), so the
#   script does not depend on anyone's personal fork.
#
# USAGE
#   ./integrate.sh run      fetch, combine, verify, write branch + manifest.
#                           Exit 0 = done, 2 = error, 3 = conflict.
#   ./integrate.sh check    report only: has upstream main or the release
#                           branch moved, has any MR head moved or changed
#                           state, since the last manifest, or is the
#                           manifest older than MAX_AGE_DAYS? Needs no
#                           clone. Exit 0 = nothing to do, 1 = re-run `run`,
#                           2 = error. Meant to be called from a scheduled
#                           job.
#   ./integrate.sh publish  push the last result to the lab repository:
#                           `main` = the upstream base (fast-forward only),
#                           `lab/main` = the combined branch (replaced, since
#                           every run rebuilds it), `integration` = this
#                           directory (script, MR list, manifests).
#
# WHAT `run` BUILDS, IN ORDER
#   1. upstream main (BASE_COMMIT, default its tip).
#   2. "Merge release-YYYY into main": the newest upstream release branch,
#      merged the way upstream does its own merge-ups, so fixes made on the
#      release branch are in the lab build before upstream merges them up.
#      Skipped when main already contains it. A conflict stops the run
#      (upstream resolves those in its own merge-up; wait for it).
#   3. "Keep main's release-only settings ...", only if step 2 changed them:
#      the release branch sets GMX_VERSION_SUFFIX (e.g. "-beta") and
#      REGRESSIONTEST_BRANCH for its own releases. Upstream keeps main's
#      values when it merges up, and so does this script.
#   4. One merge (or carried commits, see below) per MR series.
#   5. "Lab build: declare GROMACS <N> ...", only while upstream main still
#      declares an older version: once release-YYYY is branched, main is
#      the development branch for YYYY+1, but upstream changes
#      project(Gromacs VERSION ...) some time later. This commit does in
#      CMake what upstream's "Prepare main for <N> release" commit does:
#      the project version and the libgromacs SOVERSION major. The run
#      prints, and the manifest records, when upstream main declares the
#      version itself and the commit is no longer made.
#   6. "Lab integration build <version>": GMX_VERSION_STRING_OF_FORK and
#      lab-integration-manifest.txt.
#
# BRANCHES IN THE LAB REPOSITORY
#   integration  this script, mr-list.txt, manifests/, README.md. Its own
#                history, unrelated to GROMACS. The repository's default.
#   lab/main     the combined GROMACS to build and install.
#   main         a plain mirror of upstream main, never committed to.
#
#   Settings (environment variables, all optional):
#     GMX_REPO       clone to work in (default: ./gromacs next to this
#                    script; cloned on first use)
#     REFERENCE_REPO an existing GROMACS clone to borrow objects from when
#                    cloning, to save download time (the clone is still
#                    self-contained afterwards)
#     BASE_COMMIT    pin the upstream base (default: tip of upstream main)
#     RELEASE_BRANCH release branch merged in step 2 (default: the newest
#                    upstream release-YYYY; `none` to skip step 2)
#     LAB_VERSION_MAJOR  version declared in step 5 (default: the year of
#                    RELEASE_BRANCH plus one; empty to skip step 5)
#     OUTPUT_BRANCH  default lab/main
#     LAB_TAG        default smithlab
#     LAB_REPO_URL   where `publish` pushes (default below)
#     MAX_AGE_DAYS   `check` asks for a re-run after this many days (31)
#
# THE MR LIST
#   mr-list.txt, next to this script, is the only place that says which MRs
#   go in and in what order. See its header. To add an MR, add a line; to
#   drop one, delete its line. A closed (not merged) MR stops the run:
#   decide whether to drop it from the list.
#
#   Stacked MRs form a series. Only the last MR of a series is merged; the
#   earlier ones must be inside it. Each one is still recorded and checked:
#   either its head is an ancestor of the result, or (when the MR branch on
#   GitLab is older than the rebased copy inside the later MR) all of its
#   commits are present as identical patches. The manifest says which.
#
#   An MR that targets main and is merged is in the base, and is skipped.
#
#   An MR whose target is a release branch (e.g. release-2026) is carried
#   forward by cherry-picking its own commits one by one, each with a
#   "(cherry picked from commit ...)" line: while it is open, the commits of
#   its head (<target>..<head>); once merged, the commit(s) it was merged as
#   on the release branch (the squash commit, the commits of the merge, or
#   the fast-forwarded commits). Each commit already in the base (an
#   ancestor of it, or the same patch in a commit of it, e.g. after
#   upstream merged the release branch up) is not carried; when that is
#   true of all of them, the MR is skipped and the manifest says so. Their
#   changes under docs/release-notes/ are NOT carried: those notes describe
#   the release branch, the files differ between branches, and main
#   receives them when upstream merges the release branch into main. Every
#   carried commit is verified to be patch-identical to the original
#   outside that directory.
#
# ON A CONFLICT
#   The script stops, names the MR being applied, the MR(s) already applied
#   that touched the same files, and the files. It resolves nothing and
#   leaves the output branch and the previous manifest untouched. Do not
#   resolve the conflict here. Ask the owners of the two MRs to make them
#   compatible upstream (usually by rebasing one MR onto the other), then
#   re-run. If one MR must wait, comment it out in mr-list.txt. A conflict
#   in step 2 is between upstream's own branches: wait for upstream's
#   merge-up, or run with RELEASE_BRANCH=none.
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
LAB_REPO_URL=${LAB_REPO_URL:-https://github.com/Smith-Group/gromacs-smithlab.git}
MAX_AGE_DAYS=${MAX_AGE_DAYS:-31}
IN_TREE_MANIFEST=lab-integration-manifest.txt
CARRY_EXCLUDE=docs/release-notes   # not carried forward from release-branch MRs
VERSION_FILE=cmake/gmxVersionInfo.cmake
RELEASE_ONLY_SETTINGS="GMX_VERSION_SUFFIX REGRESSIONTEST_BRANCH"   # step 3

die()  { echo "ERROR: $*" >&2; exit 2; }
note() { echo "== $*"; }

# Last "key=value" value for key $1 in manifest file $2.
manifest_value() { awk -F= -v k="$1" '$1 == k { v = substr($0, length(k) + 2) } END { print v }' "$2"; }

# --- MR list -----------------------------------------------------------------

# Prints "series iid" per MR, in file order.
read_list() {
    [ -f "$MR_LIST" ] || die "missing $MR_LIST"
    sed -e 's/#.*//' "$MR_LIST" | awk 'NF >= 2 { print $1, $2 }'
}

# Prints the series names in order of first appearance.
series_order() { read_list | awk '!seen[$1]++ { print $1 }'; }

# --- upstream branches -------------------------------------------------------

# The newest upstream release-YYYY branch.
newest_release() {
    git ls-remote --heads "$UPSTREAM_URL" 'refs/heads/release-*' \
        | sed -n 's|.*refs/heads/\(release-[0-9][0-9][0-9][0-9]\)$|\1|p' | sort | tail -n 1
}

resolve_release_settings() {
    RELEASE_BRANCH=${RELEASE_BRANCH-$(newest_release)}
    [ -n "$RELEASE_BRANCH" ] || die "found no upstream release-YYYY branch; set RELEASE_BRANCH=none to skip"
    if [ -z "${LAB_VERSION_MAJOR+set}" ]; then
        case $RELEASE_BRANCH in
            release-[0-9][0-9][0-9][0-9]) LAB_VERSION_MAJOR=$(( ${RELEASE_BRANCH#release-} + 1 )) ;;
            *) LAB_VERSION_MAJOR= ;;
        esac
    fi
}

# --- GitLab API (public, read-only) ------------------------------------------

# Prints, tab-separated, for one MR:
#   state  target_branch  head_sha  content_from  content_to  title
# content_from..content_to are the commits the MR was merged as, on its
# target branch ("-" while it is not merged): the squash commit, the
# commits of the merge commit, or the fast-forwarded commits.
mr_info() {
    curl -sf "$API/merge_requests/$1" | python3 -I -c '
import json, sys
m = json.load(sys.stdin)
cfrom = cto = "-"
if m["state"] == "merged":
    if m.get("squash_commit_sha"):
        cfrom, cto = m["squash_commit_sha"] + "^", m["squash_commit_sha"]
    elif m.get("merge_commit_sha"):
        cfrom, cto = m["merge_commit_sha"] + "^1", m["merge_commit_sha"] + "^2"
    else:
        cfrom, cto = m["diff_refs"]["base_sha"], m["sha"]
print("\t".join([m["state"], m["target_branch"], m["sha"], cfrom, cto, " ".join(m["title"].split())]))
' || die "could not read !$1 from the GitLab API"
}

# --- check mode --------------------------------------------------------------

do_check() {
    local last=$MANIFEST_DIR/latest.txt
    [ -f "$last" ] || { echo "No previous manifest ($last): run './integrate.sh run' first."; exit 1; }
    local moved=0 old new
    old=$(manifest_value base_commit "$last")
    new=$(git ls-remote "$UPSTREAM_URL" "refs/heads/$BASE_BRANCH" | cut -f1)
    [ -n "$new" ] || die "cannot read upstream $BASE_BRANCH"
    if [ "$old" != "$new" ]; then
        echo "upstream $BASE_BRANCH moved: ${old:0:10} -> ${new:0:10}"; moved=1
    else
        echo "upstream $BASE_BRANCH unchanged at ${old:0:10}"
    fi
    resolve_release_settings
    local oldrb oldrc newrc
    oldrb=$(manifest_value release_branch "$last")
    oldrc=$(manifest_value release_commit "$last")
    if [ "$RELEASE_BRANCH" = none ]; then
        [ "${oldrb:-none}" = none ] || { echo "release branch merge is now switched off"; moved=1; }
    elif [ "$RELEASE_BRANCH" != "$oldrb" ]; then
        echo "release branch to merge is now $RELEASE_BRANCH (was ${oldrb:-none})"; moved=1
    else
        newrc=$(git ls-remote "$UPSTREAM_URL" "refs/heads/$RELEASE_BRANCH" | cut -f1)
        [ -n "$newrc" ] || die "cannot read upstream $RELEASE_BRANCH"
        if [ "$oldrc" != "$newrc" ]; then
            echo "upstream $RELEASE_BRANCH moved: ${oldrc:0:10} -> ${newrc:0:10}"; moved=1
        else
            echo "upstream $RELEASE_BRANCH unchanged at ${oldrc:0:10}"
        fi
    fi
    local generated age
    generated=$(manifest_value generated "$last")
    age=$(( ($(date -u +%s) - $(date -u -d "$generated" +%s)) / 86400 ))
    if [ "$age" -ge "$MAX_AGE_DAYS" ]; then
        echo "last manifest is $age days old (limit $MAX_AGE_DAYS)"; moved=1
    fi
    local series iid state target sha cfrom cto title oldsha oldstate
    while read -r series iid; do
        IFS=$'\t' read -r state target sha cfrom cto title < <(mr_info "$iid")
        oldsha=$(awk -F'\t' -v i="$iid" '$1 == "mr" && $2 == i { print $3 }' "$last")
        oldstate=$(awk -F'\t' -v i="$iid" '$1 == "mr" && $2 == i { print $4 }' "$last")
        if [ -z "$oldsha" ]; then
            echo "!$iid is new in mr-list.txt ($state)"; moved=1
        elif [ "$state" != "$oldstate" ]; then
            echo "!$iid is now $state upstream (was $oldstate)"; moved=1
        elif [ "$sha" != "$oldsha" ]; then
            echo "!$iid head moved: ${oldsha:0:10} -> ${sha:0:10}  ($title)"; moved=1
        else
            echo "!$iid unchanged ($state) at ${sha:0:10}"
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

# The first MR of a series (in list order) whose head contains commit $1.
owner_of() {
    local c=$1; shift
    local i
    for i in "$@"; do
        if G merge-base --is-ancestor "$c" "refs/lab/mr/$i" 2>/dev/null; then echo "!$i"; return; fi
    done
    echo "(unknown MR)"
}

# Reports a conflict while applying a series, puts the repository back, exits.
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
            e=$(G log -1 --format=%H "$PREPARED..${APPLIED_TIP[$s]}" -- "$f" 2>/dev/null || true)
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

# Patch-id of a commit, ignoring CARRY_EXCLUDE. Empty if it only touched that.
patch_id() { G diff "$1^" "$1" -- . ":(exclude)$CARRY_EXCLUDE" | G patch-id --stable | cut -d' ' -f1; }

# Is commit $1 already in the prepared base (steps 1-3)? Yes if it is an
# ancestor, or if a commit of the base that is not in $2 (the side the
# commit was made on) is the same patch, ignoring CARRY_EXCLUDE.
commit_in_base() {
    local c=$1 side=$2 pid d files
    G merge-base --is-ancestor "$c" "$PREPARED" && return 0
    pid=$(patch_id "$c")
    files=$(G diff --name-only "$c^" "$c" -- . ":(exclude)$CARRY_EXCLUDE")
    [ -n "$pid" ] && [ -n "$files" ] || return 1
    # shellcheck disable=SC2086
    for d in $(G rev-list --no-merges "$PREPARED" "^$side" -- $files); do
        [ "$(patch_id "$d")" = "$pid" ] && return 0
    done
    return 1
}

# The carried copy (in PREPARED..HEAD) of original commit $1, if any.
carried_copy() { G log -1 --format=%H --fixed-strings --grep="(cherry picked from commit $1)" "$PREPARED..HEAD"; }

# Counts, for release-target MR $1, its commits FROM..TO: sets N_ALL,
# N_NOTES (only touch CARRY_EXCLUDE), N_INBASE and N_PENDING (to carry).
count_carry() {
    local iid=$1 c
    N_ALL=0 N_NOTES=0 N_INBASE=0 N_PENDING=0
    for c in $(G rev-list "${FROM[$iid]}..${TO[$iid]}"); do
        N_ALL=$((N_ALL + 1))
        if [ -z "$(patch_id "$c")" ]; then N_NOTES=$((N_NOTES + 1))
        elif commit_in_base "$c" "${FROM[$iid]}"; then N_INBASE=$((N_INBASE + 1))
        else N_PENDING=$((N_PENDING + 1)); fi
    done
}

# Cherry-picks the commits FROM..TO of release-target MR $1 onto HEAD,
# leaving out CARRY_EXCLUDE and commits already in the base. Returns 1 on
# a conflict, with the cherry-pick still in progress.
carry_forward() {
    local iid=$1 c msg
    msg=$(mktemp)
    for c in $(G rev-list --reverse "${FROM[$iid]}..${TO[$iid]}"); do
        [ -n "$(patch_id "$c")" ] || continue           # only touches CARRY_EXCLUDE
        commit_in_base "$c" "${FROM[$iid]}" && continue  # already in the base
        G cherry-pick --no-commit "$c" >/dev/null 2>&1 || { rm -f "$msg"; return 1; }
        if G ls-files --error-unmatch -- "$CARRY_EXCLUDE" >/dev/null 2>&1 \
           || ! G diff --cached --quiet -- "$CARRY_EXCLUDE"; then
            G restore --staged --worktree --source=HEAD -- ":(glob)$CARRY_EXCLUDE/**"
        fi
        # Nothing left means the change is already here in another form; the
        # verification below then reports the commit as not contained.
        if G diff --cached --quiet; then continue; fi
        { G log -1 --format=%B "$c"; echo "(cherry picked from commit $c)"; } > "$msg"
        G commit -q -F "$msg" --author="$(G log -1 --format='%an <%ae>' "$c")" \
            --date="$(G log -1 --format=%aD "$c")"
    done
    rm -f "$msg"
}

# Step 2: merge RELEASE_BRANCH into the base. Sets RELEASE_COMMIT, RELEASE_HOW.
merge_release_branch() {
    RELEASE_COMMIT= RELEASE_HOW="not merged (RELEASE_BRANCH=none)"
    [ "$RELEASE_BRANCH" != none ] || { note "release branch merge switched off"; return 0; }
    RELEASE_COMMIT=$(G rev-parse "origin/$RELEASE_BRANCH")
    if G merge-base --is-ancestor "$RELEASE_COMMIT" HEAD; then
        RELEASE_HOW="already contained in upstream $BASE_BRANCH"
        note "$RELEASE_BRANCH ${RELEASE_COMMIT:0:10} is already in $BASE_BRANCH"
        return 0
    fi
    local n; n=$(G rev-list --count "HEAD..$RELEASE_COMMIT")
    note "merging $RELEASE_BRANCH ${RELEASE_COMMIT:0:10} into $BASE_BRANCH ($n commits $BASE_BRANCH lacks)"
    if ! G merge -q --no-ff --no-edit -m "Merge $RELEASE_BRANCH into $BASE_BRANCH" \
            -m "Made by integrate.sh for the lab build, the way upstream merges its release branch into $BASE_BRANCH; $RELEASE_BRANCH at $RELEASE_COMMIT." \
            "$RELEASE_COMMIT" >/dev/null 2>&1; then
        echo
        echo "CONFLICT while merging $RELEASE_BRANCH ${RELEASE_COMMIT:0:10} into $BASE_BRANCH ${BASE:0:10}."
        echo "Nothing was resolved, and $OUTPUT_BRANCH and the manifests are unchanged."
        G diff --name-only --diff-filter=U | sed 's/^/  /'
        echo "This is between upstream's own branches. Wait for upstream to merge"
        echo "$RELEASE_BRANCH into $BASE_BRANCH, or run with RELEASE_BRANCH=none."
        G merge --abort || true
        G switch -q --detach "$BASE"
        exit 3
    fi
    RELEASE_HOW="merged cleanly as $(G rev-parse --short=10 HEAD) ($n commits)"
}

# Step 3: put back main's values of RELEASE_ONLY_SETTINGS. Sets SETTINGS_KEPT.
keep_main_release_settings() {
    SETTINGS_KEPT=none
    local key ours theirs kept="" body=""
    for key in $RELEASE_ONLY_SETTINGS; do
        ours=$(G show "$BASE:$VERSION_FILE" | grep "^set($key " || true)
        theirs=$(grep "^set($key " "$GMX_REPO/$VERSION_FILE" || true)
        [ -n "$ours" ] && [ "$(grep -c . <<<"$ours")" = 1 ] && [ "$(grep -c . <<<"$theirs")" = 1 ] \
            || die "expected one 'set($key ...' line in $VERSION_FILE"
        [ "$ours" = "$theirs" ] && continue
        awk -v k="set($key " -v l="$ours" 'index($0, k) == 1 { print l; next } { print }' \
            "$GMX_REPO/$VERSION_FILE" > "$GMX_REPO/$VERSION_FILE.new"
        mv "$GMX_REPO/$VERSION_FILE.new" "$GMX_REPO/$VERSION_FILE"
        kept+=" $key"
        body+="  $theirs  ->  $ours"$'\n'
    done
    [ -n "$kept" ] || return 0
    SETTINGS_KEPT=${kept# }
    note "keeping $BASE_BRANCH's values of $SETTINGS_KEPT"
    G add "$VERSION_FILE"
    G commit -q -m "Keep main's release-only settings after merging $RELEASE_BRANCH" \
        -m "$RELEASE_BRANCH sets these for its own releases. Upstream keeps main's values when it merges a release branch into main, and so does the lab build:" \
        -m "$body"
}

# Step 5: declare LAB_VERSION_MAJOR while upstream main declares less.
# Sets VERSION_HOW.
declare_lab_version() {
    local declared so newso
    declared=$(G show "$BASE:CMakeLists.txt" | sed -n 's/^project(Gromacs VERSION \([0-9]*\)\.[0-9]*)$/\1/p')
    [ -n "$declared" ] || die "cannot read project(Gromacs VERSION ...) from CMakeLists.txt"
    if [ -z "$LAB_VERSION_MAJOR" ]; then
        VERSION_HOW="$declared (as declared by upstream $BASE_BRANCH; LAB_VERSION_MAJOR empty)"; return 0
    fi
    if [ "$declared" -ge "$LAB_VERSION_MAJOR" ]; then
        VERSION_HOW="$declared (declared by upstream $BASE_BRANCH itself; the lab no longer sets it)"
        note "upstream $BASE_BRANCH declares $declared itself: not making the lab version commit any more"
        return 0
    fi
    [ "$LAB_VERSION_MAJOR" = $((declared + 1)) ] \
        || die "upstream $BASE_BRANCH declares $declared, LAB_VERSION_MAJOR is $LAB_VERSION_MAJOR: expected $((declared + 1))"
    so=$(sed -n 's/^set(LIBRARY_SOVERSION_MAJOR \([0-9]*\))$/\1/p' "$GMX_REPO/$VERSION_FILE")
    [ -n "$so" ] || die "cannot read LIBRARY_SOVERSION_MAJOR from $VERSION_FILE"
    newso=$((so + 1))
    note "declaring GROMACS $LAB_VERSION_MAJOR (upstream $BASE_BRANCH still says $declared), SOVERSION $so -> $newso"
    sed -i "s/^project(Gromacs VERSION $declared\.0)$/project(Gromacs VERSION $LAB_VERSION_MAJOR.0)/" \
        "$GMX_REPO/CMakeLists.txt"
    sed -i -e "s/^set(LIBRARY_SOVERSION_MAJOR $so)$/set(LIBRARY_SOVERSION_MAJOR $newso)/" \
        -e "/^#  *GROMACS  *$declared  *$so\$/a #         GROMACS     $LAB_VERSION_MAJOR   $newso" \
        "$GMX_REPO/$VERSION_FILE"
    grep -q "^project(Gromacs VERSION $LAB_VERSION_MAJOR.0)$" "$GMX_REPO/CMakeLists.txt" \
        && grep -q "^set(LIBRARY_SOVERSION_MAJOR $newso)$" "$GMX_REPO/$VERSION_FILE" \
        && grep -q "^#  *GROMACS  *$LAB_VERSION_MAJOR  *$newso\$" "$GMX_REPO/$VERSION_FILE" \
        || die "could not declare version $LAB_VERSION_MAJOR"
    G add CMakeLists.txt "$VERSION_FILE"
    G commit -q -m "Lab build: declare GROMACS $LAB_VERSION_MAJOR (upstream main still says $declared)" \
        -m "Upstream main has been the development branch for $LAB_VERSION_MAJOR since release-$declared was branched, but still declares $declared. This makes the CMake changes of upstream's own \"Prepare main for <version> release\" commit: the project version, and the libgromacs SOVERSION major ($so -> $newso), so a lab build cannot be mistaken for $declared or be loaded in place of its libraries. Release-note scaffolding is not added." \
        -m "Not an upstream GROMACS commit. integrate.sh stops making it once upstream main declares $LAB_VERSION_MAJOR."
    VERSION_HOW="$LAB_VERSION_MAJOR (set by this build: upstream $BASE_BRANCH still declares $declared; SOVERSION $so -> $newso)"
}

do_run() {
    ensure_repo
    command -v curl >/dev/null || die "curl is needed"
    resolve_release_settings

    # Read the list and the state of every MR.
    declare -gA SERIES_MRS STATE TARGET SHA CFROM CTO TITLE FROM TO APPLIED_TIP HOW
    declare -ga APPLIED_SERIES=()
    local series iid
    while read -r series iid; do
        SERIES_MRS[$series]="${SERIES_MRS[$series]:-} $iid"
        IFS=$'\t' read -r "STATE[$iid]" "TARGET[$iid]" "SHA[$iid]" "CFROM[$iid]" "CTO[$iid]" "TITLE[$iid]" \
            < <(mr_info "$iid")
    done < <(read_list)

    # Fetch base branch, the release branches, and every MR head.
    local refspecs=() targets t
    targets=$( { for iid in "${!TARGET[@]}"; do echo "${TARGET[$iid]}"; done
                 [ "$RELEASE_BRANCH" = none ] || echo "$RELEASE_BRANCH"; } | sort -u)
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
            opened) FROM[$iid]=origin/${TARGET[$iid]}; TO[$iid]=refs/lab/mr/$iid ;;
            merged) FROM[$iid]=${CFROM[$iid]}; TO[$iid]=${CTO[$iid]}
                    G rev-parse -q --verify "${TO[$iid]}^{commit}" >/dev/null \
                        || die "!$iid: merged as ${TO[$iid]}, which is not on upstream ${TARGET[$iid]}" ;;
            *) die "!$iid is ${STATE[$iid]} upstream. Drop it from mr-list.txt or find out why." ;;
        esac
    done

    # Steps 1-3: the prepared base.
    G switch -q --detach "$BASE"
    merge_release_branch
    keep_main_release_settings
    PREPARED=$(G rev-parse HEAD)

    # Step 4: apply each series.
    local mrs tip
    for series in $(series_order); do
        mrs=${SERIES_MRS[$series]}
        tip=""
        for iid in $mrs; do
            if [ "${TARGET[$iid]}" = "$BASE_BRANCH" ]; then
                if [ "${STATE[$iid]}" = merged ]; then
                    G merge-base --is-ancestor "${TO[$iid]}" "$PREPARED" \
                        || die "!$iid was merged upstream after the base ${BASE:0:10}; do not pin BASE_COMMIT before it"
                    HOW[$iid]="skipped: merged upstream"; note "!$iid is merged upstream, skipped"
                else
                    tip=$iid
                fi
                continue
            fi
            count_carry "$iid"
            if [ "$N_PENDING" = 0 ]; then
                HOW[$iid]="skipped: already in the base ($N_INBASE of $N_ALL commits in $BASE_BRANCH${RELEASE_COMMIT:++$RELEASE_BRANCH}, $N_NOTES touched only $CARRY_EXCLUDE)"
                note "!$iid (${STATE[$iid]} into ${TARGET[$iid]}) is already in the base, skipped"
            else
                tip=$iid
            fi
        done
        if [ -z "$tip" ]; then note "series $series: nothing left to apply"; continue; fi

        if [ "${TARGET[$tip]}" = "$BASE_BRANCH" ]; then
            note "series $series: merging tip !$tip ${SHA[$tip]:0:10}"
            G merge -q --no-ff --no-edit -m "Merge !$tip: ${TITLE[$tip]}" \
                -m "Series '$series' of the lab integration; MR head ${SHA[$tip]}." "refs/lab/mr/$tip" \
                >/dev/null 2>&1 || conflict "merge of !$tip" "$series" "$PREPARED" "refs/lab/mr/$tip"
            HOW[$tip]="merged (series tip)"
        else
            local source
            if [ "${STATE[$tip]}" = merged ]; then
                source="its merged commits $(G rev-parse --short=10 "${FROM[$tip]}")..$(G rev-parse --short=10 "${TO[$tip]}") on ${TARGET[$tip]}"
            else
                source="its head, from ${TARGET[$tip]}"
            fi
            count_carry "$tip"
            note "series $series: carrying !$tip forward from $source ($N_PENDING of $N_ALL commits)"
            carry_forward "$tip" \
                || conflict "cherry-pick of !$tip from ${TARGET[$tip]}" "$series" "${FROM[$tip]}" "${TO[$tip]}"
            HOW[$tip]="cherry-picked from $source: $N_PENDING of $N_ALL commits carried, $N_INBASE already in the base, $N_NOTES touched only $CARRY_EXCLUDE"
        fi
        APPLIED_SERIES+=("$series")
        APPLIED_TIP[$series]=$(G rev-parse "${TO[$tip]}")
    done

    # Verify that every MR not skipped is contained in the result.
    local missing=0 from unmatched
    for series in $(series_order); do
        for iid in ${SERIES_MRS[$series]}; do
            case ${HOW[$iid]:-} in skipped:*) continue ;; esac
            if [ "${TARGET[$iid]}" != "$BASE_BRANCH" ]; then
                local c copy bad=0
                for c in $(G rev-list "${FROM[$iid]}..${TO[$iid]}"); do
                    [ -n "$(patch_id "$c")" ] || continue   # only touched CARRY_EXCLUDE
                    copy=$(carried_copy "$c")
                    if [ -n "$copy" ]; then
                        [ "$(patch_id "$c")" = "$(patch_id "$copy")" ] \
                            || { echo "DIFFERS: !$iid commit ${c:0:10} and its carried copy ${copy:0:10} are not the same patch."; bad=1; }
                    elif ! commit_in_base "$c" "${FROM[$iid]}"; then
                        echo "NOT CONTAINED: !$iid commit ${c:0:10} was not carried forward."; bad=1
                    fi
                done
                [ $bad = 0 ] || missing=1
                HOW[$iid]=${HOW[$iid]:-contained (its commits are carried copies or already in the base)}
                continue
            fi
            if G merge-base --is-ancestor "refs/lab/mr/$iid" HEAD; then
                HOW[$iid]=${HOW[$iid]:-contained (head is an ancestor)}
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

    # Step 5: the version this build declares.
    declare_lab_version

    # Step 6: version string and manifest, committed as the last commit.
    local date stamp version
    date=$(date -u +%Y%m%d)
    stamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    version="${LAB_TAG}_integration_${date}_$(G rev-parse --short=8 "$BASE")"
    sed -i "s/^set(GMX_VERSION_STRING_OF_FORK \"[^\"]*\"/set(GMX_VERSION_STRING_OF_FORK \"$version\"/" \
        "$GMX_REPO/$VERSION_FILE"
    grep -q "GMX_VERSION_STRING_OF_FORK \"$version\"" "$GMX_REPO/$VERSION_FILE" \
        || die "could not set GMX_VERSION_STRING_OF_FORK"
    {
        echo "# GROMACS lab integration manifest. Generated by integrate.sh; do not edit."
        echo "version_string_of_fork=$version"
        echo "version_major=$VERSION_HOW"
        echo "generated=$stamp"
        echo "upstream=$UPSTREAM_URL"
        echo "base_branch=$BASE_BRANCH"
        echo "base_commit=$(G rev-parse "$BASE")"
        echo "base_subject=$(G log -1 --format=%s "$BASE")"
        echo "release_branch=$RELEASE_BRANCH"
        echo "release_commit=$RELEASE_COMMIT"
        [ -z "$RELEASE_COMMIT" ] || echo "release_subject=$(G log -1 --format=%s "$RELEASE_COMMIT")"
        echo "release_merge=$RELEASE_HOW"
        echo "release_only_settings_kept_from_main=$SETTINGS_KEPT"
        echo "# mr<TAB>iid<TAB>head<TAB>state<TAB>target<TAB>series<TAB>how<TAB>title"
        while read -r series iid; do
            printf 'mr\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$iid" "${SHA[$iid]}" "${STATE[$iid]}" \
                "${TARGET[$iid]}" "$series" "${HOW[$iid]}" "${TITLE[$iid]}"
        done < <(read_list)
    } > "$GMX_REPO/$IN_TREE_MANIFEST"
    G add "$VERSION_FILE" "$IN_TREE_MANIFEST"
    G commit -q -m "Lab integration build $version" \
        -m "Not an upstream GROMACS release. Contents: $IN_TREE_MANIFEST."
    local result; result=$(G rev-parse HEAD)

    G branch -f "$OUTPUT_BRANCH" "$result"
    G switch -q "$OUTPUT_BRANCH"
    mkdir -p "$MANIFEST_DIR"
    local out=$MANIFEST_DIR/manifest-$date-${result:0:10}.txt
    { cat "$GMX_REPO/$IN_TREE_MANIFEST"; echo "result_branch=$OUTPUT_BRANCH"; echo "result_commit=$result"; } > "$out"
    cp "$out" "$MANIFEST_DIR/latest.txt"
    if git -C "$HERE" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        git -C "$HERE" add "$out" "$MANIFEST_DIR/latest.txt"
        git -C "$HERE" commit -q -m "Record the integration of $(date -u +%Y-%m-%d) ($version)" \
            -m "lab/main ${result:0:10} on upstream $BASE_BRANCH $(G rev-parse --short=10 "$BASE")${RELEASE_COMMIT:+ + $RELEASE_BRANCH ${RELEASE_COMMIT:0:10}}."
    fi

    echo
    note "done: $OUTPUT_BRANCH = ${result:0:10} in $GMX_REPO"
    note "version string: $version"
    note "version: $VERSION_HOW"
    note "$RELEASE_BRANCH: $RELEASE_HOW"
    note "manifest: $out"
    awk -F'\t' '$1 == "mr" { printf "   !%s %s  %s\n", $2, substr($3, 1, 10), $7 }' "$out"
}

# --- publish mode ------------------------------------------------------------

do_publish() {
    local last=$MANIFEST_DIR/latest.txt base result head remote_lab
    [ -f "$last" ] || die "nothing to publish: no manifest yet"
    base=$(manifest_value base_commit "$last")
    result=$(manifest_value result_commit "$last")
    head=$(G rev-parse "$OUTPUT_BRANCH")
    [ "$head" = "$result" ] || die "$OUTPUT_BRANCH is ${head:0:10}, the manifest says ${result:0:10}; re-run"
    [ -z "$(git -C "$HERE" status --porcelain --untracked-files=no)" ] \
        || die "uncommitted changes in $HERE; commit them first"
    note "publishing to $LAB_REPO_URL"
    note "main     -> ${base:0:10} (fast-forward only)"
    G push -q "$LAB_REPO_URL" "$base:refs/heads/main"
    remote_lab=$(git ls-remote "$LAB_REPO_URL" "refs/heads/$OUTPUT_BRANCH" | cut -f1)
    note "$OUTPUT_BRANCH -> ${result:0:10} (replacing ${remote_lab:0:10})"
    G push -q --force-with-lease="refs/heads/$OUTPUT_BRANCH:$remote_lab" "$LAB_REPO_URL" \
        "$result:refs/heads/$OUTPUT_BRANCH"
    note "integration -> $(git -C "$HERE" rev-parse --short=10 HEAD)"
    git -C "$HERE" push -q "$LAB_REPO_URL" "HEAD:refs/heads/integration"
}

case ${1:-} in
    run)     do_run ;;
    check)   do_check ;;
    publish) do_publish ;;
    *)     sed -n '2,/^set -euo/p' "$0" | sed -e '$d' -e 's/^# \{0,1\}//'; exit 2 ;;
esac
