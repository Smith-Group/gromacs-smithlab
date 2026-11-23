# gromacs-smithlab: unofficial GROMACS with unmerged merge requests

**This is not an official GROMACS release.** The `lab/main` branch is upstream
GROMACS `main` with a set of merge requests (MRs) added that are still under
review at <https://gitlab.com/gromacs/gromacs>. The Smith group uses it until
those MRs are merged upstream. Please report problems with it to us, not to the
GROMACS developers.

## Branches

| Branch | Contents |
|---|---|
| `integration` (this one) | The integration script, the MR list and the manifests. |
| `lab/main` | The combined GROMACS. Build and install this one. |
| `main` | A plain mirror of upstream GROMACS `main`. Never committed to. |

## What is in `lab/main`

In order:

1. upstream GROMACS `main`;
2. the newest upstream release branch (now `release-2027`) merged into it, the
   way upstream does its own merge-ups, so release-branch fixes are included
   before upstream merges them up; the release branch's own version suffix and
   regression-test branch are not taken;
3. the MRs in [`mr-list.txt`](mr-list.txt). MRs for a release branch are
   carried as cherry-picks until their content reaches the base, and are then
   skipped automatically;
4. while upstream `main` still declares the previous version, a commit
   declaring the version `main` is developing (now 2028, with the matching
   library SOVERSION). It disappears once upstream declares it;
5. the version string and the manifest.

[`manifests/latest.txt`](manifests/latest.txt) lists the upstream base commit,
the release-branch commit merged into it and, for every MR, its number, the
exact head commit used, how it went in and its title. The same file is in the
root of `lab/main` as `lab-integration-manifest.txt`.

A build reports a version like `2028.0-smithlab_integration_20261123_bc88adaf-dev`
(`gmx -version`, and the header of every `md.log`), so a result can always be
traced back to its manifest.

## Building

```
git clone -b lab/main https://github.com/Smith-Group/gromacs-smithlab.git
```

then build as described in the GROMACS installation guide.

## Updating

`integrate.sh` rebuilds `lab/main` from scratch. See its header for details.

```
./integrate.sh check     # has anything moved since the last manifest?
./integrate.sh run       # rebuild lab/main and write a new manifest
./integrate.sh publish   # push main, lab/main and integration here
```

Run `check` when an MR changes, and at least monthly to move onto newer upstream
`main` and the release branch. To add or drop an MR, edit
[`mr-list.txt`](mr-list.txt). If `run` stops on a conflict, it changes nothing.
Ask the owners of the two MRs to make them compatible upstream, or comment one
of them out, then run it again.
