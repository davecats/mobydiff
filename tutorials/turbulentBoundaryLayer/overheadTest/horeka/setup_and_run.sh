#!/bin/bash
# One-shot setup + submit for the 2:1 refinement timing/profiling campaign.
#
# RUN THIS ON A HOREKA LOGIN NODE, from inside the transferred package dir:
#     bash setup_and_run.sh
#
# It (1) builds parallel HDF5 if needed, (2) clones the solver at the EXACT
# commit into a SEPARATE code directory, (3) copies the run data and submit
# scripts into a SEPARATE run directory, and (4) submits the SLURM job.
# Code and run directories are kept apart so the clone stays pristine and
# several campaigns can share one build.
#
# Nothing here needs a restart field: the boundaryLayer case is analytic and
# cold-starts from the ini, so the whole package is a few hundred kB.
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration -- edit for your project space.
# ---------------------------------------------------------------------------
REPO_URL="${REPO_URL:-git@github.com:davecats/mobydiff.git}"
BRANCH="${BRANCH:-optimiseBlockRefinement_parentBoundaryLayer}"
# The commit under test. PIN IT: the cluster must run exactly the validated
# code, independent of any branch movement afterwards.
COMMIT="${COMMIT:-a11e355a47e1535db4f3e9bae0dcc489eaec3567}"

# Two SEPARATE directories, neither inside the other. Put them on the parallel
# WORKSPACE, not $HOME: $HOME is quota-limited and slow for MPI-IO.
#   e.g. WS=/hkfs/work/workspace/scratch/<user>-<ws>
WS="${WS:-$HOME}"
CODE_DIR="${CODE_DIR:-$WS/moby-2to1-code}"
RUN_DIR="${RUN_DIR:-$WS/moby-2to1-run}"
HDF5_DIR="${HDF5_DIR:-$HOME/hdf5}"

CASE_DIR="$(cd "$(dirname "$0")" && pwd -P)"

echo "REPO      : $REPO_URL"
echo "BRANCH    : $BRANCH"
echo "COMMIT    : $COMMIT"
echo "CODE_DIR  : $CODE_DIR"
echo "RUN_DIR   : $RUN_DIR"
echo "HDF5_DIR  : $HDF5_DIR"
echo "CASE_DIR  : $CASE_DIR"
[ "$WS" = "$HOME" ] && echo "NOTE: WS is \$HOME. Set WS=/hkfs/work/workspace/scratch/<user>-<ws> to use the workspace."
echo

# ---------------------------------------------------------------------------
# 0. Parallel HDF5 (once), with the nvhpc compilers.
# ---------------------------------------------------------------------------
module load toolkit/nvidia-hpc-sdk/25.3
bash "$CASE_DIR/build_hdf5.sh" "$HDF5_DIR"

# ---------------------------------------------------------------------------
# 1. The code, at the exact commit, in its own directory.
#
# If the commit is not reachable from the remote (e.g. the branch has not been
# pushed), this stops with the exact command needed rather than silently
# testing a different commit. To bypass git entirely, rsync your tree to
# $CODE_DIR beforehand -- an existing $CODE_DIR/.git at the right commit is
# accepted as-is.
# ---------------------------------------------------------------------------
if [ ! -d "$CODE_DIR/.git" ]; then
    echo "=== cloning $REPO_URL -> $CODE_DIR ==="
    git clone --branch "$BRANCH" "$REPO_URL" "$CODE_DIR" || {
        echo
        echo "ERROR: clone failed. Either the branch is not on the remote, or the"
        echo "cluster has no credential for it. Fixes:"
        echo "  * push the branch from the workstation:"
        echo "      git push -u origin $BRANCH"
        echo "  * or rsync the tree over and re-run:"
        echo "      rsync -avP --exclude build_ '<workstation>:<repo>/' $CODE_DIR/"
        exit 1
    }
fi
git -C "$CODE_DIR" fetch --all --tags --quiet || true
if ! git -C "$CODE_DIR" cat-file -e "${COMMIT}^{commit}" 2>/dev/null; then
    echo "ERROR: commit $COMMIT is not in $CODE_DIR."
    echo "It is probably not pushed. On the workstation:"
    echo "    git push -u origin $BRANCH"
    exit 1
fi
git -C "$CODE_DIR" checkout --quiet --detach "$COMMIT"
echo "checked out: $(git -C "$CODE_DIR" rev-parse HEAD)"
git -C "$CODE_DIR" --no-pager log --oneline -1
# A dirty clone would make the results untraceable to a commit.
if [ -n "$(git -C "$CODE_DIR" status --porcelain --untracked-files=no)" ]; then
    echo "WARNING: $CODE_DIR has uncommitted tracked changes; results will not"
    echo "         correspond to $COMMIT alone."
fi

# ---------------------------------------------------------------------------
# 2. The run data and scripts, in their own directory.
# ---------------------------------------------------------------------------
echo "=== staging run data -> $RUN_DIR ==="
mkdir -p "$RUN_DIR"
cp -v "$CASE_DIR/submit.sh" "$CASE_DIR/run_matrix.sh" \
      "$CASE_DIR/build_hdf5.sh" "$CASE_DIR/collect_profile.py" "$RUN_DIR/"
mkdir -p "$RUN_DIR/configs"
cp -v "$CASE_DIR/configs/"*.ini "$RUN_DIR/configs/"
cp -v "$CASE_DIR/README.md" "$RUN_DIR/" 2>/dev/null || true

# ---------------------------------------------------------------------------
# 3. Submit.
# ---------------------------------------------------------------------------
echo "=== submitting ==="
cd "$RUN_DIR"
sbatch --export=ALL,CODE_DIR="$CODE_DIR",RUN_DIR="$RUN_DIR",HDF5_ROOT="$HDF5_DIR" submit.sh
echo
echo "submitted. Monitor with:"
echo "    squeue --me"
echo "    tail -f $RUN_DIR/slurm-*.out"
echo "Results land in $RUN_DIR/results/ (one directory per run + summary.md)."
