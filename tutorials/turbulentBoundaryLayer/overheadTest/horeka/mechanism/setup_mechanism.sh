#!/bin/bash
# One-shot setup + submit for the MECHANISM probe on HoreKa.
#
# RUN ON A LOGIN NODE from inside the transferred package's mechanism/ dir:
#     bash setup_mechanism.sh
#
# Reuses the main campaign's code clone and HDF5 if they are already there
# (the commit is the same), so normally this only stages and submits.
set -euo pipefail

REPO_URL="${REPO_URL:-git@github.com:davecats/mobydiff.git}"
BRANCH="${BRANCH:-optimiseBlockRefinement_parentBoundaryLayer}"
# The SAME commit the 2026-09-07 matrix measured -- this probe explains that
# data, so it must not run different code.
COMMIT="${COMMIT:-a11e355a47e1535db4f3e9bae0dcc489eaec3567}"

WS="${WS:-$HOME}"
CODE_DIR="${CODE_DIR:-$WS/moby-2to1-code}"
RUN_DIR="${RUN_DIR:-$WS/moby-2to1-run}"
HDF5_DIR="${HDF5_DIR:-$HOME/hdf5}"

MECH_DIR="$(cd "$(dirname "$0")" && pwd -P)"
CASE_DIR="$(cd "$MECH_DIR/.." && pwd -P)"

echo "COMMIT    : $COMMIT"
echo "CODE_DIR  : $CODE_DIR"
echo "RUN_DIR   : $RUN_DIR"
[ "$WS" = "$HOME" ] && echo "NOTE: WS is \$HOME. Set WS=/hkfs/work/workspace/scratch/<user>-<ws>."
echo

module purge
module load toolkit/nvidia-hpc-sdk/25.3
bash "$CASE_DIR/build_hdf5.sh" "$HDF5_DIR"

if [ ! -d "$CODE_DIR/.git" ]; then
    git clone --branch "$BRANCH" "$REPO_URL" "$CODE_DIR" || {
        echo "ERROR: clone failed. Push the branch, or rsync the tree to $CODE_DIR." >&2
        exit 1; }
fi
git -C "$CODE_DIR" fetch --all --tags --quiet || true
git -C "$CODE_DIR" cat-file -e "${COMMIT}^{commit}" 2>/dev/null || {
    echo "ERROR: commit $COMMIT not in $CODE_DIR (push the branch)." >&2; exit 1; }
git -C "$CODE_DIR" checkout --quiet --detach "$COMMIT"
echo "checked out: $(git -C "$CODE_DIR" rev-parse HEAD)"

echo "=== staging -> $RUN_DIR ==="
mkdir -p "$RUN_DIR/mechanism" "$RUN_DIR/configs"
cp -v "$CASE_DIR/build_hdf5.sh" "$RUN_DIR/"
cp -v "$MECH_DIR/run_placement.sh" "$MECH_DIR/collect_mechanism.py" \
      "$MECH_DIR/submit_mechanism.sh" "$MECH_DIR/README.md" "$RUN_DIR/mechanism/"
cp -v "$CASE_DIR/configs/"*.ini "$RUN_DIR/configs/"
cp -v "$CASE_DIR/HANDOUT_cluster_session.md" "$RUN_DIR/" 2>/dev/null || true
# Staging is the step that silently cost round 1 five runs: the cluster copy of
# configs/ predated the commit that added blk8/nb16, and run_placement.sh could
# only report MISSING CONFIG. Name what is there, every time.
echo "configs staged: $(ls "$RUN_DIR/configs/" | tr '\n' ' ')"

echo "=== submitting ==="
cd "$RUN_DIR"
sbatch --export=ALL,CODE_DIR="$CODE_DIR",RUN_DIR="$RUN_DIR",HDF5_ROOT="$HDF5_DIR" \
       mechanism/submit_mechanism.sh
echo
echo "monitor : squeue --me  |  tail -f $RUN_DIR/slurm-*.out"
echo "results : $RUN_DIR/results_mechanism/mechanism.md"
