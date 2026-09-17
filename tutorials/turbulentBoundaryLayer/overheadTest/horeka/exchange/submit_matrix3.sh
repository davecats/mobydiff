#!/bin/bash
#SBATCH --job-name=moby_matrix3
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=04:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# Task 0 of the divergence-halo session: re-measure the campaign matrix AGAIN.
#
# results_horeka_2026-09-14.md was measured at 55bee89 (the map(to: c) fix).
# Since then five increments took ~18% off the step -- the two register cuts
# (jacobi_apply k2 88->64; compute_rdenom 110->80 and jacobi_compute_phi 94->80)
# and the step-work increments (ibm_mu elided on body-free ranks, the trip force
# confined to its envelope, the Laplacian coefficient pack). They landed
# UNEVENLY across buckets, so every ratio that report publishes -- the block tax,
# the strong-scaling efficiencies, the 2:1 coarse-cell-equivalent, the red-black
# erosion -- has a moved denominator.
#
# The SAME 23-run matrix (run_matrix.sh unchanged, --map-by numa --bind-to core,
# 200 steps) run THREE times in ONE allocation, so each column isolates one
# increment and the pair (ref, new) gives the total:
#   ref = 55bee89  -- the `new` column of results_horeka_2026-09-14.md, so its
#                     column here is the CONTROL: it must reproduce that report.
#   mid = 3c2903a  -- + the register cuts and the step-work increments (09-14).
#   new = 95312d7  -- + the divergence-halo exchange (09-15).
#
# All three are PINNED WORKTREES: the live tree keeps moving and a job that
# starts two days from now must not build whatever it has become.
#
# NEW then REF then MID: if the wall clock bites, the two columns that give the
# TOTAL survive, and only the per-increment split is lost.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"; MID_DIR="${MID_DIR:?}"
# What each column IS varies by run; the collector writes scaling_<pair>.md
# named from these, so a table can never be quoted against the wrong question.
LBL_REF="${LBL_REF:-ref}"; LBL_MID="${LBL_MID:-mid}"; LBL_NEW="${LBL_NEW:-new}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka"

# Stage the driver AND the configs it resolves relative to itself, then run the
# copy -- bash reads a script by byte offset and an edit mid-run splices two
# versions together (job 5139461).
STG="$RUN_DIR/matrix3_staged"
rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC/run_matrix.sh" "$SRC/mechanism/collect_scaling.py" "$STG/"
cp -r "$SRC/configs" "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

# A CMake cache records the ABSOLUTE source path, so a worktree that was
# renamed (git worktree move) carries a cache pointing at its old name and
# cmake refuses to configure -- which killed job 5145816 after it had waited a
# day in the queue. Check and wipe rather than trust.
for d in "$REF_DIR" "$MID_DIR" "$CODE_DIR"; do
    cache="$d/build_gpu/CMakeCache.txt"
    if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$d" "$cache"; then
        echo "=== $d: build_gpu cache was configured for a different source dir -- discarding"
        grep -m1 "^CMAKE_HOME_DIRECTORY:" "$cache" | sed 's/^/    /'
        rm -rf "$d/build_gpu"
    fi
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
NEW="$CODE_DIR/build_gpu/moby_solve"; REF="$REF_DIR/build_gpu/moby_solve"
MID="$MID_DIR/build_gpu/moby_solve"
[ -x "$NEW" ] && [ -x "$REF" ] && [ -x "$MID" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_matrix3}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  ($LBL_NEW)"
    echo "mid    : $(git -C "$MID_DIR" rev-parse HEAD)  ($LBL_MID)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)  ($LBL_REF)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps : ${NSTEPS:-200}"
    echo "layout : --map-by numa --bind-to core (as every earlier matrix)"
} | tee "$RES/provenance.txt"

echo "############ NEW ($LBL_NEW) ############"
NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$NEW" "$RES/new"
echo "############ REF ($LBL_REF) ############"
NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$REF" "$RES/ref"
echo "############ MID ($LBL_MID) ############"
NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$MID" "$RES/mid"

# Three pairings: the total, and each step of the way on its own.
for pair in "ref new total" "ref mid ${LBL_MID}" "mid new ${LBL_NEW}"; do
    set -- $pair
    python3 "$STG/collect_scaling.py" "$RES/$1" "$RES/$2" > "$RES/scaling_$3.md" 2>&1 \
        && { echo; echo "===== $3 ($1 -> $2)"; cat "$RES/scaling_$3.md"; } \
        || echo "=== collector failed for $3; raw logs intact ==="
done
echo "=== matrix3 job finished ==="
