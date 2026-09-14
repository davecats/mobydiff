#!/bin/bash
#SBATCH --job-name=moby_matrix3
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=03:00:00
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
# 200 steps) run twice in ONE allocation:
#   ref = 55bee89  -- the `new` column of results_horeka_2026-09-14.md, so its
#                     column here is the CONTROL: it must reproduce that report.
#   new = 3c2903a  -- the head of the register + step-work line.
#
# Both sides are PINNED WORKTREES. The live tree is about to grow the task-1
# divergence-halo exchange and a job that starts hours from now must not build it.
#
# NEW first: if the wall clock bites, the current numbers are the deliverable.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
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

for d in "$REF_DIR" "$CODE_DIR"; do
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
NEW="$CODE_DIR/build_gpu/moby_solve"; REF="$REF_DIR/build_gpu/moby_solve"
[ -x "$NEW" ] && [ -x "$REF" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_matrix3}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  (registers + step work)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)  (map(to: c); = the 2026-09-14 'new' column)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps : ${NSTEPS:-200}"
    echo "layout : --map-by numa --bind-to core (as every earlier matrix)"
} | tee "$RES/provenance.txt"

echo "############ NEW (registers + step work, 3c2903a) ############"
NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$NEW" "$RES/new"
echo "############ REF (map(to: c), 55bee89) ############"
NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$REF" "$RES/ref"

python3 "$STG/collect_scaling.py" "$RES/ref" "$RES/new" > "$RES/scaling.md" 2>&1 \
    && { echo; cat "$RES/scaling.md"; } || echo "=== collector failed; raw logs intact ==="
echo "=== matrix3 job finished ==="
