#!/bin/bash
#SBATCH --job-name=moby_scaling2
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=02:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# Step 0 after the map(to: c) fix: re-measure the whole campaign.
#
# Every strong-scaling efficiency, block tax and red-black ratio in
# results_horeka_2026-09-10.md section 9 was measured with 6-8 ms/step of
# per-launch descriptor marshalling that no longer exists (job 5142027/5142047).
# None of those tables may be quoted until this runs. The SAME 23-run matrix
# (run_matrix.sh, unchanged, --map-by numa --bind-to core, 200 steps) is run
# twice in ONE allocation -- once with a worktree at the pre-fix commit and once
# with the working tree -- so the two columns differ only in the fix.
#
# NEW first: if the wall clock bites, the corrected numbers are the deliverable
# and the reference column is the nice-to-have.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka"

# Stage the driver AND the configs it resolves relative to itself, then run the
# copy -- bash reads a script by byte offset and an edit mid-run splices two
# versions together (job 5139461).
STG="$RUN_DIR/scaling2_staged"
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

RES="${RESDIR:-$RUN_DIR/results_scaling2}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  (map(to: c))"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)  (pre-fix)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps : ${NSTEPS:-200}"
    echo "layout : --map-by numa --bind-to core (as every earlier matrix)"
} | tee "$RES/provenance.txt"

echo "############ NEW (map(to: c)) ############"
NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$NEW" "$RES/new"
echo "############ REF (pre-fix) ############"
NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$REF" "$RES/ref"

python3 "$STG/collect_scaling.py" "$RES/ref" "$RES/new" > "$RES/scaling.md" 2>&1 \
    && { echo; cat "$RES/scaling.md"; } || echo "=== collector failed; raw logs intact ==="
echo "=== scaling2 job finished ==="
