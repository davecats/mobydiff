#!/bin/bash
#SBATCH --job-name=moby_timeline
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:50:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The nsys per-kernel timeline that results_horeka_exchange_2026-09-10.md section
# 8 named as the required next step before any implementation.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
# Run a STAGED copy: bash reads a script by byte offset, so editing the source
# while a job executes it splices two versions together (job 5139461).
EXCH="$RUN_DIR/timeline_staged"
rm -rf "$EXCH"; mkdir -p "$EXCH"
cp "$SRC"/run_timeline.sh "$SRC"/collect_timeline.py "$EXCH/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

cd "$CODE_DIR" || exit 1
./compile.sh gpu || exit 1
EXE="$CODE_DIR/build_gpu/moby_solve"
[ -x "$EXE" ] || { echo "ERROR: no $EXE" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_timeline}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "commit : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "dirty  : $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps : ${NSTEPS:-30}"
} | tee "$RES/provenance.txt"

CFG_DIR="$SRC/../configs" bash "$EXCH/run_timeline.sh" "$EXE" "$RES"

python3 "$EXCH/collect_timeline.py" "$RES" > "$RES/timeline.md" 2>&1 \
    && { echo; cat "$RES/timeline.md"; } || echo "=== collector failed; CSVs intact ==="
echo "=== timeline job finished ==="
