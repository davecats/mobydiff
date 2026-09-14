#!/bin/bash
#SBATCH --job-name=moby_ncu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:40:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/ncu_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC/run_ncu.sh" "$SRC/collect_ncu.py" "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1
cd "$CODE_DIR" && ./compile.sh gpu || exit 1
EXE="$CODE_DIR/build_gpu/moby_solve"

RES="${RESDIR:-$RUN_DIR/results_ncu}"; mkdir -p "$RES"
{ echo "job    : ${SLURM_JOB_ID:-none}"; echo "date   : $(date '+%F %T %Z')"
  echo "node   : ${SLURM_JOB_NODELIST:-?}"
  echo "commit : $(git -C "$CODE_DIR" rev-parse HEAD)"
  echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

CFG_DIR="$SRC/../configs" bash "$STG/run_ncu.sh" "$EXE" "$RES"
python3 "$STG/collect_ncu.py" "$RES" > "$RES/ncu.md" 2>&1 \
    && { echo; cat "$RES/ncu.md"; } || echo "=== collector failed; CSVs intact ==="
echo "=== ncu job finished ==="
