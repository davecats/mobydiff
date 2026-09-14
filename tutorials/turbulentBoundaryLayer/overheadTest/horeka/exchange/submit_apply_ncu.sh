#!/bin/bash
#SBATCH --job-name=moby_apply_ncu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:40:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The pre-registered reading for the apply register cut: did occupancy and DRAM
# utilisation actually follow the register count? Both binaries, same node, same
# probe -- the 2026-09-14 ncu numbers were taken on another node, so the BEFORE
# side is re-measured here rather than quoted across allocations.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/apply_ncu_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC/run_ncu.sh" "$SRC/collect_ncu.py" "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1
for d in "$REF_DIR" "$CODE_DIR"; do ( cd "$d" && ./compile.sh gpu ) || exit 1; done

RES="${RESDIR:-$RUN_DIR/results_apply_ncu}"; mkdir -p "$RES"
{ echo "job    : ${SLURM_JOB_ID:-none}"; echo "date   : $(date '+%F %T %Z')"
  echo "node   : ${SLURM_JOB_NODELIST:-?}"
  echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  dirty $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
  echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
  echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

for side in ref new; do
    exe="$REF_DIR/build_gpu/moby_solve"
    [ "$side" = new ] && exe="$CODE_DIR/build_gpu/moby_solve"
    echo "############ ncu $side ############"
    CFG_DIR="$SRC/../configs" bash "$STG/run_ncu.sh" "$exe" "$RES/$side"
    python3 "$STG/collect_ncu.py" "$RES/$side" > "$RES/ncu_$side.md" 2>&1 \
        && { echo; cat "$RES/ncu_$side.md"; } || echo "=== collector failed for $side; CSVs intact ==="
done
echo "=== apply ncu job finished ==="
