#!/bin/bash
#SBATCH --job-name=moby_tune
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:40:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#
# Calibrate moby_tune on HoreKa. Two questions:
#   1. Does it find the mapping the analysis found, without being told the
#      topology?  Expect a pair drawn from {0,1} -- the NUMA-0 cards.
#   2. Does it agree with itself at 8 and at 16 ranks?  The answer is a property
#      of the node, so it must not depend on the rank count.
# A third, implicit: does the tie logic fire when it should stay quiet.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
module purge; module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1
RES="${RESDIR:-$RUN_DIR/results_tune}"; mkdir -p "$RES"; cd "$RES" || exit 1
{ echo "job   : ${SLURM_JOB_ID:-none}"; echo "nodes : ${SLURM_JOB_NODELIST:-?}"
  echo "commit: $(git -C "$CODE_DIR" rev-parse HEAD)"; } | tee provenance.txt

echo; echo "############ 8 ranks / 2 nodes ############"
RANKS=8 NODES=2 bash "$CODE_DIR/tools/moby_tune.sh" \
    "$CODE_DIR/build_gpu/moby_solve" "$RUN_DIR/configs/rect_jacobi.ini" "$RES/r8"

echo; echo "############ 16 ranks / 4 nodes ############"
RANKS=16 NODES=4 bash "$CODE_DIR/tools/moby_tune.sh" \
    "$CODE_DIR/build_gpu/moby_solve" "$RUN_DIR/configs/rect_jacobi.ini" "$RES/r16"

echo; echo "############ single node: must decline to tune ############"
RANKS=4 NODES=1 bash "$CODE_DIR/tools/moby_tune.sh" \
    "$CODE_DIR/build_gpu/moby_solve" "$RUN_DIR/configs/rect_jacobi.ini" "$RES/r4"

echo "=== tune job finished ==="
