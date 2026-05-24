#!/usr/bin/env bash
#SBATCH --job-name=fdm-nsys
#SBATCH --account=hk-project-exasim
#SBATCH --partition=accelerated
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=profiling/results/slurm-%x-%j.out
#SBATCH --error=profiling/results/slurm-%x-%j.err

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$ROOT"
mkdir -p profiling/results
source profiling/scripts/load_hpc_modules.sh

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-MANDATORY}"

if [ "${DO_BUILD:-1}" = "1" ]; then
    ./compile.sh gpu
fi

SIZE="${SIZE:-512}"
RANKS="${RANKS:-4}"
NSTEPS="${NSTEPS:-5}"
NITER="${NITER:-3}"
GPUS_PER_NODE="${GPUS_PER_NODE:-4}"
NODES="$(( (RANKS + GPUS_PER_NODE - 1) / GPUS_PER_NODE ))"
TASKS_PER_NODE="$(( (RANKS + NODES - 1) / NODES ))"
DIMS="$(python3 profiling/scripts/choose_mpi_dims.py "$SIZE" "$SIZE" "$SIZE" "$RANKS")"
RESULTS_ROOT="${RESULTS_ROOT:-$ROOT/profiling/results/${SLURM_JOB_NAME:-nsys}-${SLURM_JOB_ID:-manual}}"
CASE_DIR="$RESULTS_ROOT/nsys_${SIZE}_${RANKS}r"
mkdir -p "$CASE_DIR"

NX="$SIZE" NY="$SIZE" NZ="$SIZE" DIMS="$DIMS" NSTEPS="$NSTEPS" NITER="$NITER" \
    profiling/scripts/make_input.sh "$CASE_DIR/input.ini"

if ! command -v nsys >/dev/null 2>&1; then
    echo "nsys not found in the loaded modules" >&2
    exit 1
fi

srun --nodes "$NODES" --ntasks "$RANKS" --ntasks-per-node "$TASKS_PER_NODE" \
     --cpus-per-task "${CPUS_PER_TASK:-8}" --gpus-per-task=1 --gpu-bind="${GPU_BIND:-single:1}" \
     nsys profile --force-overwrite=true --sample=none --trace=cuda,nvtx,mpi,osrt \
     --output="$CASE_DIR/nsys_rank_%q{SLURM_PROCID}" \
     "$ROOT/build_gpu_mpi/main" "$CASE_DIR/input.ini" \
     > "$CASE_DIR/stdout.txt" 2> "$CASE_DIR/stderr.txt"

for rep in "$CASE_DIR"/*.nsys-rep; do
    [ -e "$rep" ] || continue
    nsys stats --force-export=true --report cuda_gpu_kern_sum,cuda_gpu_mem_time_sum,mpi_event_sum \
        --format csv "$rep" > "${rep%.nsys-rep}_stats.csv" 2> "${rep%.nsys-rep}_stats.err" || true
done
