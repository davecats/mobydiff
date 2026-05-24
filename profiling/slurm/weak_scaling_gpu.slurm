#!/usr/bin/env bash
#SBATCH --job-name=fdm-weak
#SBATCH --account=hk-project-exasim
#SBATCH --partition=accelerated
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
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

LOCAL_N="${LOCAL_N:-256}"
NSTEPS="${NSTEPS:-20}"
NITER="${NITER:-3}"
RANKS_LIST="${RANKS_LIST:-1 2 4 8 16}"
GPUS_PER_NODE="${GPUS_PER_NODE:-4}"
RESULTS_ROOT="${RESULTS_ROOT:-$ROOT/profiling/results/${SLURM_JOB_NAME:-weak}-${SLURM_JOB_ID:-manual}}"
mkdir -p "$RESULTS_ROOT"

for ranks in $RANKS_LIST; do
    nodes=$(( (ranks + GPUS_PER_NODE - 1) / GPUS_PER_NODE ))
    if [ -n "${SLURM_JOB_NUM_NODES:-}" ] && [ "$nodes" -gt "$SLURM_JOB_NUM_NODES" ]; then
        echo "Skipping $ranks ranks: needs $nodes nodes, allocation has $SLURM_JOB_NUM_NODES"
        continue
    fi
    tasks_per_node=$(( (ranks + nodes - 1) / nodes ))
    dims="$(python3 profiling/scripts/choose_mpi_dims.py "$ranks" "$ranks" "$ranks" "$ranks")"
    read -r dx dy dz <<< "$dims"
    nx=$(( LOCAL_N * dx ))
    ny=$(( LOCAL_N * dy ))
    nz=$(( LOCAL_N * dz ))
    echo "Running weak scaling: local=${LOCAL_N}^3 global=${nx}x${ny}x${nz} ranks=$ranks nodes=$nodes dims=$dims"
    NX="$nx" NY="$ny" NZ="$nz" DIMS="$dims" NSTEPS="$NSTEPS" NITER="$NITER" \
    RANKS="$ranks" NODES="$nodes" TASKS_PER_NODE="$tasks_per_node" \
    CASE_NAME="weak_${LOCAL_N}_${ranks}r" RESULTS_ROOT="$RESULTS_ROOT" \
    profiling/scripts/run_case.sh
done

python3 profiling/scripts/collect_timings.py "$RESULTS_ROOT" > "$RESULTS_ROOT/timings.csv"
cat "$RESULTS_ROOT/timings.csv"
