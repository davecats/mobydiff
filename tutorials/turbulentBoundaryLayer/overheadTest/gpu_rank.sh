#!/usr/bin/env bash
# Per-rank GPU pinning for multi-GPU runs. Drop-in for run_overhead.sh's BIN:
#
#   BIN=./gpu_rank.sh NRANKS=2 GPUS=0,1 ./run_overhead.sh <config>
#
# WHY THIS IS NEEDED. gpu_runtime.f90 never calls omp_set_default_device, so
# every rank offloads to OpenMP device 0. Under one mpirun that is the SAME
# physical GPU for all of them -- an unpinned 2-rank run is not a 2-GPU run, it
# is two ranks fighting over one card, and it looks like catastrophic scaling
# rather than a configuration mistake. Giving each rank its own
# CUDA_VISIBLE_DEVICES makes "device 0" resolve to a different card per rank.
# The solver prints "OpenMP target devices available: 1" per rank when this is
# working; a 2 there means the pinning did not take.
#
# GPUS: comma-separated physical devices, assigned round-robin by LOCAL rank
# (default 0,1 -- both istmcetus A6000s). Check nvidia-smi first: the machine is
# shared and GPU 0 belongs to the production campaign.
set -eu

GPUS="${GPUS:-0,1}"
IFS=',' read -r -a devs <<< "$GPUS"
lr="${OMPI_COMM_WORLD_LOCAL_RANK:-${SLURM_LOCALID:-0}}"
export CUDA_VISIBLE_DEVICES="${devs[$(( lr % ${#devs[@]} ))]}"

BIN="${SOLVER:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)/build_gpu/moby_solve}"
exec "$BIN" "$@"
