#!/usr/bin/env bash
#SBATCH --job-name=fdm-strong-700
#SBATCH --account=hk-project-exasim
#SBATCH --partition=accelerated
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=profiling/results/slurm-%x-%j.out
#SBATCH --error=profiling/results/slurm-%x-%j.err

set -euo pipefail

export SIZE="${SIZE:-700}"
export RANKS_LIST="${RANKS_LIST:-1 2 4 8}"
export NSTEPS="${NSTEPS:-20}"
export NITER="${NITER:-3}"

bash profiling/slurm/strong_scaling_gpu.slurm
