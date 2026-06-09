#!/usr/bin/env bash
#SBATCH --job-name=fdm-strong-1024
#SBATCH --account=hk-project-exasim
#SBATCH --partition=accelerated
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=profiling/results/slurm-%x-%j.out
#SBATCH --error=profiling/results/slurm-%x-%j.err

set -euo pipefail

export SIZE="${SIZE:-1024}"
export RANKS_LIST="${RANKS_LIST:-4 8 16}"
export NSTEPS="${NSTEPS:-20}"
export NITER="${NITER:-3}"

bash profiling/slurm/strong_scaling_gpu.slurm
