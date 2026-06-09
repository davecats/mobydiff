#!/usr/bin/env bash
# Source this file before building or running on the cluster.

module purge
module load toolkit/nvidia-hpc-sdk/25.3

export HDF5_ROOT="${HDF5_ROOT:-/hkfs/work/workspace/scratch/xt8786-mobydiff/mobydiff/lib}"
