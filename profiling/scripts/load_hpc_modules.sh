#!/usr/bin/env bash
# Source this file before building or running on the cluster.

set -euo pipefail

module purge
module load toolkit/nvidia-hpc-sdk/25.3

export HDF5_ROOT="${HDF5_ROOT:-/hkfs/work/workspace/scratch/xt8786-mobydiff/mobydiff/lib}"
export FC="${FC:-mpifort}"
export CC="${CC:-mpicc}"
export MPI_FC_WRAPPER="${MPI_FC_WRAPPER:-mpifort}"
export MPI_C_WRAPPER="${MPI_C_WRAPPER:-mpicc}"
export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-MANDATORY}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

# The GPU path passes OpenMP device buffers directly to MPI.
export OMPI_MCA_pml="${OMPI_MCA_pml:-ucx}"
export OMPI_MCA_coll_hcoll_enable="${OMPI_MCA_coll_hcoll_enable:-0}"
export UCX_TLS="${UCX_TLS:-rc,sm,self,cuda_copy,cuda_ipc}"
export UCX_MEMTYPE_CACHE="${UCX_MEMTYPE_CACHE:-y}"
export UCX_RNDV_THRESH="${UCX_RNDV_THRESH:-0}"
