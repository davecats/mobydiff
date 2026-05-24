#!/usr/bin/env bash
# Source this file from SLURM scripts before building or running on the HPC system.

set -euo pipefail

module purge
module load toolkit/nvidia-hpc-sdk/25.3-nompi
module load mpi/openmpi/5.0
module load lib/hdf5/1.14

export FC="${FC:-nvfortran}"
export CC="${CC:-nvc}"
export MPI_FC_WRAPPER="${MPI_FC_WRAPPER:-mpifort}"
export FDM_FAST_MPI="${FDM_FAST_MPI:-1}"

# Different sites expose HDF5 module prefixes under different variable names.
if [ -z "${HDF5_ROOT:-}" ]; then
    if [ -n "${HDF5_HOME:-}" ]; then
        export HDF5_ROOT="$HDF5_HOME"
    elif [ -n "${EBROOTHDF5:-}" ]; then
        export HDF5_ROOT="$EBROOTHDF5"
    elif [ -n "${HDF5_DIR:-}" ]; then
        export HDF5_ROOT="$HDF5_DIR"
    fi
fi
