#!/bin/bash
#SBATCH --job-name=moby_build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60gb
#SBATCH --time=01:00:00
#SBATCH --partition=dev-cpu
#SBATCH --account=hk-project-p0024333
#SBATCH --export=NONE
#
# Build mobydiff's GPU binary on HoreKa 2 (2026-11-12) on a CPU node of the
# dev-cpu partition -- the login nodes kill long compiles and the GPU queue is
# days long: the H100 target is given explicitly (-gpu=cc90, no GPU needed at
# build time; compile.sh reads OPENMP_OFFLOAD_FLAGS from the environment).
#     sbatch --export=NONE,CODE_DIR=$PROJECT/moby-junction build_hk2.sh
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"
module purge
module load NVHPC/25.9-CUDA-12.9.1 HDF5/1.14.6-NVHPC-25.9-CUDA-12.9.1 CMake/3.31.12   # the toolchain-free CMake (cmake is not on the node PATH)
export HDF5_ROOT="$EBROOTHDF5"
export OPENMP_OFFLOAD_FLAGS="-mp=gpu -gpu=cc90"
cd "$CODE_DIR" && rm -rf build_gpu
echo "=== building $(git rev-parse --short HEAD) on $(hostname) with $(nvfortran --version | sed -n 2p)"
./compile.sh gpu && ls -la build_gpu/moby_solve build_gpu/moby_prepare && module list 2>&1 | grep -E 'NVHPC|HDF5'
