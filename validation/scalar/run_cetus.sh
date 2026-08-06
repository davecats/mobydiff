#!/usr/bin/env bash
# Run one ini on istmcetus' GPU (2x A6000 cc86; the local build_gpu binary
# runs there -- shared filesystem, same paths). CUDA_VISIBLE_DEVICES is
# pinned because the machine is shared: check nvidia-smi first.
#
#   ssh istmcetus bash /home/ws/xt8786/Codes/mobydiff.scalar/validation/scalar/run_cetus.sh <ini> [device]
set -uo pipefail
cd "$(dirname "$0")"

ini=${1:?usage: run_cetus.sh <ini> [device]}
dev=${2:-0}

# istmcetus has NO nvhpc 25.9 under /opt/Nvidia and no modulefile for it --
# only the 25.9 tree on the shared filesystem, which is the same one the
# local module points at. Without this, `which mpirun` is /usr/bin/mpirun
# (system OpenMPI, gfortran) and the run dies in opal_init.
NV=/net/istmolympus/md0/software/Nvidia/nvhpc/Linux_x86_64/25.9
export PATH=$NV/compilers/bin:$NV/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
export LD_LIBRARY_PATH=$NV/compilers/lib:$NV/comm_libs/12.9/hpcx/latest/ompi/lib:${LD_LIBRARY_PATH:-}
export CUDA_VISIBLE_DEVICES=$dev
echo "host $(hostname)  device $dev  ini $ini"
nvidia-smi --query-gpu=index,name,memory.used --format=csv,noheader
exec mpirun -n 1 ../../build_gpu/moby_solve "$ini"
