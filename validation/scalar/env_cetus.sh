#!/usr/bin/env bash
# Run ANY command on istmcetus with the 25.9 toolchain in scope.
#
# istmcetus has NO nvhpc 25.9 under /opt/Nvidia and no modulefile for it --
# only the 25.9 tree on the shared filesystem, which is the same one the
# local `module load toolkits/nvhpc/25.9` points at. Without this, `which
# mpirun` is /usr/bin/mpirun (system OpenMPI, gfortran modules) and every
# run dies in opal_init with an unhelpful missing-help-file message.
#
#   ssh istmcetus bash <abs path>/env_cetus.sh <device> <command...>
#
# The machine is SHARED (2x A6000): check nvidia-smi before picking a device.
set -uo pipefail
cd "$(dirname "$0")"

dev=${1:?usage: env_cetus.sh <device> <command...>}
shift

NV=/net/istmolympus/md0/software/Nvidia/nvhpc/Linux_x86_64/25.9
export PATH=$NV/compilers/bin:$NV/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
export LD_LIBRARY_PATH=$NV/compilers/lib:$NV/comm_libs/12.9/hpcx/latest/ompi/lib:${LD_LIBRARY_PATH:-}
export CUDA_VISIBLE_DEVICES=$dev

echo "host $(hostname)  device $dev  cmd: $*"
exec "$@"
