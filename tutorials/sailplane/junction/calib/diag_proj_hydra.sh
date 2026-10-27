#!/bin/bash
# Projection convergence per step (2026-10-27): ONE step from the uniform
# cold start for several [pressure] settings, 32 istmhydra ranks each,
# sequential; snapshot_divergence.py on the proj_<v>_1.h5 files afterwards.
set -u
cd "$(dirname "$0")"
export PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.9/compilers/bin:/opt/Nvidia/nvhpc/Linux_x86_64/25.9/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
export OMPI_MCA_hwloc_base_binding_policy=none OMP_NUM_THREADS=1
echo "start $(date)" > proj.status
for v in cheb12 cheb24 cheb48 cheb96 rb6 rb20; do
    /usr/bin/time -f "wall %e s" mpirun -n 32 ../../../../build_cpu/moby_solve .proj_$v.ini > proj_$v.log 2>&1
    echo "$v exit $? $(date)" >> proj.status
done
echo "end $(date)" >> proj.status
