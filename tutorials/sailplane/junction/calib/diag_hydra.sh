#!/bin/bash
# Diagnostic of the cold-start divergence (2026-10-27): 50 steps of the
# alpha 1.5 / 3.5 c case on 32 istmhydra CPU ranks, a snapshot every 10
# steps (.diag50.ini), to see when and where the velocity leaves [-3, 3].
set -u
cd "$(dirname "$0")"
export PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.9/compilers/bin:/opt/Nvidia/nvhpc/Linux_x86_64/25.9/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
export OMPI_MCA_hwloc_base_binding_policy=none OMP_NUM_THREADS=1
echo "start $(date)" > diag50.status
mpirun -n 32 ../../../../build_cpu/moby_solve .diag50.ini > diag50.log 2>&1
echo "exit $? $(date)" >> diag50.status
