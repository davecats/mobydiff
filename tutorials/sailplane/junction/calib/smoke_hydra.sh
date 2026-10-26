#!/bin/bash
# Rank-identity smoke test of the alpha 1.5 / 3.5 c case on istmhydra's CPUs
# (the 190 M-cell case does not fit the RTX 3060): 20 steps on 16 and on 32
# ranks of the workstation build_cpu, then h5maxdiff on un vn wn pn (+ the
# RANS scalars) -> max_abs 0. Waits for prep_all.sh to finish first.
#   ssh istmhydra 'cd <this dir> && nohup ./smoke_hydra.sh > smoke_hydra.nohup 2>&1 < /dev/null &'
set -u
cd "$(dirname "$0")"
export PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.9/compilers/bin:/opt/Nvidia/nvhpc/Linux_x86_64/25.9/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
export OMPI_MCA_hwloc_base_binding_policy=none OMP_NUM_THREADS=1
until grep -q "^end" prep_all.status 2>/dev/null; do sleep 60; done
echo "start $(date) on $(hostname)" > smoke.status
for n in 16 32; do
    /usr/bin/time -f "solve wall %e s, max RSS %M kB" mpirun -n $n ../../../../build_cpu/moby_solve .smoke$n.ini > smoke$n.log 2>&1
    echo "smoke$n exit $? $(date)" >> smoke.status
done
../../../../tools/h5maxdiff smoke16_20.h5 smoke32_20.h5 > smoke_diff.txt 2>&1
echo "h5maxdiff exit $? $(date)" >> smoke.status
echo "end $(date)" >> smoke.status
