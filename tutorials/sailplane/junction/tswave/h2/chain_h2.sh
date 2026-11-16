#!/bin/bash
# The h/2 chain on istmcetus (A6000, CUDA_VISIBLE_DEVICES pinned by the
# caller): h2base.ini, h2w.ini, h2p.ini; the case file is prepared
# beforehand on the workstation CPU, and ts_analysis.py runs there too
# (istmcetus has no h5py). Status lines in status.txt, logs <stage>.log.
set -u
cd "$(dirname "$0")"
export PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.9/compilers/bin:/opt/Nvidia/nvhpc/Linux_x86_64/25.9/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
export OMPI_MCA_hwloc_base_binding_policy=none
BIN=../.bin/moby_solve_gpu_forced
for ini in h2base h2w h2p; do
    echo "start $ini $(date '+%F %T') $(hostname) gpu${CUDA_VISIBLE_DEVICES:-?} $(tail -1 ../.bin/PROVENANCE | head -c 7)" >> status.txt
    /usr/bin/time -f "$ini wall %e s, RSS %M kB" mpirun -n 1 $BIN $ini.ini > $ini.log 2>&1; rc=$?
    grep -q "main loop ended" $ini.log || rc=$(( rc == 0 ? 99 : rc ))
    echo "exit $rc $ini $(date '+%F %T')" >> status.txt
    [ $rc -eq 0 ] || { echo "$ini failed rc $rc"; exit 1; }
done
echo CHAIN DONE
