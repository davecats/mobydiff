#!/bin/bash
# Prepare the four calibration cases on an interactive ISTM host that shares
# this NFS tree (run 2026-10-26 on istmhydra: 256 EPYC 7662 cores at 2.1 GHz,
# 251 GB). The c/512 prepare is 3 h on the workstation's 4 Ryzen ranks and
# about 50 GB resident; moby_prepare is MPI-parallel only (no host OpenMP),
# and with an explicit [blocks] nb the case file does not depend on the rank
# count (checked 2026-10-26: the sailplane tutorial case prepared on 2 and
# 4 ranks is bit-identical in every dataset; the NACA nose case was prepared
# on 20 ranks and solved on 4), so: RANKS ranks per case, two cases at a time.
# No modulefile in a non-interactive ssh shell: nvhpc 25.9 by PATH.
#   ssh istmhydra 'cd <this dir> && nohup ./prep_all.sh > prep_all.nohup 2>&1 < /dev/null &'
set -u
cd "$(dirname "$0")"
export PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.9/compilers/bin:/opt/Nvidia/nvhpc/Linux_x86_64/25.9/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
export OMPI_MCA_hwloc_base_binding_policy=none
PREP=../../../../build_cpu/moby_prepare
RANKS=${RANKS:-48}
STATUS=prep_all.status
echo "start $(date) on $(hostname), $RANKS ranks per case" > $STATUS

one() {   # one case, its own log with wall time and max RSS at the end
    local c=$1
    /usr/bin/time -f "prepare wall %e s, max RSS %M kB (one rank, the time wrapper sees mpirun's children)" \
        mpirun -n $RANKS $PREP $c.ini > $c.prep.log 2>&1
    echo "$c exit $? $(date)" >> $STATUS
}

one calib_a1.5_y3.5 & one calib_a1.0_y3.5 & wait
one calib_a2.0_y3.5 & one calib_a1.5_y4.5 & wait
echo "end $(date)" >> $STATUS
