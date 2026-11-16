#!/bin/bash
# The c/512 twin on the workstation's RTX 3060: seed, development, then the
# three projection branches from the common state. Every branch's snapshots
# are reduced to span-averaged moments (sd_reduce.py) and deleted, the last
# one (the restart state) kept. Status lines in status.txt.
#     setsid nohup ./chain4.sh > chain4.out 2>&1 < /dev/null &
cd "$(dirname "$0")"
source /etc/profile.d/lmod.sh 2>/dev/null
module load toolkits/nvhpc/25.9
export OMPI_MCA_hwloc_base_binding_policy=none
BIN=.bin/moby_solve_gpu
for st in ${STAGES:-seed dev c12 c24 rb}; do
    ini=sd4_$st.ini
    grep -q "^exit 0 $ini " status.txt 2>/dev/null && continue
    echo "start $ini $(date '+%F %T') $(hostname) RTX 3060 $(head -c 7 .bin/PROVENANCE)" >> status.txt
    mpirun -n 1 $BIN $ini > sd4_$st.log 2>&1; rc=$?
    grep -q "main loop ended" sd4_$st.log || rc=$(( rc == 0 ? 99 : rc ))
    echo "exit $rc $ini $(date '+%F %T')" >> status.txt
    [ $rc -eq 0 ] || exit 1
    case $st in c12|c24|rb)
        files=$(ls -1 sd4_${st}_[0-9]*.h5 | sed 's/.*_\([0-9]*\)\.h5/\1 &/' | sort -n | cut -d' ' -f2)
        python3 sd_reduce.py $files > reduce4_$st.log 2>&1 || { echo "reduce FAILED $st" >> status.txt; exit 1; }
        for f in $(echo "$files" | head -n -1); do [ -s "${f%.h5}.red.npz" ] && rm -f "$f"; done ;;
    esac
done
echo "CHAIN DONE $(date '+%F %T')" >> status.txt
