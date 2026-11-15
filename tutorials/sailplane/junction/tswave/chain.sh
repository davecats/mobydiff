#!/bin/bash
# The chain of the TS-wave test after tsbase.ini (two flow-throughs from the
# impulsive start): tsbase2.ini continues the base to t 11 000 (the t 5000
# field was still in the start-up transient downstream), ts_analysis.py
# --base checks it against Blasius, tsw.ini forces the wave for 6400 theta/U,
# tsw_phases.ini samples one period in 8 snapshots, ts_analysis.py reads the
# growth out against Orr-Sommerfeld. Logs: <stage>.log, analysis_base.log,
# analysis.log; the status lines go to status.txt. The binary is the staged
# .bin/moby_solve_gpu_forced (8238e75; the forcing keys are bit-exact at
# their defaults, so the base stage runs on it too).
set -u
cd "$(dirname "$0")"
source /etc/profile.d/modules.sh 2>/dev/null; module load toolkits/nvhpc/25.9
export OMPI_MCA_hwloc_base_binding_policy=none
BIN=.bin/moby_solve_gpu_forced
BASE=tsbase_27500.h5
for ini in tsbase2 tsw tsw_phases; do
    echo "start $ini $(date '+%F %T') RTX 3060 $(tail -1 .bin/PROVENANCE | head -c 7)" >> status.txt
    /usr/bin/time -f "$ini wall %e s, RSS %M kB" mpirun -n 1 $BIN $ini.ini > $ini.log 2>&1; rc=$?
    grep -q "main loop ended" $ini.log || rc=$(( rc == 0 ? 99 : rc ))
    echo "exit $rc $ini $(date '+%F %T')" >> status.txt
    [ $rc -eq 0 ] || { echo "$ini failed rc $rc"; exit 1; }
    [ $ini = tsbase2 ] && python3 ts_analysis.py --base $BASE --ini tsw.ini > analysis_base.log 2>&1 &
done
wait
python3 ts_analysis.py --base $BASE --forced tswp_*.h5 --ini tsw.ini > analysis.log 2>&1
echo "exit $? analysis $(date '+%F %T')" >> status.txt
cat analysis_base.log analysis.log
echo CHAIN DONE
