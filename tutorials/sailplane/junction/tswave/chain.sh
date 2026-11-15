#!/bin/bash
# The forced chain of the TS-wave test, run after tsbase.ini has ended:
# tsw.ini (1.5 flow-throughs of harmonic forcing on the base field), then
# tsw_phases.ini (8 snapshots over one period), then ts_analysis.py.
# Invoked in the background right after tsbase was launched; it waits for
# the "exit 0 tsbase" line of status.txt. Logs: tsw.log, tsw_phases.log,
# analysis.log; the status lines go to status.txt.
set -u
cd "$(dirname "$0")"
source /etc/profile.d/modules.sh 2>/dev/null; module load toolkits/nvhpc/25.9
export OMPI_MCA_hwloc_base_binding_policy=none
SOLVER=${SOLVER:-../../../../build_gpu/moby_solve}
until grep -q "^exit . tsbase" status.txt; do sleep 60; done
grep -q "^exit 0 tsbase" status.txt || { echo "tsbase failed: $(grep '^exit' status.txt | tail -1)"; exit 1; }
cp "$SOLVER" .bin/moby_solve_gpu_forced; git -C ../../../.. rev-parse HEAD >> .bin/PROVENANCE
python3 ts_analysis.py --base tsbase_12500.h5 --ini tsw.ini > analysis_base.log 2>&1 &
for ini in tsw tsw_phases; do
    echo "start $ini $(date '+%F %T') RTX 3060 $(git -C ../../../.. rev-parse --short HEAD)" >> status.txt
    /usr/bin/time -f "$ini wall %e s, RSS %M kB" mpirun -n 1 ./.bin/moby_solve_gpu_forced $ini.ini > $ini.log 2>&1; rc=$?
    grep -q "main loop ended" $ini.log || rc=$(( rc == 0 ? 99 : rc ))
    echo "exit $rc $ini $(date '+%F %T')" >> status.txt
    [ $rc -eq 0 ] || { echo "$ini failed rc $rc"; exit 1; }
done
wait
python3 ts_analysis.py --base tsbase_12500.h5 --forced tswp_*.h5 --ini tsw.ini > analysis.log 2>&1
echo "exit $? analysis $(date '+%F %T')" >> status.txt
cat analysis_base.log analysis.log
echo CHAIN DONE
