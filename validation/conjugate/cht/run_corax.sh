#!/bin/bash
# Drive the high-resolution campaign on istmcorax (RTX 5090, cc120 build).
# Logs go in this directory because /tmp is NOT shared between the hosts.
export PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.9/compilers/bin:/opt/Nvidia/nvhpc/Linux_x86_64/25.9/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
cd "$(dirname "$0")"
export BIN=$PWD/../../../build_gpu_corax/moby_solve
export PREP=$PWD/../../../build_cpu/moby_prepare
LOG=$PWD/corax.log
case "${1:-all}" in
  clean) rm -f cht180_nb32.h5 cs_1.h5 IC_cht.h5 IC_stat.h5 dev_*.h5 stat_*.h5 \
               mint_*.h5 cht_stats*.h5 cht_vel_stats.h5; echo "[hr] cleaned" >> $LOG ;;
  time)  sed -e 's|^nsteps.*|nsteps = 100|' -e 's|^field_interval.*|field_interval = 0|' \
             .dev.ini > .timing.ini
         mpirun -n 1 "$BIN" .timing.ini 2>&1 | grep -E "seconds_per_step|conjugate interface" >> $LOG ;;
  *)     for ph in "$@"; do
           echo "[hr] === $ph  $(date +%H:%M:%S)" >> $LOG
           ./run_cht.sh $ph >> $LOG 2>&1
           rc=$?; echo "[hr] $ph rc=$rc" >> $LOG
           [ $rc -ne 0 ] && { echo "[hr] ABORT" >> $LOG; exit 1; }
         done
         echo "[hr] DONE $(date +%H:%M:%S)" >> $LOG ;;
esac
