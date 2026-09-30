#!/bin/bash
# run_leg.sh <dir> <binary> <gpu>  -- leg A (transient, stats off) then leg B (stats on), status file per leg
d=$1; bin=$2; export CUDA_VISIBLE_DEVICES=$3
source /etc/profile.d/lmod.sh; module load toolkits/nvhpc/25.9
cd "$d" && echo "start $(date '+%F %T')" > run.status
mpirun -n 1 "$bin" legA.ini > legA.log 2>&1; rc=$?; echo "legA exit=$rc $(date '+%F %T')" >> run.status; [ $rc -ne 0 ] && exit $rc
last=$(ls -t legA_*[0-9].h5 | head -1); sed -i "s|@LEGA@|$last|" legB.ini
mpirun -n 1 "$bin" legB.ini > legB.log 2>&1; rc=$?; echo "legB exit=$rc $(date '+%F %T')" >> run.status; echo done >> run.status
