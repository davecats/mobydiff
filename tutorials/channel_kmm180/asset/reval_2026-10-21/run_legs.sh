#!/usr/bin/env bash
# istmcorax chain 3 (re-validation 2026-10-21): the KMM180 two-leg statistics rerun
# (leg A transient to t = 35 with statistics off, leg B to t = 55 with statistics on),
# from the archived dyw+ 0.5 restart, red-black niter 3 / sor 1.5, cflmax 1.0.
source ~/outlet_runs/env.sh
export OMPI_MCA_hwloc_base_binding_policy=none
BIN=/home/ws/xt8786/Codes/mobydiff.cleanup/build_gpu_corax/moby_solve
ST=~/corax_chain/chain3.status
echo "start $(date '+%F %T') $(mpirun -n 1 $BIN /dev/null 2>&1 | grep -o 'commit: .*')" > $ST
cd ~/kmm_reval || exit 1
rm -f legA_*.h5 legB_*.h5
t0=$SECONDS
mpirun -n 1 $BIN legA.ini > legA.log 2>&1; rc=$?
echo "kmm legA exit=$rc wall=$((SECONDS-t0))s $(date '+%F %T')" >> $ST
if [ $rc -eq 0 ]; then
  last=$(ls -v legA_[0-9]*.h5 | tail -1); sed -i "s|^file = legA_.*|file = $last|" legB.ini
  t0=$SECONDS
  mpirun -n 1 $BIN legB.ini > legB.log 2>&1; rc=$?
  echo "kmm legB exit=$rc wall=$((SECONDS-t0))s $(date '+%F %T')" >> $ST
fi
echo "done $(date '+%F %T')" >> $ST
