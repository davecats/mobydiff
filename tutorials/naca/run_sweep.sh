#!/usr/bin/env bash
# Run one or more sweep angles SEQUENTIALLY on one host/GPU:
#   run_sweep.sh corax -2 -1 0 1     (build_gpu_corax on istmcorax)
#   run_sweep.sh cetus 2 3           (build_gpu on istmcetus GPU1)
#   run_sweep.sh local 4 5           (build_gpu on the local GPU)
# Angle tags: -2 -> m2 (file names forces_aoam2.txt etc.).
# Launch detached from the workstation, e.g.:
#   ssh istmcorax 'setsid bash -l <this file> corax -2 -1 0 1 &' (cetus
#   NEEDS the login shell -l; see memory remote-hosts).
set -uo pipefail
host=$1; shift
cd "$(dirname "$0")"
case $host in
  corax)
    export NV=/opt/Nvidia/nvhpc/Linux_x86_64/25.9
    export PATH=$NV/compilers/bin:$NV/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
    BIN=../../build_gpu_corax/main ;;
  cetus)
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load toolkits/nvhpc/25.9
    export CUDA_VISIBLE_DEVICES=1
    BIN=../../build_gpu/main ;;
  local)
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load toolkits/nvhpc/25.9
    BIN=../../build_gpu/main ;;
  *) echo "host must be corax|cetus|local" >&2; exit 1 ;;
esac

for a in "$@"; do
  tag=${a/-/m}
  sed -e "s/__AOA__/${a}.0/" -e "s/__TAG__/${tag}/" naca_base.ini > .aoa_${tag}.ini
  echo "== aoa ${a} (${host})"
  mpirun -n 1 "$BIN" .aoa_${tag}.ini > aoa_${tag}.log 2>&1
  echo "NACA_AOA_${tag}_DONE exit=$?" >> aoa_${tag}.log
done
echo "SWEEP_${host}_DONE"
