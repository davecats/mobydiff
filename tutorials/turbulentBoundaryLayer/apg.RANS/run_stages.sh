#!/bin/bash
# stage 1 (coarse) -> interpolation -> stage 2 (production) on the local GPU
module load toolkits/nvhpc/25.9 2>/dev/null
export OMPI_MCA_hwloc_base_binding_policy=none
# the solver needs the [boundary] y_max_p_profile = table feature
B=${MOBY_SOLVE:-../../../build_gpu/moby_solve}
P=../smooth.RANS/assets/postpro
cd "$(dirname "$0")"
mpirun -n 1 $B coarse.ini > coarse.log 2>&1 || { echo "coarse exit $?" >> stages.txt; exit 1; }
last=$(ls coarse_[0-9]*.h5 | sort -t_ -k2 -n | tail -1); echo "coarse done: $last" >> stages.txt
python3 $P/interp_restart.py $last ../smooth.RANS/production.case.h5 production_ic.h5 >> stages.txt 2>&1 || exit 1
mpirun -n 1 $B production.ini > production.log 2>&1; echo "production exit $?" >> stages.txt
