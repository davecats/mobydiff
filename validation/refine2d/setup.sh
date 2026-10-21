#!/usr/bin/env bash
# The zero-coefficient xz-refined body twin `ibm_ml3_xz_zero.h5` of unibody_xz.ini
# (not in git): multilevel_body/setup.sh step 2 on unibody_xz.ini -- moby_prepare
# with the cylinder STL and keep_buried, then zero_coef.py. An archived file in the
# legacy block order is refused by the solver since 2026-10-01 ("blocks table is
# not in its key order"), so regenerate rather than copy.
#   [PREP=/abs/moby_prepare] [PY=python] ./setup.sh
set -euo pipefail
cd "$(dirname "$0")"
PREP="${PREP:-../../build_cpu/moby_prepare}"; PY="${PY:-python3}"
sed -e 's|^\[ibm\]|[ibm]\nstl_file = ../cylinder/cylinder.stl|' unibody_xz.ini > .prep_unibody_xz.ini
grep -q '^keep_buried = true' .prep_unibody_xz.ini || { echo "keep_buried not set"; exit 1; }
rm -f ibm_ml3_xz_zero.h5
mpirun -n 1 "$PREP" .prep_unibody_xz.ini ibm_ml3_xz_zero.h5 > .prep_unibody_xz.log 2>&1 || { tail -5 .prep_unibody_xz.log; exit 1; }
"$PY" ../multilevel_body/zero_coef.py ibm_ml3_xz_zero.h5
grep -h "block refinement\|leaves" .prep_unibody_xz.log | tail -2
