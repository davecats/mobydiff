#!/usr/bin/env bash
# Prerequisites for the NACA 0012 SST runs: extruded NACA 0012 STL + the
# 5-level case file (coefficients + masks + dwall tiles for RANS) on the
# 512x512x8 base grid. Since the prepare/solve split (P3) the case file is
# written by moby_prepare -- no venv/mobygeom involved (only the STL
# generator still needs the geometry venv for shapely/mapbox_earcut). The
# mobygeom cross-check reference lives in
# validation/prepare/run_gates_big.sh.
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-$HOME/ibmc/bin/python}"
ROOT=../..
PREPARE="${PREPARE:-$ROOT/build_cpu/moby_prepare}"

echo "== 1. NACA 0012 STL (chord 1, LE at (4.5, 6.0), extruded past both z faces)"
$PY "$ROOT/tools/make_airfoil_stl.py" naca --code 0012 --xc 4.5 --yc 6.0 \
    --chord 1.0 --lz 0.1875 --out n0012.stl

echo "== 2. 5-level case file (moby_prepare)"
# keep_buried is LOAD-BEARING for the forces: with the buried core
# removed, its closed faces absorb the body's pressure loading outside the
# coef bookkeeping and sum(coef u dV) loses the (pressure-dominated) lift --
# measured C_L 0.018 vs the flow's actual circulation-implied 0.37 at
# aoa = 4 (the drag, friction-dominated, stayed plausible).
sed -e 's|^coeff_file = .*|stl_file = n0012.stl|' \
    -e 's|^refine_levels = 4|refine_levels = 4\nkeep_buried = true|' \
    aoa4.ini > prep_n0012.ini
mpirun -n 4 --oversubscribe --bind-to none -x OMP_NUM_THREADS=4 \
    "$PREPARE" prep_n0012.ini ibm_coeff_n0012.h5

echo "setup done"
