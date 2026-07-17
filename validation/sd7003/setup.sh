#!/usr/bin/env bash
# Prerequisites for the SD7003 transition benchmark: spline-resampled
# Selig STL + the 5-level case file (coefficients + masks + dwall tiles;
# keep_buried is load-bearing for the forces -- see
# validation/naca0012/README.md). Since the prepare/solve split (P3) the
# case file is written by moby_prepare -- no venv/mobygeom involved (only
# the STL generator still needs the geometry venv). The mobygeom
# cross-check reference lives in validation/prepare/run_gates_big.sh.
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-$HOME/ibmc/bin/python}"
ROOT=../..
PREPARE="${PREPARE:-$ROOT/build_cpu/moby_prepare}"

echo "== 1. SD7003 STL (UIUC Selig coordinates, resampled to 720 points)"
$PY "$ROOT/tools/make_airfoil_stl.py" selig --file sd7003.dat \
    --xc 4.5 --yc 6.0 --chord 1.0 --lz 0.1875 --resample 720 --out sd7003.stl

echo "== 2. 5-level case file (moby_prepare, keep_buried)"
sed -e 's|^coeff_file = .*|stl_file = sd7003.stl|' \
    -e 's|^refine_levels = 4|refine_levels = 4\nkeep_buried = true|' \
    aoa4.ini > prep_sd7003.ini
mpirun -n 4 --oversubscribe --bind-to none -x OMP_NUM_THREADS=4 \
    "$PREPARE" prep_sd7003.ini ibm_coeff_sd7003.h5

echo "setup done"
