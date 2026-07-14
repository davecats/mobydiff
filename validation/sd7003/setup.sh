#!/usr/bin/env bash
# Prerequisites for the A3 INCREMENT 3 SD7003 transition benchmark:
# spline-resampled Selig STL + the 5-level block-table coefficient file
# (dwall tiles for RANS; --keep-buried is load-bearing for the forces —
# see validation/naca0012/README.md).
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-$HOME/ibmc/bin/python}"
ROOT=../..
MOBYGRID="${MOBYGRID:-$ROOT/build_cpu/mobygrid}"

echo "== 1. grid.h5 (mobygrid on the run grid)"
mpirun -n 1 "$MOBYGRID" aoa4.ini grid.h5

echo "== 2. SD7003 STL (UIUC Selig coordinates, resampled to 720 points)"
$PY "$ROOT/tools/make_airfoil_stl.py" selig --file sd7003.dat \
    --xc 4.5 --yc 6.0 --chord 1.0 --lz 0.1875 --resample 720 --out sd7003.stl

echo "== 3. 5-level block-table coefficient file (keep-buried)"
$PY "$ROOT/tools/mobygeom.py" block-table \
    --geometry sd7003.stl --grid-file grid.h5 --re 6.0e4 \
    --block-nb 8 --levels 5 --keep-buried --output ibm_coeff_sd7003.h5

echo "setup done"
