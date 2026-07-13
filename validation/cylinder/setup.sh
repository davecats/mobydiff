#!/usr/bin/env bash
# Prerequisites for the cylinder A1/A2 gates: extruded cylinder STL + per-Re
# IBM coefficient files (coef = SOLID/re, so one file per Reynolds number).
# Requires build_cpu/mobygrid and the geometry venv (trimesh + shapely + h5py).
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-/home/davide/ibmc/bin/python}"
ROOT=../..
MOBYGRID="${MOBYGRID:-$ROOT/build_cpu/mobygrid}"

echo "== 1. grid.h5 (mobygrid on the run grid)"
mpirun -n 1 "$MOBYGRID" cyl_re40.ini grid.h5

echo "== 2. cylinder STL (D = 1 at (6.0, 8.02), extruded past both z faces)"
$PY "$ROOT/tools/make_airfoil_stl.py" cylinder --xc 6.0 --yc 8.02 --d 1.0 \
    --lz 0.25 --out cylinder.stl

echo "== 3. coefficient files"
for re in 40 100; do
    $PY "$ROOT/tools/mobygeom.py" stl-ibm-coeff \
        --geometry cylinder.stl --grid-file grid.h5 --re $re \
        --output ibm_coeff_re${re}.h5 --no-tiled-output
done
echo "setup done"
