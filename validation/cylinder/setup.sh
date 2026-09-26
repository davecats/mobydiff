#!/usr/bin/env bash
# Prerequisites for the cylinder A1/A2 gates: extruded cylinder STL + per-Re
# IBM coefficient files (coef = SOLID/re, so one file per Reynolds number).
#
# UPDATED 2026-09-26. This used to run build_cpu/mobygrid + tools/mobygeom.py
# stl-ibm-coeff (legacy global-grid files); mobygrid was DELETED in the
# prepare/solve split P3 and mobygeom's geometry subcommands retired, so the
# script could not run. The coefficient files now come from moby_prepare (the
# case ini with coeff_file swapped for stl_file). `remove_solid = false` keeps
# the blocks buried inside the cylinder, as the legacy files did (they carried
# no block_active table, so nothing was removed).
#
# Requires build_cpu/moby_prepare; the geometry venv (trimesh + shapely +
# mapbox_earcut) only if cylinder.stl is missing (it is committed).
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-/home/davide/ibmc/bin/python}"
ROOT=../..
PREP="${PREP:-$ROOT/build_cpu/moby_prepare}"

if [ ! -f cylinder.stl ]; then
    echo "== 1. cylinder STL (D = 1 at (6.0, 8.02), extruded past both z faces)"
    $PY "$ROOT/tools/make_airfoil_stl.py" cylinder --xc 6.0 --yc 8.02 --d 1.0 \
        --lz 0.25 --out cylinder.stl
fi

echo "== 2. coefficient files (moby_prepare)"
for re in 40 100; do
    sed -e 's|^coeff_file = .*|stl_file = cylinder.stl|' \
        -e 's|^nb = 8$|nb = 8\nremove_solid = false|' \
        cyl_re${re}.ini > .prep_re${re}.ini
    mpirun -n 1 "$PREP" .prep_re${re}.ini ibm_coeff_re${re}.h5
    rm -f .prep_re${re}.ini
done
echo "setup done"
