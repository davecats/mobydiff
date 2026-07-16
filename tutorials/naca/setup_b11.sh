#!/usr/bin/env bash
# B11 prerequisites: fine STL (n = 2880: nose facet sagitta 0.25 fine
# cells), grid, 12-level xz block-table WITHOUT keep-buried (interior
# removed; forces from the external control volume) and the SAME
# per-level refine boxes as b11_base.ini (the solver builder
# cross-checks the table row-by-row).
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-$HOME/ibmc/bin/python}"
ROOT=../..
MOBYGRID="${MOBYGRID:-$ROOT/build_cpu/mobygrid}"

echo "== 1. grid_b11.h5"
sed -e 's/__AOA__/0.0/' -e 's/__TAG__/setup/' b11_base.ini > .b11_setup.ini
mpirun -n 1 "$MOBYGRID" .b11_setup.ini grid_b11.h5

echo "== 2. NACA 0012 STL, span y, n = 2880, nose (50, 48)"
$PY "$ROOT/tools/make_airfoil_stl.py" naca --code 0012 --xc 50.0 --yc 48.0 \
    --chord 1.0 --lz 0.1875 --n 2880 --span y --out n0012_b11.stl

echo "== 3. 12-level xz block-table (remove-buried, dwall, per-level boxes)"
$PY "$ROOT/tools/mobygeom.py" block-table \
    --geometry n0012_b11.stl --grid-file grid_b11.h5 --re 4.0e5 \
    --block-nb 8 --levels 12 --refine-dims xz --jobs "${JOBS:-20}" \
    --refine-box 49.5 51.5 -1.0 10.0 47.44 48.56 7 \
    --refine-box 48.5 53.5 -1.0 10.0 46.44 49.56 6 \
    --refine-box 47.5 55.5 -1.0 10.0 45.44 50.56 5 \
    --refine-box 46.5 57.5 -1.0 10.0 44.44 51.56 4 \
    --refine-box 45.5 59.5 -1.0 10.0 43.44 52.56 3 \
    --refine-box 44.5 61.5 -1.0 10.0 42.44 53.56 2 \
    --refine-box 43.5 63.5 -1.0 10.0 41.44 54.56 1 \
    --output ibm_coeff_b11.h5

echo "setup done"
