#!/usr/bin/env bash
# Prerequisites for the Re 4e5 NACA 0012 polar sweep: span-y STL + grid +
# the 7-level (refine_levels = 6) xz block-table coefficient file with
# dwall tiles (RANS) and KEEP-BURIED interiors (penalization forces).
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-$HOME/ibmc/bin/python}"
ROOT=../..
MOBYGRID="${MOBYGRID:-$ROOT/build_cpu/mobygrid}"

echo "== 1. grid.h5 (mobygrid on the run grid; any angle works, grid only)"
sed -e 's/__AOA__/0.0/' -e 's/__TAG__/setup/' naca_base.ini > .grid_setup.ini
mpirun -n 1 "$MOBYGRID" .grid_setup.ini grid.h5

echo "== 2. NACA 0012 STL, span y (chord x, LIFT z; LE at (4.5, z 6.0))"
$PY "$ROOT/tools/make_airfoil_stl.py" naca --code 0012 --xc 4.5 --yc 6.0 \
    --chord 1.0 --lz 0.1875 --span y --out n0012_spany.stl

echo "== 3. 7-level xz block-table coefficient file (keep-buried, dwall)"
$PY "$ROOT/tools/mobygeom.py" block-table \
    --geometry n0012_spany.stl --grid-file grid.h5 --re 4.0e5 \
    --block-nb 8 --levels 7 --refine-dims xz --keep-buried \
    --output ibm_coeff_n0012_re4e5_xz_l7.h5

echo "setup done"
