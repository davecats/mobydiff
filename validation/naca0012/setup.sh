#!/usr/bin/env bash
# Prerequisites for the A3 INCREMENT 2 NACA 0012 SST sanity runs: extruded
# NACA 0012 STL + the 5-level block-table coefficient file (with dwall
# tiles for RANS) on the 512x512x8 base grid.
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-$HOME/ibmc/bin/python}"
ROOT=../..
MOBYGRID="${MOBYGRID:-$ROOT/build_cpu/mobygrid}"

echo "== 1. grid.h5 (mobygrid on the run grid)"
mpirun -n 1 "$MOBYGRID" aoa4.ini grid.h5

echo "== 2. NACA 0012 STL (chord 1, LE at (4.5, 6.0), extruded past both z faces)"
$PY "$ROOT/tools/make_airfoil_stl.py" naca --code 0012 --xc 4.5 --yc 6.0 \
    --chord 1.0 --lz 0.1875 --out n0012.stl

echo "== 3. 5-level block-table coefficient file"
# --keep-buried is LOAD-BEARING for the forces: with the buried core
# removed, its closed faces absorb the body's pressure loading outside the
# coef bookkeeping and sum(coef u dV) loses the (pressure-dominated) lift —
# measured C_L 0.018 vs the flow's actual circulation-implied 0.37 at
# aoa = 4 (the drag, friction-dominated, stayed plausible).
$PY "$ROOT/tools/mobygeom.py" block-table \
    --geometry n0012.stl --grid-file grid.h5 --re 1.0e5 \
    --block-nb 8 --levels 5 --keep-buried --output ibm_coeff_n0012.h5

echo "setup done"
