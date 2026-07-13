#!/usr/bin/env bash
# Prerequisites for the A3 INCREMENT 0 multi-level refine_body gates:
# 3-level block-table coefficient file for the committed cylinder STL
# (validation/cylinder/cylinder.stl, D = 1 at (6.0, 8.02)) on the small
# 128x128x8 gate grid, plus the zero-force twin for the uniform-flow gate.
# Requires build_cpu/mobygrid and the geometry venv (trimesh + h5py).
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-$HOME/ibmc/bin/python}"
ROOT=../..
MOBYGRID="${MOBYGRID:-$ROOT/build_cpu/mobygrid}"

echo "== 1. grid.h5 (mobygrid on the gate grid)"
mpirun -n 1 "$MOBYGRID" dwall.ini grid.h5

echo "== 2. 3-level block-table file (touch/buried masks, coef, dwall tiles)"
$PY "$ROOT/tools/mobygeom.py" block-table \
    --geometry "$ROOT/validation/cylinder/cylinder.stl" \
    --grid-file grid.h5 --re 40 \
    --block-nb 8 --levels 3 --output ibm_coeff_ml3.h5

echo "== 3. zero-force twin (uniform-flow gate)"
$PY make_uniform_twin.py ibm_coeff_ml3.h5 grid.h5 ibm_coeff_ml3_zero.h5

echo "setup done"
