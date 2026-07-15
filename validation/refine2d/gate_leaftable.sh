#!/usr/bin/env bash
# R2D-0 gate (b): the mobygeom Python leaf builder (--refine-dims) must
# match the solver's build_leaf_table row-by-row, xz AND xyz, on a
# non-cubic box-refined lattice. The STL sphere sits OUTSIDE the domain,
# so body classification runs (all-fluid masks) without adding
# refinement beyond the box — the solver side (leaftable_test) takes no
# geometry input.
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-$HOME/ibmc/bin/python}"
ROOT=../..
MOBYGRID="${MOBYGRID:-$ROOT/build_cpu/mobygrid}"
LEAFTEST="${LEAFTEST:-$ROOT/build_cpu/leaftable_test}"

# refine box (physical): a corner-crossing slab, 3 levels (refine_levels 2)
BOX="0.5 0.9 0.3 0.7 0.8 1.4"
NB=8
LEVELS=3

echo "== 1. grid.h5 (mobygrid on the gate grid)"
mpirun -n 1 "$MOBYGRID" leaf_grid.ini grid.h5

echo "== 2. out-of-domain sphere STL"
$PY "$ROOT/tools/mobygeom.py" make-sphere-stl --output sphere_out.stl \
    --center 0.75 -5.0 1.0 --radius 0.2

fail=0
for dims in xz xyz; do
  echo "== 3. mobygeom block-table --refine-dims $dims"
  $PY "$ROOT/tools/mobygeom.py" block-table \
      --geometry sphere_out.stl --grid-file grid.h5 \
      --block-nb $NB --levels $LEVELS --refine-box $BOX \
      --refine-dims $dims --no-dwall --output leaf_$dims.h5

  echo "== 4. solver leaf table ($dims)"
  mpirun -n 1 "$LEAFTEST" 48 32 64 1.5 1.0 2.0 $NB $((LEVELS-1)) $dims $BOX \
      > leaftable_$dims.out

  $PY compare_leaftable.py leaftable_$dims.out leaf_$dims.h5 $dims || fail=1
done

[ $fail -eq 0 ] && echo "GATE B PASS" || { echo "GATE B FAIL"; exit 1; }
