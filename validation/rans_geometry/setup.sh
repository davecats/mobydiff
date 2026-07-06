#!/usr/bin/env bash
# Generate the block-table coefficient files (with dwall_blocks) for the
# RANS T1 geometry gates, from the les_ibm wall STLs and grid.
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-/home/davide/ibmc/bin/python}"
ROOT=../..
LES_IBM=$ROOT/validation/channel_interface/les_ibm

echo "== 1. single-level block-table file (flat_l1.ini)"
$PY "$ROOT/tools/mobygeom.py" block-table \
    --geometry "$LES_IBM/wall_lo.stl" "$LES_IBM/wall_hi.stl" \
    --grid-file "$LES_IBM/grid.h5" --re 180 \
    --block-nb 8 --levels 1 --output ibm_coeff_blocks_l1.h5

echo "== 2. two-level block-table file (flat_refine.ini)"
$PY "$ROOT/tools/mobygeom.py" block-table \
    --geometry "$LES_IBM/wall_lo.stl" "$LES_IBM/wall_hi.stl" \
    --grid-file "$LES_IBM/grid.h5" --re 180 \
    --block-nb 8 --levels 2 --output ibm_coeff_blocks_l2.h5

echo "== 3. the regenerated l2 file must differ from the committed les_ibm"
echo "==    file ONLY by the added dwall_blocks dataset"
$PY - <<'EOF'
import h5py
import numpy as np

new = h5py.File("ibm_coeff_blocks_l2.h5", "r")
old = h5py.File("../channel_interface/les_ibm/ibm_coeff_blocks.h5", "r")
for name in old:
    a, b = old[name][...], new[name][...]
    assert a.shape == b.shape and np.array_equal(a, b), f"dataset {name} differs"
extra = set(new) - set(old)
assert extra == {"dwall_blocks"}, f"unexpected extra datasets: {extra}"
print("coef_blocks/masks/blocks identical; dwall_blocks added")
EOF

echo "== done"
