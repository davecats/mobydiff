#!/usr/bin/env bash
# Prerequisites for the A3 INCREMENT 0 multi-level refine_body gates:
# a 3-level block-table coefficient file for the committed cylinder STL
# (validation/cylinder/cylinder.stl, D = 1 at (6.0, 8.02)) on the 128x128x8
# gate grid, plus the zero-force twin the uniform-flow gate needs.
#
# REPAIRED 2026-09-26. The original built these with `build_cpu/mobygrid` and
# `tools/mobygeom.py block-table`. mobygrid was DELETED in the prepare/solve
# split P3 and mobygeom's geometry subcommands retired there, so this script
# could not run at all -- and because `*.h5` is gitignored, NEITHER coefficient
# file was ever committed. That left three gates unrunnable on any checkout:
# both gates in this directory and `validation/scalar/uniform3.ini`, which reads
# the twin. It is now two `moby_prepare` runs and a five-line zeroing.
#
# Needs only the CPU build (moby_prepare has its own STL reader since P1) and
# h5py. The geometry venv is NOT needed any more.
#
#   ./run_gates.sh          after this
set -euo pipefail
cd "$(dirname "$0")"

PREP="${PREP:-../../build_cpu/moby_prepare}"
PY="${PY:-python3}"
[ -x "$PREP" ] || { echo "no moby_prepare at $PREP (build with ./compile.sh cpu)" >&2; exit 1; }

# The prepare inis are DERIVED from the gate inis rather than written out, so
# the grid, nb and refine_levels cannot drift out of sync with the cases that
# consume the files -- the failure the original script was one edit away from.
# `stl_file` replaces `coeff_file`: prepare reads the geometry, the solver reads
# the file prepare writes.
stl_ini() {   # stl_ini <src.ini> <dst.ini> [extra blocks key]
    sed -e 's|^coeff_file = .*|stl_file = ../cylinder/cylinder.stl|' "$1" > "$2"
    # NOT `[ -n "$3" ] && sed ...`: under `set -e` that construct returns 1
    # from the whole function when the optional argument is absent, and the
    # script dies after step 1 with no message. It did exactly that once.
    if [ -n "${3:-}" ]; then
        sed -i "s|^refine_body = true|refine_body = true\n$3|" "$2"
    fi
}

echo "== 1. the real 3-level file (buried leaves removed, graded coefficients)"
stl_ini dwall.ini .prep_ml3.ini
mpirun -n 1 "$PREP" .prep_ml3.ini ibm_coeff_ml3.h5

# The twin keeps the SAME touch-driven refinement -- the interfaces under test --
# but must not break a constant field, which needs two things:
#   * buried leaves KEPT, so no closed faces appear in the flow. `[blocks]
#     keep_buried` zeroes the buried masks in classify_refinement_masks, which
#     is exactly what mobygeom's --keep-buried wrote by hand;
#   * coef = 0, so the immersed boundary exerts no force. Nothing in the config
#     can express that, hence zero_coef.py.
# The original twin had to rebuild the leaf table in Python (importing
# mobygeom's build_leaf_table_py) because it zeroed the masks AFTER the table
# was written. Here prepare writes a table that already matches the zeroed
# masks, so the Python step is reduced to zeroing one dataset -- which also
# removes the last mobygeom import from this directory.
echo "== 2. the zero-force twin (buried leaves kept, then coefficients zeroed)"
stl_ini dwall.ini .prep_ml3_zero.ini "keep_buried = true"
mpirun -n 1 "$PREP" .prep_ml3_zero.ini ibm_coeff_ml3_zero.h5
$PY zero_coef.py ibm_coeff_ml3_zero.h5

rm -f .prep_ml3.ini .prep_ml3_zero.ini
echo "setup done: ibm_coeff_ml3.h5 (real), ibm_coeff_ml3_zero.h5 (zero-force twin)"
