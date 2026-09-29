#!/usr/bin/env bash
# P1b gates (docs/prepare_solve_strategy.md): moby_prepare vs mobygeom on
# the BIG committed geometries. Only the ASCII+transform sailplane remains:
# the NACA 0012 / SD7003 5-level airfoil legs were removed with
# validation/naca0012 and validation/sd7003 (2026-09-26; their P1b results
# stay recorded in README.md, the drivers are in git history). Each case: prepare the case file, generate a mobygeom
# block-table reference FROM THE CASE FILE's grid datasets (P3: the case
# file replaced the mobygrid handshake), compare (blocks + masks
# identical, zero solid-classification flips, graded coef to a per-case
# tolerance, interior dwall to round-off class), then a short solve.
#
# EXPENSIVE (mobygeom references take tens of minutes; prepares ~4-6 min
# each) and needs the ibmc venv. Run pieces by hand when iterating -- this
# script records the exact procedure.
#
# Usage: ./run_gates_big.sh [build_dir]       (default ../../build_cpu)
set -uo pipefail

cd "$(dirname "$0")"
BUILD="$(realpath "${1:-../../build_cpu}")"
HERE="$(realpath .)"
PY="${PY:-python3}"
PYG="${PYG:-$HOME/ibmc/bin/python}"
RUN="mpirun -n 4 --oversubscribe --bind-to none -x OMP_NUM_THREADS=4"
pass=0; fail=0

check() {
    local name="$1"; shift
    if "$@" >/dev/null 2>&1; then
        echo "PASS  $name"; pass=$((pass+1))
    else
        echo "FAIL  $name"; fail=$((fail+1))
    fi
}

# ---- sailplane: ASCII STL + scale/translate, spaces in the path ----------
# Since numerics review step 7 the tutorial ini itself carries [blocks] nb =
# 10 and the STL declaration ([case] file = sailplane_case.h5): the gate
# prepares it into its own file, compares against the mobygeom reference,
# and checks the solve is rank-count independent (the legacy file the old
# solve leg compared against was retired with the legacy reader).
check "sailplane: prepare (ASCII + transform)" bash -c \
    "cd ../../tutorials/sailplane && mpirun -n 2 --oversubscribe $BUILD/moby_prepare input.ini prep_sailplane_case.h5"
check "sailplane: mobygeom reference (grid from the case file)" bash -c \
    "cd ../../tutorials/sailplane && $PYG $HERE/../../tools/mobygeom.py block-table \
        --geometry 'FRUE V0 ohneRundung.stl' --grid-file prep_sailplane_case.h5 \
        --re 1.0e5 --scale 0.001 --translate 17.288649999032064 0.0 4.150549 \
        --block-nb 10 --levels 1 --no-dwall \
        --output sailplane_ref_blocks.h5 --jobs 4"
check "sailplane: vs mobygeom (2 grazing outliers of 93M)" bash -c \
    "cd ../../tutorials/sailplane && $PY $HERE/compare_case.py prep_sailplane_case.h5 sailplane_ref_blocks.h5 --coef-tol 2e-3"
check "sailplane: solve from the prepared file, 1 == 2 ranks (1 step)" bash -c '
    cd ../../tutorials/sailplane
    for r in 1 2; do
        sed -e "s|^file = .*|file = prep_sailplane_case.h5|" \
            -e "s|^field_prefix = .*|field_prefix = sp_r$r|" input.ini > solve_r$r.ini
        mpirun -n $r --oversubscribe '"$BUILD"'/moby_solve solve_r$r.ini >/dev/null 2>&1 || exit 1
    done
    '"$PY"' ../../tools/compare_fields.py --tolerance 0 sp_r1_1.h5 sp_r2_1.h5 un vn wn pn'

echo
echo "passed: $pass  failed: $fail"
[ "$fail" -eq 0 ]
