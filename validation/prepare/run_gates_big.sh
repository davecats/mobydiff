#!/usr/bin/env bash
# P1b gates (docs/prepare_solve_strategy.md): moby_prepare vs mobygeom on
# the BIG committed geometries -- the NACA 0012 / SD7003 5-level airfoil
# cases (keep_buried, dwall, deep lattices) and the ASCII+transform
# sailplane. Each case: prepare the case file, generate a mobygeom
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

airfoil_case() { # <dir> <label> <stl> <re>
    local dir="$1" label="$2" stl="$3" re="$4"
    (
        cd "$dir"
        [ -f "$stl" ] || PY="$PYG" ./setup.sh
        sed -e "s|^coeff_file = .*|stl_file = $stl|" \
            -e 's|^refine_levels = 4|refine_levels = 4\nkeep_buried = true|' \
            -e "s|^runtime_file = .*|runtime_file = forces_prep_unused.txt|" \
            aoa4.ini > "prep_${label}.ini"
    )
    check "$label: prepare (L5, keep_buried)" bash -c \
        "cd $dir && $RUN $BUILD/moby_prepare prep_${label}.ini prep_${label}_case.h5"
    check "$label: mobygeom reference (grid from the case file)" bash -c \
        "cd $dir && $PYG $HERE/../../tools/mobygeom.py block-table \
            --geometry $stl --grid-file prep_${label}_case.h5 --re $re \
            --block-nb 8 --levels 5 --keep-buried \
            --output ${label}_mobygeom_ref.h5 --jobs 16"
    check "$label: vs mobygeom" bash -c \
        "cd $dir && $PY $HERE/compare_case.py prep_${label}_case.h5 ${label}_mobygeom_ref.h5 --coef-tol 2e-3"
}

airfoil_case ../naca0012 n0012 n0012.stl 1.0e5
airfoil_case ../sd7003 sd7003 sd7003.stl 6.0e4

# ---- sailplane: ASCII STL + scale/translate, spaces in the path ----------
(
    cd ../../tutorials/sailplane
    sed -e 's|^coeff_file = .*|stl_file = "FRUE V0 ohneRundung.stl"\nstl_scale = 0.001\nstl_translate = 17.288649999032064 0.0 4.150549|' \
        -e 's|^\[ibm\]|[blocks]\nnb = 10\n\n[ibm]|' input.ini > prep_blocks.ini
)
check "sailplane: prepare (ASCII + transform)" bash -c \
    "cd ../../tutorials/sailplane && mpirun -n 2 --oversubscribe $BUILD/moby_prepare prep_blocks.ini prep_sailplane_case.h5"
check "sailplane: mobygeom reference (grid from the case file)" bash -c \
    "cd ../../tutorials/sailplane && $PYG $HERE/../../tools/mobygeom.py block-table \
        --geometry 'FRUE V0 ohneRundung.stl' --grid-file prep_sailplane_case.h5 \
        --re 1.0e5 --scale 0.001 --translate 17.288649999032064 0.0 4.150549 \
        --block-nb 10 --levels 1 --no-dwall \
        --output sailplane_ref_blocks.h5 --jobs 4"
check "sailplane: vs mobygeom (2 grazing outliers of 93M)" bash -c \
    "cd ../../tutorials/sailplane && $PY $HERE/compare_case.py prep_sailplane_case.h5 sailplane_ref_blocks.h5 --coef-tol 2e-3"
check "sailplane: solve prep == committed legacy (1 step)" bash -c '
    cd ../../tutorials/sailplane
    sed -e "s|^coeff_file = .*|coeff_file = prep_sailplane_case.h5|" \
        -e "s|^\[ibm\]|[blocks]\nnb = 10\n\n[ibm]|" \
        -e "s|^field_prefix = .*|field_prefix = sp_prep|" input.ini > solve_prep.ini
    sed -e "s|^\[ibm\]|[blocks]\nnb = 10\n\n[ibm]|" \
        -e "s|^field_prefix = .*|field_prefix = sp_ref|" input.ini > solve_ref.ini
    mpirun -n 2 --oversubscribe '"$BUILD"'/moby_solve solve_ref.ini >/dev/null 2>&1 &&
    mpirun -n 2 --oversubscribe '"$BUILD"'/moby_solve solve_prep.ini >/dev/null 2>&1 &&
    '"$PY"' ../../tools/compare_fields.py --tolerance 0 sp_ref_1.h5 sp_prep_1.h5 un vn wn pn'

echo
echo "passed: $pass  failed: $fail"
echo "(short solve-equivalence gates for naca/sd7003 run separately on the"
echo " GPU -- see validation/prepare/README.md P1b status)"
[ "$fail" -eq 0 ]
