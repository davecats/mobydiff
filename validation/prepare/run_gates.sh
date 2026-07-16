#!/usr/bin/env bash
# P0 gates for moby_prepare (docs/prepare_solve_strategy.md): the solver run
# from a prepared case file must be BIT-EXACT vs its inline analytic run,
# and the case file must be independent of the prepare rank count.
#
# Usage:  ./run_gates.sh [build_dir]        (default ../../build_cpu_nofma)
# Needs the nofma CPU build (main + moby_prepare) and python3+h5py.
set -uo pipefail

cd "$(dirname "$0")"
BUILD="${1:-../../build_cpu_nofma}"
PY="${PY:-python3}"
CMP="$PY ../../tools/compare_fields.py --tolerance 0"
pass=0; fail=0

check() { # <name> <cmd...>
    local name="$1"; shift
    if "$@" >/dev/null 2>&1; then
        echo "PASS  $name"; pass=$((pass+1))
    else
        echo "FAIL  $name"; fail=$((fail+1))
    fi
}

twin() { # <base.ini> <case.h5> <new_prefix> -> writes <base>_file.ini
    local base="$1" case="$2" prefix="$3"
    sed -e "s|^enabled = true|enabled = true\ncoeff_file = ${case}|" \
        -e "s|^field_prefix = .*|field_prefix = ${prefix}|" \
        "${base}.ini" > "${base}_file.ini"
}

run_case() { # <base> <steps-suffix-of-final-field> <ransgeom:0|1>
    local base="$1" step="$2" rans="$3"
    local case="${base}_case.h5"
    local pfx pfxf
    pfx="$(grep '^field_prefix' "${base}.ini" | awk '{print $3}')"
    pfxf="${pfx}f"

    rm -f "$case" "${base}_case_np4.h5" "${pfx}_${step}.h5" "${pfxf}_${step}.h5" \
        "${pfx}_ransgeom.h5" "${pfxf}_ransgeom.h5" "${base}_file.ini"

    # prepare (1 rank) + the rank-count-independence twin (4 ranks)
    check "${base}: prepare (1 rank)" mpirun -n 1 "$BUILD/moby_prepare" "${base}.ini" "$case"
    check "${base}: prepare (4 ranks)" mpirun -n 4 --oversubscribe "$BUILD/moby_prepare" "${base}.ini" "${base}_case_np4.h5"
    check "${base}: case file 1==4 ranks" $PY h5same.py "$case" "${base}_case_np4.h5"

    # inline analytic reference vs the solve from the prepared file
    twin "$base" "$case" "$pfxf"
    check "${base}: inline solve" mpirun -n 1 "$BUILD/main" "${base}.ini"
    check "${base}: solve from case file" mpirun -n 1 "$BUILD/main" "${base}_file.ini"
    check "${base}: fields bit-exact" $CMP "${pfx}_${step}.h5" "${pfxf}_${step}.h5" un vn wn pn
    if [ "$rans" = 1 ]; then
        check "${base}: ransgeom identical" $PY h5same.py "${pfx}_ransgeom.h5" "${pfxf}_ransgeom.h5"
    fi
}

run_case wavy 1 1
run_case wavy_refine 1 1
run_case wavysolid 5 0

# Optional GPU cross-check (GPU_BUILD=../../build_gpu_nofma ./run_gates.sh):
# the GPU solve from the CPU-prepared case file vs the CPU solve from the
# same file must match at tolerance 0.
if [ -n "${GPU_BUILD:-}" ]; then
    sed -e "s|^field_prefix = wavyrf|field_prefix = wavyrG|" wavy_refine_file.ini \
        > wavy_refine_file_gpu.ini
    check "wavy_refine: GPU solve from case file" \
        mpirun -n 1 "$GPU_BUILD/main" wavy_refine_file_gpu.ini
    check "wavy_refine: GPU==CPU from case file" \
        $CMP wavyrf_1.h5 wavyrG_1.h5 un vn wn pn
fi

echo
echo "passed: $pass  failed: $fail"
[ "$fail" -eq 0 ]
