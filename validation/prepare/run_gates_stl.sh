#!/usr/bin/env bash
# P1 gates for moby_prepare's STL geometry (docs/prepare_solve_strategy.md):
# the STL indicator (geometry_stl.f90) drives the SAME machinery as the
# analytic path, so the gates compare against mobygeom's committed/generated
# block-table references and an exact shift-invariance twin.
#
#   flat / flat_refine : les_ibm wall slabs vs the COMMITTED mobygeom files
#                        (validation/rans_geometry/ibm_coeff_blocks_l{1,2}.h5)
#                        -- blocks + masks identical, coef <= 1e-6 rel,
#                        interior dwall <= 2e-9 -- plus a 1-step solve from
#                        the prepared file vs the committed file.
#   sphere             : curved body vs a freshly generated mobygeom
#                        reference (needs the ibmc venv: PYG, default
#                        ~/ibmc/bin/python; skipped if absent).
#   sphere_shift       : the SAME mesh exactly float32-translated onto the
#                        x-periodic boundary; masks/blocks/coef/dwall must
#                        be the exactly rolled copy (gates minimum-image).
#
# Usage: ./run_gates_stl.sh [build_dir]      (default ../../build_cpu_nofma)
set -uo pipefail

cd "$(dirname "$0")"
BUILD="${1:-../../build_cpu_nofma}"
PY="${PY:-python3}"
PYG="${PYG:-$HOME/ibmc/bin/python}"
pass=0; fail=0

check() {
    local name="$1"; shift
    if "$@" >/dev/null 2>&1; then
        echo "PASS  $name"; pass=$((pass+1))
    else
        echo "FAIL  $name"; fail=$((fail+1))
    fi
}

# ---- flat slabs vs the committed mobygeom references --------------------
for base in flat flat_refine; do
    ref=../rans_geometry/ibm_coeff_blocks_l1.h5
    [ "$base" = flat_refine ] && ref=../rans_geometry/ibm_coeff_blocks_l2.h5
    rm -f "${base}_case.h5"
    check "$base: prepare (4 ranks)" mpirun -n 4 --oversubscribe "$BUILD/moby_prepare" "$base.ini" "${base}_case.h5"
    check "$base: vs committed mobygeom" $PY compare_case.py "${base}_case.h5" "$ref"
done

# 1-step solve: prepared file vs the committed mobygeom file (coef differ
# at the bisection tolerance, so fields/ransgeom match tightly, not to 0).
pfx_ref=flatS; pfx_new=flatP
sed -e "s|^stl_file = .*|coeff_file = ../rans_geometry/ibm_coeff_blocks_l1.h5|" \
    -e "s|^field_prefix = .*|field_prefix = ${pfx_ref}|" flat.ini > flat_solve_ref.ini
sed -e "s|^stl_file = \(.*\)|stl_file = \1\ncoeff_file = flat_case.h5|" \
    -e "s|^field_prefix = .*|field_prefix = ${pfx_new}|" flat.ini > flat_solve_new.ini
check "flat: solve from committed file" mpirun -n 1 "$BUILD/main" flat_solve_ref.ini
check "flat: solve from prepared file" mpirun -n 1 "$BUILD/main" flat_solve_new.ini
check "flat: 1-step fields match (1e-10)" $PY ../../tools/compare_fields.py \
    --tolerance 1e-10 "${pfx_ref}_1.h5" "${pfx_new}_1.h5" un vn wn pn
check "flat: ransgeom dwall match (1e-9)" $PY - <<EOF
import h5py, numpy as np, sys
a = h5py.File("${pfx_ref}_ransgeom.h5","r"); b = h5py.File("${pfx_new}_ransgeom.h5","r")
dd = max(float(np.max(np.abs(a[k][...]-b[k][...]))) for k in ("dwall","yeff"))
wc = int(np.sum(a["wallcell"][...] != b["wallcell"][...]))
print("dwall/yeff max", dd, "wallcell diffs", wc)
sys.exit(0 if dd <= 1e-9 and wc == 0 else 1)
EOF

# ---- sphere vs mobygeom + exact shift-invariance ------------------------
if [ -x "$PYG" ] && "$PYG" -c "import igl, trimesh" >/dev/null 2>&1; then
    rm -f sphere.stl sphere_shift.stl sphere_grid.h5 sphere_ref.h5 sphere_case.h5 spheres_case.h5
    check "sphere: STL" "$PYG" ../../tools/mobygeom.py make-sphere-stl \
        --output sphere.stl --center 0.5 0.5 0.5 --radius 0.2
    check "sphere: shifted STL (exact)" $PY - <<'EOF'
import struct
import numpy as np
with open('sphere.stl','rb') as f:
    data = bytearray(f.read())
n = struct.unpack_from('<I', data, 80)[0]
for t in range(n):
    base = 84 + 50*t
    vals = list(struct.unpack_from('<12f', data, base))
    for vi in (3, 6, 9):
        vals[vi] = float(np.float32(vals[vi]) - np.float32(0.5))
    struct.pack_into('<12f', data, base, *vals)
with open('sphere_shift.stl','wb') as f:
    f.write(bytes(data))
EOF
    check "sphere: mobygrid" mpirun -n 1 "$BUILD/mobygrid" sphere.ini sphere_grid.h5
    check "sphere: mobygeom reference" "$PYG" ../../tools/mobygeom.py block-table \
        --geometry sphere.stl --grid-file sphere_grid.h5 --re 100 \
        --block-nb 8 --levels 2 --output sphere_ref.h5 --jobs 8
    check "sphere: prepare" mpirun -n 4 --oversubscribe "$BUILD/moby_prepare" sphere.ini sphere_case.h5
    check "sphere: vs mobygeom" $PY compare_case.py sphere_case.h5 sphere_ref.h5
    check "sphere_shift: prepare" mpirun -n 4 --oversubscribe "$BUILD/moby_prepare" sphere_shift.ini spheres_case.h5
    check "sphere_shift: exact roll of centred case" $PY shift_check.py sphere_case.h5 spheres_case.h5 32
else
    echo "SKIP  sphere gates (no ibmc venv at $PYG)"
fi

echo
echo "passed: $pass  failed: $fail"
[ "$fail" -eq 0 ]
