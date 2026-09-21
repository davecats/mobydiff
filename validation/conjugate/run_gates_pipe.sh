#!/usr/bin/env bash
# Pipe-campaign prerequisites F1 and F5 (docs/next_session_pipe_cht.md).
#
#   F1  [scalar.N] source_dir -- the Kasagi source on a chosen direction.
#   F2  [scalar.N] solid_thickness -- the solid's second material band, i.e.
#       the insulating jacket that makes the pipe wall a uniform-thickness
#       shell on a Cartesian grid.
#   F5  make_geometry_stl.py annulus -- the pipe's immersed body.
#
#   ./run_gates_pipe.sh [source_dir|guard|annulus|ranks|
#                        band|insulate|bandguard|bandannulus|banddet|all]
#
# Environment: BIN   (default ../../build_cpu/moby_solve)
#              PREP  (default the moby_prepare next to BIN)
#              NBIN/NPREP (nofma CPU pair, banddet group)
#              GBIN/GPREP (nofma GPU pair, banddet group)
#              RANKS (default 4, prepare only)
#              PY    (default python3; needs h5py, and ~/ibmc/bin/python here)
#
# The source_dir group runs in seconds; the annulus group prepares two 64^2x8
# cases and takes a couple of minutes (the 16384-facet BVH is the cost). The
# F2 groups run on a 4x16x4 slab and take seconds each.
#
# NOT covered here, deliberately: that source_dir OFF is bit-exact. That is
# the standard suite's job -- validation/scalar/run_bitexact{,_s3}.sh against
# a reference binary -- and it is where the "inert by construction" claim is
# actually settled, on 16 cases, CPU and GPU.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

BIN=${BIN:-$ROOT/build_cpu/moby_solve}
PREP=${PREP:-$(dirname "$BIN")/moby_prepare}
NBIN=${NBIN:-$ROOT/build_cpu_nofma/moby_solve}
NPREP=${NPREP:-$ROOT/build_cpu_nofma/moby_prepare}
GBIN=${GBIN:-$ROOT/build_gpu_nofma/moby_solve}
GPREP=${GPREP:-$ROOT/build_gpu_nofma/moby_prepare}
RANKS=${RANKS:-4}
PY=${PY:-python3}
sel=${1:-all}
status=0

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
run() { echo "   \$ $*"; "$@"; }
report() { if [ "$1" -eq 0 ]; then echo "   PASS"; else echo "   FAIL"; status=1; fi }

solve() {  # <template> <prefix> <steps> <sed-expr...>
    local tpl=$1 pre=$2 steps=$3; shift 3
    local ini=".$pre.ini"
    sed -e "s|@PREFIX@|$pre|" -e "s|@STEPS@|$steps|" "$@" "$tpl" > "$ini"
    mpirun -n 1 "$BIN" "$ini" > "$pre.log" 2>&1
}

# --- F1 (1): the closed form, all three branches ---------------------------
# Uniform (u,v,w) = (1,2,3) is an exact steady solution and a uniform scalar
# has no gradient, so theta = source * u_dir * t EXACTLY -- three different
# numbers for the three directions, which is what makes it a direction test.
if want source_dir; then
    echo "== F1 (1) source_dir vs the closed form theta = source * u_dir * t"
    for d in x y z; do
        solve kasagi_uniform.ini "ku_$d" 20 -e "s|@SDIR@|$d|" || { report 1; continue; }
    done
    run $PY ./check_source_dir.py closed-form ku_x_20.h5 ku_y_20.h5 ku_z_20.h5 \
        --components 1.0 2.0 3.0 --time 0.04
    report $?

    # --- F1 (2): the same arithmetic in x and in z, with structure ---------
    echo "== F1 (2) diagonal channel: source_dir x and z are BIT-IDENTICAL"
    for d in x z; do
        solve kasagi_channel.ini "kc_$d" 50 -e "s|@SDIR@|$d|" || { report 1; continue; }
    done
    run $PY ./check_source_dir.py identical kc_x_50.h5 kc_z_50.h5
    report $?
fi

# --- F1 (3): the config guards --------------------------------------------
if want guard; then
    echo "== F1 (3) config guards"
    # NOTE the substitution order: the direction placeholder is filled with
    # the PROBE's value, not filled with a good one and patched afterwards
    # (which silently leaves a valid ini and "passes" by accepting nothing).
    probe_ini() {  # <name> <sdir> [post-sed]
        sed -e "s|@SDIR@|$2|" -e "s|@PREFIX@|kg|" -e "s|@STEPS@|1|" kasagi_uniform.ini \
            | sed "${3:-}" > ".kg_$1.ini"
    }
    probe_ini uniform-source z 's|^source_type = velocity|source_type = uniform|'
    probe_ini bad-direction q
    for name in uniform-source bad-direction; do
        if mpirun -n 1 "$BIN" ".kg_$name.ini" > ".kg_$name.log" 2>&1; then
            echo "   $name: ACCEPTED -- it must be a hard config error"; report 1
        else
            echo "   $name: rejected -- $(grep -m1 'ERROR STOP' ".kg_$name.log")"; report 0
        fi
    done
fi

# --- F5: the annulus, through moby_prepare --------------------------------
# pipe_case <tag> <axis> <outer-flag> <outer-value>
# The grid follows the axis: 8 cells along it over 0.2 (periodic), 64 across
# the 1.3-wide cross-section. Built by substituting the FINAL values per
# direction rather than by permuting placeholders -- sed applies expressions
# in order, so a later expression matching "@N@" never fires once an earlier
# one has already written the number. That is not hypothetical: it is how the
# first version of this function ran the x case on the z case's grid.
pipe_case() {
    local tag=$1 axis=$2 oflag=$3 oval=$4
    local n=(64 64 64) l=(1.3 1.3 1.3) per=(false false false)
    local d; d=$(case $axis in x) echo 0;; y) echo 1;; z) echo 2;; esac)
    n[$d]=8; l[$d]=0.2; per[$d]=true
    $PY "$ROOT/tools/make_geometry_stl.py" annulus ".$tag.stl" --axis "$axis" --centre 0.65 0.65 \
        --r-inner 0.5 "$oflag" "$oval" --facets 16384 --a0 -0.8 --a1 1.0 \
        --domain-half 0.65 || return 1
    sed -e "s|@STL@|.$tag.stl|" -e "s|@CASE@|$tag.h5|" -e "s|@NB@|8|" \
        -e "s|^nx = @N@|nx = ${n[0]}|"  -e "s|^ny = @N@|ny = ${n[1]}|" \
        -e "s|^nz = @NZ@|nz = ${n[2]}|" -e "s|^lx = @L@|lx = ${l[0]}|" \
        -e "s|^ly = @L@|ly = ${l[1]}|"  -e "s|^lz = @LZ@|lz = ${l[2]}|" \
        -e "s|^periodic_x = .*|periodic_x = ${per[0]}|" \
        -e "s|^periodic_y = .*|periodic_y = ${per[1]}|" \
        -e "s|^periodic_z = .*|periodic_z = ${per[2]}|" \
        pipe_geom.ini | sed '/^coeff_file/d' > ".$tag.prep.ini"
    if grep -v '^;' ".$tag.prep.ini" | grep -q '@'; then
        echo "   unsubstituted placeholder in .$tag.prep.ini"; return 1
    fi
    mpirun -n "$RANKS" "$PREP" ".$tag.prep.ini" "$tag.h5" > "$tag.prep.log" 2>&1
}

if want annulus; then
    echo "== F5 (1) annulus about z, square outer surface"
    pipe_case pipe64 z --box-half 1.25 || report 1
    run $PY "$ROOT/tools/check_annulus.py" pipe64.h5 --centre 0.65 0.65 --r-inner 0.5 \
        --facets 16384 --box-half 1.25 --domain-half 0.65
    report $?

    echo "== F5 (2) annulus about x, cylindrical outer surface"
    pipe_case pipex x --r-outer 1.6 || report 1
    run $PY "$ROOT/tools/check_annulus.py" pipex.h5 --axis x --centre 0.65 0.65 --r-inner 0.5 \
        --facets 16384 --r-outer 1.6 --domain-half 0.65
    report $?
fi

if want ranks; then
    echo "== F5 (3) prepared case file is independent of the rank count"
    [ -f pipe64.h5 ] || pipe_case pipe64 z --box-half 1.25 || report 1
    mpirun -n 1 "$PREP" .pipe64.prep.ini pipe64_np1.h5 > pipe64_np1.prep.log 2>&1 || report 1
    run $PY ../prepare/h5same.py pipe64.h5 pipe64_np1.h5
    report $?
fi

# ==========================================================================
# F2: the solid's second material band
# ==========================================================================

# One banded-slab case, mirroring run_gates_c1.sh's slab_case.
#   band_case <y_wall> <kappa_s> <cap> <depth> <kappa_o> <cap_o> <init>
#             <prefix> <nsteps> <write> [bin] [prep] [ranks]
# BAND_EXTRA_SED, if set, is one more sed expression applied to the template
# (the refinement probe uses it to add a [blocks] refine box to BOTH the
# prepare and the solve input, which must describe the same leaf table).
band_case() {
    local yw=$1 ka=$2 cap=$3 dep=$4 ko=$5 co=$6 ini=$7 pre=$8 ns=$9 wr=${10}
    local bin=${11:-$BIN} prep=${12:-$PREP} nr=${13:-1}
    $PY ./make_slab_stl.py "$pre.stl" --y-top "$yw" > /dev/null || return 1
    sed -e "s|@STL@|$pre.stl|" -e "s|@CASE@|$pre.h5|" \
        -e "s|@KAPPA@|$ka|" -e "s|@CAP@|$cap|" -e "s|@DEPTH@|$dep|" \
        -e "s|@OUTERK@|$ko|" -e "s|@OUTERC@|$co|" -e "s|@INIT@|$ini|" \
        -e "s|@PREFIX@|$pre|" -e "s|@NSTEPS@|$ns|" -e "s|@WRITE@|$wr|" \
        -e "${BAND_EXTRA_SED:-s|^\\[case\\]|[case]|}" \
        band_slab.ini > ".$pre.full.ini"
    # moby_prepare COMPUTES the coefficients, so its input must not name a
    # coefficient file; the solve input takes the file and drops the STL.
    sed '/^coeff_file/d' ".$pre.full.ini" > ".$pre.prep.ini"
    sed '/^stl_file/d'   ".$pre.full.ini" > ".$pre.ini"
    mpirun -n "$nr" "$prep" ".$pre.prep.ini" "$pre.h5" > "$pre.prep.log" 2>&1 || {
        tail -5 "$pre.prep.log"; return 1; }
    mpirun -n "$nr" "$bin" ".$pre.ini" > "$pre.log" 2>&1 || { tail -20 "$pre.log"; return 1; }
    return 0
}

# --- F2 (1): the three-material slab, band boundary swept through a cell ---
# The C1 slab gate one layer deeper. The piecewise-linear three-layer profile
# is an exact fixed point ONLY if the band face carries the true series
# resistance, i.e. only if the level-set weight on psi = phi + d is the true
# cut fraction of THAT iso-surface. The gate starts at the fixed point and
# checks it does not move; y_wall is held at a half-cell cut (C1 already
# swept it) and the BAND boundary is what moves here.
if want band; then
    echo "== F2 (1) three-material slab: the exact profile is a fixed point"
    yw=0.25
    for f in 0.05 0.35 0.65 0.95; do
        dep=$($PY -c "print(repr(0.25 - (0.09375 + $f/16.0)))")
        for ko in 0.01 1.0 100.0; do
            tag="bnd_f$(echo $f | tr . p)_k$(echo $ko | tr . p)"
            band_case "$yw" 2.0 2.0 "$dep" "$ko" 1.0 0.0 "$tag" 1 1 \
                || { report 1; continue; }
            $PY ./seed_slab_ic.py "${tag}_1.h5" "${tag}_ic.h5" \
                --y-wall "$yw" --kappa 2.0 --depth "$dep" --outer-kappa "$ko" \
                > /dev/null || { report 1; continue; }
            sed -e "s|^nsteps.*|nsteps = 500|" -e "s|^field_interval.*|field_interval = 500|" \
                -e "s|^field_prefix = $tag|field_prefix = ${tag}_fp|" \
                -e "s|^\\[output\\]|[restart]\\nfile = ${tag}_ic.h5\\n\\n[output]|" \
                ".$tag.ini" > ".${tag}_fp.ini"
            mpirun -n 1 "$BIN" ".${tag}_fp.ini" > "${tag}_fp.log" 2>&1 \
                || { tail -20 "${tag}_fp.log"; report 1; continue; }
            run $PY ./check_conjugate.py slab "${tag}_fp_501.h5" \
                --y-wall "$yw" --kappa 2.0 --depth "$dep" --outer-kappa "$ko" \
                --prev "${tag}_ic.h5"
            [ $? -eq 0 ] || status=1
        done
    done
fi

# --- F2 (2): the INSULATED band, cold start ------------------------------
# kappa_outer = 0 is what the campaign uses, and it is a different statement:
# the band face coefficient must be EXACTLY zero, not small. Then no heat
# crosses the band boundary, the shell and the fluid come to the y = L
# Dirichlet value 1, and the outer band is inert at solid_init = 0.5 -- a
# value NEITHER domain face carries, so a leak in either direction shows.
# The y = 0 Dirichlet face is disconnected and must stay disconnected.
if want insulate; then
    echo "== F2 (2) insulated outer band: inert at solid_init, no leak"
    yw=0.25
    dep=$($PY -c "print(repr(0.25 - 0.125))")
    band_case "$yw" 1.0 1.0 "$dep" 0.0 1.0 0.5 ins 40000 39000 || report 1
    run $PY ./check_conjugate.py slab ins_40000.h5 \
        --y-wall "$yw" --kappa 1.0 --depth "$dep" --outer-kappa 0.0 \
        --outer-init 0.5 --prev ins_39000.h5 --tolerance 1e-12
    report $?
fi

# --- F2 (3): the config guards --------------------------------------------
if want bandguard; then
    echo "== F2 (3) config guards"
    # One prepared case file to run the probes against; the geometry is
    # irrelevant, only the config parsing and the init-time band check are.
    [ -f bg.h5 ] || band_case 0.25 1.0 1.0 0.125 0.0 1.0 0.0 bg 1 0 || report 1
    band_probe() {  # <name> <sed-expr>
        sed "$2" .bg.ini > ".bg_$1.ini"
    }
    # thinner than the grid (h = 1/16): the shell then has holes, so a face
    # joins the fluid straight to the outer band. This one is caught at INIT,
    # from the real phi field, not by a config rule.
    band_probe thin 's|^solid_thickness = .*|solid_thickness = 0.02|'
    band_probe negative 's|^solid_thickness = .*|solid_thickness = -0.1|'
    band_probe non-conjugate 's|^ibm_wall = conjugate|ibm_wall = dirichlet|'
    band_probe with-tangential \
        's|^solid_thickness = .*|solid_thickness = 0.125\ntangential_correction = true|'
    for name in thin negative non-conjugate with-tangential; do
        if mpirun -n 1 "$BIN" ".bg_$name.ini" > ".bg_$name.log" 2>&1; then
            echo "   $name: ACCEPTED -- it must be a hard config error"; report 1
        else
            echo "   $name: rejected -- $(grep -m1 'ERROR STOP' ".bg_$name.log")"; report 0
        fi
    done

    # --- the 2:1 precondition applies to the BAND boundary too -------------
    # The conjugate face coefficient is a same-level arm, so it cannot read
    # across a refinement interface -- at the body surface (C1 gate 3c) and
    # equally at a solid_thickness band boundary, which C1's check could not
    # see. Built as a PAIR so it cannot pass for the wrong reason: the body
    # surface sits at y = 0.40625, strictly inside a block row, and only the
    # BAND boundary (0.40625 - 0.15625 = 0.25) lands on the 2:1 face the
    # refine box creates. Same geometry, same box, band off -> must RUN.
    echo "   -- the 2:1 precondition at a band boundary"
    export BAND_EXTRA_SED='s|^nb = 4|nb = 4\nrefine = 0.0 0.25 0.0 0.25 0.0 0.25\nrefine_levels = 1|'
    if band_case 0.40625 1.0 1.0 0.15625 0.0 1.0 0.0 bgr 1 0 > /dev/null 2>&1; then
        echo "   band on a 2:1 face: ACCEPTED -- it must be a hard error"; report 1
    else
        echo "   band on a 2:1 face: rejected -- $(grep -m1 'ERROR STOP' bgr.log)"; report 0
    fi
    if band_case 0.40625 1.0 1.0 0.0 0.0 1.0 0.0 bgrc 1 0 > /dev/null 2>&1; then
        echo "   control, same box with the band off: runs"; report 0
    else
        echo "   control, same box with the band off: REJECTED -- the probe above"
        echo "   would then be testing the body surface, not the band"
        tail -5 bgrc.log; report 1
    fi
    unset BAND_EXTRA_SED
fi

# --- F2 (4): the band on the REAL annulus ---------------------------------
# The premise the whole feature rests on: inside the domain the level set
# -d < phi < 0 IS the annulus r_inner < r < r_inner + d. Checked against the
# ANALYTIC polygon distance, so it is a statement about the geometry rather
# than a restatement of what the solver stored.
if want bandannulus; then
    echo "== F2 (4) the level-set band IS the annulus"
    [ -f pipe64.h5 ] || pipe_case pipe64 z --box-half 1.25 || report 1
    run $PY "$ROOT/tools/check_annulus.py" pipe64.h5 --centre 0.65 0.65 --r-inner 0.5 \
        --facets 16384 --box-half 1.25 --domain-half 0.65 --band-depth 0.1
    report $?
fi

# --- F2 (5): determinism ---------------------------------------------------
# A banded run must not depend on the decomposition or the device. Both
# statements are max_abs 0, on nofma binaries.
if want banddet; then
    echo "== F2 (5) banded run: 1 == 4 ranks == GPU, exactly"
    yw=0.25
    dep=$($PY -c "print(repr(0.25 - 0.125))")
    band_case "$yw" 2.0 0.5 "$dep" 0.01 1.0 0.5 bdet_r1 200 200 "$NBIN" "$NPREP" 1 \
        || report 1
    band_case "$yw" 2.0 0.5 "$dep" 0.01 1.0 0.5 bdet_r4 200 200 "$NBIN" "$NPREP" 4 \
        || report 1
    run $PY $ROOT/tools/compare_fields.py bdet_r1_200.h5 bdet_r4_200.h5 --tolerance 0
    report $?
    if [ -x "$GBIN" ]; then
        band_case "$yw" 2.0 0.5 "$dep" 0.01 1.0 0.5 bdet_gpu 200 200 "$GBIN" "$NPREP" 1 \
            || report 1
        run $PY $ROOT/tools/compare_fields.py bdet_r1_200.h5 bdet_gpu_200.h5 --tolerance 0
        report $?
    else
        echo "   GPU binary $GBIN not found -- SKIPPED"; status=1
    fi
fi

echo
if [ $status -eq 0 ]; then echo "pipe prerequisite gates (F1, F2, F5): ALL PASS"
else echo "pipe prerequisite gates (F1, F2, F5): FAILURES"; fi
exit $status
