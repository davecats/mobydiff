#!/usr/bin/env bash
# C1 gates: CONJUGATE HEAT TRANSFER at the immersed interface
# (docs/next_session_conjugate.md Section 10, increment C1).
#
#   ./run_gates_c1.sh [slab|converge|weight|capacity|contact|peclet|limits|conserve|guard|refine|stats|det|all]
#
# Environment: BIN   (default ../../build_cpu/moby_solve)
#              PREP  (default the moby_prepare next to BIN)
#              NBIN/NPREP (nofma CPU pair, det group; default build_cpu_nofma)
#              GBIN  (nofma GPU binary, det group; default build_gpu_nofma)
#              RANKS (default 1)
#
# Everything here runs in minutes on one core.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

BIN=${BIN:-$ROOT/build_cpu/moby_solve}
PREP=${PREP:-$(dirname "$BIN")/moby_prepare}
NBIN=${NBIN:-$ROOT/build_cpu_nofma/moby_solve}
NPREP=${NPREP:-$ROOT/build_cpu_nofma/moby_prepare}
GBIN=${GBIN:-$ROOT/build_gpu_nofma/moby_solve}
RANKS=${RANKS:-1}
PY=${PY:-python3}
CMP="$PY $ROOT/tools/compare_fields.py --tolerance 0"
sel=${1:-all}
status=0

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
run() { echo "   \$ $*"; "$@"; }
report() { if [ "$1" -eq 0 ]; then echo "   PASS"; else echo "   FAIL"; status=1; fi }

# One slab case: <y_wall> <kappa> <capacity> <prefix> <nsteps> <write> [bin] [prep] [ranks]
slab_case() {
    local yw=$1 ka=$2 cap=$3 pre=$4 ns=$5 wr=$6
    local bin=${7:-$BIN} prep=${8:-$PREP} nr=${9:-$RANKS}
    $PY ./make_slab_stl.py "$pre.stl" --y-top "$yw" > /dev/null || return 1
    sed -e "s|@STL@|$pre.stl|" -e "s|@CASE@|$pre.h5|" \
        -e "s|@KAPPA@|$ka|" -e "s|@CAP@|$cap|" -e "s|@PREFIX@|$pre|" \
        -e "s|@NSTEPS@|$ns|" -e "s|@WRITE@|$wr|" slab.ini > ".$pre.full.ini"
    # moby_prepare COMPUTES the coefficient file, so its input must not name
    # one; the solve input takes the file and drops the STL.
    sed '/^coeff_file/d' ".$pre.full.ini" > ".$pre.prep.ini"
    sed '/^stl_file/d'   ".$pre.full.ini" > ".$pre.ini"
    mpirun -n "$nr" "$prep" ".$pre.prep.ini" "$pre.h5" > "$pre.prep.log" 2>&1 || {
        tail -5 "$pre.prep.log"; return 1; }
    mpirun -n "$nr" "$bin" ".$pre.ini" > "$pre.log" 2>&1 || { tail -20 "$pre.log"; return 1; }
    return 0
}

# --- (1) the 1D two-material slab: EXACT -----------------------------------
# The cut position is swept through a full cell (the last solid cell centre
# sits at y = 3.5 h = 0.21875, h = 1/16, so the cut arm here also happens to
# sit ON a block boundary -- the ghost-inclusive-tile case) and kappa_s over
# five decades. The steady solution is piecewise linear, which is an exact
# FIXED POINT of the discrete operator ONLY if the cut face's series
# resistance is exact -- i.e. only if the level-set weight w is the true cut
# fraction. Nothing else in this gate can absorb an error in w.
#
# The gate STARTS at the fixed point and checks it does not move: any error
# in w changes k_face and the profile leaves at once, at a rate set by the
# error rather than by the slowest eigenmode. The `converge` group is the
# companion statement that the solver also REACHES this profile from a cold
# start (which takes O(10^5) steps at the extreme contrasts, hence the split).
if want slab; then
    echo "== (1) 1D two-material slab: the exact profile is a fixed point"
    for w in 0.05 0.20 0.35 0.50 0.65 0.80 0.95; do
        yw=$($PY -c "print(repr(0.21875 + $w/16.0))")
        for ka in 0.01 1.0 10.0 1000.0; do
            tag=$($PY -c "print('slab_w%s_k%s' % ('$w'.replace('.','p'), '$ka'.replace('.','p')))")
            slab_case "$yw" "$ka" "$ka" "$tag" 1 1 || { report 1; continue; }
            $PY ./seed_slab_ic.py "${tag}_1.h5" "${tag}_ic.h5" \
                --y-wall "$yw" --kappa "$ka" > /dev/null || { report 1; continue; }
            sed -e "s|^nsteps.*|nsteps = 500|" -e "s|^field_interval.*|field_interval = 500|" \
                -e "s|^field_prefix = $tag|field_prefix = ${tag}_fp|" \
                -e "s|^\\[output\\]|[restart]\\nfile = ${tag}_ic.h5\\n\\n[output]|" \
                ".$tag.ini" > ".${tag}_fp.ini"
            mpirun -n "$RANKS" "$BIN" ".${tag}_fp.ini" > "${tag}_fp.log" 2>&1 \
                || { tail -20 "${tag}_fp.log"; report 1; continue; }
            run $PY ./check_conjugate.py slab "${tag}_fp_501.h5" \
                --y-wall "$yw" --kappa "$ka" --prev "${tag}_ic.h5"
            [ $? -eq 0 ] || status=1
        done
    done
fi

# --- (1a) ...and the solver REACHES it from a cold start -------------------
# One cold-start run per kappa decade, from a uniform field, long enough that
# the residual between the last two writes is at round-off. This is the
# statement the fixed-point gate cannot make.
if want converge; then
    echo "== (1a) cold start converges to the same exact profile"
    yw=$($PY -c "print(repr(0.21875 + 0.5/16.0))")
    for ka in 0.01 1.0 10.0; do
        tag="cvg_$(echo $ka | tr . p)"
        slab_case "$yw" "$ka" "$ka" "$tag" 200000 190000 || { report 1; continue; }
        run $PY ./check_conjugate.py slab "${tag}_200000.h5" \
            --y-wall "$yw" --kappa "$ka" --prev "${tag}_190000.h5" --tolerance 1e-11
        [ $? -eq 0 ] || status=1
    done
fi

# --- (1b) the level-set weight, straight out of the case file --------------
# Independent of the solver: phi = +-dwall_blocks signed by coef_p_blocks, so
# w is rebuilt from the two datasets the scheme actually reads and compared
# with the analytic cut position of the STL plane.
if want weight; then
    echo "== (1b) level-set weight w vs the analytic cut position"
    for w in 0.05 0.50 0.95; do
        yw=$($PY -c "print(repr(0.21875 + $w/16.0))")
        tag="wt_$(echo $w | tr . p)"
        slab_case "$yw" 1.0 1.0 "$tag" 1 0 || { report 1; continue; }
        run $PY ./check_conjugate.py weight "$tag.h5" --y-wall "$yw"
        [ $? -eq 0 ] || status=1
    done
fi

# --- (1c) the capacity is irrelevant at steady state -----------------------
# Same kappa_s, three different C_s: the steady field must be the SAME. It
# also proves the capacity division is where it belongs (on the flux
# divergence) rather than folded into the conductivity.
if want capacity; then
    echo "== (1c) steady state independent of the solid capacity"
    yw=$($PY -c "print(repr(0.21875 + 0.35/16.0))")
    for cap in 0.5 1.0 8.0; do
        tag="cap_$(echo $cap | tr . p)"
        slab_case "$yw" 10.0 "$cap" "$tag" 150000 150000 || { report 1; continue; }
        run $PY ./check_conjugate.py slab "${tag}_150000.h5" --y-wall "$yw" --kappa 10.0 --tolerance 1e-11
        [ $? -eq 0 ] || status=1
    done
fi

# --- (1d) the material-max Peclet limiter ----------------------------------
# alpha_s/alpha_f = kappa_s/C_s = 200 with the limiter ON: dt must shrink by
# that factor and the run must stay bounded. The control is the SAME case
# with pecletmax = 0, which is expected to blow up -- if it does not, the
# gate is not testing anything.
if want peclet; then
    echo "== (1d) Peclet limiter over materials (alpha_s/alpha_f = 200)"
    yw=$($PY -c "print(repr(0.21875 + 0.5/16.0))")
    slab_case "$yw" 200.0 1.0 "pec_on" 4000 4000
    rc=$?
    if [ $rc -eq 0 ]; then

        mx=$($PY -c "
import h5py,numpy as np
print('%.6e' % np.abs(h5py.File('pec_on_4000.h5')['theta'][...]).max())")
        echo "   limiter on: bounded, max|theta| = $mx"
        rc=$($PY -c "print(0 if float('$mx') < 10.0 else 1)")
    fi
    report $rc
    sed -e 's|^pecletmax.*|pecletmax = 0.0|' -e 's|^dtmax.*|dtmax = 1.0e-2|' \
        -e 's|^dt = .*|dt = 1.0e-2|' .pec_on.ini > .pec_off.ini
    sed -i 's|field_prefix = pec_on|field_prefix = pec_off|' .pec_off.ini
    mpirun -n 1 "$BIN" .pec_off.ini > pec_off.log 2>&1
    mxo=$($PY -c "
import h5py,numpy as np,sys
try: print('%.6e' % np.abs(h5py.File('pec_off_4000.h5')['theta'][...]).max())
except Exception: print('inf')" 2>/dev/null)
    echo "   control (pecletmax = 0, dt = 1e-2): max|theta| = $mxo  (must be unbounded)"
    $PY -c "
v = '$mxo'
import math
bad = (v == 'inf') or not math.isfinite(float(v)) or float(v) > 10.0
raise SystemExit(0 if bad else 1)"
    report $?
fi

# --- (1e) contact resistance -----------------------------------------------
# [scalar.N] contact_resistance adds R_c to the series resistance of the cut
# face. The steady solution stays PIECEWISE LINEAR -- with a jump q*R_c at
# the interface -- so it stays an exact fixed point and the gate is the same
# equality as (1). R_c = 0 must reproduce the plain harmonic mean exactly,
# which is the arithmetic identity the key is worth having.
if want contact; then
    echo "== (1e) contact resistance: the jump q*R_c, exactly"
    yw=$($PY -c "print(repr(0.21875 + 0.35/16.0))")
    for rc in 0.0 0.25 4.0; do
        tag="rc_$(echo $rc | tr . p)"
        slab_case "$yw" 10.0 10.0 "$tag" 1 1 || { report 1; continue; }
        sed -i "s|^solid_init.*|solid_init = 0.0\ncontact_resistance = $rc|" ".$tag.ini"
        $PY ./seed_slab_ic.py "${tag}_1.h5" "${tag}_ic.h5" \
            --y-wall "$yw" --kappa 10.0 --contact "$rc" > /dev/null || { report 1; continue; }
        sed -e "s|^nsteps.*|nsteps = 500|" -e "s|^field_interval.*|field_interval = 500|" \
            -e "s|^field_prefix = $tag|field_prefix = ${tag}_fp|" \
            -e "s|^\\[output\\]|[restart]\\nfile = ${tag}_ic.h5\\n\\n[output]|" \
            ".$tag.ini" > ".${tag}_fp.ini"
        mpirun -n "$RANKS" "$BIN" ".${tag}_fp.ini" > "${tag}_fp.log" 2>&1 \
            || { tail -20 "${tag}_fp.log"; report 1; continue; }
        run $PY ./check_conjugate.py slab "${tag}_fp_501.h5" \
            --y-wall "$yw" --kappa 10.0 --contact "$rc" --prev "${tag}_ic.h5"
        [ $? -eq 0 ] || status=1
    done
fi

# --- (3d) the statistics branch runs on a conjugate case -------------------
# scalar_stats.f90's y-face diffusivity and convective mask were extended
# with the SAME helper and the SAME branch as the transport kernel, so the
# rows keep reporting the flux the kernel applied (that invariant is the
# whole reason S4's accumulators are built from the transport expression).
# C1 SMOKE ONLY: the case runs, on CPU and on GPU, and every row is finite.
# The QUANTITATIVE gate on the conjugate flux columns belongs to C3, the
# Nusselt increment, which is where the strategy doc puts it.
if want stats; then
    echo "== (3d) scalar_stats on a conjugate case (smoke; C3 gates the numbers)"
    sed -e 's|^\[scalar\]|[scalar]\nstats_sample_interval = 20\nstats_write_interval = 20\nstats_file = wavy_stats.h5|' \
        -e 's|^nsteps.*|nsteps = 20|' -e 's|^field_interval.*|field_interval = 20|' \
        -e 's|field_prefix = wavy|field_prefix = wstat|' wavy.ini > .wstat.ini
    for tag in cpu gpu; do
        bin=$BIN; [ "$tag" = gpu ] && bin=$GBIN
        sed "s|stats_file = wavy_stats.h5|stats_file = wavy_stats_$tag.h5|" .wstat.ini > ".wstat_$tag.ini"
        mpirun -n 1 "$bin" ".wstat_$tag.ini" > "wstat_$tag.log" 2>&1
        if [ $? -ne 0 ]; then tail -12 "wstat_$tag.log"; report 1; continue; fi
        $PY -c "
import h5py, numpy as np, sys
f = h5py.File('wavy_stats_$tag.h5')
bad = 0
for k in f:
    d = f[k][...]
    if np.issubdtype(d.dtype, np.floating) and not np.isfinite(d).all(): bad += 1
print('   $tag: datasets', len(list(f)), ' non-finite', bad)
sys.exit(1 if bad else 0)"
        report $?
    done
fi

# --- (2) the two limits ----------------------------------------------------
# kappa_s -> infinity must reproduce the `dirichlet` mode (the solid becomes
# isothermal at the y_min wall value, so the fluid sees a Dirichlet interface)
# and kappa_s -> 0 the `adiabatic` one. Both are LIMITS, not identities: the
# S3 modes place their effective boundary on the STAIRCASE (a penalised solid
# cell centre / a masked staggered face) while the conjugate interface sits at
# its true position, so the residual is the O(h) discretisation difference and
# the gate is that it behaves like one. The kappa dependence is the sharp
# half: the interface flux must go to its limit like kappa_s.
if want limits; then
    echo "== (2) kappa_s -> infinity == dirichlet, kappa_s -> 0 == adiabatic"
    yw=$($PY -c "print(repr(0.21875 + 0.5/16.0))")
    for ka in 1.0e3 1.0e5 1.0e7; do
        tag="lim_hi_$(echo $ka | tr .+ pp)"
        slab_case "$yw" "$ka" "$ka" "$tag" 30000 30000 || { report 1; continue; }
        # The dirichlet twin: the body held at the y_min wall value.
        sed -e "s|ibm_wall = conjugate|ibm_wall = dirichlet|" \
            -e "/^solid_k/d" -e "/^solid_rhocp/d" -e "s|^solid_init.*|ibm_value = 0.0|" \
            -e "s|field_prefix = $tag|field_prefix = ${tag}_ref|" ".$tag.ini" > ".${tag}_ref.ini"
        mpirun -n 1 "$BIN" ".${tag}_ref.ini" > "${tag}_ref.log" 2>&1
        run $PY ./check_conjugate.py limit "${tag}_30000.h5" "${tag}_ref_30000.h5" \
            --y-wall "$yw" --tolerance 0.05
        [ $? -eq 0 ] || status=1
    done
    # C_s = kappa_s here (alpha_s = 1) for the same reason the high group
    # does it, and it is LOAD-BEARING at low kappa_s: with C_s = 1 the solid
    # has alpha_s = kappa_s and equilibrates over t ~ L^2/alpha_s = 6e3, so a
    # 30k-step run measures the solid's CHARGING flux, not the steady one --
    # measured |q| came out 7x the closed-form value and the kappa rate read
    # 0.58 instead of 1.
    for ka in 1.0e-3 1.0e-5 1.0e-7; do
        tag="lim_lo_$(echo $ka | tr .+- ppm)"
        slab_case "$yw" "$ka" "$ka" "$tag" 30000 30000 || { report 1; continue; }
        sed -e "s|ibm_wall = conjugate|ibm_wall = adiabatic|" \
            -e "/^solid_k/d" -e "/^solid_rhocp/d" -e "/^solid_init/d" \
            -e "s|field_prefix = $tag|field_prefix = ${tag}_ref|" ".$tag.ini" > ".${tag}_ref.ini"
        mpirun -n 1 "$BIN" ".${tag}_ref.ini" > "${tag}_ref.log" 2>&1
        run $PY ./check_conjugate.py limit "${tag}_30000.h5" "${tag}_ref_30000.h5" \
            --y-wall "$yw" --tolerance 0.05
        [ $? -eq 0 ] || status=1
    done
    echo "   -- interface flux vs kappa_s (the rate is the sharp statement):"
    $PY ./flux_limit.py --y-wall "$yw" \
        --high lim_hi_1p0e3_30000.h5 lim_hi_1p0e5_30000.h5 lim_hi_1p0e7_30000.h5 \
        --high-kappa 1e3 1e5 1e7 \
        --low lim_lo_1p0em3_30000.h5 lim_lo_1p0em5_30000.h5 lim_lo_1p0em7_30000.h5 \
        --low-kappa 1e-3 1e-5 1e-7
    report $?
fi

# --- (3) conservation ------------------------------------------------------
# Insulated composite box on the OBLIQUE analytic wavy wall, with the flow on
# so the masked convective flux is part of the statement. sum(C theta dV)
# must be conserved to round-off; the reference capacity map is the analytic
# wavy wall, so it is independent of the solver's own marker.
if want conserve; then
    echo "== (3) sum(C theta dV) conserved in an insulated composite box"
    run mpirun -n "$RANKS" "$BIN" wavy.ini > wavy.log 2>&1
    if [ $? -ne 0 ]; then tail -20 wavy.log; report 1; else
        run $PY ./check_conjugate.py conserve wavy_100.h5 wavy_200.h5 \
            --capacity 2.0 --wavy
        [ $? -eq 0 ] || status=1
    fi
fi

# --- (3b) the config guards ------------------------------------------------
# Every one of these silently produces a wrong answer if it is allowed
# through, so each must be a HARD ERROR (strategy doc Section 12).
if want guard; then
    echo "== (3b) config guards"
    guard() {  # <label> <sed expr...>
        local label=$1; shift
        local ini=".guard.ini"
        cp wavy.ini "$ini"
        for e in "$@"; do sed -i "$e" "$ini"; done
        mpirun -n 1 "$BIN" "$ini" > guard.log 2>&1
        if [ $? -ne 0 ]; then echo "   $label: rejected  PASS"; else
            echo "   $label: ACCEPTED  FAIL"; status=1; fi
    }
    guard "remove_solid = true"     's|^remove_solid.*|remove_solid = true|'
    guard "ibm_value with conjugate" 's|^solid_init.*|ibm_value = 1.0|'
    guard "solid_k without conjugate" 's|^ibm_wall.*|ibm_wall = dirichlet|' '/^solid_init/d'
    guard "no immersed body"        's|^enabled = true|enabled = false|'
    guard "solid_k <= 0"            's|^solid_k.*|solid_k = 0.0|'
fi

# --- (3c) the 2:1 precondition, both ways ----------------------------------
# The cut-face coefficient is a SAME-LEVEL two-point arm: on a coarse/fine
# block face the halo value is a restriction or a prolong of the other level
# and neither phi nor the flux means what the scheme assumes. refine_body's
# one-block 26-neighbour buffer is what guarantees this never happens, and
# the solver CHECKS it rather than assuming it (strategy doc Section 12).
#   negative: a hand-placed refinement box whose 2:1 face cuts the wall must
#             be a hard error naming the fix;
#   positive: refine_body + keep_buried puts the whole surface inside the
#             finest level, so the same body runs -- and that is also the
#             multi-level conjugate path's smoke test.
if want refine; then
    echo "== (3c) cut faces must not sit on a 2:1 block face"
    # The box must put a 2:1 face THROUGH the surface, and the surface here
    # is nearly horizontal -- an x-plane interface carries no cut arm however
    # far it cuts across the domain (measured: 0 bad faces). nb = 4 makes the
    # y block rows 0.03125 tall, and refining only the bottom row puts a 2:1
    # y-face at y = 0.03125, which the 0.010 .. 0.035 wall crosses.
    sed -e 's|^nb = 8|nb = 4\nrefine = 0.0 1.0 0.0 0.03125 0.0 0.25\nrefine_levels = 1|' \
        -e 's|^nsteps.*|nsteps = 5|' -e 's|^field_interval.*|field_interval = 5|' \
        -e 's|field_prefix = wavy|field_prefix = ref_bad|' wavy.ini > .ref_bad.ini
    mpirun -n 1 "$BIN" .ref_bad.ini > ref_bad.log 2>&1
    if [ $? -ne 0 ] && grep -q "2:1 block face" ref_bad.log; then
        echo "   hand-placed refinement box across the wall: rejected  PASS"
    else
        echo "   hand-placed refinement box across the wall: ACCEPTED  FAIL"; status=1
    fi
    sed -e 's|^nb = 8|nb = 8\nrefine_body = true\nkeep_buried = true\nrefine_levels = 1|' \
        -e 's|^nsteps.*|nsteps = 20|' -e 's|^field_interval.*|field_interval = 20|' \
        -e 's|field_prefix = wavy|field_prefix = ref_ok|' wavy.ini > .ref_ok.ini
    run mpirun -n 1 "$BIN" .ref_ok.ini > ref_ok.log 2>&1
    if [ $? -eq 0 ]; then
        grep -E "conjugate interface:" ref_ok.log | sed 's/^/  /'
        echo "   refine_body + keep_buried: runs  PASS"
    else
        tail -12 ref_ok.log; echo "   refine_body + keep_buried: FAILED"; status=1
    fi
    # ...and the same case WITHOUT keep_buried must be a hard config error.
    sed 's|^keep_buried = true|keep_buried = false|' .ref_ok.ini > .ref_nb.ini
    mpirun -n 1 "$BIN" .ref_nb.ini > ref_nb.log 2>&1
    if [ $? -ne 0 ]; then echo "   refine_body without keep_buried: rejected  PASS"
    else echo "   refine_body without keep_buried: ACCEPTED  FAIL"; status=1; fi
fi

# --- (4) determinism -------------------------------------------------------
# 1 rank == 4 ranks and CPU == GPU at TOLERANCE 0, on both geometry paths.
if want det; then
    echo "== (4) determinism: 1 == 4 ranks, CPU == GPU (tolerance 0)"
    sed -e 's|^nsteps.*|nsteps = 20|' -e 's|^field_interval.*|field_interval = 20|' \
        -e 's|^field_prefix = wavy|field_prefix = det_r1|' wavy.ini > .det_r1.ini
    sed 's|field_prefix = det_r1|field_prefix = det_r4|' .det_r1.ini > .det_r4.ini
    sed 's|field_prefix = det_r1|field_prefix = det_gpu|' .det_r1.ini > .det_gpu.ini
    run mpirun -n 1 "$NBIN" .det_r1.ini  > det_r1.log  2>&1
    run mpirun -n 4 "$NBIN" .det_r4.ini  > det_r4.log  2>&1
    run mpirun -n 1 "$GBIN" .det_gpu.ini > det_gpu.log 2>&1
    echo "   analytic path, 1 vs 4 ranks:"
    run $CMP det_r1_20.h5 det_r4_20.h5; [ $? -eq 0 ] || status=1
    echo "   analytic path, CPU vs GPU:"
    run $CMP det_r1_20.h5 det_gpu_20.h5; [ $? -eq 0 ] || status=1

    # File path: prepared once with the CPU nofma binary (the canonical
    # preprocessor), solved on 1 rank, 4 ranks and the GPU from that one file.
    yw=$($PY -c "print(repr(0.21875 + 0.5/16.0))")
    slab_case "$yw" 7.0 3.0 "detf_r1" 20 20 "$NBIN" "$NPREP" 1
    if [ $? -ne 0 ]; then report 1; else
        sed 's|field_prefix = detf_r1|field_prefix = detf_r4|' .detf_r1.ini > .detf_r4.ini
        sed 's|field_prefix = detf_r1|field_prefix = detf_gpu|' .detf_r1.ini > .detf_gpu.ini
        run mpirun -n 4 "$NBIN" .detf_r4.ini  > detf_r4.log  2>&1
        run mpirun -n 1 "$GBIN" .detf_gpu.ini > detf_gpu.log 2>&1
        echo "   file path, 1 vs 4 ranks:"
        run $CMP detf_r1_20.h5 detf_r4_20.h5; [ $? -eq 0 ] || status=1
        echo "   file path, CPU vs GPU:"
        run $CMP detf_r1_20.h5 detf_gpu_20.h5; [ $? -eq 0 ] || status=1
    fi
fi

echo
if [ $status -eq 0 ]; then echo "ALL C1 GATES PASS"; else echo "SOME C1 GATES FAILED"; fi
exit $status
