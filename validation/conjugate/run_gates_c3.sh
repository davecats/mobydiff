#!/usr/bin/env bash
# C3 gates: the FLUID-FRACTION-WEIGHTED CAPACITY, the conjugate NUSSELT
# diagnostic and the TIME-STEP decision
# (docs/next_session_conjugate.md Section 10, increment C3).
#
#   ./run_gates_c3.sh [fraction|transient|nusselt|budget|channel|guard|all]
#
# Environment: BIN   (default ../../build_cpu/moby_solve)
#              PREP  (the moby_prepare next to BIN)
#              NBIN/NPREP/GBIN (the nofma triple, determinism)
#              CBIN/CPREP (the archived C2 pair, ~/c2_ref_binaries -- the
#                          POINTWISE-capacity control of the transient gate)
#              RANKS (default 1)
#
# The three items share one cost: the capacity and the time-step convention
# both move dt at cut cells, so C1's and C2's suites are re-run after them
# (that is not done here -- run ./run_gates_c1.sh and ./run_gates_c2.sh).
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

BIN=${BIN:-$ROOT/build_cpu/moby_solve}
PREP=${PREP:-$(dirname "$BIN")/moby_prepare}
NBIN=${NBIN:-$ROOT/build_cpu_nofma/moby_solve}
NPREP=${NPREP:-$ROOT/build_cpu_nofma/moby_prepare}
GBIN=${GBIN:-$ROOT/build_gpu_nofma/moby_solve}
CBIN=${CBIN:-$HOME/c2_ref_binaries/moby_solve_cpu}
CPREP=${CPREP:-$HOME/c2_ref_binaries/moby_prepare_cpu}
RANKS=${RANKS:-1}
PY=${PY:-python3}
CMP="$PY $ROOT/tools/compare_fields.py --tolerance 0"
sel=${1:-all}
status=0

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
run() { echo "   \$ $*"; "$@"; }
report() { if [ "$1" -eq 0 ]; then echo "   PASS"; else echo "   FAIL"; status=1; fi }

# --- (0) the closed form itself, against brute force -----------------------
# plane_box_fraction is the one new piece of arithmetic in the increment, and
# it is checked where it cannot share a mistake with its own derivation: the
# unit test compares it with a 300^3 midpoint count, and with hand-derivable
# exact values in each of the three degeneracy regimes (the plane parallel to
# two axes -- a grid-aligned wall, the common case -- to one, and to none).
if want fraction; then
    echo "== (0) the plane-in-box fluid fraction vs brute force (unit test)"
    run mpirun -n 1 "$(dirname "$BIN")/scalar_test"
    report $?
fi

# --- (0b) the one new config guard -----------------------------------------
# The interface-heat diagnostic reports the BASELINE cut-face flux. With the
# C2 tangential correction on, the kernel applies a different flux at those
# faces, and a diagnostic that does not report the flux the kernel applied is
# worse than none (the invariant S4 exists for). The correction ships disabled
# by measurement, so the combination is rejected rather than carrying a
# second, ungated copy of its six-face stencil in the statistics.
if want guard; then
    echo "== (0b) heat_interval with tangential_correction is a hard error"
    sed -e 's|^ibm_wall = conjugate|ibm_wall = conjugate\ntangential_correction = true|' \
        -e 's|^\[time\]|[scalar]\nheat_interval = 10\n\n[time]|' wavy.ini > .g6.ini
    mpirun -n 1 "$BIN" .g6.ini > g6.log 2>&1
    if [ $? -ne 0 ] && grep -q "heat_interval with tangential_correction" g6.log; then
        echo "   rejected  PASS"
    else
        echo "   ACCEPTED  FAIL"; status=1
    fi
fi

# --- (1) the transient two-material slab -----------------------------------
# One decaying eigenmode of the two-material slab, seeded exactly and compared
# at t_end: field error, interface flux (through the solver's own heat
# diagnostic) and the realised DECAY RATE. Refined 16 -> 32 -> 64 with dt
# FIXED at ~h^2/100 so the RK3 temporal error (order 3, i.e. h^6 here) cannot
# be mistaken for the spatial one.
#
# The control is the archived C2 binary on the SAME inis: pointwise capacity,
# everything else identical. That is the measurement the increment turns on.
transient_case() {   # <ny> <kappa> <cap> <tag> <bin> <prep>
    local ny=$1 ka=$2 cap=$3 pre=$4 bin=$5 prep=$6
    local yw=0.31 tend=0.05 dt ns wr
    dt=$($PY -c "print(repr(0.01/$ny**2))")
    ns=$($PY -c "print(int(round($tend/(0.01/$ny**2))))")
    dt=$($PY -c "print(repr($tend/$ns))")
    wr=$((ns / 2))
    $PY ./make_slab_stl.py "$pre.stl" --y-top "$yw" > /dev/null || return 1
    sed -e "s|@STL@|$pre.stl|" -e "s|@CASE@|$pre.h5|" -e "s|@NY@|$ny|" \
        -e "s|@KAPPA@|$ka|" -e "s|@CAP@|$cap|" -e "s|@PREFIX@|$pre|" \
        -e "s|@NSTEPS@|1|" -e "s|@WRITE@|1|" -e "s|@DT@|$dt|" \
        -e "s|@HEAT@|-1|" -e "s|@HEATFILE@|$pre.heat.txt|" transient.ini > ".$pre.full.ini"
    sed '/^coeff_file/d' ".$pre.full.ini" > ".$pre.prep.ini"
    sed '/^stl_file/d'   ".$pre.full.ini" > ".$pre.ini"
    mpirun -n "$RANKS" "$prep" ".$pre.prep.ini" "$pre.h5" > "$pre.prep.log" 2>&1 || {
        tail -5 "$pre.prep.log"; return 1; }
    # step 1 only, to mint a snapshot with the right layout to seed into
    mpirun -n 1 "$bin" ".$pre.ini" > "$pre.log" 2>&1 || { tail -20 "$pre.log"; return 1; }
    $PY ./check_transient.py seed "${pre}_1.h5" "${pre}_ic.h5" \
        --y-wall "$yw" --kappa "$ka" --cap "$cap" > /dev/null || return 1
    # ...and the run proper, from the seeded mode
    sed -e "s|^nsteps.*|nsteps = $ns|" -e "s|^field_interval.*|field_interval = $ns|" \
        -e "s|^heat_interval.*|heat_interval = $wr|" \
        -e "s|^field_prefix = $pre|field_prefix = ${pre}_r|" \
        -e "s|^\\[output\\]|[restart]\\nfile = ${pre}_ic.h5\\n\\n[output]|" \
        ".$pre.ini" > ".${pre}_r.ini"
    rm -f "$pre.heat.txt"
    mpirun -n "$RANKS" "$bin" ".${pre}_r.ini" > "${pre}_r.log" 2>&1 || {
        tail -20 "${pre}_r.log"; return 1; }
    return 0
}

if want transient; then
    echo "== (1) transient two-material slab: the capacity's own gate"
    ka=10.0; cap=4.0
    $PY ./check_transient.py mode --y-wall 0.31 --kappa $ka --cap $cap
    : > c3_transient.dat
    for ny in 16 32 64; do
        for leg in c3 c2; do
            bin=$BIN; prep=$PREP
            [ "$leg" = c2 ] && { bin=$CBIN; prep=$CPREP; }
            [ -x "$bin" ] || { echo "   (skipping the $leg leg: no binary $bin)"; continue; }
            tag="tr_${leg}_$ny"
            transient_case "$ny" "$ka" "$cap" "$tag" "$bin" "$prep" || { report 1; continue; }
            ns=$($PY -c "print(int(round(0.05/(0.01/$ny**2))))")
            echo "   -- ny = $ny, $leg ($( [ $leg = c3 ] && echo 'fraction-weighted' || echo 'pointwise' ) capacity)"
            run $PY ./check_transient.py check "${tag}_r_$((ns + 1)).h5" \
                --y-wall 0.31 --kappa $ka --cap $cap --heat "$tag.heat.txt" \
                --emit c3_transient.dat --tag "$leg $ny"
            [ $? -eq 0 ] || status=1
        done
    done
    echo "   -- observed order (L2, Linf, interface flux, DECAY RATE):"
    $PY - <<'EOF'
import math
rows = {}
for line in open("c3_transient.dat"):
    leg, ny, l2, linf, q, mu = line.split()
    rows.setdefault(leg, []).append((int(ny), float(l2), float(linf),
                                     float(q), float(mu)))
for leg in sorted(rows):
    v = sorted(rows[leg])
    print(f"   {leg}:")
    for i, (ny, l2, linf, q, mu) in enumerate(v):
        o = ""
        if i:
            p = v[i - 1]
            o = "   order %5.2f %5.2f %5.2f %5.2f" % tuple(
                math.log(p[j] / c) / math.log(2.0) if c > 0 else float("nan")
                for j, c in ((1, l2), (2, linf), (3, q), (4, mu)))
        print("     ny %3d   L2 %.4e  Linf %.4e  q %.4e  mu %.4e%s"
              % (ny, l2, linf, q, mu, o))
EOF
fi

# --- (2) the interface-heat (Nusselt) diagnostic ---------------------------
# C1 left the conjugate flux columns SMOKE-gated. Two statements here, and a
# third in the `budget` group:
#   (a) a CLOSED FORM. A slab whose outer face is insulated and whose solid
#       carries a volumetric source must deliver every watt of it across the
#       interface at steady state: H = C_s S y_w A. No reference run -- and
#       exact rather than approximate, because the fluid fraction of a cell
#       cut by a PLANE is exact, so the fraction-weighted source integrates to
#       the true solid volume.
#   (b) the solver's sum against an INDEPENDENT Python transcription of the
#       same discrete sum (the S4 idiom), on the same snapshot.
if want nusselt; then
    echo "== (2) the interface-heat diagnostic: closed form, and vs Python"
    yw=0.3125; ka=4.0; cap=1.0; src=2.0; ny=16
    dt=$($PY -c "print(repr(0.01/$ny**2))")
    ns=400000; wr=200000
    $PY ./make_slab_stl.py nus.stl --y-top "$yw" > /dev/null
    sed -e "s|@STL@|nus.stl|" -e "s|@CASE@|nus.h5|" -e "s|@NY@|$ny|" \
        -e "s|@KAPPA@|$ka|" -e "s|@CAP@|$cap|" -e "s|@PREFIX@|nus|" \
        -e "s|@NSTEPS@|$ns|" -e "s|@WRITE@|$wr|" -e "s|@DT@|$dt|" \
        -e "s|@HEAT@|$wr|" -e "s|@HEATFILE@|nus.heat.txt|" transient.ini > .nus.full.ini
    # the solid generates; the y_min end is INSULATED so all of it must leave
    # through the interface, and the fluid end is held at 0 to absorb it.
    sed -i -e "s|^solid_init.*|solid_init = 0.0\nsolid_source = $src|" \
        -e "s|^y_min_type.*|y_min_type = neumann|" \
        -e "s|^y_min_value.*|y_min_value = 0.0|" .nus.full.ini
    sed '/^coeff_file/d' .nus.full.ini > .nus.prep.ini
    sed '/^stl_file/d'   .nus.full.ini > .nus.ini
    rm -f nus.heat.txt
    mpirun -n "$RANKS" "$PREP" .nus.prep.ini nus.h5 > nus.prep.log 2>&1 || {
        tail -5 nus.prep.log; report 1; }
    mpirun -n "$RANKS" "$BIN" .nus.ini > nus.log 2>&1 || { tail -20 nus.log; report 1; }
    run $PY ./check_nusselt.py slab --heat nus.heat.txt --y-wall "$yw" \
        --kappa "$ka" --cap "$cap" --source "$src" --area 0.0625 --tolerance 1e-11
    [ $? -eq 0 ] || status=1
    run $PY ./check_nusselt.py sum "nus_$ns.h5" --case nus.h5 --kappa "$ka" \
        --heat nus.heat.txt
    [ $? -eq 0 ] || status=1
fi

# --- (2b) the control-volume cross-check, with the flow on -----------------
# A2 validated the penalization force against a Gauss/CV border flux; the same
# move here, but stated as an instantaneous balance rather than a steady one.
# Over the cells the solver classifies FLUID the discrete flux form
# telescopes, so d/dt sum C theta dV is EXACTLY the interface heat minus what
# leaves through a chosen plane. The case is C1's insulated wavy box with the
# flow on, so the geometry is oblique, the analytic dwall path is the one
# exercised, and no steady state has to be paid for.
#
# The plane flux reads the VELOCITY, so it is only as good as the projection:
# the A2 landmine. Both niter are run, and the residual must not depend on it.
if want budget; then
    echo "== (2b) fluid-side CV budget: dE/dt == interface heat - plane flux"
    # niter 40 vs 200: the plane flux reads the velocity, so a residual that
    # moved with the projection's convergence would be the A2 landmine.
    # dt vs dt/2: the residual that is left must fall by 4, or it is not the
    # centred difference's truncation but the diagnostic's error.
    : > c3_budget.dat
    for tag in 40 200 40h; do
        ni=$tag; dtl=2.0e-4; nst=202; wanted=200
        [ "$tag" = 40h ] && { dtl=1.0e-4; nst=404; wanted=400; ni=40; }
        sed -e "s|^nsteps.*|nsteps = $nst|" -e 's|^field_interval.*|field_interval = 1|' \
            -e "s|^dt = .*|dt = $dtl|" -e "s|^dtmax.*|dtmax = $dtl|" \
            -e "s|^niter.*|niter = $ni|" \
            -e "s|^\\[time\\]|[scalar]\\nheat_interval = 1\\nheat_file = wbud_$tag.heat.txt\\n\\n[time]|" \
            -e "s|field_prefix = wavy|field_prefix = wbud_$tag|" wavy.ini > ".wbud_$tag.ini"
        rm -f "wbud_$tag.heat.txt" wbud_${tag}_*.h5
        mpirun -n "$RANKS" "$BIN" ".wbud_$tag.ini" > "wbud_$tag.log" 2>&1 || {
            tail -20 "wbud_$tag.log"; report 1; continue; }
        echo "   -- niter = $ni, dt = $dtl, whole fluid (the interface is the only border)"
        run $PY ./check_nusselt.py budget "wbud_${tag}_$((wanted)).h5" \
            "wbud_${tag}_$((wanted + 1)).h5" "wbud_${tag}_$((wanted + 2)).h5" \
            --heat "wbud_$tag.heat.txt" --cap 2.0 --wavy --tolerance 1e-5 \
            --emit c3_budget.dat --tag "whole $tag"
        [ $? -eq 0 ] || status=1
        echo "   -- niter = $ni, dt = $dtl, a CV bounded by a plane in the fluid (A2 style)"
        run $PY ./check_nusselt.py budget "wbud_${tag}_$((wanted)).h5" \
            "wbud_${tag}_$((wanted + 1)).h5" "wbud_${tag}_$((wanted + 2)).h5" \
            --heat "wbud_$tag.heat.txt" --cap 2.0 --wavy --plane 0.125 \
            --tolerance 1e-5 --emit c3_budget.dat --tag "plane $tag"
        [ $? -eq 0 ] || status=1
    done
    echo "   -- the residual against dt (it must fall like dt^2):"
    sed 's/^/     /' c3_budget.dat
fi

# --- (3) the conducting channel wall ---------------------------------------
# The scheme on the FULL production stack: WALE LES + file-based IBM + (in the
# `refine` leg) 2:1 block refinement, on validation/channel_interface/les_ibm's
# geometry -- flat walls sitting MID-CELL, so the cut fractions are
# non-trivial while the interface stays grid-aligned, i.e. a case where the
# cut-face coefficient is exact and anything the gate sees is the capacity,
# the diagnostic or the stack.
#
# Three statements: it runs and is DETERMINISTIC (1 == 4 ranks, CPU == GPU, at
# tolerance 0, on the nofma pair); the fluid-cell energy budget identifies the
# interface heat with nu_t active; and on the single-level leg the diagnostic
# equals an independent Python transcription.
#
# The coefficient file is prepared here rather than reused: the committed
# ibm_coeff.h5 predates S3 and carries neither coef_p_blocks nor dwall_blocks.
LES=../channel_interface/les_ibm
chan_case() {   # <tag> <refine-lines> <case.h5> <restart|->
    local tag=$1 ref=$2 case=$3 rst=$4
    sed -e "s|@CASE@|$case|" -e "s|@PREFIX@|$tag|" -e "s|@NSTEPS@|1|" \
        -e "s|@WRITE@|1|" -e "s|@HEAT@|-1|" -e "s|@HEATFILE@|$tag.heat.txt|" \
        -e "s|@REFINE@|$ref|" -e "s|@RESTART@|$rst|" chan_conj.ini > ".$tag.full.ini"
    sed -e '/^coeff_file/d' -e '/^\[restart\]/,$d' \
        -e "s|^\\[ibm\\]|[ibm]\\nstl_file = $LES/wall_lo.stl\\nstl_file = $LES/wall_hi.stl|" \
        ".$tag.full.ini" > ".$tag.prep.ini"
    if [ "$rst" = "-" ]; then
        # COLD START. The refined leg has to: refine_body + keep_buried keeps
        # 3328 leaves where the committed IC_refine.h5 (prepared without
        # keep_buried, which a conjugate run may not do) has 2560, so their
        # block tables are different files' worth of layout. The budget
        # identity does not care what the flow is.
        sed -e '/^stl_file/d' -e '/^\[restart\]/,$d' ".$tag.full.ini" > ".$tag.ini"
    else
        sed '/^stl_file/d' ".$tag.full.ini" > ".$tag.ini"
    fi
    [ -f "$case" ] || mpirun -n 4 "$NPREP" ".$tag.prep.ini" "$case" > "$tag.prep.log" 2>&1 || {
        tail -8 "$tag.prep.log"; return 1; }
    return 0
}

if want channel; then
    echo "== (3) conducting channel wall (LES + file IBM + 2:1), on les_ibm's geometry"
    for leg in flat refine; do
        ref=""; case=chan_conj.h5; rst=$LES/IC.h5
        [ "$leg" = refine ] && { ref='refine_body = true\nrefine_levels = 1'
                                 case=chan_conj_ref.h5; rst=- ; }
        chan_case "cc_$leg" "$ref" "$case" "$rst" || { report 1; continue; }
        # TWO STAGES. The scalar starts as a step across the interface, so the
        # first steps carry a violent transient whose third time derivative --
        # which is exactly what the centred dE/dt truncates against -- is
        # enormous: measured, the budget residual is 15 % at step 3 and 1e-7 by
        # step 200. So relax first (one snapshot, 18 MB each), then take two
        # more steps writing every one; the restart file is the third
        # snapshot the centred difference needs.
        sed -e 's|^nsteps.*|nsteps = 200|' -e 's|^field_interval.*|field_interval = 200|' \
            -e 's|^dtmax.*|dtmax = 1.0e-4|' ".cc_$leg.ini" > ".cc_${leg}_a.ini"
        # The restart line is REPLACED, not prepended: the template already
        # carries a [restart] section (empty for the cold-start leg), and two
        # of them means the last one wins -- which silently restarted stage B
        # from the ORIGINAL IC and reset the step counter.
        sed -e 's|^nsteps.*|nsteps = 2|' -e 's|^field_interval.*|field_interval = 1|' \
            -e 's|^heat_interval.*|heat_interval = 1|' \
            -e '/^\[restart\]/,$d' ".cc_${leg}_a.ini" > ".cc_${leg}_b.ini"
        printf '[restart]\nfile = cc_%s_200.h5\n' "$leg" >> ".cc_${leg}_b.ini"
        # The budget legs run on the GPU: they are a physics statement, not a
        # bit-exactness one, and the refined case is 3328 leaves -- 50 minutes
        # of one CPU rank against under a minute of device time. The
        # determinism trio below is the one that needs the nofma pair.
        rm -f "cc_$leg.heat.txt" cc_${leg}_*.h5
        run mpirun -n 1 "$GBIN" ".cc_${leg}_a.ini" > "cc_${leg}_a.log" 2>&1
        if [ $? -ne 0 ]; then tail -15 "cc_${leg}_a.log"; report 1; continue; fi
        grep -E "conjugate interface:" "cc_${leg}_a.log" | sed 's/^/  /'
        run mpirun -n 1 "$GBIN" ".cc_${leg}_b.ini" > "cc_${leg}_b.log" 2>&1
        if [ $? -ne 0 ]; then tail -15 "cc_${leg}_b.log"; report 1; continue; fi
        echo "   -- fluid-cell energy budget ($leg), nu_t active"
        run $PY ./check_nusselt.py budget "cc_${leg}_200.h5" "cc_${leg}_201.h5" \
            "cc_${leg}_202.h5" --heat "cc_$leg.heat.txt" --cap 2.0 --case "$case" \
            --tolerance 1e-4
        [ $? -eq 0 ] || status=1
        if [ "$leg" = flat ]; then
            echo "   -- the diagnostic vs an independent Python transcription"
            run $PY ./check_nusselt.py sum "cc_${leg}_202.h5" --case "$case" --kappa 5.0 \
                --re 180.0 --pr 0.71 --heat "cc_$leg.heat.txt"
            [ $? -eq 0 ] || status=1
        fi
        # determinism on the full stack, at TOLERANCE 0
        sed -e 's|^nsteps.*|nsteps = 20|' -e 's|^field_interval.*|field_interval = 20|' \
            -e 's|^heat_interval.*|heat_interval = 20|' \
            -e "s|field_prefix = cc_$leg|field_prefix = cd_${leg}_r1|" ".cc_$leg.ini" > ".cd_${leg}_r1.ini"
        sed -e "s|field_prefix = cd_${leg}_r1|field_prefix = cd_${leg}_r4|" \
            -e 's|^dims = 1 1 1|dims = 0 0 0|' ".cd_${leg}_r1.ini" > ".cd_${leg}_r4.ini"
        sed "s|field_prefix = cd_${leg}_r1|field_prefix = cd_${leg}_gpu|" \
            ".cd_${leg}_r1.ini" > ".cd_${leg}_gpu.ini"
        mpirun -n 1 "$NBIN" ".cd_${leg}_r1.ini"  > "cd_${leg}_r1.log"  2>&1
        mpirun -n 4 "$NBIN" ".cd_${leg}_r4.ini"  > "cd_${leg}_r4.log"  2>&1
        mpirun -n 1 "$GBIN" ".cd_${leg}_gpu.ini" > "cd_${leg}_gpu.log" 2>&1
        echo "   -- determinism ($leg): 1 vs 4 ranks"
        run $CMP "cd_${leg}_r1_20.h5" "cd_${leg}_r4_20.h5"; [ $? -eq 0 ] || status=1
        echo "   -- determinism ($leg): CPU vs GPU"
        run $CMP "cd_${leg}_r1_20.h5" "cd_${leg}_gpu_20.h5"; [ $? -eq 0 ] || status=1
    done
fi

echo
if [ $status -eq 0 ]; then echo "ALL C3 GATES PASS"; else echo "SOME C3 GATES FAILED"; fi
exit $status
