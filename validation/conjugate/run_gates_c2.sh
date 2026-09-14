#!/usr/bin/env bash
# C2 gates: THE TANGENTIAL TERM the C1 baseline drops -- its indicator
# e_face, and the correction behind [scalar.N] tangential_correction
# (docs/next_session_conjugate.md Section 10, increment C2).
#
#   ./run_gates_c2.sh [flux|indicator|bvp|cylinder|stzero|residual|dt|c1|all]
#
# Environment: BIN   (default ../../build_cpu/moby_solve)
#              PREP  (the moby_prepare next to BIN)
#              NBIN/NPREP/GBIN (the nofma triple, passed through to the C1
#                               suite's det group)
#              RANKS (default 1)
#
# `flux` and `cylinder` run NO solver at all: the field is the manufactured
# solution and phi comes from the case file, so they measure the SCHEME
# rather than a transient. `indicator` is the cross-check that the Fortran
# computes what the checker does; `bvp` is the one global statement the
# manufactured family allows; `dt` is the time-step penalty; `c1` re-runs
# the entire C1 suite with the correction ON.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

BIN=${BIN:-$ROOT/build_cpu/moby_solve}
PREP=${PREP:-$(dirname "$BIN")/moby_prepare}
RANKS=${RANKS:-1}
PY=${PY:-python3}
sel=${1:-all}
status=0

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
run() { echo "   \$ $*"; "$@"; }
report() { if [ "$1" -eq 0 ]; then echo "   PASS"; else echo "   FAIL"; status=1; fi }

# One oblique-plane case file: <theta> <n> <tag> [extra sed...]
# The domain is h-cubic and 4 cells deep in the periodic z; at 45 degrees it
# is half as wide, so the plane still crosses the whole x range and leaves a
# quarter-domain margin at both ends (a plane through a domain CORNER would
# put cut faces on the stale double-physical corner ghosts).
oblique_case() {
    local th=$1 n=$2 tag=$3; shift 3
    local nx ny lx ly lz x0
    if [ "$th" = 45 ]; then nx=$((n / 2)); lx=0.5; x0=0.25; else nx=$n; lx=1.0; x0=0.5; fi
    ny=$n; ly=1.0
    lz=$($PY -c "print(repr(4.0/$n))")
    # y0 is deliberately NOT on the cell lattice: at 45 degrees with
    # y0 = 0.5 a whole diagonal of cell centres lands EXACTLY on the plane,
    # w collapses to 0 or 1 and the measured flux is off by the full
    # conductivity contrast. That is a degenerate geometry, not a scheme
    # error, and 0.5117 keeps every centre off the surface at every h here.
    $PY ./make_geometry_stl.py plane "$tag.stl" --theta "$th" --x0 "$x0" --y0 0.5117 \
        --span 1.0 --depth 1.0 --z0 -0.25 --z1 "$($PY -c "print(repr(4.0/$n + 0.25))")" \
        > /dev/null || return 1
    # The case-specific substitutions run FIRST, so they claim the
    # placeholders they name and the defaults below fill in the rest.
    cp oblique.ini ".$tag.full.ini"
    for e in "$@"; do sed -i "$e" ".$tag.full.ini"; done
    sed -i -e "s|@STL@|$tag.stl|" -e "s|@CASE@|$tag.h5|" -e "s|@PREFIX@|$tag|" \
        -e "s|@NX@|$nx|" -e "s|@NY@|$ny|" -e "s|@NZ@|4|" \
        -e "s|@LX@|$lx|" -e "s|@LY@|$ly|" -e "s|@LZ@|$lz|" -e "s|@NB@|4|" \
        -e "s|@KAPPA@|1.0|" -e "s|@CAP@|1.0|" -e "s|@NSTEPS@|1|" -e "s|@WRITE@|1|" \
        -e "s|@GX@|0.0|" -e "s|@GY@|0.0|" -e "s|@DT@|1.0e-5|" -e "s|@IND@|0|" \
        -e "s|@TANG@|false|" ".$tag.full.ini"
    sed '/^coeff_file/d' ".$tag.full.ini" > ".$tag.prep.ini"
    sed '/^stl_file/d'   ".$tag.full.ini" > ".$tag.ini"
    mpirun -n "$RANKS" "$PREP" ".$tag.prep.ini" "$tag.h5" > "$tag.prep.log" 2>&1 || {
        tail -5 "$tag.prep.log"; return 1; }
    return 0
}

# --- (1) the oblique plane: the measurement the increment exists for -------
# The exact solution is T = q_n xi/k + A (t.x), so the ratio r = A/q_n is a
# DIAL. For each (theta, kappa_s, r, h) the checker reports, per cut face:
# the indicator e_face; the baseline's relative error in the interface-NORMAL
# flux, against Section 3's closed form; and the corrected one. The case
# file depends only on the GEOMETRY, so one prepare serves every (kappa, r).
if want flux; then
    echo "== (1) oblique plane: interface-flux error vs the tangential ratio, vs h"
    : > c2_flux.dat
    for th in 30 45; do
        for n in 64 128 256; do
            tag="obl${th}_$n"
            oblique_case "$th" "$n" "$tag" || { report 1; continue; }
            x0=0.5; [ "$th" = 45 ] && x0=0.25
            for ka in 10.0 1000.0; do
                for r in 0.0 0.01 0.1 1.0; do
                    run $PY ./check_oblique.py flux "$tag.h5" --theta "$th" --x0 "$x0" --y0 0.5117 \
                        --kappa "$ka" --q-n 1.0 --amp "$r" --emit c2_flux.dat
                    [ $? -eq 0 ] || status=1
                done
                # r = infinity: pure tangential, the worst case there is --
                # and the one the corrected flux must reproduce EXACTLY,
                # since s_t = grad_d T makes the straddling correction inert.
                # The tolerance is the STL distance floor divided by h (the
                # phi column of the table), not a scheme tolerance: the
                # baseline's error here is 0.8, six orders above it.
                run $PY ./check_oblique.py flux "$tag.h5" --theta "$th" --x0 "$x0" --y0 0.5117 \
                    --kappa "$ka" --q-n 0.0 --amp 1.0 --tolerance 1e-6
                [ $? -eq 0 ] || status=1
            done
        done
    done
    echo "   (table: theta kappa ratio h nfaces e_max e_rms qbase_max qbase_rms"
    echo "           qcorr_max qcorr_rms fbase_max fbase_rms phi_err -> c2_flux.dat)"
fi

# --- (1a) the SOLVER's own e_face against the checker's --------------------
# The indicator is a deliverable, so it has to be shown to compute what it
# claims. The r = infinity field is the one whose Neumann data is exactly
# consistent on every face, so the solver's ghosts are the analytic field
# and the two implementations must agree to round-off. e_face = 1 exactly
# there, by hand: s_t = grad_d T when grad T is tangential.
if want indicator; then
    echo "== (1a) the solver's e_face vs the checker's (r = infinity)"
    gx=$($PY -c "import math;print(repr(math.cos(math.radians(30))))")
    gy=$($PY -c "import math;print(repr(math.sin(math.radians(30))))")
    oblique_case 30 64 ind30 "s|@IND@|1|" "s|@GX@|$gx|" "s|@GY@|$gy|" \
        "s|@KAPPA@|10.0|" "s|@DT@|1.0e-8|" || report 1
    mpirun -n 1 "$BIN" .ind30.ini > ind30.log 2>&1
    $PY ./seed_manufactured.py plane ind30_1.h5 ind30_ic.h5 --theta 30 \
        --kappa 10 --q-n 0 --amp 1 > /dev/null
    sed -e 's|^nsteps.*|nsteps = 1|' \
        -e "s|^\\[output\\]|[restart]\\nfile = ind30_ic.h5\\n\\n[output]|" \
        -e 's|field_prefix = ind30|field_prefix = ind30_s|' .ind30.ini > .ind30_s.ini
    mpirun -n 1 "$BIN" .ind30_s.ini > ind30_s.log 2>&1
    grep "conjugate indicator" ind30_s.log | head -1 | sed 's/^/  /'
    run $PY ./check_oblique.py flux ind30.h5 --theta 30 --y0 0.5117 --kappa 10 --q-n 0 --amp 1
    [ $? -eq 0 ] || status=1
fi

# --- (1b) the one global statement: the r = infinity BVP -------------------
# Every finite ratio has a JUMP in grad T, so its exact boundary data differs
# between the two materials on the two domain faces the plane cuts -- and the
# solver has only constant per-face rows. r = infinity has no jump, so all
# six faces carry the exact constant Neumann value and the BVP is well posed
# (up to the pure-Neumann constant, which the checker removes). It is
# therefore the only ratio at which the local flux error can be watched
# turning into a SOLUTION error.
if want bvp; then
    echo "== (1b) r = infinity BVP: field error vs h, correction off and on"
    : > c2_bvp.dat
    gx=$($PY -c "import math;print(repr(math.cos(math.radians(30))))")
    gy=$($PY -c "import math;print(repr(math.sin(math.radians(30))))")
    # ONE grid. What this gate states is a COMPARISON -- the same problem,
    # the same time, the correction off and on -- not an order, which gate 1
    # already carries. Finer grids are priced out rather than uninteresting:
    # the explicit step at kappa_s = 10 is 0.4 h^2/kappa_s while the response
    # has to spread over the fluid diffusion time L^2/alpha_f = 1, so n = 64
    # costs 1e5 steps per leg and n = 128 costs 4e5.
    for n in 32; do
        # The run starts AT the exact solution and relaxes to the DISCRETE
        # steady state; that response spreads over the whole domain on the
        # fluid diffusion time L^2/alpha_f = 1, so t_end = 1 with the
        # explicit step 0.4 h^2/kappa_s. Two writes, so the residual between
        # them makes "converged" a measurement rather than a claim.
        ns=$($PY -c "print(int(round(1.0*10.0/(0.4*(1.0/$n)**2))))")
        hf=$((ns / 2))
        for t in false true; do
            tag="bvp_${n}_$t"
            oblique_case 30 "$n" "$tag" "s|@GX@|$gx|" "s|@GY@|$gy|" \
                "s|@KAPPA@|10.0|" "s|@TANG@|$t|" "s|@NSTEPS@|1|" \
                "s|@WRITE@|1|" "s|@DT@|1.0e-4|" || { report 1; continue; }
            mpirun -n 1 "$BIN" ".$tag.ini" > "$tag.log" 2>&1
            $PY ./seed_manufactured.py plane "${tag}_1.h5" "${tag}_ic.h5" \
                --theta 30 --x0 0.5 --y0 0.5117 --kappa 10 --q-n 0 --amp 1 > /dev/null
            sed -e "s|^nsteps.*|nsteps = $ns|" -e "s|^field_interval.*|field_interval = $hf|" \
                -e "s|^\\[output\\]|[restart]\\nfile = ${tag}_ic.h5\\n\\n[output]|" \
                -e "s|field_prefix = $tag|field_prefix = ${tag}_r|" ".$tag.ini" > ".${tag}_r.ini"
            mpirun -n "$RANKS" "$BIN" ".${tag}_r.ini" > "${tag}_r.log" 2>&1 || {
                tail -12 "${tag}_r.log"; report 1; continue; }
            run $PY ./check_oblique.py field "${tag}_r_$((ns + 1)).h5" --theta 30 --x0 0.5 \
                --y0 0.5117 --kappa 10 --q-n 0 --amp 1 --tag "n=$n tang=$t" --emit c2_bvp.dat
            [ $? -eq 0 ] || status=1
            # Convergence is a MEASUREMENT: the same error at half the
            # time. (Do not read the raw snapshot difference instead -- the
            # pure-Neumann problem has a uniform O(h) drift, because the
            # discrete boundary fluxes cancel only to the accuracy of the
            # staircase solid/fluid split of the boundary faces. The
            # checker removes the mean, so the error is unaffected.)
            echo -n "   the same error at half the time:"
            $PY ./check_oblique.py field "${tag}_r_$hf.h5" --theta 30 --x0 0.5 \
                --y0 0.5117 --kappa 10 --q-n 0 --amp 1 | tail -1
        done
    done
    cat c2_bvp.dat | sed 's/^/   /'
    # ...and WHY: the local truncation error of the two schemes on the exact
    # field. The correction makes the pointwise face flux exact and the cell
    # balance worse -- a finite-difference divergence wants the face AVERAGE.
    echo "   local truncation error on the exact field (r = infinity):"
    run $PY ./check_oblique.py residual obl30_64.h5 --theta 30 --x0 0.5 --y0 0.5117 \
        --kappa 10 --q-n 0 --amp 1
fi

# --- (2) the cylindrical shell: the O(kappa h/a) error in w ----------------
# A radial field has NO tangential gradient at all, so s_t = 0 exactly and
# the correction must be inert. What is left in the cut-face flux is the
# level-set weight's curvature error, which this isolates.
if want cylinder; then
    echo "== (2) cylindrical shell: the curvature error in w (s_t = 0 by symmetry)"
    : > c2_cyl.dat
    : > c2_dip.dat
    for n in 64 128 256; do
        tag="cyl_$n"
        lz=$($PY -c "print(repr(4.0/$n))")
        $PY ./make_geometry_stl.py cylinder "$tag.stl" --centre 0.5 0.5 --radius 0.25 \
            --facets 4096 --z0 -0.25 --z1 "$($PY -c "print(repr(4.0/$n + 0.25))")" > /dev/null
        sed -e "s|@STL@|$tag.stl|" -e "s|@CASE@|$tag.h5|" -e "s|@PREFIX@|$tag|" \
            -e "s|@NX@|$n|" -e "s|@NY@|$n|" -e "s|@NZ@|4|" \
            -e "s|@LX@|1.0|" -e "s|@LY@|1.0|" -e "s|@LZ@|$lz|" -e "s|@NB@|4|" \
            -e "s|@KAPPA@|10.0|" -e "s|@CAP@|1.0|" -e "s|@NSTEPS@|1|" -e "s|@WRITE@|1|" \
            -e "s|@GX@|0.0|" -e "s|@GY@|0.0|" -e "s|@DT@|1.0e-5|" -e "s|@IND@|0|" \
            -e "s|@TANG@|false|" oblique.ini > ".$tag.full.ini"
        sed '/^coeff_file/d' ".$tag.full.ini" > ".$tag.prep.ini"
        sed '/^stl_file/d'   ".$tag.full.ini" > ".$tag.ini"
        mpirun -n "$RANKS" "$PREP" ".$tag.prep.ini" "$tag.h5" > "$tag.prep.log" 2>&1 || {
            tail -5 "$tag.prep.log"; report 1; continue; }
        for ka in 10.0 1000.0; do
            run $PY ./check_cylinder.py flux "$tag.h5" --radius 0.25 --kappa "$ka" \
                --source 1.0 --emit c2_cyl.dat
            [ $? -eq 0 ] || status=1
            # ...and the cylinder in a UNIFORM GRADIENT, which is the same
            # curved interface WITH a tangential gradient: harmonic on both
            # sides, so the exact divergence is zero and the truncation test
            # applies, and its local ratio sweeps 0 -> infinity around the
            # body. This is where the s_t estimate is measured against a
            # known exact value on a curved interface.
            run $PY ./check_cylinder.py dipole "$tag.h5" --radius 0.25 \
                --kappa "$ka" --emit c2_dip.dat
            [ $? -eq 0 ] || status=1
        done
    done
fi

# --- (3b) STAGE 0: does a ONE-SIDED s_t converge? --------------------------
# KILL GATE 1 of docs/next_session_tangential.md. The shipped de-bias models
# the straddle of the arm difference and is derived for a PLANE, so on a
# curved interface it does not converge (43/54/46 %, gate (2) above). A
# same-side stencil never straddles, and s_t is CONTINUOUS across the
# interface so either side estimates the same number. This measures whether
# that converges -- everything downstream depends on it.
#
# The plane is the unit test (piecewise-linear field + exact phi => a
# one-sided estimator must return the geometry floor); the CURVED cases are
# the gate. `multipole --mode 2` exists because the dipole's interior is
# exactly LINEAR, which would flatter any solid-side score.
if want stzero; then
    echo "== (3b) stage 0: one-sided s_t, the convergence gate"
    run $PY ./check_st.py plane --tag obl30 --theta 30 --kappa 10 \
        --q-n 1.0 --amp 1.0 --tolerance 1.0e-8
    [ $? -eq 0 ] || status=1
    for ka in 10.0 1000.0; do
        run $PY ./check_st.py cylinder  --tag cyl --radius 0.25 --kappa "$ka"
        [ $? -eq 0 ] || status=1
        run $PY ./check_st.py multipole --tag cyl --radius 0.25 --kappa "$ka" --mode 2
        [ $? -eq 0 ] || status=1
    done
fi

# --- (3c) STAGE 1b: the residual table, which is THE test ------------------
# Everything upstream measures an s_t error. This measures what the cell
# balance sees -- the cut-cell truncation residual on the exact field, whose
# exact value is zero -- and it is the measurement that decided C2.
#
# Read the table as: C1 baseline | k_area + the SHIPPED s_t (C2's "way out")
# | + the stage-0 one-sided s_t | + stage 1b's split and extension. Planes
# should be EXACT; curved cases should be BOUNDED where the baseline grows
# like 1/h, and never worse at any contrast.
if want residual; then
    echo "== (3c) stage 1b: the cut-cell truncation residual"
    # THREE angles, because 45 degrees is degenerate: there the C1 baseline is
    # accidentally EXACT (1e-09 where 30 degrees reads 5.18), so it
    # discriminates nothing. 20 and 35 are generated here; 30 comes from the
    # flux gate above.
    for th in 20 35; do
        [ -f "obl${th}_256.h5" ] && continue
        oblique_case "$th" 256 "obl${th}_256" || { report 1; continue; }
    done
    for th in 20 30 35; do
        for ka in 10.0 1000.0; do
            for r in 0.01 1.0; do
                run $PY ./check_oblique.py residual "obl${th}_256.h5" --theta "$th" \
                    --x0 0.5 --y0 0.5117 --kappa "$ka" --q-n 1.0 --amp "$r"
            done
            run $PY ./check_oblique.py residual "obl${th}_256.h5" --theta "$th" \
                --x0 0.5 --y0 0.5117 --kappa "$ka" --q-n 0.0 --amp 1.0
        done
    done
    for ka in 0.01 10.0 1000.0; do
        for n in 64 128 256; do
            for md in 2 1; do
                run $PY ./check_cylinder.py dipole "cyl_$n.h5" --radius 0.25 \
                    --kappa "$ka" --mode "$md"
            done
        done
    done
fi

# --- (4) the time-step penalty of the correction ---------------------------
# The correction is an EXPLICIT spatial operator at the cut faces whose size
# grows with the conductivity contrast, so it enters the explicit limit --
# and C1 already showed that limit is subtler than the strategy doc assumed.
# scalar_conjugate_peclet_rate adds its Gershgorin row sum; this measures
# what that costs at kappa_s = 1e3, and that the run is actually stable.
if want dt; then
    echo "== (4) time-step penalty at kappa_s = 1e3"
    gx=$($PY -c "import math;print(repr(math.cos(math.radians(30))))")
    gy=$($PY -c "import math;print(repr(math.sin(math.radians(30))))")
    # BOTH values of pecletmax, because the run at the shipped default is
    # itself a measurement: at kappa_s = 1e3 on an OBLIQUE interface the
    # C1 baseline is marginal there (see the README), so 0.2 is where the
    # correction's own penalty can be read off a stable pair.
    for pm in 0.4 0.2; do
    for t in false true; do
        tag="dtk_${pm}_$t"
        oblique_case 30 64 "$tag" "s|@KAPPA@|1000.0|" "s|@TANG@|$t|" \
            "s|@GX@|$gx|" "s|@GY@|$gy|" "s|@DT@|1.0e-3|" || { report 1; continue; }
        mpirun -n 1 "$BIN" ".$tag.ini" > "$tag.log" 2>&1
        $PY ./seed_manufactured.py plane "${tag}_1.h5" "${tag}_ic.h5" --theta 30 --y0 0.5117 \
            --kappa 1000 --q-n 0 --amp 1 > /dev/null
        sed -e 's|^nsteps.*|nsteps = 2000|' -e 's|^field_interval.*|field_interval = 2000|' \
            -e "s|^\\[output\\]|[restart]\\nfile = ${tag}_ic.h5\\n\\n[output]|" \
            -e "s|^pecletmax.*|pecletmax = $pm|" \
            -e "s|field_prefix = $tag|field_prefix = ${tag}_r|" ".$tag.ini" > ".${tag}_r.ini"
        mpirun -n 1 "$BIN" ".${tag}_r.ini" > "${tag}_r.log" 2>&1
        $PY -c "
import h5py, numpy as np
f = h5py.File('${tag}_r_2001.h5')
th = f['theta'][...]
print('   pecletmax = $pm  tangential_correction = %-5s  dt = %.6e  max|theta| = %.6e  finite = %s'
      % ('$t', f.attrs['dt'], np.abs(th).max(), bool(np.isfinite(th).all())))"
        [ $? -eq 0 ] || status=1
    done
    done
fi

# --- (3) every C1 gate, with the correction ON -----------------------------
# The grid-aligned C1 cases have s_t = 0 BY CONSTRUCTION, so they must
# reproduce their recorded numbers to round-off; the wavy one is oblique, so
# it is a genuine re-measurement of conservation, the guards, the 2:1
# precondition and determinism with the correction active.
if want c1; then
    echo "== (3) the whole C1 suite with tangential_correction = true"
    TANG=true ./run_gates_c1.sh all
    [ $? -eq 0 ] || status=1
fi

echo
if [ $status -eq 0 ]; then echo "ALL C2 GATES PASS"; else echo "SOME C2 GATES FAILED"; fi
exit $status
