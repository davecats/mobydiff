#!/usr/bin/env bash
# Passive-scalar S2 gates (docs/next_session_scalar.md, increment S2 --
# the turbulent closure D_face = 1/(Re Pr) + nut/Pr_t).
#
#   ./run_gates_s2.sh [group]   groups: kays wferr sst les band det
#
# Environment: BIN   solver binary   (default ../../build_gpu/moby_solve)
#              CBIN  CPU binary      (default ../../build_cpu/moby_solve)
#              UNIT  unit-test binary(default ../../build_cpu/scalar_test)
#              TRELAX end time of the LES relaxation leg (default 30.0)
#              TEND   end time of the LES statistics leg (default 34.0)
#              SNAP   snapshot interval of the statistics leg, in steps (400)
#              SKIP   snapshots of the statistics leg dropped (default 2)
#
# The LES legs start from the ../channel_interface/les campaign's DEVELOPED
# velocity field at t = 24.8 with theta seeded by make_theta_ic.py, so only
# the SCALAR spins up. Leg 1 relaxes it; retarget_theta.py then replaces the
# MEAN profile by the one the run's own measured total diffusivity implies
# (keeping every fluctuation), which removes the remaining slow-mode wait;
# leg 2 measures. Whether that is really stationary is checked, not assumed:
# the analysis reports the constant-flux residual and a first-half /
# second-half comparison of theta_tau.
#
# The `det` group's CPU == GPU / 1 == 4 rank comparisons MUST use the nofma
# builds (compile_nofma.sh): default FMA contraction differs between host and
# device and between vectorisation widths.
set -uo pipefail
cd "$(dirname "$0")"

BIN=${BIN:-../../build_gpu/moby_solve}
CBIN=${CBIN:-../../build_cpu/moby_solve}
UNIT=${UNIT:-../../build_cpu/scalar_test}
NBIN=${NBIN:-../../build_cpu_nofma/moby_solve}
NGBIN=${NGBIN:-../../build_gpu_nofma/moby_solve}
TRELAX=${TRELAX:-30.0}
TEND=${TEND:-34.0}
SNAP=${SNAP:-400}
SKIP=${SKIP:-2}
sel=${1:-all}
status=0

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
run() {  # ini ranks logname [binary]
    local ini=$1 ranks=$2 log=$3 bin=${4:-$BIN}
    echo "== $log (ranks $ranks)"
    if ! mpirun -n "$ranks" "$bin" "$ini" > "$log.log" 2>&1; then
        echo "   RUN FAILED -- see $log.log"; status=1; return 1
    fi
    grep -E "seconds_per_step" "$log.log" | head -1 | sed 's/^/   /'
}
newest() { ls -t "$1"_*.h5 2>/dev/null | head -1; }

# --- (m) the Kays-Crawford correlation unit test ---------------------------
if want kays; then
    echo "== Kays-Crawford correlation (src/test_scalar.f90)"
    mpirun -n 1 "$UNIT" | sed 's/^/   /' || status=1
fi

# --- (k) wall functions + scalars: REJECTED in S2, ACCEPTED since S5a ------
# The S2 gate was that the combination error-stopped ("a thermal wall
# function is increment S5"). S5a implemented it, so the same ini is now the
# positive control: the solver must START and report the thermal wall
# function's per-scalar constants. The physics gates live in
# run_gates_s5.sh.
if want wferr; then
    echo "== [scalar] + [rans] wall_treatment = wall_function must be ACCEPTED (S5a)"
    sed -e "s/^nsteps.*/nsteps = 5/" -e "s/^t_final.*/t_final = 0.0/" \
        wferr.ini > .wferr_short.ini
    if mpirun -n 1 "$CBIN" .wferr_short.ini > wferr.log 2>&1; then
        grep -m1 "thermal wall function" wferr.log | sed 's/^/   /'
        grep -q "thermal wall function" wferr.log \
            && echo "   PASS" || { echo "   FAIL: no thermal wall function reported"; status=1; }
    else
        echo "   FAIL: the solver stopped -- see wferr.log"; status=1
    fi
fi

# --- (i) steady k-omega SST, resolved walls: theta+ log slope --------------
if want sst; then
    rm -f turbsst_*.h5
    run turbsst.ini 1 turbsst && {
        # theta_kc runs prt_model = kays IN THE SOLVER, so the checker must
        # be told: without --prt-model its nut-integral prediction is built
        # with a constant Pr_t and the comparison is against the wrong model
        # (it read 0.68 % / 1.9e-01 instead of 0.094 % / 2.6e-08, and
        # reported "Pr_t 0.8500 .. 0.8500", which is its own assumption and
        # not a measurement). The S2 README numbers were produced with the
        # flag; this makes the group reproduce them by itself.
        for s in theta theta_kc; do
            model=constant; [ "$s" = theta_kc ] && model=kays
            python3 check_scalar_turb.py rans "$(newest turbsst)" --scalar $s \
                --prt-model $model --dump "turbsst_$s.dat" | sed 's/^/   /' || status=1
        done
    }
fi

# The last snapshot of a leg is the POST-LOOP write_field, and a run that
# stops on t_final takes ONE EXTRA step whose dt is the accumulated round-off
# in t_current -- here 10401 steps instead of 10400, t_current =
# 29.999999999975433 against a 1e-12 stopping tolerance, so dt = 2.46e-11.
# The projection's pressure on that step is amplified by 1/dt: measured
# 2026-08-07, |pn| = 1.5e6 at step 60001 against 9.1 at 60000, with the
# VELOCITY identical to 8.5e-6. Restarting the statistics leg from it NaNs
# the run in 54 steps (dt -> 0, cfl NaN). Retarget from the last PERIODIC
# snapshot instead -- physically the same field, one step earlier.
# (docs/next_session_verification.md B1; NOT the niter = 6 pn-drift mode,
# which is what this was blamed on before the mechanism was measured.)
last_periodic() {  # prefix interval
    ls "$1"_*.h5 2>/dev/null | sed -E 's/.*_([0-9]+)\.h5$/\1 &/' | sort -n -k1,1 \
        | awk -v iv="$2" '$1 % iv == 0 { f = $2 } END { print f }'
}

# --- LES legs: relax the scalar, re-target the mean, then measure ---------
les_legs() {   # case  (turbles | turbslab)
    local case=$1
    rm -f "${case}_"*.h5 "RT_${case}.h5"
    sed -e "s/^t_final = .*/t_final = $TRELAX/" -e "s/^field_interval = .*/field_interval = 800/" \
        "$case.ini" > ".${case}_relax.ini"
    run ".${case}_relax.ini" 1 "${case}_relax" || return 1
    python3 retarget_theta.py "$(last_periodic $case 800)" "RT_${case}.h5" \
        --stats "${case}_*.h5" --skip 2 | sed 's/^/   /'
    mkdir -p relax && mv "${case}_"*.h5 relax/
    sed -e "s/^t_final = .*/t_final = $TEND/" -e "s/^field_interval = .*/field_interval = $SNAP/" \
        -e "s|^file = IC_${case}.h5|file = RT_${case}.h5|" "$case.ini" > ".${case}_stats.ini"
    run ".${case}_stats.ini" 1 "$case" || return 1
    rm -f ".${case}_relax.ini" ".${case}_stats.ini"
}

# --- (h) developed turbulent channel under LES: theta+ vs Kader -----------
if want les; then
    les_legs turbles && \
        python3 check_scalar_turb.py channel 'turbles_*.h5' --skip "$SKIP" \
            --dump turbles_profile.dat | sed 's/^/   /' || status=1
fi

# --- (j) 2:1 wall-band refined channel: NO spurious scalar band ------------
if want band; then
    les_legs turbslab && \
        python3 check_scalar_turb.py band 'turbslab_*.h5' --reference 'turbles_*.h5' \
            --skip "$SKIP" --dump turbslab_band.dat | sed 's/^/   /' || status=1
fi

# --- (l) determinism on a scalar + turbulence case ------------------------
if want det; then
    rm -f dles_*.h5
    for tag in r1 r4 gpu; do
        sed -e "s/^field_prefix = .*/field_prefix = dles_$tag/" detles.ini > ".dles_$tag.ini"
        rm -f "dles_${tag}_"*.h5
        case $tag in
            r1)  run ".dles_$tag.ini" 1 "dles_$tag" "$NBIN" ;;
            r4)  run ".dles_$tag.ini" 4 "dles_$tag" "$NBIN" ;;
            gpu) run ".dles_$tag.ini" 1 "dles_$tag" "$NGBIN" ;;
        esac
        rm -f ".dles_$tag.ini"
    done
    echo "   1 rank vs 4 ranks:"
    python3 ../../tools/compare_fields.py "$(newest dles_r1)" "$(newest dles_r4)" \
        un vn wn pn nut theta --tolerance 0 | sed 's/^/     /' || status=1
    echo "   CPU vs GPU:"
    python3 ../../tools/compare_fields.py "$(newest dles_r1)" "$(newest dles_gpu)" \
        un vn wn pn nut theta --tolerance 0 | sed 's/^/     /' || status=1
fi

echo
[ $status -eq 0 ] && echo "scalar S2 gates: no failures" || echo "scalar S2 gates: FAILURES"
exit $status
