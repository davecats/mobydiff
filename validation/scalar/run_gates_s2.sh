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

# --- (k) wall functions + scalars must be a hard config error --------------
if want wferr; then
    echo "== [scalar] + [rans] wall_treatment = wall_function must be rejected"
    if mpirun -n 1 "$CBIN" wferr.ini > wferr.log 2>&1; then
        echo "   FAIL: the solver started"; status=1
    else
        grep -m1 "ERROR STOP" wferr.log | sed 's/^/   /'
        grep -q "wall_function is not implemented" wferr.log \
            && echo "   PASS" || { echo "   FAIL: wrong message"; status=1; }
    fi
fi

# --- (i) steady k-omega SST, resolved walls: theta+ log slope --------------
if want sst; then
    rm -f turbsst_*.h5
    run turbsst.ini 1 turbsst && {
        for s in theta theta_kc; do
            python3 check_scalar_turb.py rans "$(newest turbsst)" --scalar $s \
                --dump "turbsst_$s.dat" | sed 's/^/   /' || status=1
        done
    }
fi

# --- LES legs: relax the scalar, re-target the mean, then measure ---------
les_legs() {   # case  (turbles | turbslab)
    local case=$1
    rm -f "${case}_"*.h5 "RT_${case}.h5"
    sed -e "s/^t_final = .*/t_final = $TRELAX/" -e "s/^field_interval = .*/field_interval = 800/" \
        "$case.ini" > ".${case}_relax.ini"
    run ".${case}_relax.ini" 1 "${case}_relax" || return 1
    python3 retarget_theta.py "$(newest $case)" "RT_${case}.h5" \
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
