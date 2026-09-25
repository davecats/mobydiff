#!/usr/bin/env bash
# Passive-scalar S5a gates (docs/next_session_scalar.md, increment S5a --
# the THERMAL WALL FUNCTION: Kader/Jayatilleke under [rans] wall_treatment =
# wall_function, which S2 rejected as a hard config error).
#
#   ./run_gates_s5.sh [group]   groups: unit ref sweep det
#
# Environment: BIN   solver binary   (default ../../build_cpu/moby_solve)
#              GBIN  GPU binary      (default ../../build_gpu/moby_solve)
#              UNIT  unit-test binary(default ../../build_cpu/scalar_test)
#              NBIN/NGBIN  the nofma CPU/GPU binaries (the `det` group)
#              RANKS ranks for the channel legs (default 4)
#
# Shape of the validation, after T3 (validation/rans_sst): every case is the
# T3 wall-function channel with the turbsst scalar pair bolted on, and the
# reference is the RESOLVED scalar channel (turbsst.ini) restarted with the
# S4 statistics on -- so theta_tau and the wall flux on both sides come out
# of the same instrument, the solver's own statistics, and not out of a
# hand-rolled post-processing of the profile.
#
# Every run is step-bound (nsteps), so the sample steps are known in advance:
# the statistics are switched on with the interval EQUAL to the run length,
# which takes exactly ONE sample, at the last step, of the very field the
# snapshot carries. That matters twice -- the accumulator is cumulative, so
# sampling from step 0 would average in the transient, and the wall gate
# below is an IDENTITY between the statistics' flux and the snapshot's own
# k/theta, which needs both to be the same instant. `stats_on` edits the
# [scalar] section only: [case.channel] carries keys of the same name for
# the channel statistics, which must stay off.
set -uo pipefail
cd "$(dirname "$0")"

BIN=${BIN:-../../build_cpu/moby_solve}
GBIN=${GBIN:-../../build_gpu/moby_solve}
UNIT=${UNIT:-../../build_cpu/scalar_test}
NBIN=${NBIN:-../../build_cpu_nofma/moby_solve}
NGBIN=${NGBIN:-../../build_gpu_nofma/moby_solve}
RANKS=${RANKS:-4}
sel=${1:-all}
status=0

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
newest() { ls -t "$1"_*.h5 2>/dev/null | head -1; }

# Switch the [scalar] statistics on in a copy of an ini: sample/write every
# $2 steps (the section-range address keeps [case.channel] untouched).
stats_on() {  # in out interval
    sed -e "/^\[scalar\]/,/^\[scalar\./ s/^stats_sample_interval.*/stats_sample_interval = $3/" \
        -e "/^\[scalar\]/,/^\[scalar\./ s/^stats_write_interval.*/stats_write_interval = $3/" \
        "$1" > "$2"
}

run() {  # ini ranks logname [binary]
    local ini=$1 ranks=$2 log=$3 bin=${4:-$BIN}
    echo "== $log (ranks $ranks)"
    if ! mpirun -n "$ranks" "$bin" "$ini" > "$log.log" 2>&1; then
        echo "   RUN FAILED -- see $log.log"; status=1; return 1
    fi
    grep -E "seconds_per_step" "$log.log" | head -1 | sed 's/^/   /'
}

# --- (m2) the correlations, host-side against mpmath ----------------------
if want unit; then
    echo "== thermal wall function correlations (src/test_scalar.f90)"
    mpirun -n 1 "$UNIT" | sed 's/^/   /' || status=1
fi

# --- (w0) the RESOLVED reference: turbsst + the S4 statistics -------------
if want ref; then
    rm -f wfsst_*.h5 wfsst_stats*.h5
    stats_on wfsst.ini .wfsst_st.ini 200
    run .wfsst_st.ini "$RANKS" wfsst && {
        python3 ../../tools/scalar_stats.py profile wfsst_stats.h5 \
            --pr 0.71 --walls 1.0 -1.0 | sed 's/^/   /' || status=1
    }
fi

# --- (w)/(x) the wall-function sweep --------------------------------------
if want sweep; then
    for tag in y05 y15 y30 y45; do
        rm -f wfs180_${tag}_*.h5 wfs180_${tag}_stats*.h5
        stats_on wfs180_$tag.ini .wfs180_${tag}_st.ini 20000
        run .wfs180_${tag}_st.ini "$RANKS" wfs180_$tag || continue
        python3 check_scalar_wf.py wall wfs180_${tag}_stats.h5 \
            "$(newest wfs180_$tag)" --scalar theta --kader-tol 0.25 \
            | sed 's/^/   /' || status=1
    done
    echo "== the sweep against the RESOLVED reference"
    python3 check_scalar_wf.py compare wfsst_stats.h5 \
        wfs180_y05_stats.h5 wfs180_y15_stats.h5 wfs180_y30_stats.h5 \
        wfs180_y45_stats.h5 --tolerance 0.12 | sed 's/^/   /' || status=1
    # theta_kc (index 2) carries a LOOSER tolerance BY DESIGN, and the last
    # column says why: the thermal wall function is defined with a CONSTANT
    # Pr_t, so a wall-function run's wall cells treat theta and theta_kc
    # identically and their theta_tau ratio is ~1, while the resolved
    # reference -- whose wall cells are exactly where Kays-Crawford damps
    # the eddy diffusivity hardest (Pe_t -> 0) -- separates the pair by 7%.
    # The wall-function theta_kc is therefore EXPECTED to sit further above
    # its resolved counterpart than theta does.
    echo "== the Kays-Crawford scalar of the same runs (theta_kc, index 2)"
    python3 check_scalar_wf.py compare wfsst_stats.h5 \
        wfs180_y05_stats.h5 wfs180_y15_stats.h5 wfs180_y30_stats.h5 \
        wfs180_y45_stats.h5 --index 2 --tolerance 0.16 --pair \
        | sed 's/^/   /' || status=1
fi

# --- (z) determinism of the wall-function scalar path ---------------------
# nofma builds: the wall function's log() differs by an ulp between host and
# device libm (the T3 finding), so CPU vs GPU is a tolerance, 1 == 4 ranks an
# equality.
#
# The rank comparison PINS [mpi] dims to a z split. This channel is NOT
# rank-independent under an x split -- and that is PRE-EXISTING, nothing to
# do with the scalars: wf180_y30.ini WITHOUT any [scalar] section, run with
# the S4 reference binary, reproduces the same deviation to the last bit
# (un 8.800384e-03 after 20 steps, already 7.6e-04 after ONE step, so it is
# the initial state or the first substage on an x-split rank box, not an
# accumulation). Splitting z is exact there, so pinning it makes this gate
# a statement about the scalar wall-function path rather than a re-run of
# someone else's bug.
if want det; then
    sed -e "s/^nsteps.*/nsteps = 20/" -e "s/^field_interval.*/field_interval = 20/" \
        wfs180_y30.ini > .wfdet.ini
    for r in 1 4; do
        sed -e "s/^field_prefix.*/field_prefix = wfdet_r$r/" .wfdet.ini > .wfdet_r$r.ini
        printf '\n[mpi]\ndims = 1 1 %s\n' "$r" >> .wfdet_r$r.ini
        rm -f wfdet_r${r}_*.h5
        run .wfdet_r$r.ini "$r" wfdet_r$r "$NBIN" || true
    done
    sed -e "s/^field_prefix.*/field_prefix = wfdet_gpu/" .wfdet.ini > .wfdet_gpu.ini
    rm -f wfdet_gpu_*.h5
    run .wfdet_gpu.ini 1 wfdet_gpu "$NGBIN" || true
    ds="un vn wn pn nut k omega theta theta_kc"
    echo "== 1 rank vs 4 ranks (tolerance 0)"
    python3 ../../tools/compare_fields.py wfdet_r1_20.h5 wfdet_r4_20.h5 $ds \
        --tolerance 0 | sed 's/^/   /' || status=1
    echo "== CPU vs GPU (the wall-function log() ulp class, tolerance 1e-11)"
    python3 ../../tools/compare_fields.py wfdet_r1_20.h5 wfdet_gpu_20.h5 $ds \
        --tolerance 1e-11 | sed 's/^/   /' || status=1
fi

echo
[ $status -eq 0 ] && echo "S5a gates: ALL PASS" || echo "S5a gates: FAILURES"
exit $status
