#!/usr/bin/env bash
# The conjugate turbulent-channel validation (Re_tau = 180, Pr = 0.71).
#
#   ./run_cht.sh [prepare|ic|develop|reseed|stats|all]
#
# Environment: BIN (default ../../../build_gpu/moby_solve), PREP, RANKS.
#
# THREE PHASES, because two of them are transients nobody measures:
#   prepare  moby_prepare on the two wall STLs -> cht180_nb32.h5 (single
#            level: 240 leaves, 0 refinement levels -- there are no 2:1
#            interfaces in this case, deliberately);
#   ic       map the developed KMM180 DNS field into the fluid gap and seed
#            each scalar with its series-resistance profile;
#   develop  t = 0 .. 10, statistics OFF: the mapped velocity regenerates the
#            small scales the re-interpolation lost, and the scalars find the
#            near-wall profile the analogy IC gets deliberately wrong;
#   stats    t = 10 .. 31, statistics ON: the measurement.
#
# nb = 32 rather than the default 8 is a pure decomposition choice -- measured
# 4.3x cheaper projection exchange at identical physics -- and it is REQUIRED
# to be set at all, because a prepared case file's block table must not depend
# on the rank count (moby_prepare errors out without it).
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../../.. && pwd)

BIN=${BIN:-$ROOT/build_gpu/moby_solve}
PREP=${PREP:-$ROOT/build_cpu/moby_prepare}
RANKS=${RANKS:-1}
PY=${PY:-python3}
sel=${1:-all}
CASE=cht180_nb32.h5

# t per step is ~1.5e-3 (dt = 1.44e-3, set by the cut-cell Peclet limit).
DEV_STEPS=${DEV_STEPS:-6800}      # t -> ~10.2
STAT_STEPS=${STAT_STEPS:-14000}   # t -> ~31.3
SAMPLE=${SAMPLE:-25}
FLUSH=${FLUSH:-2000}

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }

ini() {  # <prefix> <nsteps> <write> <restart> <sample> <flush> <statsfile>
    sed -e "s|@CASE@|$CASE|" -e "s|@PREFIX@|$1|" -e "s|@NSTEPS@|$2|" \
        -e "s|@WRITE@|$3|" -e "s|@RESTART@|$4|" -e "s|@SAMPLE@|$5|" \
        -e "s|@STATS@|$6|" -e "s|@STATSFILE@|$7|" \
        cht180.ini > ".$1.ini"
}

if want prepare; then
    echo "== prepare (single level, no refinement)"
    $PY ../make_slab_stl.py wall_lo.stl --y-bottom -0.3 --y-top 1.0 \
        --x-range 0 12.566370614359172 --z-range 0 6.283185307179586 --pad 0.3
    $PY ../make_slab_stl.py wall_hi.stl --y-bottom 3.0 --y-top 4.3 \
        --x-range 0 12.566370614359172 --z-range 0 6.283185307179586 --pad 0.3
    ini prep 1 1 - 0 0 x
    sed -e '/^coeff_file/d' -e '/^\[restart\]/,$d' \
        -e 's|^\[ibm\]|[ibm]\nstl_file = wall_lo.stl\nstl_file = wall_hi.stl|' \
        .prep.ini > .prep_in.ini
    mpirun -n 8 "$PREP" .prep_in.ini "$CASE" 2>&1 | tail -3
fi

if want ic; then
    echo "== initial condition (KMM180 -> the gap, + the scalar profiles)"
    ini mint 1 1 - 0 0 x
    sed '/^\[restart\]/,$d' .mint.ini > .mint_cs.ini
    mpirun -n 1 "$BIN" .mint_cs.ini > mint.log 2>&1 || { tail -5 mint.log; exit 1; }
    mv -f mint_1.h5 cs_1.h5 2>/dev/null || true
    $PY ./make_cht_ic.py cs_1.h5 IC_cht.h5
fi

if want develop; then
    echo "== develop: t = 0 .. ~10 (statistics off)"
    ini dev "$DEV_STEPS" "$DEV_STEPS" IC_cht.h5 0 0 x
    time mpirun -n "$RANKS" "$BIN" .dev.ini > dev.log 2>&1 || { tail -15 dev.log; exit 1; }
    grep -E "seconds_per_step" dev.log | tail -1
fi

if want reseed; then
    # The solid mean is SET, not waited for: with a wall one half-height thick
    # it relaxes on d^2/alpha_s, which is t ~ 8e5 for the C_s = 1e4 scalar --
    # and the capacity cannot change the steady state at all (C1 gate 1c). The
    # fluid side has equilibrated by now, so the run knows its own fluid
    # resistance and the exact steady solid profile follows in closed form.
    # Fluctuations are untouched; they equilibrate on the penetration time.
    echo "== re-seed the solid mean from the developed fluid"
    $PY ./reseed_solid.py "dev_${DEV_STEPS}.h5" IC_stat.h5
fi

if want stats; then
    echo "== statistics: t = ~10 .. ~31"
    rm -f cht_stats.h5 cht_vel_stats.h5
    ini stat "$STAT_STEPS" "$FLUSH" IC_stat.h5 "$SAMPLE" "$FLUSH" cht_stats.h5
    time mpirun -n "$RANKS" "$BIN" .stat.ini > stat.log 2>&1 || { tail -15 stat.log; exit 1; }
    grep -E "seconds_per_step" stat.log | tail -1
fi
