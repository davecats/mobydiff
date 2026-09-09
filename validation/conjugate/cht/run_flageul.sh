#!/usr/bin/env bash
# THE FLAGEUL-MATCHED conjugate channel: their grid resolution and their
# outer-wall boundary condition, on top of the Kasagi source and Re_tau = 149.
#
#   dy+ 0.490 at the walls -> 4.87 at the centreline  (theirs: 0.49 -> 4.8)
#   dx+ 8.36, dz+ 4.18                                 (theirs: 14.8, 5.1)
#   d+  149, solid thickness = one half-height          (theirs: 149)
#   outer solid face: IMPOSED HEAT FLUX                 (theirs: imposed flux)
#
# The y line cannot come from any built-in distribution -- they all cluster at
# the DOMAIN ENDS, and this case needs the two interior fluid/solid interfaces
# resolved AND landing exactly on cell faces. make_ynodes.py writes it and
# [grid.y] nodes_file reads it.
#
# The outer-face BC is a NONZERO NEUMANN, which the solver already supported:
# bcValue is d(theta)/dy, so an imposed flux q is q/(kappa_s D_f), per scalar.
# The problem is then pure Neumann and the temperature level is free, so the
# seed puts theta = 0 AT THE INTERFACE -- which keeps the kappa_s = 0.1 scalar
# O(10) in the fluid instead of the O(600) the Dirichlet version carried.
#
# This removes BOTH structural differences from Flageul et al. at once:
#   * the SOURCE is Kasagi's beta*u_x ([scalar.N] source_type = velocity), so
#     the mean flux follows the CUMULATIVE FLOW RATE -- their flux profile
#     exactly, where run_bulk.sh's uniform source only approximates it;
#   * Re_tau is 149, not 180, matching theirs.
# The case file is re-prepared at Re = 149 (the IBM coefficients carry the
# 1/Re scaling; the geometry -- dwall, masks, leaf table -- comes out
# IDENTICAL, verified dataset by dataset).
#
# `source` is scaled so theta_tau = J_wall = source*int(u dy) ~ 1; the scalar
# problem is LINEAR, so that is numerical conditioning, not physics.
#
#   ./run_kasagi.sh [ic|develop|stats|all]
#
# Same grid, same case file and the same six conjugate scalars as run_cht.sh;
# the only difference is the thermal problem. Here a uniform volumetric source
# S = 1 heats the FLUID and both outer solid faces are held at zero, so the
# steady mean balance is dJ/dy = S and the wall-normal flux falls to zero at
# the centreline -- the structure of Flageul's case, whose source is Kasagi's
# f_T u_x. Ours is uniform rather than proportional to u(y); that is a <= 7 %
# difference in the flux profile and 2-4 % in the near-wall comparison region,
# against the 100 % difference at the centreline that run_cht.sh's constant
# flux gives.
#
# TWO PHASES, not four. There is no prepare (the case file is run_cht.sh's,
# unchanged -- the geometry is identical) and no reseed:
#   ic       take a DEVELOPED snapshot of the antisymmetric campaign, keep its
#            velocity, and overwrite the scalars with the bulk-heating profile.
#            J_wall = S h holds exactly and independently of kappa_s, so the
#            solid mean theta_i = J R_s is known in closed form and the seed is
#            already the answer -- which is why the reseed phase is gone;
#   develop  statistics off: the scalars shed the Reynolds-analogy seed's
#            wrong sublayer offset. The velocity is developed already.
#   stats    the measurement.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../../.. && pwd)

BIN=${BIN:-$ROOT/build_gpu/moby_solve}
RANKS=${RANKS:-1}
PY=${PY:-python3}
sel=${1:-all}
CASE=cht149_flageul.h5
SEED=${SEED:-kstat_80000.h5}       # the developed Re_tau 149 Kasagi field

DEV_STEPS=${DEV_STEPS:-80000}     # t -> ~24.5 (dt ~ 3.07e-4)
SETTLE_STEPS=${SETTLE_STEPS:-65000}    # t -> 20, at the correct capacities
CLEAN_STEPS=${CLEAN_STEPS:-650000}    # t -> 220: a 200-t.u. window = 29800
                                      # wall units, matching their 29000
EXT_STEPS=${EXT_STEPS:-311000}    # `extend`: doubles the window to ~28400 wall
                                  # units, matching Flageul's 29000
EXT_FROM=${EXT_FROM:-Fstat_391000.h5}
STAT_STEPS=${STAT_STEPS:-311000}  # t -> 95.5, the SAME averaging window
                                  # as the Kasagi run (14200 wall units)
SAMPLE=${SAMPLE:-25}
FLUSH=${FLUSH:-2000}
FIELD=${FIELD:-20000}

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }

# A run is SUCCESSFUL if the solver says so in its own log, not if mpirun
# returns 0. Teardown can fail long after the physics is finished -- this run
# lost a whole chained leg to PMIX "NO-PERMISSIONS" errors raised AFTER
# "main loop ended", with the snapshot already on disk -- so the exit code
# alone is the wrong test.
finished() {   # <logfile>
    grep -q "main loop ended" "$1"
}

ini() {  # <prefix> <nsteps> <write> <restart> <sample> <flush> <statsfile>
    sed -e "s|@CASE@|$CASE|" -e "s|@PREFIX@|$1|" -e "s|@NSTEPS@|$2|" \
        -e "s|@WRITE@|$3|" -e "s|@RESTART@|$4|" -e "s|@SAMPLE@|$5|" \
        -e "s|@STATS@|$6|" -e "s|@STATSFILE@|$7|" \
        cht149_flageul.ini > ".$1.ini"
}

if want ic; then
    echo "== initial condition: bulk-heating scalars on a developed velocity"
    $PY ./make_flageul_ic.py "$SEED" fmint_1.h5 IC_flageul.h5 || exit 1
fi

if want develop; then
    echo "== develop: t = 0 .. ~10 (statistics off)"
    ini Fdev "$DEV_STEPS" "$DEV_STEPS" IC_flageul.h5 0 0 x
    time mpirun -n "$RANKS" "$BIN" .Fdev.ini > Fdev.log 2>&1
    finished Fdev.log || { echo "DEVELOP FAILED"; tail -15 Fdev.log; exit 1; }
    grep -E "seconds_per_step" Fdev.log | tail -1
fi

# ---------------------------------------------------------------------------
# THE THREE-LEG PLAN. The first statistics campaign was contaminated: the solid
# diffusion time is d^2/alpha_s ~ 106 time units and the window was t = 24..120,
# i.e. LESS THAN ONE such time, started from a 24-time-unit develop leg. The
# damage is measurable -- decomposing the deep-solid variance into spatial and
# temporal parts, the temporal term falls by a factor 45 when the transient is
# dropped, and the outer-face value goes from 4.6x the reference to 1.2x.
#   settle  : run at the CORRECT capacities with no statistics, after
#             reseed_solid_exact.py has put every solid on its exact steady
#             mean and the source has been corrected to balance the outflux
#   clean   : statistics FROM ZERO over one uncontaminated window
if want settle; then
    echo "== settle: $SETTLE_STEPS steps at the correct capacities, statistics OFF"
    ini Cdev "$SETTLE_STEPS" "$SETTLE_STEPS" IC_clean.h5 0 0 x
    time mpirun -n "$RANKS" "$BIN" .Cdev.ini > Cdev.log 2>&1
    finished Cdev.log || { echo "SETTLE FAILED"; tail -15 Cdev.log; exit 1; }
    grep -E "seconds_per_step" Cdev.log | tail -1
fi

if want clean; then
    echo "== clean statistics FROM ZERO ($CLEAN_STEPS steps)"
    rm -f flageul_clean.h5 Cstat_vel.h5
    ini Cstat "$CLEAN_STEPS" "$FIELD" "Cdev_${SETTLE_STEPS}.h5" "$SAMPLE" "$FLUSH" flageul_clean.h5
    time mpirun -n "$RANKS" "$BIN" .Cstat.ini > Cstat.log 2>&1
    finished Cstat.log || { echo "CLEAN FAILED"; tail -15 Cstat.log; exit 1; }
    grep -E "seconds_per_step" Cstat.log | tail -1
fi

if want extend; then
    # CONTINUE the accumulated statistics rather than starting a new window:
    # read_stats_restart reads the stats file itself, so as long as it is NOT
    # deleted the sums carry over. The deep-solid variance is limited by the
    # window length, not by resolution (measured: the instantaneous
    # within-plane variance there is ~0.012 against an accumulated 0.059, so
    # the excess is temporal wander the window has not averaged out), which is
    # why this is the fix and a finer solid grid is not.
    echo "== extend: +$EXT_STEPS steps, statistics CONTINUED from $EXT_FROM"
    [ -f flageul_stats.h5 ] || { echo "no flageul_stats.h5 to extend"; exit 1; }
    ini Fstat "$EXT_STEPS" "$FIELD" "$EXT_FROM" "$SAMPLE" "$FLUSH" flageul_stats.h5
    time mpirun -n "$RANKS" "$BIN" .Fstat.ini > Fext.log 2>&1
    finished Fext.log || { echo "EXTEND FAILED"; tail -15 Fext.log; exit 1; }
    grep -E "seconds_per_step" Fext.log | tail -1
fi

if want stats; then
    echo "== statistics"
    rm -f flageul_stats.h5 Fstat_vel.h5
    ini Fstat "$STAT_STEPS" "$FIELD" "Fdev_${DEV_STEPS}.h5" "$SAMPLE" "$FLUSH" flageul_stats.h5
    time mpirun -n "$RANKS" "$BIN" .Fstat.ini > Fstat.log 2>&1
    finished Fstat.log || { echo "STATS FAILED"; tail -15 Fstat.log; exit 1; }
    grep -E "seconds_per_step" Fstat.log | tail -1
fi
