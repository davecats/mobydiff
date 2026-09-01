#!/usr/bin/env bash
# The BULK-HEATING conjugate turbulent channel (Re_tau = 180, Pr = 0.71).
#
#   ./run_bulk.sh [ic|develop|stats|all]
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
CASE=cht180_nb32.h5
SEED=${SEED:-stat_80000.h5}       # a developed field of the antisymmetric run

DEV_STEPS=${DEV_STEPS:-7000}      # t -> ~10
STAT_STEPS=${STAT_STEPS:-70000}   # t -> ~101, i.e. ~18000 wall units
SAMPLE=${SAMPLE:-25}
FLUSH=${FLUSH:-2000}
FIELD=${FIELD:-20000}

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }

ini() {  # <prefix> <nsteps> <write> <restart> <sample> <flush> <statsfile>
    sed -e "s|@CASE@|$CASE|" -e "s|@PREFIX@|$1|" -e "s|@NSTEPS@|$2|" \
        -e "s|@WRITE@|$3|" -e "s|@RESTART@|$4|" -e "s|@SAMPLE@|$5|" \
        -e "s|@STATS@|$6|" -e "s|@STATSFILE@|$7|" \
        cht180_bulk.ini > ".$1.ini"
}

if want ic; then
    echo "== initial condition: bulk-heating scalars on a developed velocity"
    $PY ./make_bulk_ic.py "$SEED" IC_bulk.h5 || exit 1
fi

if want develop; then
    echo "== develop: t = 0 .. ~10 (statistics off)"
    ini bdev "$DEV_STEPS" "$DEV_STEPS" IC_bulk.h5 0 0 x
    time mpirun -n "$RANKS" "$BIN" .bdev.ini > bdev.log 2>&1 || { tail -15 bdev.log; exit 1; }
    grep -E "seconds_per_step" bdev.log | tail -1
fi

if want stats; then
    echo "== statistics"
    rm -f bulk_stats.h5 bstat_vel.h5
    ini bstat "$STAT_STEPS" "$FIELD" "bdev_${DEV_STEPS}.h5" "$SAMPLE" "$FLUSH" bulk_stats.h5
    time mpirun -n "$RANKS" "$BIN" .bstat.ini > bstat.log 2>&1 || { tail -15 bstat.log; exit 1; }
    grep -E "seconds_per_step" bstat.log | tail -1
fi
