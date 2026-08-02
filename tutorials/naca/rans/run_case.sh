#!/usr/bin/env bash
# NACA 0012 / alpha 5 / Re 4e5 validation case — RE-RUN the flow.
#
#   ./run_case.sh restart   continue the production (nose-refined) case from
#                           the shipped converged state c11_nose_640000.h5
#                           (t = 34.75) — minutes to hours, for verification.
#   ./run_case.sh scratch   full reproduction: staged L10 -> L11 -> nose band.
#                           Days on one GPU; see the README timings.
#   ./run_case.sh base      the level-11-everywhere BASELINE (c11_aoa5.ini),
#                           the "before" of the nose-refinement comparison.
#                           Needs its own converged state c11_aoa5_450013.h5.
#
# Post-processing (forces / Cp / Cf + the OpenFOAM overlay): ./postprocess.sh
#
# Prerequisites: the solver built (../../..; compile.sh cpu && compile.sh gpu).
# Prepared geometry is regenerated from assets/geometry/n0012_b11.stl if absent.
# NOTE moby_prepare is MPI-parallel only (the CPU build has no OpenMP, so the
# classify pragmas are inert) — give it many ranks, not many threads.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$HERE/../../..
PY=${PY:-python3}
MPIRUN=${MPIRUN:-mpirun}
NRANK=${NRANK:-20}                       # prepare ranks; the classify stage scales with this
cd "$HERE"

COEF_BASE=assets/geometry/ibm_coeff_c11.h5
COEF_NOSE=assets/geometry/ibm_coeff_c11_nose.h5
STATE_NOSE=c11_nose_640000.h5
STATE_BASE=c11_aoa5_450013.h5

prepare() {  # $1 = prep ini, $2 = output case file
    [ -f "$2" ] && { echo "== $2 exists, skipping prepare"; return; }
    echo "== moby_prepare $1 -> $2  (~40 min at $NRANK ranks for the nose case)"
    $MPIRUN -n "$NRANK" --oversubscribe --bind-to none \
        "$ROOT/build_cpu/moby_prepare" "$1" "$2"
}

case "${1:-restart}" in
restart)
    prepare .prep_c11_nose.ini "$COEF_NOSE"
    [ -f "$STATE_NOSE" ] || { echo "no $STATE_NOSE; use '$0 scratch'"; exit 1; }
    echo "== restarting the nose-refined case from $STATE_NOSE"
    sed "s|^file = .*|file = $STATE_NOSE|" c11_aoa5_nose.ini > .run.ini
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .run.ini
    ;;
base)
    prepare .prep_c11.ini "$COEF_BASE"
    [ -f "$STATE_BASE" ] || { echo "no $STATE_BASE (not shipped); use '$0 scratch'"; exit 1; }
    echo "== restarting the level-11 baseline from $STATE_BASE"
    sed "s|^file = .*|file = $STATE_BASE|" c11_aoa5.ini > .run.ini
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .run.ini
    ;;
scratch)
    # Stage 1: the SAME physics on the L10 twin (2x coarser wall band,
    # dt 1e-4) carries the whole transient at ~5x lower cost.
    sed 's/refine_levels = 11/refine_levels = 10/' .prep_c11.ini > .prep_l10.ini
    prepare .prep_l10.ini assets/geometry/.ibm_coeff_l10.h5
    prepare .prep_c11.ini "$COEF_BASE"
    sed -e 's/refine_levels = 11/refine_levels = 10/' \
        -e "s|$COEF_BASE|assets/geometry/.ibm_coeff_l10.h5|" \
        -e 's/dt = 5.0e-5/dt = 1.0e-4/' -e 's/dtmax = 5.0e-5/dtmax = 1.0e-4/' \
        -e 's/^t_final = .*/t_final = 12.0/' \
        -e 's/field_prefix = c11_aoa5/field_prefix = .l10_stage/' \
        -e '/^\[restart\]/,+1d' c11_aoa5.ini > .stage1.ini
    echo "== stage 1: L10 transient to t = 12 (~120k steps)"
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .stage1.ini
    LAST=$(ls -t .l10_stage_*.h5 | head -1)
    echo "== interpolating $LAST onto the L11 layout"
    $PY interp_restart.py "$LAST" "$COEF_BASE" .restart_l11.h5
    sed -e 's|^file = .*|file = .restart_l11.h5|' -e 's/^t_final = .*/t_final = 30.0/' \
        c11_aoa5.ini > .stage2.ini
    echo "== stage 2: L11 to convergence (t = 30; extend if C_L still drifts)"
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .stage2.ini
    # Stage 3: add the level-12 nose band and re-converge. Cf settles within
    # ~1 t.u.; Cp needs ~3 (it follows the global circulation) -> t = 35.
    prepare .prep_c11_nose.ini "$COEF_NOSE"
    LAST=$(ls -t c11_aoa5_*.h5 | head -1)
    echo "== interpolating $LAST onto the nose-refined layout"
    $PY interp_restart.py "$LAST" "$COEF_NOSE" .restart_nose.h5
    echo "== stage 3: nose band to t = 35 (~42 h on an RTX 3060)"
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" c11_aoa5_nose.ini
    ;;
*)
    echo "usage: $0 [restart|base|scratch]  (post-processing: ./postprocess.sh)"; exit 1;;
esac
