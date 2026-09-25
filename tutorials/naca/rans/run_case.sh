#!/usr/bin/env bash
# NACA 0012 / alpha 5 / Re 4e5 validation case — RE-RUN the flow.
#
#   ./run_case.sh restart   continue the production (nose-refined) case from
#                           the shipped converged state c11_nose_660000.h5
#                           (t = 34.0) — minutes to hours, for verification.
#   ./run_case.sh scratch   full reproduction: staged L10 -> L11 -> nose band.
#                           Days on one GPU; see the README timings. Each
#                           stage stops on its own steadiness measure
#                           ([case.airfoil] steady_tol, tuned per stage
#                           below and overridable in .steady.conf), with
#                           t_final as the safety net; per-stage output
#                           goes to .stage{1,2,3}.log.
#
# c11_aoa5.ini / .prep_c11.ini are the level-11 stage of that protocol; they
# are inputs to `scratch`, not a case to run on their own.
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
STATE_NOSE=c11_nose_660000.h5

prepare() {  # $1 = prep ini, $2 = output case file
    [ -f "$2" ] && { echo "== $2 exists, skipping prepare"; return; }
    echo "== moby_prepare $1 -> $2  (~40 min at $NRANK ranks for the nose case)"
    $MPIRUN -n "$NRANK" --oversubscribe --bind-to none \
        "$ROOT/build_cpu/moby_prepare" "$1" "$2"
}

# --- steady-state stop, per stage ---------------------------------------
# Each `scratch` stage ends on the control-volume budget's OWN unsteady
# term ([case.airfoil] steady_tol: |2 dmom/dt|/qref, in coefficient units),
# not on a hand-tuned t_final — that measure is what the forces trace's
# `dmomdt` column reports, so a stage documents its own convergence.
# t_final stays as a SAFETY NET: a stage that never reaches the tolerance
# still terminates, at the time the previous campaign used.
#
# TWO things decide these values, both MEASURED in the 2026-08 from-scratch
# reproduction (README "Convergence and stopping"):
#
# (1) The SAMPLING INTERVAL, not just the tolerance. The control volume's
#     momentum carries a coherent f ~ 35 U/c oscillation that cancels in the
#     reported coefficients but not in a one-interval difference, so the
#     measure floors at ~2A/dt_sample: at the inis' 20-step sampling that is
#     0.19 during the transient and ~8e-3 even when converged — unfireable.
#     Each stage therefore samples ~0.25 t.u. apart (FSI below), which drops
#     the floor two decades and lets the measure read the physical drift.
# (2) The WINDOW: steady_samples x FSI x dt = ONE time unit per stage (the
#     built-in 3 samples would be 1.5e-3 t.u., which any plateau trips).
#
# The SPANS are not margin. After the L10 -> L11 interpolation C_L relaxes
# with tau ~ 9 t.u., so stage 2's 18 t.u. is ~2 time constants; shortening
# it to 6 was tried and under-reads lift by 2 %. Stage 3 keeps the published
# 5 t.u. and should RUN IT OUT: its tolerance certifies the momentum budget,
# which the converged circulation already satisfies, while Cp needs ~3 t.u.
#
# The keys live HERE, not in the committed inis, so c11_aoa5.ini and
# c11_aoa5_nose.ini keep reproducing the documented t_final behaviour.
# Exercised in the reference run: stage 2 stopped on 2e-3 at t = 29, stage 3
# on 1e-3 at t = 31 (and was continued to 34 by design); stage 1 ran to its
# net. FSI is per stage because dt is.
S1_TOL=${S1_TOL:-2.0e-3};  S1_SAMPLES=${S1_SAMPLES:-4};  S1_FSI=${S1_FSI:-2500}
S2_TOL=${S2_TOL:-2.0e-3};  S2_SAMPLES=${S2_SAMPLES:-4};  S2_FSI=${S2_FSI:-5000}
S3_TOL=${S3_TOL:-1.0e-3};  S3_SAMPLES=${S3_SAMPLES:-4};  S3_FSI=${S3_FSI:-10000}
S1_TEND=${S1_TEND:-12.0};  S2_SPAN=${S2_SPAN:-18.0};     S3_SPAN=${S3_SPAN:-5.0}

# Re-read before every stage, so a multi-day detached run can be retuned
# between stages (editing THIS script while bash is executing it is not
# safe; a sourced file is).
steady_conf() { [ -f "$HERE/.steady.conf" ] && . "$HERE/.steady.conf" || true; }

# sed program setting the stopping keys for a stage: tolerance, window and
# the sampling interval the measure is differenced over ($1 $2 $3).
steady_sed() {
    printf '/^\\[case\\.airfoil\\]/a steady_tol = %s\n/^\\[case\\.airfoil\\]/a steady_samples = %s\ns/^force_sample_interval = .*/force_sample_interval = %s/\n' \
        "$1" "$2" "$3"
}

# t_current of a snapshot: the stages hand over at whatever time the
# tolerance stopped them, so the downstream safety nets are RELATIVE.
t_of() { $PY -c "import h5py,sys;print('%.4f'%h5py.File(sys.argv[1],'r').attrs['t_current'])" "$1"; }

case "${1:-restart}" in
restart)
    prepare .prep_c11_nose.ini "$COEF_NOSE"
    [ -f "$STATE_NOSE" ] || { echo "no $STATE_NOSE; use '$0 scratch'"; exit 1; }
    echo "== restarting the nose-refined case from $STATE_NOSE"
    sed "s|^file = .*|file = $STATE_NOSE|" c11_aoa5_nose.ini > .run.ini
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .run.ini
    ;;
scratch)
    # Stage 1: the SAME physics on the L10 twin (2x coarser wall band,
    # dt 1e-4) carries the whole transient at ~5x lower cost.
    sed 's/refine_levels = 11/refine_levels = 10/' .prep_c11.ini > .prep_l10.ini
    prepare .prep_l10.ini assets/geometry/.ibm_coeff_l10.h5
    prepare .prep_c11.ini "$COEF_BASE"
    steady_conf
    sed -e 's/refine_levels = 11/refine_levels = 10/' \
        -e "s|$COEF_BASE|assets/geometry/.ibm_coeff_l10.h5|" \
        -e 's/dt = 5.0e-5/dt = 1.0e-4/' -e 's/dtmax = 5.0e-5/dtmax = 1.0e-4/' \
        -e "s/^t_final = .*/t_final = $S1_TEND/" \
        -e 's/^runtime_file = .*/runtime_file = forces_stage1_l10.txt/' \
        -e 's/field_prefix = c11_aoa5/field_prefix = .l10_stage/' \
        -e '/^\[restart\]/,+1d' \
        -e "$(steady_sed "$S1_TOL" "$S1_SAMPLES" "$S1_FSI")" c11_aoa5.ini > .stage1.ini
    echo "== $(date +%F\ %T) stage 1: L10 transient, steady_tol $S1_TOL over" \
         "$S1_SAMPLES samples (net t = $S1_TEND) -> .stage1.log"
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .stage1.ini > .stage1.log 2>&1
    LAST=$(ls -t .l10_stage_*.h5 | head -1)
    echo "== interpolating $LAST (t = $(t_of "$LAST")) onto the L11 layout"
    $PY interp_restart.py "$LAST" "$COEF_BASE" .restart_l11.h5

    steady_conf
    T2=$($PY -c "print('%.4f'%($(t_of "$LAST") + $S2_SPAN))")
    sed -e 's|^file = .*|file = .restart_l11.h5|' -e "s/^t_final = .*/t_final = $T2/" \
        -e 's/^runtime_file = .*/runtime_file = forces_stage2_c11.txt/' \
        -e "$(steady_sed "$S2_TOL" "$S2_SAMPLES" "$S2_FSI")" c11_aoa5.ini > .stage2.ini
    echo "== $(date +%F\ %T) stage 2: L11, steady_tol $S2_TOL over $S2_SAMPLES" \
         "samples (net t = $T2) -> .stage2.log"
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .stage2.ini > .stage2.log 2>&1

    # Stage 3: add the level-12 nose band and re-converge. Cf settles within
    # ~1 t.u.; Cp needs ~3 (it follows the global circulation).
    prepare .prep_c11_nose.ini "$COEF_NOSE"
    LAST=$(ls -t c11_aoa5_*.h5 | head -1)
    echo "== interpolating $LAST (t = $(t_of "$LAST")) onto the nose-refined layout"
    $PY interp_restart.py "$LAST" "$COEF_NOSE" .restart_nose.h5

    steady_conf
    T3=$($PY -c "print('%.4f'%($(t_of "$LAST") + $S3_SPAN))")
    sed -e 's|^file = .*|file = .restart_nose.h5|' -e "s/^t_final = .*/t_final = $T3/" \
        -e "$(steady_sed "$S3_TOL" "$S3_SAMPLES" "$S3_FSI")" c11_aoa5_nose.ini > .stage3.ini
    echo "== $(date +%F\ %T) stage 3: nose band, steady_tol $S3_TOL over $S3_SAMPLES" \
         "samples (net t = $T3; ~12.7 h per time unit on an RTX 3060) -> .stage3.log"
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .stage3.ini > .stage3.log 2>&1
    echo "== $(date +%F\ %T) converged state: $(ls -t c11_nose_*.h5 | head -1)"
    ;;
*)
    echo "usage: $0 [restart|scratch]  (post-processing: ./postprocess.sh)"; exit 1;;
esac
