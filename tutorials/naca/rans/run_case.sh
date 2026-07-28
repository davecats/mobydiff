#!/usr/bin/env bash
# NACA 0012 / alpha 5 / Re 4e5 validation case — run + compare vs OpenFOAM.
#
#   ./run_case.sh restart   continue from a converged/near-converged state
#                           (c11_aoa5_450013.h5 or c11_aoa5ab_370013.h5) —
#                           minutes to hours, for verification.
#   ./run_case.sh scratch   full reproduction: staged L10 -> L11 protocol
#                           (~10-14 h on a modern GPU).
#   ./run_case.sh post      post-process the newest snapshot + OF overlay.
#
# Prerequisites: the solver built (../../..; compile.sh cpu && compile.sh gpu),
# python with numpy/h5py/scipy/matplotlib for the post-processing.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$HERE/../../..
PY=${PY:-python3}
MPIRUN=${MPIRUN:-mpirun}
cd "$HERE"

prepare() {  # $1 = prep ini, $2 = output case file
    [ -f "$2" ] && { echo "== $2 exists, skipping prepare"; return; }
    echo "== moby_prepare $1 -> $2"
    $MPIRUN -n 4 --oversubscribe "$ROOT/build_cpu/moby_prepare" "$1" "$2"
}

case "${1:-post}" in
restart)
    prepare .prep_c11.ini ibm_coeff_c11.h5
    RESTART=""
    for f in c11_aoa5_450013.h5 c11_aoa5ab_370013.h5; do
        [ -f "$f" ] && RESTART=$f && break
    done
    [ -n "$RESTART" ] || { echo "no restart state found; use '$0 scratch'"; exit 1; }
    echo "== restarting from $RESTART"
    sed "s|^file = .*|file = $RESTART|" c11_aoa5.ini > .run.ini
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .run.ini
    ;;
scratch)
    # Stage 1: the SAME physics on the L10 twin (2x coarser wall band,
    # dt 1e-4) carries the whole transient at ~5x lower cost.
    sed 's/refine_levels = 11/refine_levels = 10/' .prep_c11.ini > .prep_l10.ini
    prepare .prep_l10.ini .ibm_coeff_l10.h5
    prepare .prep_c11.ini ibm_coeff_c11.h5
    sed -e 's/refine_levels = 11/refine_levels = 10/' \
        -e 's/ibm_coeff_c11.h5/.ibm_coeff_l10.h5/' \
        -e 's/dt = 5.0e-5/dt = 1.0e-4/' -e 's/dtmax = 5.0e-5/dtmax = 1.0e-4/' \
        -e 's/^t_final = .*/t_final = 12.0/' \
        -e 's/field_prefix = c11_aoa5/field_prefix = .l10_stage/' \
        -e '/^\[restart\]/,+1d' c11_aoa5.ini > .stage1.ini
    echo "== stage 1: L10 transient to t = 12 (~120k steps)"
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .stage1.ini
    LAST=$(ls -t .l10_stage_*.h5 | head -1)
    echo "== interpolating $LAST onto the L11 layout"
    $PY interp_restart.py "$LAST" ibm_coeff_c11.h5 .restart_l11.h5
    sed -e 's|^file = .*|file = .restart_l11.h5|' -e 's/^t_final = .*/t_final = 18.0/' \
        c11_aoa5.ini > .stage2.ini
    echo "== stage 2: L11 to convergence (t = 18; extend if C_L still drifts)"
    $MPIRUN -n 1 "$ROOT/build_gpu/moby_solve" .stage2.ini
    ;;
post)
    SNAP=$(ls -t c11_aoa5_*.h5 2>/dev/null | head -1)
    [ -n "$SNAP" ] || { echo "no c11_aoa5_*.h5 snapshot found"; exit 1; }
    echo "== post-processing $SNAP"
    $PY cv_forces.py "$SNAP" --boxes 1.5 2.5 --aoa 5
    $PY surface_cp_cf.py "$SNAP" --out cpcf_c11_aoa5_final.npz --plot cpcf_c11_aoa5_final.png
    $PY compare_openfoam.py --npz cpcf_c11_aoa5_final.npz --out cpcf_c11_final_vs_openfoam.png
    echo "== targets: C_L 0.514 (OF 0.5142), C_D 0.013 (OF 0.0134), Cp_min ~-1.79 (OF -1.780)"
    ;;
*)
    echo "usage: $0 [restart|scratch|post]"; exit 1;;
esac
