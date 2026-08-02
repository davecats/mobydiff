#!/usr/bin/env bash
# NACA 0012 / alpha 5 / Re 4e5 — regenerate the validation statistics and
# figures from a converged snapshot, and overlay against OpenFOAM.
#
#   ./postprocess.sh [snapshot.h5]
#
# Default: the shipped converged state of the production (nose-refined)
# case. The case is detected from the snapshot name, which selects the
# prepared case file and how many levels carry the surface:
#
#   c11_nose_*  -> nose-refined: levels 12 (nose band) + 11 elsewhere
#   c11_aoa5_*  -> level-11 baseline: one level
#
# ALWAYS post-process a REGULAR-CADENCE snapshot. The final one a run
# writes sits on a dt-clipped micro-step, which inflates the stored
# incremental pn and corrupts every pressure-based quantity.
#
# Outputs:
#   assets/referenceStats/cpcf_<case>_final{.npz,_cp.dat,_cf.dat}
#   assets/figures/cpcf_<case>_final.png
#   assets/figures/cpcf_<case>_vs_openfoam.png
#   assets/figures/cf_crosscheck.png          (baseline only, see below)
#
# Prerequisites: python with numpy/h5py/scipy/matplotlib.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
PY=${PY:-python3}
cd "$HERE"

SNAP=${1:-}
if [ -z "$SNAP" ]; then
    SNAP=$(ls -t c11_nose_*.h5 c11_aoa5_*.h5 2>/dev/null | head -1)
fi
[ -n "$SNAP" ] && [ -f "$SNAP" ] || { echo "no snapshot found (pass one explicitly)"; exit 1; }

case "$(basename "$SNAP")" in
c11_nose_*) COEF=assets/geometry/ibm_coeff_c11_nose.h5; LEVELS=2; TAG=c11_nose ;;
*)          COEF=assets/geometry/ibm_coeff_c11.h5;      LEVELS=1; TAG=c11_aoa5 ;;
esac
# per-case figure names, so a baseline run never overwrites the production ones
case "$TAG" in
c11_nose)   FIG=assets/figures/cpcf_c11_nose_final.png
            FIGOF=assets/figures/cpcf_c11_nose_vs_openfoam.png ;;
*)          FIG=assets/figures/cpcf_c11_aoa5_final.png
            FIGOF=assets/figures/cpcf_c11_final_vs_openfoam.png ;;
esac
NPZ=assets/referenceStats/cpcf_${TAG}_final.npz
echo "== post-processing $SNAP  (coef $COEF)"

$PY postProcess/cv_forces.py "$SNAP" --boxes 1.5 2.5 --aoa 5
# --coef: the IBM coefficient tiles gate the penalization band out of the
# wall-gradient fit (see surface_cp_cf.py); without it Cf under-reads badly.
$PY postProcess/surface_cp_cf.py "$SNAP" --coef "$COEF" \
    --out "$NPZ" --plot "$FIG"
$PY postProcess/compare_openfoam.py --npz "$NPZ" \
    --out "$FIGOF"

# Two independent Cf estimators, as a check that the wall-gradient measurement
# is measuring the field and not itself. Single-level surfaces only -- it is a
# validation tool, and the three-estimator agreement was established on the
# level-11 baseline.
if [ "$LEVELS" = "1" ]; then
    $PY postProcess/cf_crosscheck.py "$SNAP" --coef "$COEF" --npz "$NPZ" \
        --plot assets/figures/cf_crosscheck.png
    mv -f cf_crosscheck_*.npz assets/referenceStats/
else
    echo "== skipping cf_crosscheck.py (single-level surfaces only)"
fi

if [ "$LEVELS" = "2" ]; then
    $PY postProcess/plot_nose_comparison.py
    echo "== targets (nose-refined, t = 34.75): C_L 0.520 +- 0.005 (OF 0.5142),"
    echo "   C_D 0.0128 +- 0.0008 (OF 0.0134), Cp_min -1.7686 (OF -1.7797)"
else
    echo "== targets (level-11 baseline): C_L 0.514 (OF 0.5142), C_D 0.0130 (OF 0.0134),"
    echo "   Cp_min -1.775 (OF -1.780)"
fi
