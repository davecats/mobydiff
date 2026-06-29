#!/usr/bin/env bash
# Reflux ON vs OFF developed-stats study for the 2:1 wall-band channel.
#
# WHY: the u'/v' interface BANDS are a MOMENTUM-REFLUX artifact, not an interface
# transfer/ownership or energy defect. The reflux replaces the coarse interface
# cell's flux F_coarse by avg(F_fine); for the normal flux (q)^2,
# avg(of squares) >> (avg)^2 (Jensen), so the reflux injects the fine-side resolved
# Reynolds-stress flux into the under-resolved coarse cell -> the band (Cevheri &
# Stoesser energy accumulation). On the benchmark (250 steps) reflux OFF erases the
# u' spike (excess 1.56->1.00) and v' step (kink 0.31->0.04) at zero stability /
# divergence cost. See docs/next_session_edges_les.md "MECHANISM
# CONFIRMED". The reflux exists to conserve the MEAN interface momentum flux
# (-<u'v'>), so this study checks the TRADE on the converged statistics: does
# reflux OFF degrade the mean profile / Reynolds shear, or is it simply better?
#
# This runs TWO developed two-leg cases (transient t=0..5 discarded, then stats
# t=5..25): reflux ON (the current default) and reflux OFF. Each is ~7 h on 2 GPUs
# (~14 h on 1). The solver stays GENERIC -- this only toggles the existing
# [blocks] momentum_reflux config flag; nothing channel-specific is added.
#
# Usage (on the run machine, after building):
#   module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
#   ./compile.sh gpu
#   cd validation/channel_interface/developed
#   ./run_reflux_study.sh [gpu|cpu] [nranks]      # default: gpu 2
#
# Needs tutorials/channel_kmm180/channel_kmm180_restart.h5 (the IC source); the
# refined IC.h5 is generated automatically on first use.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ARCH="${1:-gpu}"
RANKS="${2:-2}"
ROOT="$(cd "$HERE/../../.." && pwd)"
PY="${PYTHON:-python3}"

echo "==== reflux study: arch=$ARCH ranks=$RANKS ===="

# reflux ON (current default) and reflux OFF, each a full transient+stats run.
"$PY" "$HERE/run_developed.py" --arch "$ARCH" --ranks "$RANKS" --name reflux_on
"$PY" "$HERE/run_developed.py" --arch "$ARCH" --ranks "$RANKS" --no-reflux --name reflux_off

ON="$HERE/runs/reflux_on/stats/channel_stats.h5"
OFF="$HERE/runs/reflux_off/stats/channel_stats.h5"
echo
echo "==== done. Overlay the converged statistics: ===="
echo "  $PY $ROOT/tools/plot_channel_stats.py reflux_study.png \\"
echo "      $ON:reflux-on \\"
echo "      $OFF:reflux-off"
echo
echo "  To add the uniform-fine ground truth, append e.g.:"
echo "      <path-to-reference>/channel_stats.h5:uniform-fine"
echo
echo "  PASS for the fix = reflux-off REMOVES the u'/v' interface band in u'/v' rms"
echo "  WITHOUT degrading the mean U(y+) or the Reynolds shear -<u'v'> at the"
echo "  interface vs the uniform-fine reference. If reflux-off DOES degrade the"
echo "  mean flux, the reflux is needed but must be made band-free in a GENERIC"
echo "  (non-channel-specific) way -- see the handoff doc."
