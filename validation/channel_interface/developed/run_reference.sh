#!/usr/bin/env bash
# Uniform-fine REFERENCE run (256x128x256, single level) for the 2:1-interface
# study: the ground truth the refined reflux_on / reflux_off runs are compared
# against. Two-leg, matching those runs (transient t=0..5 discarded, then stats
# t=5..25), same dtmax / niter / sor / Chebyshev (all already in reference.ini).
#
# The refined fine level shares this uniform grid bitwise (reference.ini grids are
# `subdivided = true`), so the interface bands sit on cells this run also has --
# the apples-to-apples comparison for the residual interface wiggles.
#
# Usage (on the run machine, after building):
#   module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
#   ./compile.sh gpu
#   cd validation/channel_interface/developed
#   ./run_reference.sh [gpu|cpu] [nranks]          # default: gpu 1
#
# Needs tutorials/channel_kmm180/channel_kmm180_restart.h5 (the IC source); the
# uniform REF_IC.h5 is generated automatically on first use. Env overrides
# T_TRANSIENT / T_AVERAGE (defaults 5 / 20) -- used tiny for smoke tests.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
ARCH="${1:-gpu}"
RANKS="${2:-1}"
T_TRANSIENT="${T_TRANSIENT:-5.0}"
T_AVERAGE="${T_AVERAGE:-20.0}"
FIELD1="${FIELD1:-2000}"   # leg-1 field write interval (restart picked from the last one)
PY="${PYTHON:-python3}"

BIN="$ROOT/build_$ARCH/main"
REFINI="$ROOT/validation/channel_interface/reference.ini"
SRC="$ROOT/tutorials/channel_kmm180/channel_kmm180_restart.h5"
RUN="$HERE/runs/reference"
IC="$RUN/REF_IC.h5"
T_END=$(awk "BEGIN{print $T_TRANSIENT + $T_AVERAGE}")
mkdir -p "$RUN/transient" "$RUN/stats"

[ -x "$BIN" ] || { echo "binary $BIN not found -- build with ./compile.sh $ARCH"; exit 1; }

# 0. uniform-fine initial condition
if [ ! -f "$IC" ]; then
    [ -f "$SRC" ] || { echo "source restart $SRC not found (see tutorials/channel_kmm180)"; exit 1; }
    echo "== generating uniform REF_IC.h5"
    "$PY" "$ROOT/tools/make_channel_restart.py" --mode reference --source "$SRC" --out "$IC"
fi

# Helper: write an input.ini from reference.ini with the given key overrides, and
# append an [mpi] dims section when running multi-rank (reference.ini has none).
gen_ini() {  # $1 dest  $2 restart  $3 t_final  $4 stats_sample  $5 stats_write  $6 field_interval
    sed -e "s#RESTART_PLACEHOLDER#$2#" \
        -e "s/^t_final = .*/t_final = $3/" \
        -e "s/^stats_sample_interval = .*/stats_sample_interval = $4/" \
        -e "s/^stats_write_interval = .*/stats_write_interval = $5/" \
        -e "s/^field_interval = .*/field_interval = $6/" \
        "$REFINI" > "$1"
    if [ "$RANKS" -gt 1 ]; then printf '\n[mpi]\ndims = %d 1 1\n' "$RANKS" >> "$1"; fi
}

echo "== LEG 1 transient t=0..$T_TRANSIENT (stats off), arch=$ARCH ranks=$RANKS =="
gen_ini "$RUN/transient/input.ini" "$IC" "$T_TRANSIENT" -1 -1 "$FIELD1"
( cd "$RUN/transient" && mpirun -n "$RANKS" "$BIN" input.ini > run.log 2>&1 )
LAST=$(ls -t "$RUN"/transient/channel_field_*.h5 | head -1)
echo "   transient final field: $LAST"

echo "== LEG 2 stats t=$T_TRANSIENT..$T_END (stats on, fresh) =="
rm -f "$RUN"/stats/channel_stats.h5 "$RUN"/stats/channel_stats_l*.h5
gen_ini "$RUN/stats/input.ini" "$LAST" "$T_END" 50 20000 0
( cd "$RUN/stats" && mpirun -n "$RANKS" "$BIN" input.ini > run.log 2>&1 )

STATS="$RUN/stats/channel_stats.h5"
echo
echo "== DONE. reference stats: $STATS =="
echo "   overlay against the refined runs:"
echo "   $PY $ROOT/tools/plot_channel_stats.py reflux_study.png \\"
echo "       $HERE/runs/reflux_on/stats/channel_stats.h5:reflux-on \\"
echo "       $HERE/runs/reflux_off/stats/channel_stats.h5:reflux-off \\"
echo "       $STATS:uniform"
