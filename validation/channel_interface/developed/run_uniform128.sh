#!/usr/bin/env bash
# Uniform BASE-resolution reference (128x64x128, no refinement), two-leg,
# matching the reflux / reference runs (transient t=0..5 discarded, stats
# t=5..25). Its grid lines ARE the refined case's level-0 (coarse core) lines
# bitwise, so refined-core vs this isolates the 2:1-interface effect from the
# coarse-resolution deficit. Cheap (~1M cells) -- runs fast even on one laptop GPU.
#
# Usage (after building):
#   module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
#   ./compile.sh gpu
#   cd validation/channel_interface/developed
#   ./run_uniform128.sh [gpu|cpu] [nranks]          # default: gpu 1
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
ARCH="${1:-gpu}"; RANKS="${2:-1}"
T_TRANSIENT="${T_TRANSIENT:-5.0}"; T_AVERAGE="${T_AVERAGE:-20.0}"; FIELD1="${FIELD1:-2000}"
PY="${PYTHON:-python3}"
BIN="$ROOT/build_$ARCH/main"
INI="$ROOT/validation/channel_interface/uniform128.ini"
SRC="$ROOT/tutorials/channel_kmm180/channel_kmm180_restart.h5"
RUN="$HERE/runs/uniform128"; IC="$RUN/BASE_IC.h5"
T_END=$(awk "BEGIN{print $T_TRANSIENT + $T_AVERAGE}")
mkdir -p "$RUN/transient" "$RUN/stats"
[ -x "$BIN" ] || { echo "binary $BIN not found -- build with ./compile.sh $ARCH"; exit 1; }

if [ ! -f "$IC" ]; then
    [ -f "$SRC" ] || { echo "source restart $SRC not found"; exit 1; }
    echo "== generating base-128 BASE_IC.h5"
    "$PY" "$ROOT/tools/make_channel_restart.py" --mode base --source "$SRC" --out "$IC"
fi

gen_ini() {  # dest restart t_final stats_sample stats_write field_interval
    sed -e "s#RESTART_PLACEHOLDER#$2#" \
        -e "s/^t_final = .*/t_final = $3/" \
        -e "s/^stats_sample_interval = .*/stats_sample_interval = $4/" \
        -e "s/^stats_write_interval = .*/stats_write_interval = $5/" \
        -e "s/^field_interval = .*/field_interval = $6/" "$INI" > "$1"
    if [ "$RANKS" -gt 1 ]; then printf '\n[mpi]\ndims = %d 1 1\n' "$RANKS" >> "$1"; fi
}

echo "== LEG 1 transient t=0..$T_TRANSIENT (stats off), arch=$ARCH ranks=$RANKS =="
gen_ini "$RUN/transient/input.ini" "$IC" "$T_TRANSIENT" -1 -1 "$FIELD1"
( cd "$RUN/transient" && MOBY_STEPDIV=1 mpirun -n "$RANKS" "$BIN" input.ini > run.log 2>&1 )
LAST=$(ls -t "$RUN"/transient/channel_field_*.h5 | head -1)
echo "   transient final field: $LAST"

echo "== LEG 2 stats t=$T_TRANSIENT..$T_END (stats on, fresh) =="
rm -f "$RUN"/stats/channel_stats.h5 "$RUN"/stats/channel_stats_l*.h5
gen_ini "$RUN/stats/input.ini" "$LAST" "$T_END" 50 20000 0
( cd "$RUN/stats" && MOBY_STEPDIV=1 mpirun -n "$RANKS" "$BIN" input.ini > run.log 2>&1 )
echo "== DONE. uniform-128 stats: $RUN/stats/channel_stats.h5 =="
