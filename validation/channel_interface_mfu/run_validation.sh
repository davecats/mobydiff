#!/usr/bin/env bash
# Minimal-flow-unit (MFU) variant of the 2:1-interface channel validation.
#
# Same Re_tau = 180 setup as ../channel_interface, but in a small 2x2x2 box at
# (close to) the same near-wall resolution: base 24x64x40, fine level 48x128x80.
# The y-line is identical to the full case, so the interfaces sit at the same
# y+ = 112 / 55. ICs are interpolated from the full channel restart by sampling
# its 2x2 corner window; the transient leg lets the MFU re-equilibrate.
#
# ~17x fewer cells than the full case -> fast turnaround. Note: the MFU has far
# fewer homogeneous samples per plane, so converged turbulence statistics need
# longer time-averaging than the full case; this case is meant for quick
# stability / qualitative checks and fast iteration, not converged statistics.
#
#   ./run_validation.sh [gpu|cpu] [nranks]
#
# Post-process with:
#   python3 ../../tools/channel_interface_validation.py \
#       --reference runs/reference/stats --refined runs/refined_y110/stats --out plots_y110
# Visualise a cross-section of any field with:
#   python3 ../../tools/plot_field_section.py runs/refined_y110/stats/channel_field_XXXX.h5

set -euo pipefail
cd "$(dirname "$0")"

ARCH="${1:-gpu}"
NRANKS="${2:-1}"
BIN="../../build_${ARCH}/main"
TRANSIENT_T=5.0
AVERAGE_T=20.0
# Minimal flow unit box and base resolution.
BOX="--lx 2 --lz 2 --nx 24 --nz 40"

[ -x "$BIN" ] && BIN="$(cd "$(dirname "$BIN")" && pwd)/main" || {
    echo "binary $BIN not found (build with ./compile.sh $ARCH)"; exit 1; }

echo "== quick interface-decay gate"
( cd ../../tutorials/interface_decay
  rm -f decay_*.h5
  mpirun -n 1 "$BIN" input.ini > run.log 2>&1 )
python3 ../../tools/check_interface_decay.py ../../tutorials/interface_decay

echo "== generating MFU initial conditions from channel_kmm180_restart.h5"
mkdir -p ic
python3 ../../tools/make_channel_restart.py --mode reference $BOX \
    --source ../../tutorials/channel_kmm180/channel_kmm180_restart.h5 \
    --out ic/reference_ic.h5
python3 ../../tools/make_channel_restart.py --mode refined --band-cells 24 $BOX \
    --source ../../tutorials/channel_kmm180/channel_kmm180_restart.h5 \
    --out ic/refined_y110_ic.h5
python3 ../../tools/make_channel_restart.py --mode refined --band-cells 16 $BOX \
    --source ../../tutorials/channel_kmm180/channel_kmm180_restart.h5 \
    --out ic/refined_y55_ic.h5

run_case () {
    local name="$1" ic="$2"
    local dirA="runs/$name/transient" dirB="runs/$name/stats"
    echo "== $name: transient leg (t = 0 .. $TRANSIENT_T)"
    mkdir -p "$dirA" "$dirB"
    sed -e "s|^t_final = .*|t_final = $TRANSIENT_T|" \
        -e "s|^file = RESTART_PLACEHOLDER|file = ../../../ic/$ic|" \
        -e "s|^stats_sample_interval = .*|stats_sample_interval = 0|" \
        -e "s|^stats_write_interval = .*|stats_write_interval = 0|" \
        -e "s|^field_interval = .*|field_interval = 0|" \
        "$name.ini" > "$dirA/input.ini"
    ( cd "$dirA" && mpirun -n "$NRANKS" "$BIN" input.ini > run.log 2>&1 )
    local final
    final=$(ls -1 "$dirA"/channel_field_*.h5 | sort -t_ -k3 -n | tail -1)
    echo "   transient final field: $final"

    echo "== $name: statistics leg (t = $TRANSIENT_T .. $(python3 -c "print($TRANSIENT_T+$AVERAGE_T)"))"
    sed -e "s|^t_final = .*|t_final = $(python3 -c "print($TRANSIENT_T+$AVERAGE_T)")|" \
        -e "s|^file = RESTART_PLACEHOLDER|file = ../../../$final|" \
        "$name.ini" > "$dirB/input.ini"
    ( cd "$dirB" && mpirun -n "$NRANKS" "$BIN" input.ini > run.log 2>&1 )
    echo "   done: stats in $dirB"
}

run_case reference reference_ic.h5
run_case refined_y110 refined_y110_ic.h5
run_case refined_y55 refined_y55_ic.h5

echo "== all runs complete; post-process with tools/channel_interface_validation.py"
