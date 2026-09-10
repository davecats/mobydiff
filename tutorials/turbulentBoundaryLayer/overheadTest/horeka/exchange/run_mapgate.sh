#!/bin/bash
# The bit-exactness gate that `select_target_device` (commit 059248e) still owed.
#
#   run_mapgate.sh <new_exe> <ref_exe> <results_dir>
#
# The debt was recorded as "the WORKSTATION gate (nofma, 7-case suite)". The
# workstation is not reachable from HoreKa, so the suite is run here instead,
# against the pre-change binary, on the SAME node in the SAME allocation.
#
# NOFMA IS NOT USED, DELIBERATELY. -Mnofma exists to stop the compiler inventing
# 1-2 ulp differences between two textually different sources that are
# arithmetically identical. Here the arithmetic source is BYTE-IDENTICAL between
# the two binaries -- only the integer that picks a physical GPU changes -- so
# both sides contract identically and a production-flag comparison is a strictly
# tighter test than a nofma one. Any nonzero max_abs is a real difference.
#
# Every case is a single-node run, which is exactly where the rule reduces to the
# identity map; the multi-node fields were already gated max_abs 0 at 8x2 in
# results_horeka_2026-09-10.md section 7.
set -uo pipefail

EXE="${1:?usage: run_mapgate.sh <new_exe> <ref_exe> <results_dir>}"
REF="${2:?usage: run_mapgate.sh <new_exe> <ref_exe> <results_dir>}"
RES="${3:?usage: run_mapgate.sh <new_exe> <ref_exe> <results_dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../../../.." && pwd)"
H5MAXDIFF="${H5MAXDIFF:-$ROOT/tools/h5maxdiff}"
NSTEPS="${NSTEPS:-200}"

# name : ini (relative to the repo root) : ranks. The standard 7-case list.
CASES="${CASES:-
min_channel:tutorials/min_channel/input_gpu.ini:4
min_channel1:tutorials/min_channel/input_gpu.ini:1
beltrami_slaby:validation/beltrami/slab_y.ini:1
les_ibm:validation/channel_interface/les_ibm/channel_ibm.ini:1
les_ibm_refine:validation/channel_interface/les_ibm/channel_ibm_refine.ini:1
turb180:validation/rans_sst/turb180.ini:1
wf180_y30:validation/rans_sst/wf180_y30.ini:1
lam30t:validation/rans_sst/lam30t.ini:1
}"

mkdir -p "$RES"
echo "=== map gate start $(date '+%F %T')  nsteps=$NSTEPS"
for entry in $CASES; do
    name="${entry%%:*}"; rest="${entry#*:}"
    ini="${rest%%:*}"; ranks="${rest##*:}"
    src="$ROOT/$ini"
    [ -f "$src" ] || { echo "MISSING CASE: $ini"; continue; }
    for side in ref new; do
        run="$RES/${name}_${side}"
        [ -f "$run/run.log" ] && continue
        rm -rf "$run"
        # The inis name their coefficient/IC/grid files by RELATIVE path, so the
        # run has to happen inside a copy of the case directory -- and a COPY,
        # not the repo directory itself, so the gate never leaves snapshots or
        # runtime files in the working tree.
        cp -a "$(dirname "$src")" "$run"
        cp "$src" "$run/config.ini"
        sed -i -e "s/^nsteps *=.*/nsteps = $NSTEPS/" -e "s/^t_final *=.*/t_final = 0.0/" \
               "$run/config.ini"
        ( cd "$run" && mpirun -n "$ranks" --bind-to core --map-by numa \
              "$( [ "$side" = ref ] && echo "$REF" || echo "$EXE" )" \
              config.ini > run.log 2>&1 )
        rc=$?
        [ $rc -ne 0 ] && { echo "  $name $side FAILED (exit $rc) -- see $run/run.log"; continue; }
    done
    pfx=$(awk -F= '/^field_prefix/ {gsub(/ /,"",$2); print $2}' "$RES/${name}_new/config.ini" 2>/dev/null)
    a=$(ls "$RES/${name}_ref/${pfx}"_*.h5 2>/dev/null | tail -1)
    b=$(ls "$RES/${name}_new/${pfx}"_*.h5 2>/dev/null | tail -1)
    if [ -n "$a" ] && [ -n "$b" ]; then
        echo "--- $name  $(basename "$a")"
        "$H5MAXDIFF" "$a" "$b" | tee "$RES/${name}.txt" | sed 's/^/    /'
    else
        echo "--- $name: NO SNAPSHOT PAIR" | tee "$RES/${name}.txt"
    fi
done
echo "=== map gate done $(date '+%F %T')"
