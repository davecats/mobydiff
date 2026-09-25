#!/usr/bin/env bash
# Passive-scalar increments S0/S1: the "[scalar] count = 0 is bit-exact BY
# CONSTRUCTION" gate. Runs the standard 7-case suite (min_channel, les_ibm
# +- refine_body, Beltrami y-slab, turb180, wf180_y30, lam30t) with a
# REFERENCE binary and the NEW binary and compares every field dataset with
# tools/compare_fields.py at tolerance 0.
#
# BOTH binaries must be nofma builds (validation/scalar/compile_nofma.sh):
# default FMA contraction introduces 1-2 ulp differences for arithmetically
# identical source.
#
#   REF=/path/to/base_cpu_nofma NEW=../../build_cpu_nofma/moby_solve \
#       RANKS=4 ./run_bitexact.sh
#
# Environment: REF (required), NEW (default ../../build_<mode>_nofma/moby_solve),
#              RANKS (default 1), STEPS (default 20), MODE (cpu|gpu, labels only)
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

REF=${REF:?set REF to the reference (pre-change) nofma binary}
MODE=${MODE:-cpu}
NEW=${NEW:-$ROOT/build_${MODE}_nofma/moby_solve}
RANKS=${RANKS:-1}
STEPS=${STEPS:-20}
sel=${1:-all}
status=0
tag="bx_${MODE}"

# case | dir | ini | steps | ranks (0 = $RANKS) | datasets
# The les_ibm inis pin [mpi] dims = 1 1 1, so they run on one rank.
CASES=(
  "min_channel|$ROOT/tutorials/min_channel|input.ini|$STEPS|0|un vn wn pn"
  "les_ibm|$ROOT/validation/channel_interface/les_ibm|channel_ibm.ini|$STEPS|1|un vn wn pn nut"
  "les_ibm_refine|$ROOT/validation/channel_interface/les_ibm|channel_ibm_refine.ini|$STEPS|1|un vn wn pn nut"
  "beltrami_yslab|$ROOT/validation/beltrami|slab_y.ini|5|0|un vn wn pn"
  "turb180|$ROOT/validation/rans_sst|turb180.ini|$STEPS|0|un vn wn pn nut k omega"
  "wf180_y30|$ROOT/validation/rans_sst|wf180_y30.ini|$STEPS|0|un vn wn pn nut k omega"
  "lam30t|$ROOT/validation/rans_sst|lam30t.ini|$STEPS|0|un vn wn pn nut k omega gamma rethetat"
)

# Rewrite a case ini for a short fixed-step run with its own output prefix.
short_ini() {  # in out steps prefix
    sed -e "s/^nsteps.*/nsteps = $3/" \
        -e "s/^t_final.*/t_final = 0.0/" \
        -e "s/^field_interval.*/field_interval = $3/" \
        -e "s/^field_prefix.*/field_prefix = $4/" \
        -e "s/^stats_sample_interval.*/stats_sample_interval = 0/" \
        -e "s/^stats_write_interval.*/stats_write_interval = 0/" "$1" > "$2"
}

for entry in "${CASES[@]}"; do
    IFS='|' read -r name dir ini steps cranks datasets <<< "$entry"
    [ "$sel" = all ] || [ "$sel" = "$name" ] || continue
    [ "$cranks" = 0 ] && cranks=$RANKS
    if [ ! -f "$dir/$ini" ]; then
        echo "== $name SKIPPED (missing $dir/$ini)"; status=1; continue
    fi
    echo "== $name ($MODE, $cranks rank(s), $steps steps)"
    ( cd "$dir" || exit 1
      for side in ref new; do
          bin=$REF; [ "$side" = new ] && bin=$NEW
          pfx="${tag}_${name}_${side}"
          short_ini "$ini" ".${pfx}.ini" "$steps" "$pfx"
          rm -f "${pfx}_"*.h5
          if ! mpirun -n "$cranks" "$bin" ".${pfx}.ini" > ".${pfx}.log" 2>&1; then
              echo "   RUN FAILED ($side) -- see $dir/.${pfx}.log"; exit 2
          fi
      done ) || { status=1; continue; }
    refh5=$(ls -t "$dir/${tag}_${name}_ref_"*.h5 2>/dev/null | head -1)
    newh5=$(ls -t "$dir/${tag}_${name}_new_"*.h5 2>/dev/null | head -1)
    if [ -z "$refh5" ] || [ -z "$newh5" ]; then
        echo "   NO OUTPUT"; status=1; continue
    fi
    if python3 "$ROOT/tools/compare_fields.py" "$refh5" "$newh5" $datasets \
            --tolerance 0 | sed 's/^/   /'; then
        echo "   PASS (max_abs 0)"
    else
        echo "   FAIL"; status=1
    fi
done

echo
[ $status -eq 0 ] && echo "bit-exactness suite ($MODE): ALL PASS" \
                  || echo "bit-exactness suite ($MODE): FAILURES"
exit $status
