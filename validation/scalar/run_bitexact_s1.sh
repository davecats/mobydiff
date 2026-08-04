#!/usr/bin/env bash
# S2 companion to run_bitexact.sh: a scalar run with turbulence OFF must be
# bit-exact against the S1 binaries. run_bitexact.sh proves `[scalar] count =
# 0` (no scalar allocated at all); this proves the OTHER half -- the transport
# kernel's arithmetic is unchanged when no eddy viscosity exists, which is a
# by-construction property of the S2 kernel (every face diffusivity keeps the
# molecular value exactly; the eddy term sits behind `if (useNut)`), not a
# cancellation.
#
# BOTH binaries must be nofma builds (compile_nofma.sh).
#
#   REF=/path/to/s1_cpu_nofma NEW=../../build_cpu_nofma/moby_solve ./run_bitexact_s1.sh
#
# Environment: REF (required), NEW, RANKS (default 1), MODE (cpu|gpu, labels)
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

REF=${REF:?set REF to the reference (S1) nofma binary}
MODE=${MODE:-cpu}
NEW=${NEW:-$ROOT/build_${MODE}_nofma/moby_solve}
RANKS=${RANKS:-1}
status=0
tag="sx_${MODE}"

# ini | steps | ranks (0 = $RANKS) | datasets
CASES=(
  "uniform3.ini|20|0|un vn wn pn theta phi"
  "det.ini|20|0|un vn wn pn s1"
  "conduction.ini|50|0|un vn wn pn s1"
  "prsweep.ini|50|0|un vn wn pn s1"
  "conserve.ini|20|0|un vn wn pn s1"
)

for entry in "${CASES[@]}"; do
    IFS='|' read -r ini steps cranks datasets <<< "$entry"
    name=${ini%.ini}
    [ "$cranks" = 0 ] && cranks=$RANKS
    echo "== $name ($MODE, $cranks rank(s), $steps steps)"
    for side in ref new; do
        bin=$REF; [ "$side" = new ] && bin=$NEW
        pfx="${tag}_${name}_${side}"
        sed -e "s/^nsteps.*/nsteps = $steps/" -e "s/^t_final.*/t_final = 0.0/" \
            -e "s/^field_interval.*/field_interval = $steps/" \
            -e "s/^field_prefix.*/field_prefix = $pfx/" "$ini" > ".${pfx}.ini"
        rm -f "${pfx}_"*.h5
        if ! mpirun -n "$cranks" "$bin" ".${pfx}.ini" > ".${pfx}.log" 2>&1; then
            echo "   RUN FAILED ($side) -- see .${pfx}.log"; status=1; continue 2
        fi
        rm -f ".${pfx}.ini"
    done
    refh5=$(ls -t "${tag}_${name}_ref_"*.h5 2>/dev/null | head -1)
    newh5=$(ls -t "${tag}_${name}_new_"*.h5 2>/dev/null | head -1)
    if python3 "$ROOT/tools/compare_fields.py" "$refh5" "$newh5" $datasets \
            --tolerance 0 | sed 's/^/   /'; then
        echo "   PASS (max_abs 0)"
    else
        echo "   FAIL"; status=1
    fi
done

echo
[ $status -eq 0 ] && echo "scalar-with-turbulence-off suite ($MODE): ALL PASS" \
                  || echo "scalar-with-turbulence-off suite ($MODE): FAILURES"
exit $status
