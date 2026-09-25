#!/usr/bin/env bash
# S3 companion to run_bitexact.sh / run_bitexact_s1.sh: a scalar run with NO
# IMMERSED BODY must be bit-exact against the S2 binaries.
#
# run_bitexact.sh proves `[scalar] count = 0` (no scalar allocated at all);
# run_bitexact_s1.sh proves the eddy term is dormant without a turbulence
# model; this proves the third by-construction property of the S3 kernel --
# without a body the penalization branch is not entered and the six face
# masks collapse to the FACE_CLOSED flags, so the S2 arithmetic survives
# byte for byte (with `useIbm` off, `adiab` is .false. for every scalar and
# the masks reduce to `clw`... exactly).
#
# It doubles as the cheap, sharp form of "every S1/S2 gate still reads the
# same number": the S2 gate cases are turbulent-channel CAMPAIGNS costing
# hours to days (turbsst alone is 186711 steps), but if the S3 binary
# reproduces the S2 binary's fields at TOLERANCE 0 on those very inis, no
# statistic derived from them can have moved.
#
# BOTH binaries must be nofma builds (compile_nofma.sh).
#
#   REF=~/s2_ref/moby_solve_cpu_nofma NEW=../../build_cpu_nofma/moby_solve ./run_bitexact_s3.sh
#
# Environment: REF (required), NEW, RANKS (default 1), MODE (cpu|gpu, labels)
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

REF=${REF:?set REF to the reference (S2) nofma binary}
MODE=${MODE:-cpu}
NEW=${NEW:-$ROOT/build_${MODE}_nofma/moby_solve}
RANKS=${RANKS:-1}
sel=${1:-all}
status=0
tag="s3x_${MODE}"

# ini | steps | ranks (0 = $RANKS) | datasets
# The first five are the S1 gate cases (no turbulence model), the last four
# the S2 ones (LES / SST / the 2:1 band / the determinism twin).
CASES=(
  "uniform3.ini|20|0|un vn wn pn theta phi"
  "det.ini|20|0|un vn wn pn s1"
  "conduction.ini|50|0|un vn wn pn s1"
  "prsweep.ini|50|0|un vn wn pn s1"
  "conserve.ini|20|0|un vn wn pn s1"
  "turbles.ini|20|0|un vn wn pn nut theta"
  "turbslab.ini|20|0|un vn wn pn nut theta"
  "turbsst.ini|20|0|un vn wn pn nut k omega theta theta_kc"
  "detles.ini|20|0|un vn wn pn nut theta"
)

for entry in "${CASES[@]}"; do
    IFS='|' read -r ini steps cranks datasets <<< "$entry"
    name=${ini%.ini}
    [ "$sel" = all ] || [ "$sel" = "$name" ] || continue
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
[ $status -eq 0 ] && echo "scalar-without-a-body suite ($MODE): ALL PASS" \
                  || echo "scalar-without-a-body suite ($MODE): FAILURES"
exit $status
