#!/usr/bin/env bash
# Bit-exactness on the OUTLET cases: every case below has a declared outlet
# face, which the 7-case and 9-case suites (validation/scalar/run_bitexact*.sh)
# do not. Runs each with a REFERENCE binary and the NEW binary and compares
# un/vn/wn/pn with tools/compare_fields.py at tolerance 0.
#
# BOTH binaries must be nofma builds (./compile.sh cpu_nofma | gpu_nofma).
#
#   REF=/abs/ref/build_cpu_nofma/moby_solve NEW=/abs/build_cpu_nofma/moby_solve \
#       [RANKS=1] [MODE=cpu] ./run_bitexact_outlet.sh [case]
#
# Paths must be ABSOLUTE (the driver changes directory per case). The output
# prefix is keyed on MODE: give concurrent runs their own label.
#
# A REF binary older than the outlet planes (un_xmax ..., 2026-10-01) does
# not read them: the cases that restart from an IC carrying a plane
# (pois_io, lamboseen, blasius2d) then differ by construction. And never run
# this next to ./run_gates.sh: that driver regenerates IC_pois.h5 and
# IC_vortex.h5 under the two sides.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

REF=${REF:?set REF to the reference nofma binary (absolute path)}
MODE=${MODE:-cpu}
NEW=${NEW:-$ROOT/build_${MODE}_nofma/moby_solve}
RANKS=${RANKS:-1}
sel=${1:-all}
status=0
tag="obx_${MODE}"

# case | dir | ini | steps | ranks (0 = $RANKS)
CASES=(
  "oblique|$ROOT/validation/freestream|oblique.ini|20|0"
  "pois_io|$ROOT/validation/freestream|pois_io.ini|50|0"
  "lamboseen|$ROOT/validation/freestream|lamboseen.ini|100|0"
  "outlet_high|$ROOT/validation/freestream|outlet_high.ini|100|0"
  "outlet_low|$ROOT/validation/freestream|outlet_low.ini|100|0"
  "blasius2d|$ROOT/validation/blasius|blasius2d.ini|50|1"
  "cyl_re100|$ROOT/validation/cylinder|cyl_re100.ini|20|1"
)

# A short fixed-step variant with its own output prefix. A pinned
# `[case] file` is dropped, so each side prepares and reads its OWN case file
# (<prefix>.case.h5): a stretched node line differs by ulps between builds
# and the solver refuses a case file whose lines are not its own.
short_ini() {  # in out steps prefix
    awk '/^[ \t]*\[/ {sec = $0; gsub(/[ \t]/, "", sec)}
         sec == "[case]" && /^[ \t]*file[ \t]*=/ {next}
         {print}' "$1" |
    sed -e "s/^nsteps.*/nsteps = $3/" \
        -e "s/^t_final.*/t_final = 0.0/" \
        -e "s/^field_interval.*/field_interval = 0/" \
        -e "s/^field_prefix.*/field_prefix = $4/" \
        -e "s/^runtime_file.*/runtime_file = $4.forces.txt/" > "$2"
}

for entry in "${CASES[@]}"; do
    IFS='|' read -r name dir ini steps cranks <<< "$entry"
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
          rm -f "${pfx}_"*.h5 "${pfx}.case.h5"
          "$ROOT/tools/prepare_if_missing.sh" "$cranks" "$bin" ".${pfx}.ini" ".${pfx}.prep.log" \
              || { echo "   PREPARE FAILED ($side) -- see $dir/.${pfx}.prep.log"; exit 2; }
          if ! mpirun -n "$cranks" "$bin" ".${pfx}.ini" > ".${pfx}.log" 2>&1; then
              echo "   RUN FAILED ($side) -- see $dir/.${pfx}.log"; exit 2
          fi
      done ) || { status=1; continue; }
    refh5=$(ls -t "$dir/${tag}_${name}_ref_"*.h5 2>/dev/null | head -1)
    newh5=$(ls -t "$dir/${tag}_${name}_new_"*.h5 2>/dev/null | head -1)
    if [ -z "$refh5" ] || [ -z "$newh5" ]; then
        echo "   NO OUTPUT"; status=1; continue
    fi
    if python3 "$ROOT/tools/compare_fields.py" "$refh5" "$newh5" un vn wn pn \
            --tolerance 0 | sed 's/^/   /'; then
        echo "   PASS (max_abs 0)"
    else
        echo "   FAIL"; status=1
    fi
done

echo
[ $status -eq 0 ] && echo "outlet bit-exactness ($MODE): ALL PASS" \
                  || echo "outlet bit-exactness ($MODE): FAILURES"
exit $status
