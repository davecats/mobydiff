#!/usr/bin/env bash
# A3 INCREMENT 1 gate (b): RANS scalar inlet values (docs/next_session_airfoil.md).
#
# Environment:  BIN=<solver>   (default ../../build_cpu/main)
set -o pipefail
cd "$(dirname "$0")"
set -u

BIN=${BIN:-../../build_cpu/main}
CMP="python3 ../../tools/compare_fields.py"
status=0

run() { # ini ranks log
    local ini=$1 ranks=$2 log=$3
    echo "== $log (ranks $ranks) =="
    if ! mpirun -n "$ranks" "$BIN" "$ini" > "$log.log" 2>&1; then
        echo "   FAILED — see $log.log"; status=1; return 1
    fi
    tail -n 2 "$log.log" | sed 's/^/   /'
}

for r in 1 4; do
    sed -e "s/^field_prefix.*/field_prefix = inl_r${r}/" inlet_channel.ini > .r.ini
    rm -f inl_r${r}_*.h5
    run .r.ini $r inl_r$r
    rm -f .r.ini
done

echo "-- declared patch types (from the 1-rank log):"
grep "RANS domain face" inl_r1.log | sed 's/^/   /'
grep -q "x_min: inlet (declared)" inl_r1.log && \
grep -q "x_max: outlet (declared)" inl_r1.log && \
grep -q "y_min: wall (declared)" inl_r1.log \
    && echo "patch-type report: PASS" || { echo "patch-type report: FAIL"; status=1; }

python3 check_inlet.py "$(ls -t inl_r1_*.h5 | head -1)" \
    --tu 5.0 --nut-ratio 10.0 --re 100.0 || status=1

$CMP "$(ls -t inl_r1_*.h5 | head -1)" "$(ls -t inl_r4_*.h5 | head -1)" \
    un vn wn pn k omega nut --tolerance 0 \
    && echo "inlet 1==4 ranks: EXACT" || { echo "inlet ranks MISMATCH"; status=1; }

echo
[ $status -eq 0 ] && echo "ALL GATES PASS" || echo "GATE FAILURES (see above)"
exit $status
