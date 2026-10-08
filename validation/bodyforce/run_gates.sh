#!/usr/bin/env bash
# Body-force gates (2026-10-06): the correction kernel visits the force's
# SUPPORT only (bodyforce_type sup_blocks / sup_var_lo:hi; the trip = v on its
# block list), which the profile of job 5184720 showed is 8.5 % of the
# boundary-layer step when done over the whole domain.
#
#   REF=/abs/ref/build_cpu_nofma/moby_solve NEW=/abs/build_cpu_nofma/moby_solve \
#       [RANKS=1] ./run_gates.sh
#
# 1. trip_small.ini (trip ON, 7 of 32 blocks) REF vs NEW: un/vn/wn/pn max_abs 0.
# 2. its twin with trip_amp = 0 vs NEW: vn MUST differ (the gate is not vacuous).
# 3. force_sine.ini (generic [force] profile on three components) REF vs NEW:
#    max_abs 0 -- the whole-domain support of every non-trip source.
# Both binaries nofma; each side prepares its own case file.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
REF=${REF:?set REF to the reference nofma binary (absolute path)}
NEW=${NEW:-$ROOT/build_cpu_nofma/moby_solve}
RANKS=${RANKS:-1}
CMP="python3 $ROOT/tools/compare_fields.py --tolerance 0"
status=0

run() { # <label> <bin> <ini>
    local lbl="$1" bin="$2" ini="$3" pfx
    pfx=$(grep '^field_prefix' "$ini" | awk '{print $3}')
    mkdir -p "run_$lbl"; cp "$ini" "run_$lbl/"
    ( cd "run_$lbl" && rm -f "${pfx}"_*.h5 "${pfx}.case.h5" \
        && mpirun -n "$RANKS" "$(dirname "$bin")/moby_prepare" "$ini" > prepare.log 2>&1 \
        && mpirun -n "$RANKS" "$bin" "$ini" > run.log 2>&1 ) || { echo "   $lbl FAILED (see run_$lbl/run.log)"; return 1; }
}
last() { ls "run_$1"/$2_*.h5 | sort -V | tail -1; }

echo "== trip ON: REF vs NEW (nofma, $RANKS ranks)"
run trip_ref "$REF" trip_small.ini && run trip_new "$NEW" trip_small.ini \
    && $CMP "$(last trip_ref trip)" "$(last trip_new trip)" un vn wn pn || status=1
echo "== trip OFF twin vs NEW: vn must differ"
sed 's/^trip_amp = 0.03/trip_amp = 0.0/' trip_small.ini > trip_off.ini
run trip_off "$NEW" trip_off.ini || status=1
if $CMP "$(last trip_off trip)" "$(last trip_new trip)" vn > /dev/null 2>&1; then
    echo "   FAIL: the trip did not move vn (vacuous gate)"; status=1
else
    echo "   PASS: the trip moves vn"
fi
echo "== generic [force] sine profile: REF vs NEW"
run force_ref "$REF" force_sine.ini && run force_new "$NEW" force_sine.ini \
    && $CMP "$(last force_ref force)" "$(last force_new force)" un vn wn pn || status=1
[ $status -eq 0 ] && echo "bodyforce gates: ALL PASS" || echo "bodyforce gates: FAILURES"
exit $status
