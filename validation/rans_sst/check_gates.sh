#!/usr/bin/env bash
# RANS T2 gate checker — run AFTER run_gates.sh (any machine with
# python3 + h5py + numpy; PY overrides the interpreter).
set -uo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
status=0

latest() { ls -v "$1"_*.h5 2>/dev/null | tail -1; }

check() {
    local name=$1; shift
    local f
    f=$(latest "$name")
    if [ -z "$f" ]; then
        echo "== $name: NO FIELD FILE — skipped"; status=1; return
    fi
    echo "== $name ($f)"
    if ! $PY rans_channel_check.py "$f" "$@"; then status=1; fi
    echo
}

check laminar  --mode laminar
# The pure kappa/B log line deviates from real profiles by several % over
# the overlap region; the DNS centreline U+ anchors (18.20 / 20.13, Moser
# et al.) are the sharper RANS-quality measure.
check turb180  --mode loglaw --tolerance 0.06 --uplus-center 18.20
check turb395  --mode loglaw --tolerance 0.08 --uplus-center 20.13
check ibm180   --mode loglaw --wall-lo 0.259375 --wall-hi 2.259375 --tolerance 0.10

# base180u (uniform y+_1 ~ 2.8) is INFORMATIONAL only: its under-resolved
# near-wall region feeds a spurious core-k plateau, which is exactly why
# the refined case exists. Not a pass/fail gate.
f=$(latest base180u)
if [ -n "$f" ]; then
    echo "== base180u ($f) — informational (coarse-wall control)"
    $PY rans_channel_check.py "$f" --mode loglaw --tolerance 0.10 || true
    echo
fi

# Gate (d): no interface band, and the refined core must agree with the
# RESOLVED single-level reference (turb180's natural grid), not with the
# coarse control.
REF=$(latest turb180)
check refine180 --mode band ${REF:+--reference "$REF"}

echo "overall: $([ $status -eq 0 ] && echo ALL PASS || echo FAILURES PRESENT)"
exit $status
