#!/usr/bin/env bash
# Red-black SOR across a 2:1 interface (R1, docs/next_session_redblack_interface.md).
#
#   ./run_gates.sh [cpu|gpu]
#
# Env: BIN (solver), PY (python with h5py), RANKS (ranks for the CPU gates).
#
# READ THIS BEFORE ADDING A GATE HERE. The uniform-flow-through-a-patch test,
# which the plan called the decisive transfer gate, is BLIND to the interface
# patch: uniform flow makes the divergence and therefore phi exactly zero, so
# interface_correct adds nothing and the test passes with the patch disabled
# (measured: max deviation 0.000e+00 either way). It still gates the halo
# TRANSFER operators, which is what it was written for, and it is kept for that.
# The gate with power over the patch is the refined CHANNEL below: with
# interface_correct disabled it blows up to |u| ~ 4e4 / |p| ~ 2e11 within 20
# steps, against a stable run with it enabled.
set -uo pipefail

MODE="${1:-cpu}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
BIN="${BIN:-$ROOT/build_${MODE}/moby_solve}"
PY="${PY:-/home/ws/xt8786/ibmc/bin/python}"
[ -x "$PY" ] || PY=python3
CMP="$ROOT/tools/compare_fields.py"
RANKS="${RANKS:-4}"
[ "$MODE" = gpu ] && RANKS=1

cd "$HERE"
fail=0

run() {  # run <tag> <ini> <ranks>
    local tag="$1" ini="$2" ranks="$3"
    sed "s/^field_prefix = .*/field_prefix = gate_${tag}/" "$ini" > ".gate_${tag}.ini"
    if ! mpirun -n "$ranks" "$BIN" ".gate_${tag}.ini" > ".gate_${tag}.log" 2>&1; then
        echo "RUN FAILED $tag"; tail -4 ".gate_${tag}.log"; fail=1; return 1
    fi
    grep -h "^ block refinement:" ".gate_${tag}.log" | tail -1
}

echo "== red-black + 2:1 runs, and stays finite ($MODE)"
run rb1 refined_channel.ini 1 || true
if [ -f gate_rb1_20.h5 ]; then
    "$PY" - gate_rb1_20.h5 <<'EOF' || fail=1
import sys, h5py, numpy as np
h = h5py.File(sys.argv[1], "r")
bad = {v: float(np.abs(np.asarray(h[v])).max()) for v in ("un","vn","wn","pn")}
ok = all(np.isfinite(m) for m in bad.values()) and bad["un"] < 1e3
print(("PASS" if ok else "FAIL") + " finite: " +
      "  ".join(f"max|{v}|={m:.3e}" for v, m in bad.items()))
sys.exit(0 if ok else 1)
EOF
else
    fail=1
fi

if [ "$MODE" = cpu ]; then
    echo
    echo "== rank independence: red-black + 2:1 on 1 rank == $RANKS ranks"
    run rbN refined_channel.ini "$RANKS" || true
    if [ -f gate_rb1_20.h5 ] && [ -f gate_rbN_20.h5 ]; then
        if "$PY" "$CMP" gate_rb1_20.h5 gate_rbN_20.h5 un vn wn pn --tolerance 0 > .cmp_rank.txt 2>&1
        then echo "PASS 1 rank == $RANKS ranks (EXACT)"
        else echo "FAIL 1 rank == $RANKS ranks"; head -6 .cmp_rank.txt; fail=1; fi
    else
        fail=1
    fi
fi

echo
echo "== halo transfer: uniform oblique flow through a 3-level patch, EXACT ($MODE)"
echo "   (gates the TRANSFER operators only -- see the header: blind to the patch)"
sed -e 's/^accel = chebyshev/solver = redblack/' -e 's/^sor = 0.8/sor = 1.5/' \
    -e 's/^field_prefix = .*/field_prefix = gate_uni/' \
    ../block_nb/uniform_rect.ini > .gate_uni.ini
if mpirun -n 1 "$BIN" .gate_uni.ini > .gate_uni.log 2>&1; then
    grep -h "^ block refinement:" .gate_uni.log | tail -1
    "$PY" - gate_uni_50.h5 <<'EOF' || fail=1
import sys, h5py, numpy as np
h = h5py.File(sys.argv[1], "r")
exp = (0.9396926207859084, 0.3420201433256687, 0.2)
dev = max(float(np.max(np.abs(h[v][...] - e))) for v, e in zip(("un","vn","wn"), exp))
spread = float(np.ptp(h["pn"][...]))
lv = np.asarray(h["blocks"])[:, 3]
levels = {int(l): int((lv == l).sum()) for l in np.unique(lv)}
ok = dev == 0.0 and spread == 0.0 and len(levels) >= 3
print(f"{'PASS' if ok else 'FAIL'} uniform: max dev {dev:.3e}, pn spread "
      f"{spread:.3e}, levels {levels}")
sys.exit(0 if ok else 1)
EOF
else
    echo "RUN FAILED uniform"; tail -4 .gate_uni.log; fail=1
fi

echo
[ "$fail" = 0 ] && echo "ALL PASS ($MODE)" || echo "FAILURES ($MODE)"
exit $fail
