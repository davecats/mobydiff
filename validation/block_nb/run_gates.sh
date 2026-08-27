#!/usr/bin/env bash
# Per-direction [blocks] nb gates (Phase 1).
#
# The designed property being tested is nb-INDEPENDENCE: the block decomposition
# is bookkeeping, so tiling the same grid differently must give bit-identical
# fields. That makes non-cubic nb self-gating -- no reference implementation is
# needed, the cubic run IS the reference.
#
#   ./run_gates.sh [cpu|gpu]
#
# Env: BIN (solver), PY (python with h5py), RANKS (ranks for the CPU gates).
set -uo pipefail

MODE="${1:-cpu}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
BIN="${BIN:-$ROOT/build_${MODE}/moby_solve}"
PY="${PY:-/home/ws/xt8786/ibmc/bin/python}"
[ -x "$PY" ] || PY=python3
CMP="$ROOT/tools/compare_fields.py"
RANKS="${RANKS:-1}"
[ "$MODE" = gpu ] && RANKS=1

cd "$HERE"
fail=0

# run <tag> <body.ini> <ranks> <blocks-section...>
run() {
    local tag="$1" body="$2" ranks="$3"; shift 3
    local ini=".gate_${tag}.ini"
    { cat "$body"; printf '\n[blocks]\n'; printf '%s\n' "$@"; } \
        | sed "s/^field_prefix = gate/field_prefix = gate_${tag}/" > "$ini"
    if ! mpirun -n "$ranks" "$BIN" "$ini" > ".gate_${tag}.log" 2>&1; then
        echo "RUN FAILED $tag"; tail -4 ".gate_${tag}.log"; fail=1; return 1
    fi
    grep -h "^ block refinement:" ".gate_${tag}.log" 2>/dev/null | tail -1
}

same() {
    local what="$1" a="gate_$2_20.h5" b="gate_$3_20.h5"
    if [ ! -f "$a" ] || [ ! -f "$b" ]; then echo "FAIL $what (missing snapshot)"; fail=1; return; fi
    local ds
    ds=$("$PY" -c "
import h5py
f=set(h5py.File('$a','r')); g=set(h5py.File('$b','r'))
print(' '.join(sorted((f&g)-{'x','y','z','blocks'})))")
    if "$PY" "$CMP" "$a" "$b" $ds --tolerance 0 > ".cmp_$2_$3.txt" 2>&1; then
        echo "PASS $what"
    else
        echo "FAIL $what"; head -6 ".cmp_$2_$3.txt"; fail=1
    fi
}

echo "== single level: the block layout must not change the answer ($MODE)"
run unset  base.ini "$RANKS" "; nb unset -> one block per rank box"
run cubic  base.ini "$RANKS" "nb = 8"
run rect1  base.ini "$RANKS" "nb = 32 16 8"
run rect2  base.ini "$RANKS" "nb = 16 8 4"
same "non-cubic 32 16 8 == cubic 8"      rect1 cubic
same "non-cubic 16 8 4  == cubic 8"      rect2 cubic
same "cubic 8           == nb unset"     cubic unset

if [ "$MODE" = cpu ]; then
    echo
    echo "== rank independence: non-cubic on 1 rank == 4 ranks"
    run rect1r4 base.ini 4 "nb = 32 16 8"
    same "non-cubic 1 rank == 4 ranks" rect1r4 rect1
fi

# The 2:1 interface is NOT nb-independent (README: measured, and it predates
# per-direction nb), so a refined layout cannot be gated against a differently
# tiled one. The interface gate is uniform-flow preservation instead: a constant
# field survives any CONSISTENT set of transfer operators exactly, so a
# per-direction indexing mistake anywhere in the entry generation, the gather
# maps or the ghost blend shows up as a nonzero.
echo
echo "== 2:1 interface with non-cubic nb: uniform oblique flow must be EXACT ($MODE)"
for f in uniform_rect uniform_rect_xyz; do
    pfx=$(grep '^field_prefix' "$f.ini" | awk '{print $3}')
    if mpirun -n "$RANKS" "$BIN" "$f.ini" > ".$f.log" 2>&1; then
        grep -h "block refinement" ".$f.log" | tail -1
        "$PY" - "$pfx"_50.h5 "$f" <<'EOF' || fail=1
import sys, h5py, numpy as np
h = h5py.File(sys.argv[1], "r")
exp = (0.9396926207859084, 0.3420201433256687, 0.2)
dev = max(float(np.max(np.abs(h[v][...] - e))) for v, e in zip(("un","vn","wn"), exp))
p = h["pn"][...]
spread = float(p.max() - p.min())
lv = np.asarray(h["blocks"])[:, 3]
levels = {int(l): int((lv == l).sum()) for l in np.unique(lv)}
ok = (dev == 0.0 and spread == 0.0)
# a patch that swallowed level 0 would only carry l1-l2 interfaces
ok = ok and len(levels) >= 3
print(f"{'PASS' if ok else 'FAIL'} {sys.argv[2]}: max dev {dev:.3e}, pn spread "
      f"{spread:.3e}, levels {levels}")
sys.exit(0 if ok else 1)
EOF
    else
        echo "RUN FAILED $f"; tail -4 ".$f.log"; fail=1
    fi
done

echo
[ "$fail" = 0 ] && echo "ALL PASS ($MODE)" || echo "FAILURES ($MODE)"
exit $fail
