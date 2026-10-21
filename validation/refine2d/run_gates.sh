#!/usr/bin/env bash
# refine2d (refine_dims = xz quadtree) gates, from the README's R2D-0/R2D-2 record:
#
#   [BIN=/abs/moby_solve] [RANKS=4] ./run_gates.sh [uniform|twin|order|rans|all]
#
#   uniform  uniform_xz (xz patch, 3 levels) and unibody_xz (refine_body in xz with a
#            zero-coefficient body, ibm_ml3_xz_zero.h5): uniform oblique flow EXACT
#   twin     gate (c): all-refined xz box == the doubled-x,z single-level twin at
#            tolerance 0, and allref_xz 1 rank == RANKS ranks
#   order    the Beltrami xz patch pair bp_xz_32 / bp_xz_64: L2 velocity error order
#            >= 2.5 (2.83 recorded 2026-07-15)
#   rans     rans_xz smoke: the turb180 column with xz-refined wall bands, 50 steps
#            finite, 1 rank == RANKS ranks EXACT incl. k omega nut
# Env: BIN (solver), PY (python with h5py), RANKS (default 4).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
BIN="${BIN:-$ROOT/build_cpu/moby_solve}"; PY="${PY:-$HOME/ibmc/bin/python}"; [ -x "$PY" ] || PY=python3
RANKS="${RANKS:-4}"; CMP="$ROOT/tools/compare_fields.py"; PIM="$ROOT/tools/prepare_if_missing.sh"
sel="${1:-all}"; want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
cd "$HERE"; fail=0
run() {  # <tag> <ini> <ranks>  (field_prefix -> <tag>)
    sed "s/^field_prefix = .*/field_prefix = $1/" "$2" > ".$1.ini"
    "$PIM" "$3" "$BIN" ".$1.ini" > /dev/null || { echo "PREPARE FAILED $1"; fail=1; return 1; }
    mpirun -n "$3" "$BIN" ".$1.ini" > ".$1.log" 2>&1 || { echo "RUN FAILED $1"; tail -3 ".$1.log"; fail=1; return 1; }
    grep -h "^ block refinement:" ".$1.log" | tail -1; return 0
}
last() { ls -1 "$1"_[0-9]*.h5 2>/dev/null | sort -t_ -k2 -n | tail -1; }
uniform_check() {  # <field> : the oblique freestream of the inis, exact
    "$PY" - "$1" <<'PYEOF' || fail=1
import sys, h5py, numpy as np
h = h5py.File(sys.argv[1], "r"); exp = (0.9396926207859084, 0.3420201433256687, 0.2)
dev = max(float(np.max(np.abs(h[v][...] - e))) for v, e in zip(("un","vn","wn"), exp)); spread = float(np.ptp(h["pn"][...]))
lv = np.asarray(h["blocks"])[:, 3]; levels = {int(l): int((lv == l).sum()) for l in np.unique(lv)}
ok = dev == 0.0 and spread == 0.0
print(f"{'PASS' if ok else 'FAIL'} uniform: max dev {dev:.3e}, pn spread {spread:.3e}, levels {levels}"); sys.exit(0 if ok else 1)
PYEOF
}
if want uniform; then
    echo "== uniform oblique flow through xz patches / an xz-refined zero-coefficient body, EXACT"
    for c in uniform_xz unibody_xz; do rm -f "g_${c}_"*.h5; run "g_$c" "$c.ini" 1 && uniform_check "$(last "g_$c")"; done
fi
if want twin; then
    echo "== gate (c): all-refined xz box == the doubled-x,z single-level twin (tolerance 0), and 1 == $RANKS ranks"
    rm -f g_ar1_*.h5 g_arN_*.h5 g_twin_*.h5
    run g_ar1 allref_xz.ini 1 && run g_twin twin_xz.ini 1 \
        && { "$PY" "$CMP" "$(last g_twin)" "$(last g_ar1)" un vn wn pn --tolerance 0 > .cmp_twin.txt 2>&1 \
             && echo "PASS allref_xz == twin_xz (EXACT)" || { echo "FAIL allref_xz vs twin_xz"; head -6 .cmp_twin.txt; fail=1; }; }
    run g_arN allref_xz.ini "$RANKS" \
        && { "$PY" "$CMP" "$(last g_ar1)" "$(last g_arN)" un vn wn pn --tolerance 0 > .cmp_rank.txt 2>&1 \
             && echo "PASS allref_xz 1 rank == $RANKS ranks (EXACT)" || { echo "FAIL allref_xz rank independence"; head -6 .cmp_rank.txt; fail=1; }; }
fi
if want order; then
    echo "== Beltrami xz patch pair: L2 velocity error at t = 0.25, order"
    e=()
    for n in 32 64; do rm -f "g_bp${n}_"*.h5; run "g_bp$n" "bp_xz_$n.ini" 1 || continue
        l2=$("$PY" "$ROOT/tools/check_beltrami.py" "$(last "g_bp$n")" | grep -o 'vel=[0-9.e+-]*' | cut -d= -f2); echo "   base $n: L2 vel = $l2"; e+=("$l2"); done
    [ ${#e[@]} = 2 ] && "$PY" - "${e[0]}" "${e[1]}" <<'PYEOF' || fail=1
import sys, math; a, b = float(sys.argv[1]), float(sys.argv[2]); o = math.log2(a / b)
print(f"   order {o:.2f}  " + ("PASS (>= 2.5)" if o >= 2.5 else "FAIL (< 2.5)")); sys.exit(0 if o >= 2.5 else 1)
PYEOF
fi
if want rans; then
    echo "== rans_xz smoke: 50 steps finite, 1 rank == $RANKS ranks EXACT"
    rm -f g_rx1_*.h5 g_rxN_*.h5
    run g_rx1 rans_xz.ini 1 && run g_rxN rans_xz.ini "$RANKS" && {
        "$PY" - "$(last g_rx1)" <<'PYEOF' || fail=1
import sys, h5py, numpy as np
h = h5py.File(sys.argv[1], "r"); m = {v: float(np.abs(np.asarray(h[v])).max()) for v in ("un","vn","wn","pn","k","omega","nut")}
ok = all(np.isfinite(x) for x in m.values()) and m["un"] < 1e3
print(("PASS" if ok else "FAIL") + " finite: " + "  ".join(f"max|{v}|={x:.3e}" for v, x in m.items())); sys.exit(0 if ok else 1)
PYEOF
        "$PY" "$CMP" "$(last g_rx1)" "$(last g_rxN)" un vn wn pn k omega nut --tolerance 0 > .cmp_rans.txt 2>&1 \
            && echo "PASS rans_xz 1 rank == $RANKS ranks (EXACT)" || { echo "FAIL rans_xz rank independence"; head -8 .cmp_rans.txt; fail=1; }
    }
fi
echo; [ "$fail" = 0 ] && echo "refine2d gates ($sel): ALL PASS" || echo "refine2d gates ($sel): FAILURES"
exit $fail
