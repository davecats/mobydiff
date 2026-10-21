#!/usr/bin/env bash
# The R1 "NOT yet done" items of README.md, run cheaply (decision 6 of the
# 2026-10-21 re-validation handout): red-black SOR across a 2:1 interface is
# admitted as a production projection on refined cases only if these hold.
#
#   [BIN=/abs/moby_solve] [RANKS=4] ./run_r1_items.sh cpu    # (a) the omega sweep
#   [BIN=/abs/moby_solve]           ./run_r1_items.sh gpu    # (b) Beltrami y-slab, (c) refine_body 2000 steps
#
# (a) omega stability with an interface present: refined_channel.ini for 200
#     steps at sor 0.8 .. 1.7, finite or not (the gate runs at 1.5 only).
# (b) the Beltrami y-slab (../beltrami/slab_y.ini, exact decaying solution
#     through a y-normal 2:1 slab) with red-black niter 6 / sor 1.5 against
#     the same case with Chebyshev niter 12: the L2 velocity error at t = 2
#     must not grow by more than 1.5x -- a per-colour transfer error at the
#     interface shows up here where uniform flow is blind to it.
# (c) refine_body stability: the les_ibm channel_ibm_refine.ini (2:1 around
#     an immersed wall) with red-black niter 6 for 2000 steps, finite.
# NOT here: the global mass residual with a refined patch (the divsum tool
# left with the MOBY_* hooks; the body_outlet checker's interior rms
# divergence on the pocket case is the surviving mass statement) and the
# performance item (the performance handout). The developed Re_tau 180 band
# channel against Chebyshev is ../channel_interface/developed/README.md.
set -uo pipefail
MODE="${1:-cpu}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
BIN="${BIN:-$ROOT/build_${MODE}/moby_solve}"
PY="${PY:-$HOME/ibmc/bin/python}"; [ -x "$PY" ] || PY=python3
RANKS="${RANKS:-4}"; [ "$MODE" = gpu ] && RANKS=1
PIM="$ROOT/tools/prepare_if_missing.sh"
cd "$HERE"; fail=0

finite() {  # <field.h5> <bound on max|un|>
    "$PY" - "$1" "$2" <<'EOF'
import sys, h5py, numpy as np
h = h5py.File(sys.argv[1], "r"); m = {v: float(np.abs(np.asarray(h[v])).max()) for v in ("un","vn","wn","pn")}
ok = all(np.isfinite(x) for x in m.values()) and m["un"] < float(sys.argv[2])
print(("PASS" if ok else "FAIL") + "  " + "  ".join(f"max|{v}|={x:.3e}" for v, x in m.items())); sys.exit(0 if ok else 1)
EOF
}

if [ "$MODE" = cpu ]; then
    echo "== (a) omega sweep, red-black niter 3 across the 2:1 wall bands, 200 steps, $RANKS ranks"
    echo "   (the item IS the sweep: it passes when the production sor = 1.5 is stable WITH"
    echo "    margin, i.e. 1.6 is stable too; the first unstable value is the measured limit)"
    stable=""; limit=""
    for w in 0.8 1.0 1.2 1.4 1.5 1.6 1.7; do
        sed -e "s/^sor = .*/sor = $w/" -e 's/^nsteps = .*/nsteps = 200/' -e 's/^field_interval = .*/field_interval = 200/' \
            -e "s/^field_prefix = .*/field_prefix = r1w_$w/" refined_channel.ini > ".r1w_$w.ini"
        rm -f "r1w_${w}_"*.h5
        "$PIM" "$RANKS" "$BIN" ".r1w_$w.ini" > /dev/null || { echo "PREPARE FAILED"; fail=1; continue; }
        if mpirun -n "$RANKS" "$BIN" ".r1w_$w.ini" > ".r1w_$w.log" 2>&1 && [ -f "r1w_${w}_200.h5" ] \
            && out=$(finite "r1w_${w}_200.h5" 1e3); then
            echo "   sor $w: stable   ${out#PASS }"; stable="$stable $w"
        else
            echo "   sor $w: UNSTABLE ${out#FAIL }"; [ -z "$limit" ] && limit=$w
        fi
    done
    case "$stable" in *" 1.5"*" 1.6"*) echo "   PASS: stable through sor 1.6, first unstable value ${limit:-none above 1.7}";;
        *) echo "   FAIL: sor 1.5 or 1.6 unstable (stable:$stable)"; fail=1;; esac
fi

if [ "$MODE" = gpu ]; then
    echo "== (b) Beltrami y-slab: red-black niter 6 / sor 1.5 vs Chebyshev niter 12 (L2 velocity error at t = 2)"
    mkdir -p runs
    for v in cheb12 rb6; do
        d="runs/slab_y_$v"; mkdir -p "$d"; rm -f "$d"/beltrami_*.h5
        if [ $v = cheb12 ]; then subs='s/^niter = .*/niter = 12\naccel = chebyshev/'; else subs='s/^niter = .*/niter = 6\nsolver = redblack/; s/^sor = .*/sor = 1.5/'; fi
        sed -e "$subs" ../beltrami/slab_y.ini > "$d/input.ini"
        ( cd "$d" && "$PIM" 1 "$BIN" input.ini > /dev/null && mpirun -n 1 "$BIN" input.ini > run.log 2>&1 ) || { echo "RUN FAILED $v"; fail=1; continue; }
        f=$(ls -1 "$d"/beltrami_*.h5 | sort -t_ -k2 -n | tail -1)
        "$PY" "$ROOT/tools/check_beltrami.py" "$f" | grep -E 'L2|Linf' | sed "s/^/   $v: /"
        eval "l2_$v=$("$PY" "$ROOT/tools/check_beltrami.py" "$f" | grep -o 'vel=[0-9.e+-]*' | cut -d= -f2)"
    done
    "$PY" - "$l2_cheb12" "$l2_rb6" <<'EOF' || fail=1
import sys; a, b = float(sys.argv[1]), float(sys.argv[2]); r = b / a
print(f"   red-black / Chebyshev L2 velocity error = {r:.3f}  " + ("PASS (<= 1.5)" if r <= 1.5 else "FAIL (> 1.5)")); sys.exit(0 if r <= 1.5 else 1)
EOF

    echo "== (c) refine_body stability: les_ibm channel_ibm_refine.ini, red-black niter 6 / sor 1.5, 2000 steps"
    sed -e 's/^niter = .*/niter = 6\nsolver = redblack\nsor = 1.5/' -e '/^accel = /d' -e 's/^nsteps = .*/nsteps = 2000/' \
        -e 's/^t_final = .*/t_final = 0.0/' -e 's/^field_interval = .*/field_interval = 2000/' \
        -e 's/^field_prefix = .*/field_prefix = r1body/' ../channel_interface/les_ibm/channel_ibm_refine.ini > .r1body.ini
    rm -f r1body_*.h5
    ( cd ../channel_interface/les_ibm && "$PIM" 1 "$BIN" "$HERE/.r1body.ini" > /dev/null \
        && mpirun -n 1 "$BIN" "$HERE/.r1body.ini" > "$HERE/.r1body.log" 2>&1 && mv r1body_*.h5 "$HERE"/ ) \
        || { echo "   RUN FAILED: $(grep -m1 -i 'error\|nan\|stop' .r1body.log | cut -c1-80)"; fail=1; }
    [ -f r1body_2000.h5 ] && { printf "   2000 steps: "; finite r1body_2000.h5 1e3 || fail=1; }
fi

echo; [ "$fail" = 0 ] && echo "R1 items ($MODE): ALL PASS" || echo "R1 items ($MODE): FAILURES"
exit $fail
