#!/usr/bin/env bash
# The penalization-factor gate: du/dt = f - lambda u (velocity) and
# ds/dt = S - lambda (s - s_body) (a Dirichlet-penalized scalar) in a uniform
# periodic box (see decay.ini). T = 0.2, dt = 0.04 / 2^k, lambda = 20 and 200.
#
#     [SOLVER=../../build_cpu/moby_solve] [REF=/abs/older/moby_solve] ./run_decay.sh
#
# GATE, for the velocity AND the scalar:
#   lambda = 20   THIRD-order convergence: every observed order >= 2.7 and
#                 |err| <= 2e-4 at dt = 0.04;
#   lambda = 200  the steady fixed point (lambda T = 40): |err| <= 1e-13 at
#                 every dt -- the factors satisfy state + x*incr = 1;
#   the fields uniform to round-off, v = w = 0.
# With REF set, the same runs on that binary are printed with their observed
# order (implicit Euler reads ~1.0; the exact exponential of step 9 reads
# round-off for the velocity).
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
SOLVER=${SOLVER:-$ROOT/build_cpu/moby_solve}
REF=${REF:-}
rm -f decay.case.h5
"$ROOT/tools/prepare_if_missing.sh" 1 "$SOLVER" decay.ini prepare.log || { echo "PREPARE FAILED"; exit 1; }

status=0
for lam in 20 200; do
    python3 - "$lam" <<'PY'
import sys, h5py
lam, pr = float(sys.argv[1]), 0.5                  # pr: [scalar.1] in decay.ini
with h5py.File("decay.case.h5", "r+") as f:
    f["coef_blocks"][...] = lam                    # every DOF, ghosts included
    f["coef_p_blocks"][...] = lam*pr               # the scalar's rate is coef_p/Pr
PY
    for side in new ref; do
        bin=$SOLVER
        if [ $side = ref ]; then [ -n "$REF" ] || continue; bin=$REF; fi
        for k in 0 1 2 3; do
            n=$(( 5 << k )); dt=$(python3 -c "print(0.04/2**$k)")
            sed -e "s/^dt *=.*/dt = $dt/" -e "s/^dtmax *=.*/dtmax = $dt/" -e "s/^nsteps *=.*/nsteps = $n/" \
                -e "s/^field_interval *=.*/field_interval = $n/" \
                -e "s/^field_prefix *=.*/field_prefix = d_${side}_${lam}_$k/" decay.ini > .d.ini
            rm -f d_${side}_${lam}_${k}_*.h5
            mpirun -n 1 "$bin" .d.ini > .d_${side}_${lam}_$k.log 2>&1 || { echo "RUN FAILED: .d_${side}_${lam}_$k.log"; exit 1; }
        done
        python3 - "$side" "$lam" <<'PY' || status=1
import glob, math, sys
import h5py, numpy as np
side, lam = sys.argv[1], float(sys.argv[2])
T = 0.2
decay = math.exp(-lam*T)
# velocity: u0 = 1, f = 1; scalar: s0 = 1, S = 1, s_body = 0.25 (decay.ini)
exact = {"un": 1.0*decay + (1.0/lam)*(-math.expm1(-lam*T)),
         "theta": (0.25 + 1.0/lam) + (1.0 - 0.25 - 1.0/lam)*decay}
bad = False
for name in ("un", "theta"):
    errs, orders = [], []
    print(f"  {side}  lambda = {lam:g}   {name}: exact value at T = {exact[name]:.16f}")
    for k in range(4):
        p = sorted(glob.glob(f"d_{side}_{int(lam)}_{k}_*.h5"))[-1]
        with h5py.File(p) as h:
            q, v, w, t = h[name][...], h["vn"][...], h["wn"][...], float(h.attrs["t_current"])
        assert abs(t - T) < 1e-12, f"{p} ends at t = {t}"
        e = abs(q.mean() - exact[name]); errs.append(e)
        note = ""
        if k > 0 and e > 0 and errs[k-1] > 0:
            orders.append(math.log2(errs[k-1]/e)); note = f"   order {orders[-1]:5.2f}"
        spread, vw = np.ptp(q), max(abs(v).max(), abs(w).max())
        print(f"     dt = {0.04/2**k:8.5f}  x_max = {lam*0.04/2**k*8/15:6.3f}  value = {q.mean():.16f}  "
              f"|err| = {e:.2e}  spread {spread:.1e}  max|v,w| {vw:.1e}{note}")
        if side == "new" and (spread > 1e-14 or vw > 1e-14):
            print("     FAIL: the field is not uniform / v, w not zero"); bad = True
    if side != "new":
        continue
    if lam == 20 and (errs[0] > 2e-4 or len(orders) < 3 or min(orders) < 2.7):
        print("     FAIL: not third order (order >= 2.7, |err| <= 2e-4 at dt = 0.04)"); bad = True
    if lam == 200 and max(errs) > 1e-13:
        print("     FAIL: the steady fixed point is not reached to round-off"); bad = True
sys.exit(1 if bad else 0)
PY
    done
done
rm -f .d.ini
[ $status -eq 0 ] && echo "penalization decay gate: PASS" || echo "penalization decay gate: FAIL"
exit $status
