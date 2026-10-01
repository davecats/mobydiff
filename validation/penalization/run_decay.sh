#!/usr/bin/env bash
# The exact-penalization gate: du/dt = f - lambda u in a uniform periodic box
# (see decay.ini). T = 0.2, dt = 0.04 / 2^k, lambda = 20 and 200.
#
#     [SOLVER=../../build_cpu/moby_solve] [REF=/abs/older/moby_solve] ./run_decay.sh
#
# GATE: |u(T) - exact| <= 1e-13 at every dt and both lambda, the field
# uniform to round-off, v = w = 0. With REF set, the same runs on that binary
# are printed with their observed order (a pre-F7 binary reads ~1.0).
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
with h5py.File("decay.case.h5", "r+") as f:
    f["coef_blocks"][...] = float(sys.argv[1])     # every DOF, ghosts included
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
u0, f, T = 1.0, 1.0, 0.2
exact = u0*math.exp(-lam*T) + (f/lam)*(-math.expm1(-lam*T))
errs = []
print(f"  {side}  lambda = {lam:g}   exact u(T) = {exact:.16f}")
for k in range(4):
    p = sorted(glob.glob(f"d_{side}_{int(lam)}_{k}_*.h5"))[-1]
    with h5py.File(p) as h:
        u, v, w, t = h["un"][...], h["vn"][...], h["wn"][...], float(h.attrs["t_current"])
    assert abs(t - T) < 1e-12, f"{p} ends at t = {t}"
    e = abs(u.mean() - exact); errs.append(e)
    order = "" if k == 0 or e == 0 or errs[k-1] == 0 else f"   order {math.log2(errs[k-1]/e):5.2f}"
    print(f"     dt = {0.04/2**k:8.5f}  x_max = {lam*0.04/2**k*8/15:6.3f}  u = {u.mean():.16f}  |err| = {e:.2e}  "
          f"spread {np.ptp(u):.1e}  max|v,w| {max(abs(v).max(), abs(w).max()):.1e}{order}")
if side == "new" and max(errs) > 1e-13:
    print("     FAIL: the exact factor does not integrate the linear problem to round-off"); sys.exit(1)
PY
    done
done
rm -f .d.ini
[ $status -eq 0 ] && echo "penalization decay gate: PASS" || echo "penalization decay gate: FAIL"
exit $status
