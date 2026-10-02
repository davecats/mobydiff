#!/usr/bin/env bash
# The convective limit of the RK3 step, bracketed (box.ini, seed_ic.py).
#
#     [SOLVER=../../build_cpu/moby_solve] ./run_gate.sh
#
# Central convection linearised about a uniform flow has the imaginary
# spectrum i sum_d u_d sin(k_d h_d)/h_d, and RK3 is stable on the imaginary
# axis up to sqrt(3), so the limit is
#
#     dt x sum_d |u_d|/h_d <= 1.732 ,
#
# and `cflmax` bounds exactly that sum (since 2026-10-02; it bounded the
# largest single component before). In terms of the Courant number C of one
# component the limit depends on the direction of the flow:
#
#     flow       sum_d |u_d|/h_d    limit on C
#     (1,0,0)    1 |u|/h            1.732
#     (1,1,0)    2 |u|/h            0.866
#     (1,1,1)    3 |u|/h            0.577
#
# GATE (the perturbation of the uniform flow, max over the three components,
# seeded at 5e-7)
#   fixed step, one leg below (sum 1.65) and one above (sum 1.80) per flow,
#   600 steps: below it never exceeds twice its seed; above it grows by more
#   than 1e3 (or the run ends in NaN).
#   adaptive step (dtmax large), the SAME two values for the three flows:
#   cflmax = 1.60 stays bounded, cflmax = 1.85 grows and the solver announces
#   at init that the value is above the limit.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
SOLVER=${SOLVER:-$ROOT/build_cpu/moby_solve}
status=0

# leg <tag> <u> <v> <w> <sum Courant number | cflmax> <fixed|adaptive>
leg() {
    local tag=$1 u=$2 v=$3 w=$4 C=$5 mode=$6 dt dtmax cfl
    dt=$(python3 -c "import math; h=2*math.pi/16; print(repr($C*h/(abs($u)+abs($v)+abs($w))))")
    if [ "$mode" = fixed ]; then dtmax=$dt; cfl=0.0
    else dt=1.0e-4; dtmax=10.0; cfl=$C; fi
    sed -e "s/^initial_u = .*/initial_u = $u/" -e "s/^initial_v = .*/initial_v = $v/" -e "s/^initial_w = .*/initial_w = $w/" \
        -e "s/^dt = .*/dt = $dt/" -e "s/^dtmax = .*/dtmax = $dtmax/" -e "s/^cflmax = .*/cflmax = $cfl/" \
        -e "s/^field_prefix = .*/field_prefix = cl_$tag/" box.ini > .cl_$tag.ini
    rm -f cl_${tag}_*.h5 cl_$tag.case.h5 cl_${tag}_ic.h5
    "$ROOT/tools/prepare_if_missing.sh" 1 "$SOLVER" .cl_$tag.ini .cl_$tag.prep.log \
        || { echo "PREPARE FAILED $tag"; status=1; return; }
    # The template: one step of the same case; then the seeded restart.
    sed -e "s/^nsteps = .*/nsteps = 1/" -e "s/^field_interval = .*/field_interval = 0/" .cl_$tag.ini > .cl_${tag}_t.ini
    mpirun -n 1 "$SOLVER" .cl_${tag}_t.ini > .cl_${tag}_t.log 2>&1 || { echo "TEMPLATE FAILED $tag"; status=1; return; }
    python3 seed_ic.py cl_${tag}_1.h5 cl_${tag}_ic.h5 $u $v $w && rm -f cl_${tag}_1.h5
    printf '\n[restart]\nfile = cl_%s_ic.h5\n' "$tag" >> .cl_$tag.ini
    mpirun -n 1 "$SOLVER" .cl_$tag.ini > .cl_$tag.log 2>&1   # a blow-up may exit non-zero
}

# check <tag> <u> <v> <w> <stable|unstable> [announced]
check() {
    python3 - "$@" <<'PY' || status=1
import glob, sys
import h5py, numpy as np
tag, vel, want = sys.argv[1], [float(v) for v in sys.argv[2:5]], sys.argv[5]
log = open(f".cl_{tag}.log").read()
announced = "is above the RK3 limit of central convection" in log
amp = []
for name in sorted(glob.glob(f"cl_{tag}_[0-9]*.h5"), key=lambda n: int(n.split("_")[-1][:-3])):
    with h5py.File(name, "r") as f:
        amp.append(max(float(np.abs(f[v][...] - u).max()) for v, u in zip(("un", "vn", "wn"), vel)))
seed = 5.0e-7
worst = max(amp) if amp and all(np.isfinite(amp)) else float("inf")
ok = (want == "stable" and worst <= 2.0*seed) or (want == "unstable" and worst >= 1.0e3*seed)
if len(sys.argv) > 6:
    ok = ok and announced
print(f"{'PASS' if ok else 'FAIL'} {tag:10s} expected {want:8s} largest perturbation {worst:.2e} (seed {seed:.0e})"
      f"{'  announced at init' if announced else ''}")
sys.exit(0 if ok else 1)
PY
}

#    tag      u v w   sum   mode
leg  x_below  1 0 0   1.65  fixed;    check x_below  1 0 0 stable
leg  x_above  1 0 0   1.80  fixed;    check x_above  1 0 0 unstable
leg  xy_below 1 1 0   1.65  fixed;    check xy_below 1 1 0 stable
leg  xy_above 1 1 0   1.80  fixed;    check xy_above 1 1 0 unstable
leg  d_below  1 1 1   1.65  fixed;    check d_below  1 1 1 stable
leg  d_above  1 1 1   1.80  fixed;    check d_above  1 1 1 unstable
leg  ax_160   1 0 0   1.60  adaptive; check ax_160   1 0 0 stable
leg  ax_185   1 0 0   1.85  adaptive; check ax_185   1 0 0 unstable announced
leg  axy_160  1 1 0   1.60  adaptive; check axy_160  1 1 0 stable
leg  axy_185  1 1 0   1.85  adaptive; check axy_185  1 1 0 unstable announced
leg  ad_160   1 1 1   1.60  adaptive; check ad_160   1 1 1 stable
leg  ad_185   1 1 1   1.85  adaptive; check ad_185   1 1 1 unstable announced

[ $status = 0 ] && echo "courant-limit gate: PASS" || echo "courant-limit gate: FAIL"
exit $status
