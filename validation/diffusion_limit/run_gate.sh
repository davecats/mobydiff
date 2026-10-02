#!/usr/bin/env bash
# The explicit-diffusion limit of the RK3 step, bracketed (box.ini).
#
#     [SOLVER=../../build_cpu/moby_solve] ./run_gate.sh
#
# RK3 is stable on the negative real axis down to z = -2.5127 and the extreme
# eigenvalue of the discrete Laplacian is -4 sum_d nu/h_d^2, so the limit is
#
#     dt x sum_d nu/h_d^2 <= 0.628 .
#
# and `pecletmax` bounds exactly that sum (since 2026-10-02; it bounded the
# finest single direction's dt x nu/h^2 before). In terms of P = dt x nu/h^2
# of the finest direction the limit depends on how many directions are that
# fine:
#
#     grid        sum_d nu/h_d^2      limit on P
#     16x16x16    3    nu/h^2         0.2094
#     16x16x4     2.06 nu/h^2         0.3046
#     16x4x4      1.13 nu/h^2         0.5584
#
# GATE
#   fixed step (cflmax = pecletmax = 0), one leg just below and one just above
#   each limit, 600 steps:
#     below: max|u| decays monotonically over the snapshots and no warning;
#     above: the run blows up (NaN or max|u| > 1) and the solver prints its
#            one-time "exceeds the RK3 limit of explicit diffusion" warning.
#   adaptive step (cflmax 0.8, dtmax large), the SAME two values on all three
#   grids -- one number is safe or not whatever the grid:
#     pecletmax = 0.60: no warning, the flow decays;
#     pecletmax = 0.66: the warning, and the run leaves the decaying solution
#                       (NaN, or a state held by the Courant limiter).
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
SOLVER=${SOLVER:-$ROOT/build_cpu/moby_solve}
status=0

# leg <tag> <nx> <ny> <nz> <P> <fixed|adaptive>: P is the fixed step's
# nu dt/h^2 (h the finest spacing), or pecletmax of an adaptive leg.
leg() {
    local tag=$1 nx=$2 ny=$3 nz=$4 P=$5 mode=$6 dt dtmax pec cfl
    dt=$(python3 -c "import math; h=2*math.pi/max($nx,$ny,$nz); print(repr($P*h*h))")
    if [ "$mode" = fixed ]; then dtmax=$dt; pec=0.0; cfl=0.0
    else dt=1.0e-4; dtmax=1.0; pec=$P; cfl=0.8; fi
    sed -e "s/^nx = .*/nx = $nx/" -e "s/^ny = .*/ny = $ny/" -e "s/^nz = .*/nz = $nz/" \
        -e "s/^dt = .*/dt = $dt/" -e "s/^dtmax = .*/dtmax = $dtmax/" \
        -e "s/^pecletmax = .*/pecletmax = $pec/" -e "s/^cflmax = .*/cflmax = $cfl/" \
        -e "s/^field_prefix = .*/field_prefix = dl_$tag/" box.ini > .dl_$tag.ini
    rm -f dl_${tag}_*.h5 dl_$tag.case.h5
    "$ROOT/tools/prepare_if_missing.sh" 1 "$SOLVER" .dl_$tag.ini .dl_$tag.prep.log \
        || { echo "PREPARE FAILED $tag"; status=1; return; }
    mpirun -n 1 "$SOLVER" .dl_$tag.ini > .dl_$tag.log 2>&1   # a blow-up may exit non-zero
}

# check <tag> <stable|unstable>
check() {
    python3 - "$@" <<'PY' || status=1
import sys
import h5py, numpy as np
tag, want = sys.argv[1], sys.argv[2]
log = open(f".dl_{tag}.log").read()
warned = "exceeds the RK3 limit of explicit diffusion" in log
amp = []
for step in (100, 200, 300, 400, 500, 600):
    try:
        with h5py.File(f"dl_{tag}_{step}.h5", "r") as f:
            amp.append(max(float(np.abs(f[v][...]).max()) for v in ("un", "vn", "wn")))
    except OSError:
        amp.append(float("nan"))
blown = any(not np.isfinite(a) or a > 1.0 for a in amp)
decays = all(np.isfinite(amp)) and all(b < a for a, b in zip(amp, amp[1:]))
ok = (want == "stable" and decays and not warned) or (want == "unstable" and blown and warned)
note = ""
last = "nan" if not np.isfinite(amp[-1]) else f"{amp[-1]:.2e}"
print(f"{'PASS' if ok else 'FAIL'} {tag:10s} expected {want:8s} max|u| at 100/600: {amp[0]:.2e} / {last}"
      f"  warning {'yes' if warned else 'no'}{note}")
sys.exit(0 if ok else 1)
PY
}

#    tag     nx ny nz   P      mode
leg  a_below 16 16 16  0.205  fixed;    check a_below stable
leg  a_above 16 16 16  0.215  fixed;    check a_above unstable
leg  b_below 16 16  4  0.295  fixed;    check b_below stable
leg  b_above 16 16  4  0.315  fixed;    check b_above unstable
leg  c_below 16  4  4  0.545  fixed;    check c_below stable
leg  c_above 16  4  4  0.570  fixed;    check c_above unstable
leg  pa_060  16 16 16  0.60   adaptive; check pa_060 stable
leg  pa_066  16 16 16  0.66   adaptive; check pa_066 unstable
leg  pb_060  16 16  4  0.60   adaptive; check pb_060 stable
leg  pb_066  16 16  4  0.66   adaptive; check pb_066 unstable
leg  pc_060  16  4  4  0.60   adaptive; check pc_060 stable
leg  pc_066  16  4  4  0.66   adaptive; check pc_066 unstable

[ $status = 0 ] && echo "diffusion-limit gate: PASS" || echo "diffusion-limit gate: FAIL"
exit $status
