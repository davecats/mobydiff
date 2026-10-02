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
# `pecletmax` bounds P = dt x nu/h^2 of the finest SINGLE direction, so the
# largest stable P depends on how many directions are that fine:
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
#   adaptive step on 16^3 (cflmax 0.8, dtmax large):
#     pecletmax = 0.5: the init line reports the effective sum bound 1.50, the
#                      warning fires, the run blows up;
#     pecletmax = 0.2: effective sum bound 0.60, no warning, the flow decays.
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

# check <tag> <stable|unstable> [expected effective sum bound]
check() {
    python3 - "$@" <<'PY' || status=1
import glob, re, sys
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
if len(sys.argv) > 3:
    m = re.search(r"pecletmax x ratio =\s*([0-9.]+)", log)
    got = float(m.group(1)) if m else float("nan")
    ok = ok and abs(got - float(sys.argv[3])) < 0.006
    note = f"  effective sum bound {got:.2f}"
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
leg  p_050   16 16 16  0.5    adaptive; check p_050 unstable 1.50
leg  p_020   16 16 16  0.2    adaptive; check p_020 stable 0.60

[ $status = 0 ] && echo "diffusion-limit gate: PASS" || echo "diffusion-limit gate: FAIL"
exit $status
