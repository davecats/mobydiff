#!/usr/bin/env bash
# Taylor-Green vortex test: spatial convergence (uniform) + interface-artifact
# check (refined patch), both against the exact decaying solution.
#
#   ./run_tgv.sh [cpu|gpu] [niter]
#
# Convergence: runs uniform nx=ny=32/64/128 and fits the L2 velocity-error order
# (expect ~2). Artifact: runs the central 2:1 refine patch, reports its L2 error
# vs the uniform nx=64 run and writes an error map.
# [niter] (optional) overrides the projection iteration count [pressure] niter.

set -uo pipefail
cd "$(dirname "$0")"
ARCH="${1:-cpu}"
NITER="${2:-}"
BIN="../../build_${ARCH}/main"
[ -x "$BIN" ] && BIN="$(cd "$(dirname "$BIN")" && pwd)/main" || {
    echo "binary $BIN not found (build with ./compile.sh $ARCH)"; exit 1; }
source /etc/profile.d/modules.sh 2>/dev/null
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3 2>/dev/null || true
export OMP_TARGET_OFFLOAD=MANDATORY
export MOBY_TGV=1
CHK=../../tools/check_tgv.py
mkdir -p runs

echo "== uniform convergence sweep"
for N in 32 64 128; do
    d=runs/uniform_$N; mkdir -p "$d"; rm -f "$d"/tgv_*.h5   # avoid stale files
    sed -e "s|^nx = .*|nx = $N|" -e "s|^ny = .*|ny = $N|" uniform.ini > "$d/input.ini"
    [ -n "$NITER" ] && sed -i "s|^niter = .*|niter = $NITER|" "$d/input.ini"
    ( cd "$d" && mpirun -x MOBY_TGV -x OMP_TARGET_OFFLOAD -n 1 "$BIN" input.ini > run.log 2>&1 )
    f=$(ls -1 "$d"/tgv_*.h5 | sort -t_ -k2 -n | tail -1)
    python3 "$CHK" "$f"
done

echo "== refined central patch (nx=64 base)"
d=runs/refined; mkdir -p "$d"; rm -f "$d"/tgv_*.h5   # avoid stale files
cp refined.ini "$d/input.ini"
[ -n "$NITER" ] && sed -i "s|^niter = .*|niter = $NITER|" "$d/input.ini"
( cd "$d" && mpirun -x MOBY_TGV -x OMP_TARGET_OFFLOAD -n 1 "$BIN" input.ini > run.log 2>&1 )
f=$(ls -1 "$d"/tgv_*.h5 | sort -t_ -k2 -n | tail -1)
python3 "$CHK" "$f" --error-map runs/refined/error_map.png

echo "== done. L2(uniform 32/64/128) should drop ~4x per refinement (order 2)."
