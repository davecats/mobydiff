#!/usr/bin/env bash
# 3D Beltrami / ABC flow test (fully 3D exact decaying NS solution): spatial
# convergence (uniform) + 2:1-interface artifact (refined 3D patch), both against
# the exact self-similar decay. Exercises all three directions (unlike the 2D TGV).
#
#   ./run_beltrami.sh [cpu|gpu] [niter]
#
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
export MOBY_BELTRAMI=1
CHK=../../tools/check_beltrami.py
mkdir -p runs

run() {  # run() <dir> <ini>
    local d="runs/$1"; mkdir -p "$d"; rm -f "$d"/beltrami_*.h5
    cp "$2" "$d/input.ini"
    [ -n "$NITER" ] && sed -i "s/^niter = .*/niter = $NITER/" "$d/input.ini"
    ( cd "$d" && mpirun -x MOBY_BELTRAMI -x OMP_TARGET_OFFLOAD -n 1 "$BIN" input.ini > run.log 2>&1 )
    python3 "$CHK" "$(ls -1 "$d"/beltrami_*.h5 | sort -t_ -k2 -n | tail -1)"
}

echo "== uniform convergence sweep (cube)"
for N in 32 64; do
    sed -e "s|^nx = .*|nx = $N|" -e "s|^ny = .*|ny = $N|" -e "s|^nz = .*|nz = $N|" \
        uniform.ini > "runs/uniform_$N.ini"
    run "uniform_$N" "runs/uniform_$N.ini"
done

echo "== refined central 3D patch (nx=64 base)"
run "refined" refined.ini

echo "== done. uniform L2 should drop ~4x (order 2); refined shows the 3D interface artifact."
