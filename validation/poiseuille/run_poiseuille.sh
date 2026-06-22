#!/usr/bin/env bash
# Poiseuille + wall blowing/suction test (exact steady NS solution). Stresses
# momentum convection + diffusion AND the pressure projection.
#
#   ./run_poiseuille.sh [cpu|gpu] [niter]
#
# [niter] (optional) overrides the projection iteration count [pressure] niter in
# every case -- the steady smooth Poiseuille converges fine with niter=3 (much
# faster than the default). e.g.  ./run_poiseuille.sh gpu 3
#
# A. Convergence: walls in y, sweep the wall-normal resolution ny=32/64/128 (the
#    homogeneous x,z stay 4 cells; lx=lz=4/ny keeps the cells cubic). L2 ~ h^2.
# B. Orientation: the same physics with walls in x / y / z must give an IDENTICAL
#    error at a fixed resolution (the scheme is isotropic) -- a cheap, sufficient
#    proof that the x and z paths work once y is shown to converge.

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
export MOBY_POISEUILLE=1
CHK=../../tools/check_poiseuille.py
mkdir -p runs

run() {  # run() <dir> <ini>
    local d="runs/$1"; mkdir -p "$d"; rm -f "$d"/poiseuille_*.h5
    cp "$2" "$d/input.ini"
    [ -n "$NITER" ] && sed -i "s/^niter = .*/niter = $NITER/" "$d/input.ini"
    ( cd "$d" && mpirun -x MOBY_POISEUILLE -x OMP_TARGET_OFFLOAD -n 1 "$BIN" input.ini > run.log 2>&1 )
    python3 "$CHK" "$(ls -1 "$d"/poiseuille_*.h5 | sort -t_ -k2 -n | tail -1)"
}

echo "== A. wall-normal convergence (walls in y)"
for N in 32 64 128; do
    LX=$(python3 -c "print(4.0/$N)")
    sed -e "s|^ny = .*|ny = $N|" -e "s|^lx = .*|lx = $LX|" -e "s|^lz = .*|lz = $LX|" \
        walls_y.ini > "runs/wy_$N.ini"
    run "conv_$N" "runs/wy_$N.ini"
done

echo "== B. orientation equivalence (walls in x / y / z, same resolution)"
for o in x y z; do run "orient_$o" "walls_$o.ini"; done

echo "== done. A: L2 should drop ~4x per refinement (order 2)."
echo "         B: the three L2 values must match (isotropic scheme)."
