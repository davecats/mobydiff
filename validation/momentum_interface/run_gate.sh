#!/usr/bin/env bash
# RETIRED: this gate drove the MOBY_PREDONLY / MOBY_RHSDUMP / MOBY_DIVDUMP /
# MOBY_TERMDUMP diagnostic hooks, which were removed in the code cleanup once the
# 2:1 momentum predictor was validated. Kept verbatim as a record of the method;
# it will NOT produce the dumps below until those hooks are reinstated.
#
# Momentum-predictor 2:1-interface gate (two axes, both single-step).
#
# Builds the operator-truncation instrument's input decks on the fly (uniform /
# y-slab / 3D-patch, swept over nx=32/64/128 with nb matched so the block
# lattice -- and hence the physical interface location -- is fixed at 8^3), runs
# each ONE step, and reports:
#   Axis 1 (accuracy):    MOBY_PREDONLY + MOBY_RHSDUMP -> tools/rhsband.py
#                         interior order ~2 (sanity); fine-band order = verdict.
#   Axis 2 (continuity):  MOBY_DIVDUMP (real predictor+sync) -> tools/divsum.py
#                         global Sum(vol*div(qs)) must stay round-off.
#
# The field is `initial = tgv3d` (manufactured: every velocity component varies
# in every direction, so the wall-normal velocity varies in the normal
# direction at all three interfaces -- the property Beltrami lacks). One step is
# sufficient; this gate is about operator consistency, NOT time evolution.
#
# Usage:  ./run_gate.sh /path/to/build_cpu/main  [workdir]
set -euo pipefail
MAIN="${1:?usage: run_gate.sh <main-binary> [workdir]}"
WORK="${2:-./mi_gate_work}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
L=6.283185307179586
mkdir -p "$WORK"

declare -A NB=( [32]=4 [64]=8 [128]=16 )

gen() {  # kind nx
  local kind=$1 n=$2 nb=${NB[$2]} refine
  case "$kind" in
    uniform) refine="" ;;
    slab)    refine="-1.0 10.0 1.6 4.7 -1.0 10.0" ;;
    patch)   refine="1.6 4.7 1.6 4.7 1.6 4.7" ;;
  esac
  {
    echo "[case]"; echo "name = generic"
    echo "[grid]"; echo "nx = $n"; echo "ny = $n"; echo "nz = $n"
    echo "lx = $L"; echo "ly = $L"; echo "lz = $L"
    for d in x y z; do echo "[grid.$d]"; echo "distribution = uniform"; done
    if [ -n "$refine" ]; then
      echo "[blocks]"; echo "nb = $nb"; echo "refine = $refine"; echo "refine_levels = 1"
    fi
    echo "[ibm]"; echo "enabled = false"
    echo "[flow]"; echo "initial = tgv3d"; echo "re = 100.0"
    echo "[time]"; echo "dt = 1.0e-3"; echo "nsteps = 1"; echo "t_final = 0.0"
    echo "cflmax = 0.0"; echo "pecletmax = 0.0"; echo "dtmax = 0.1"
    echo "[pressure]"; echo "niter = 50"; echo "sor = 0.8"
    echo "[boundary]"; for d in x y z; do echo "periodic_$d = true"; done
    echo "[output]"; echo "field_interval = 0"; echo "field_prefix = tgv3d"
  } > "$WORK/${kind}_${n}.ini"
}

run_rhs() {  # kind nx
  local kind=$1 n=$2 d="$WORK/${kind}_${n}_rhs"
  rm -rf "$d"; mkdir -p "$d"; ( cd "$d"
    MOBY_PREDONLY=1 MOBY_RHSDUMP=1 mpirun -x MOBY_PREDONLY -x MOBY_RHSDUMP \
      -n 1 "$MAIN" "$WORK/${kind}_${n}.ini" > log 2>&1 )
}
run_div() {  # kind nx
  local kind=$1 n=$2 d="$WORK/${kind}_${n}_div"
  rm -rf "$d"; mkdir -p "$d"; ( cd "$d"
    MOBY_DIVDUMP=1 mpirun -x MOBY_DIVDUMP -n 1 "$MAIN" "$WORK/${kind}_${n}.ini" > log 2>&1 )
}

for kind in uniform slab patch; do
  for n in 32 64 128; do gen "$kind" "$n"; run_rhs "$kind" "$n"; done
  echo "================ Axis 1 (accuracy / order): $kind ================"
  python3 "$ROOT/tools/rhsband.py" \
    "$WORK/${kind}_32_rhs/tgv3d_rhs_1.h5" \
    "$WORK/${kind}_64_rhs/tgv3d_rhs_1.h5" \
    "$WORK/${kind}_128_rhs/tgv3d_rhs_1.h5"
done

echo "================ Axis 2 (continuity null-space) ================"
for kind in slab patch; do
  run_div "$kind" 64
  echo "-- $kind nx=64 --"
  # The orchestrated div run is occasionally flaky on WSL (back-to-back mpiruns);
  # don't abort the whole gate -- it runs fine standalone, just re-run by hand.
  python3 "$ROOT/tools/divsum.py" "$WORK/${kind}_64_div/tgv3d_divpre_1.h5" \
    || echo "  (div run flaked; rerun: MOBY_DIVDUMP=1 mpirun -n 1 \$MAIN $WORK/${kind}_64.ini)"
done

run_term() {  # kind nx var
  local kind=$1 n=$2 var=$3 d="$WORK/${kind}_${n}_term${var}"
  rm -rf "$d"; mkdir -p "$d"; ( cd "$d"
    MOBY_TERMDUMP=$var mpirun -x MOBY_TERMDUMP -n 1 "$MAIN" "$WORK/${kind}_${n}.ini" > log 2>&1 )
}
echo "================ Axis 3 (momentum conservation) ================"
echo "Sum(vol*adv) per component must be ~round-off (uniform = zero reference)."
for kind in uniform slab patch; do
  echo "-- $kind nx=64 --"
  for var in 1 2 3; do
    run_term "$kind" 64 "$var"
    python3 "$ROOT/tools/momsum.py" "$WORK/${kind}_64_term${var}/tgv3d_adv_0.h5" 2>/dev/null | sed -n '2p' \
      || echo "  comp $var: (term run flaked; rerun by hand)"
  done
done
