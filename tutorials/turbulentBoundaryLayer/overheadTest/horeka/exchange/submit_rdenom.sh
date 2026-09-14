#!/bin/bash
#SBATCH --job-name=moby_rdenom
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:55:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# A/B for "hoist the denominator and divergence metrics too": the same lever that
# took jacobi_apply k2 from 88 to 64 registers, applied to the projection's other
# two volume kernels (compute_rdenom 110 -> 80, jacobi_compute_phi 94 -> 80).
# Readings pre-registered in PREREGISTERED_rdenom.md -- including why this is a
# weaker case than k2 was.
#
#   REF_DIR is a worktree at the APPLY-FIX commit (6708193), not at the older
#   pre-apply one: this job must isolate THIS increment.
#
# Four things in one allocation:
#   0  registers -- cuobjdump on both binaries
#   1  ncu, both binaries, now including compute_rdenom (never profiled before)
#   2  timeline pass A at 8 and 4 ranks, BEFORE and AFTER
#   3  Pass G + the 7-case suite -- max_abs 0
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/rdenom_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC"/run_timeline.sh "$SRC"/run_exchange.sh "$SRC"/run_mapgate.sh \
   "$SRC"/run_ncu.sh "$SRC"/collect_ncu.py "$SRC"/collect_exchange.py "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

for d in "$REF_DIR" "$CODE_DIR"; do
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
NEW="$CODE_DIR/build_gpu/moby_solve"; REF="$REF_DIR/build_gpu/moby_solve"
[ -x "$NEW" ] && [ -x "$REF" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_rdenom}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  dirty $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

echo "############ 0: registers (cuobjdump, both binaries) ############"
for side in ref new; do
    exe="$REF"; [ "$side" = new ] && exe="$NEW"
    echo "--- $side"
    cuobjdump -res-usage "$exe" | awk '
      /^ *Function nvkernel_/ { name=$2; sub(/:$/,"",name); next }
      /REG:/ && name != "" { print name"  "$1" "$2" "$3" "$4; name="" }' |
      grep -E 'pressure_solver|step_momentum' | sort
done 2>&1 | tee "$RES/registers.txt"

# 57 matching launches per step once compute_rdenom is in the set (one per
# substage ahead of that substage's 6 x 3), so skip two steps and profile
# exactly one substage.
echo "############ 1: ncu, both binaries ############"
for side in ref new; do
    exe="$REF"; [ "$side" = new ] && exe="$NEW"
    CFG_DIR="$SRC/../configs" \
        KERNEL_RE='jacobi_apply|jacobi_compute_phi|compute_rdenom' \
        LAUNCH_SKIP=114 LAUNCH_COUNT=19 \
        bash "$STG/run_ncu.sh" "$exe" "$RES/ncu_$side"
    python3 "$STG/collect_ncu.py" "$RES/ncu_$side" > "$RES/ncu_$side.md" 2>&1 \
        && { echo; cat "$RES/ncu_$side.md"; } || echo "=== collector failed for $side ==="
done

SPECS="refined_yp82_rect_jacobi:8:2 rect_jacobi:8:2 refined_yp82_rect_jacobi:4:1 rect_jacobi:4:1"
for side in before after; do
    exe="$REF"; [ "$side" = after ] && exe="$NEW"
    echo "############ 2$side: timeline ############"
    CFG_DIR="$SRC/../configs" TL_PASSES="A" SPECS="$SPECS" NSTEPS="${NSTEPS:-100}" \
        bash "$STG/run_timeline.sh" "$exe" "$RES/$side"
done

echo "############ 3: Pass G -- production cases, bit-exact ############"
REF="$REF" H5MAXDIFF="$CODE_DIR/tools/h5maxdiff" CFG_DIR="$SRC/../configs" \
    PASSES="G" NSTEPS_GATE=20 bash "$STG/run_exchange.sh" "$NEW" "$RES/gate"

echo "############ 4: the 7-case suite ############"
MOBY_ROOT="$CODE_DIR" NSTEPS=200 \
    bash "$STG/run_mapgate.sh" "$NEW" "$REF" "$RES/suite"

echo "=== summary: s/step and the projection buckets ==="
for side in before after; do
  for d in "$RES/$side"/host_*; do
    [ -d "$d" ] || continue
    printf "%-7s %-34s " "$side" "$(basename "$d" | sed 's/^host_//')"
    grep -E "^timing: nsteps" "$d/run.log" | awk '{printf "s/step %s  ", $NF}'
    grep -E "proj_timing: (sweep|apply|setup) " "$d/run.log" \
      | awk '{printf "%s=%s ", $2, $(NF-2)}'
    echo
  done
done
echo "=== rdenom A/B job finished ==="
