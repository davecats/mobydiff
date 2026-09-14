#!/bin/bash
#SBATCH --job-name=moby_apply
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:55:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# A/B for "hoist the face-correction metric out of jacobi_apply's face kernel":
# does k2's register count fall below the occupancy thresholds, does the apply
# bucket shrink, and are the fields bit-exact?
#
#   REF_DIR is a worktree at the PRE-change commit, CODE_DIR the working tree.
#   Both are built here, on this node, with the same module and CMake cache.
#
# Three things in one allocation:
#   0  registers -- cuobjdump on both binaries, the primary evidence
#   1  timeline pass A (untraced brackets + s/step), BEFORE and AFTER, 4 ranks
#   2  Pass G   -- production cases, ref vs new, max_abs 0
#   3  mapgate  -- the 7-case suite against the SAME reference binary
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/apply_staged"
rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC"/run_timeline.sh "$SRC"/run_exchange.sh "$SRC"/run_mapgate.sh \
   "$SRC"/collect_timeline.py "$SRC"/collect_exchange.py "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

for d in "$REF_DIR" "$CODE_DIR"; do
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
NEW="$CODE_DIR/build_gpu/moby_solve"
REF="$REF_DIR/build_gpu/moby_solve"
[ -x "$NEW" ] && [ -x "$REF" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_apply}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  dirty $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

echo "############ 0: registers (cuobjdump, both binaries) ############"
# Built here on an A100 node, so this is the arch the runs below actually use.
# STACK/LOCAL must stay 0: a register cut that spills is a loss.
for side in ref new; do
    exe="$REF"; [ "$side" = new ] && exe="$NEW"
    echo "--- $side"
    cuobjdump -res-usage "$exe" | awk '
      /^ *Function nvkernel_/ { name=$2; sub(/:$/,"",name); next }
      /REG:/ && name != "" { print name"  "$1" "$2" "$3" "$4; name="" }' |
      grep -E 'pressure_solver|step_momentum' | sort
done 2>&1 | tee "$RES/registers.txt"

echo "############ 1a: timeline BEFORE ############"
CFG_DIR="$SRC/../configs" TL_PASSES="A" NSTEPS="${NSTEPS:-100}" \
    bash "$STG/run_timeline.sh" "$REF" "$RES/before"
echo "############ 1b: timeline AFTER ############"
CFG_DIR="$SRC/../configs" TL_PASSES="A" NSTEPS="${NSTEPS:-100}" \
    bash "$STG/run_timeline.sh" "$NEW" "$RES/after"

echo "############ 2: Pass G -- production cases, bit-exact ############"
REF="$REF" H5MAXDIFF="$CODE_DIR/tools/h5maxdiff" CFG_DIR="$SRC/../configs" \
    PASSES="G" NSTEPS_GATE=20 bash "$STG/run_exchange.sh" "$NEW" "$RES/gate"

echo "############ 3: the 7-case suite ############"
MOBY_ROOT="$CODE_DIR" NSTEPS=200 \
    bash "$STG/run_mapgate.sh" "$NEW" "$REF" "$RES/suite"

echo "=== summary: s/step and the projection buckets ==="
for side in before after; do
  for d in "$RES/$side"/host_*; do
    [ -d "$d" ] || continue
    printf "%-7s %-34s " "$side" "$(basename "$d" | sed 's/^host_//')"
    grep -E "^timing: nsteps" "$d/run.log" | awk '{printf "s/step %s  ", $NF}'
    grep -E "proj_timing: (sweep|apply) " "$d/run.log" \
      | awk '{printf "%s=%s ", $2, $(NF-2)}'
    echo
  done
done
echo "=== apply A/B job finished ==="
