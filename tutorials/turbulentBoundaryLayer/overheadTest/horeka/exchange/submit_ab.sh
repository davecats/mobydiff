#!/bin/bash
#SBATCH --job-name=moby_ab
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
# A/B for "map the comm_type parent object once": does the 7952-byte per-launch
# descriptor blob disappear, does the per-launch bracket fall, and are the fields
# bit-exact?
#
#   REF_DIR is a worktree at the PRE-change commit, CODE_DIR the working tree.
#   Both are built here, on this node, with the same module and CMake cache.
#
# Three things in one allocation:
#   1  timeline, BEFORE and AFTER -- untraced brackets + nsys launch traffic
#   2  Pass G   -- production cases, ref vs new, max_abs 0
#   3  mapgate  -- the 7-case suite against the SAME reference binary
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
EXCH="$RUN_DIR/ab_staged"
rm -rf "$EXCH"; mkdir -p "$EXCH"
cp "$SRC"/run_timeline.sh "$SRC"/run_exchange.sh "$SRC"/run_mapgate.sh \
   "$SRC"/collect_timeline.py "$SRC"/collect_exchange.py \
   "$SRC"/analyse_launch_traffic.py "$EXCH/"

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

RES="${RESDIR:-$RUN_DIR/results_ab}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  dirty $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

echo "############ 1a: timeline BEFORE ############"
CFG_DIR="$SRC/../configs" NSTEPS="${NSTEPS:-30}" \
    bash "$EXCH/run_timeline.sh" "$REF" "$RES/before"
echo "############ 1b: timeline AFTER ############"
CFG_DIR="$SRC/../configs" NSTEPS="${NSTEPS:-30}" \
    bash "$EXCH/run_timeline.sh" "$NEW" "$RES/after"

echo "############ 2: Pass G -- production cases, bit-exact ############"
REF="$REF" H5MAXDIFF="$CODE_DIR/tools/h5maxdiff" CFG_DIR="$SRC/../configs" \
    PASSES="G" NSTEPS_GATE=20 bash "$EXCH/run_exchange.sh" "$NEW" "$RES/gate"

echo "############ 3: the 7-case suite ############"
MOBY_ROOT="$CODE_DIR" NSTEPS=200 \
    bash "$EXCH/run_mapgate.sh" "$NEW" "$REF" "$RES/suite"

for side in before after; do
    python3 "$EXCH/collect_timeline.py" "$RES/$side" > "$RES/timeline_$side.md" 2>&1 || true
    python3 "$EXCH/analyse_launch_traffic.py" "$RES/$side"/nsys_*/rep_0.sqlite \
        > "$RES/traffic_$side.md" 2>&1 || true
done
echo "=== A/B job finished ==="
