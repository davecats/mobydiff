#!/bin/bash
#SBATCH --job-name=moby_lap
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
# Does packing the nine Laplacian coefficient arrays into three buy anything?
# It did NOT move step_momentum's register count (128 either way, checked off the
# queue), so the only hypothesis left is contiguity: the three stencil
# coefficients for one index now share a sector. This job decides whether to keep
# the change or revert it -- a change to the solver's most bit-exactness-critical
# kernel with no measured benefit should not stay.
#
#   REF_DIR is a worktree at 0b36b36 (the nine-array version).
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/lap_staged"; rm -rf "$STG"; mkdir -p "$STG"
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

RES="${RESDIR:-$RUN_DIR/results_lap}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  dirty $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

SPECS="refined_yp82_rect_jacobi:8:2 rect_jacobi:8:2 refined_yp82_rect_jacobi:4:1 rect_jacobi:4:1"
for side in before after; do
    exe="$REF"; [ "$side" = after ] && exe="$NEW"
    echo "############ 1$side: timeline ############"
    CFG_DIR="$SRC/../configs" TL_PASSES="A" SPECS="$SPECS" NSTEPS="${NSTEPS:-100}" \
        bash "$STG/run_timeline.sh" "$exe" "$RES/$side"
done

echo "############ 2: Pass G -- production cases, bit-exact ############"
REF="$REF" H5MAXDIFF="$CODE_DIR/tools/h5maxdiff" CFG_DIR="$SRC/../configs" \
    PASSES="G" NSTEPS_GATE=20 bash "$STG/run_exchange.sh" "$NEW" "$RES/gate"

echo "############ 3: the 7-case suite ############"
MOBY_ROOT="$CODE_DIR" NSTEPS=200 \
    bash "$STG/run_mapgate.sh" "$NEW" "$REF" "$RES/suite"

echo "############ 4: ncu -- did the contiguity help? ############"
CFG_DIR="$SRC/../configs" KERNEL_RE='step_momentum' \
    LAUNCH_SKIP=12 LAUNCH_COUNT=6 \
    bash "$STG/run_ncu.sh" "$NEW" "$RES/ncu_after"
python3 "$STG/collect_ncu.py" "$RES/ncu_after" > "$RES/ncu_after.md" 2>&1 \
    && { echo; cat "$RES/ncu_after.md"; } || echo "=== collector failed ==="

echo "=== summary: s/step and the projection buckets ==="
for side in before after; do
  for d in "$RES/$side"/host_*; do
    [ -d "$d" ] || continue
    printf "%-7s %-34s " "$side" "$(basename "$d" | sed 's/^host_//')"
    grep -E "^timing: nsteps" "$d/run.log" | awk '{printf "s/step %s  ", $NF}'
    grep -E "step_timing: (bodyforce|ibm_mu|momentum) " "$d/run.log" \
      | awk '{printf "%s=%s ", $2, $(NF-2)}'
    echo
  done
done
echo "=== lap A/B job finished ==="
