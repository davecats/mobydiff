#!/bin/bash
#SBATCH --job-name=moby_profile2
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=01:00:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The fresh profile of the production step for the second cleaning and
# profiling pass (docs/next_session_cleanup_pass2.md section 2). ONE column,
# the current head: every published ratio has a moved denominator since
# 2026-09-14, so nothing is compared here -- this is the baseline the pass
# works from. The last two baselines are results_horeka_2026-09-25.md (three
# columns, six configs) and results_horeka_2026-10-01.md (the block order).
#
# Three measurements, one allocation, in this order (the matrix first: if the
# wall clock bites, the resumable part survives):
#   1. the campaign matrix, six configs at 1/2/4/8 ranks (rough_jacobi runs
#      again since 2026-10-05), plus tutorials/min_channel at 1 and 4 ranks
#      (MATRIX_EXTRA), 200 steps, [output] profile = true -> the phase tables;
#   2. the nsys per-kernel timeline at 1 rank on four configs (30 steps):
#      device time per kernel per step, collect_kernels.py;
#   3. ncu on EVERY kernel of one step, two configs (3 steps, the first 700
#      matching launches): registers, occupancy, DRAM/SM utilisation --
#      launch__registers_per_thread and sm__warps_active, never a step time.
# 16 ranks: RANKS_SMALL=16 RANKS_BIG=16 PASSES=matrix with --nodes=4
# --partition=accelerated into the same RESDIR (run_matrix.sh skips done runs).
set -uo pipefail
RANKS_SMALL="${RANKS_SMALL:-1,2,4,8}"
RANKS_BIG="${RANKS_BIG:-4,8}"
PASSES="${PASSES:-matrix nsys ncu}"
export RANKS_SMALL RANKS_BIG

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka"

# Staged copies: bash reads a script by byte offset, an edit mid-run splices
# two versions together (job 5139461).
STG="$RUN_DIR/profile2_staged"
rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC/run_matrix.sh" "$SRC/mechanism/collect_scaling.py" "$SRC/collect_profile.py" \
   "$SRC/exchange/run_timeline.sh" "$SRC/exchange/collect_kernels.py" \
   "$SRC/exchange/run_ncu.sh" "$SRC/exchange/collect_ncu.py" "$STG/"
cp -r "$SRC/configs" "$STG/configs"
# The tutorial channel as a seventh config: no snapshots (the matrix only
# needs timings), the rest as shipped.
sed -e 's/^field_interval *=.*/field_interval = 0/' \
    "$CODE_DIR/tutorials/min_channel/input_gpu.ini" > "$STG/configs/min_channel.ini"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

cache="$CODE_DIR/build_gpu/CMakeCache.txt"
if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$CODE_DIR" "$cache"; then
    echo "=== $CODE_DIR: build_gpu cache is for a different source dir -- discarding"
    rm -rf "$CODE_DIR/build_gpu"
fi
echo "=== building $CODE_DIR ($(git -C "$CODE_DIR" rev-parse --short HEAD))"
( cd "$CODE_DIR" && ./compile.sh gpu ) || exit 1
NEW="$CODE_DIR/build_gpu/moby_solve"
[ -x "$NEW" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_profile2}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "head   : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps : ${NSTEPS:-200} (matrix), 30 (nsys), 3 (ncu)"
    echo "layout : --map-by numa --bind-to core (as every matrix)"
    echo "passes : $PASSES"
} | tee -a "$RES/provenance.txt"

case " $PASSES " in *" matrix "*)
    echo "############ matrix: six configs + min_channel ############"
    CFG_DIR="$STG/configs" NSTEPS="${NSTEPS:-200}" MATRIX_EXTRA="min_channel:1,4" \
        bash "$STG/run_matrix.sh" "$NEW" "$RES/new"
    python3 "$STG/collect_scaling.py" "$RES/new" > "$RES/scaling.md" 2>&1 \
        && { echo; cat "$RES/scaling.md"; } || echo "=== scaling collector failed ==="
    python3 "$STG/collect_profile.py" "$RES/new" > "$RES/phases.md" 2>&1 \
        && { echo; cat "$RES/phases.md"; } || echo "=== phase collector failed (see phases.md) ==="
    echo "=== s/step, every run (from the timing: line) ==="
    for d in "$RES/new"/*/; do
        [ -f "$d/run.log" ] || continue
        printf "  %-46s %s\n" "$(basename "$d")" "$(awk '/^timing:/ {v=$NF} END {print v}' "$d/run.log")"
    done
;; esac

case " $PASSES " in *" nsys "*)
    echo "############ nsys per-kernel timeline, 1 rank ############"
    CFG_DIR="$STG/configs" NSTEPS=30 TL_PASSES="B" \
        SPECS="refined_yp82_rect_jacobi:1:1 rough_jacobi:1:1 base_jacobi:1:1 min_channel:1:1 refined_yp82_rect_redblack:1:1" \
        bash "$STG/run_timeline.sh" "$NEW" "$RES/nsys"
    python3 "$STG/collect_kernels.py" "$RES/nsys" > "$RES/kernels.md" 2>&1 \
        && { echo; cat "$RES/kernels.md"; } || echo "=== kernel collector failed ==="
;; esac

case " $PASSES " in *" ncu "*)
    echo "############ ncu, every kernel of one step ############"
    # Every nvfortran kernel matches `nvkernel`; skip the init launches by
    # count (the first step holds them), profile the following 700 launches
    # (a refined Chebyshev step is ~450 launches: 3 substages x (6 x 3
    # projection kernels + exchange rounds), a red-black one more).
    CFG_DIR="$STG/configs" NSTEPS=3 KERNEL_RE="nvkernel" LAUNCH_SKIP=450 LAUNCH_COUNT=700 \
        CONFIGS="refined_yp82_rect_jacobi rough_jacobi" \
        bash "$STG/run_ncu.sh" "$NEW" "$RES/ncu"
    python3 "$STG/collect_ncu.py" "$RES/ncu" > "$RES/ncu.md" 2>&1 \
        && { echo; cat "$RES/ncu.md"; } || echo "=== ncu collector failed; CSVs intact ==="
;; esac
echo "=== profile2 finished $(date '+%F %T') ==="
