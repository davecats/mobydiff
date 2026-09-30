#!/bin/bash
#SBATCH --job-name=moby_predictor
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:50:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The predictor register fix (PREREGISTERED_predictor.md): the consolidated
# head's fused momentum kernel lost 1.3-2.0 % of the step when the S3 lockdown
# removed the `if (skew)` guards and the compiler cut it from 128 to 100
# registers without gaining occupancy (results_horeka_2026-09-26.md, ncu job
# 5163917). The guards are back around a runtime-true logical.
#
#   ref  REF_DIR   the head without the guard (0811811 or later on main)
#   new  CODE_DIR  the guard
#
# Two measurements, one allocation. (1) ncu on the predictor kernel, both
# binaries on ONE node, so the register count, occupancy and the two
# utilisation axes are read on the same silicon as the 09-26 table. (2) The
# step at 4 and 8 ranks on the three configs the 09-26 phase table placed the
# regression on, both binaries, shipped configs on both sides (both heads
# reject `[flow] convection`; nothing is injected).
#
# NOT here: bit-exactness. A guard moves no expression but may move FMA
# contraction, so the equivalence is the nofma suite, run on the workstation
# before this is submitted (the commit message records the result).
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka"
STG="$RUN_DIR/predictor_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC/run_matrix.sh" "$SRC/mechanism/collect_scaling.py" \
   "$SRC/exchange/run_ncu.sh" "$SRC/exchange/collect_ncu.py" "$STG/"
cp -r "$SRC/configs" "$STG/configs"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

for d in "$REF_DIR" "$CODE_DIR"; do
    cache="$d/build_gpu/CMakeCache.txt"
    if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$d" "$cache"; then
        echo "=== $d: build_gpu cache is for a different source dir -- discarding"
        rm -rf "$d/build_gpu"
    fi
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
NEW="$CODE_DIR/build_gpu/moby_solve"; REF="$REF_DIR/build_gpu/moby_solve"
[ -x "$NEW" ] && [ -x "$REF" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_predictor}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  (guard)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)  (no guard)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps : ${NSTEPS:-200}"
    echo "layout : --map-by numa --bind-to core (as every matrix)"
} | tee "$RES/provenance.txt"

# (1) ncu, the predictor's two kernels: 6 matching launches per step (2 kernels
# x 3 substages), so skip two steps and profile exactly one.
for side in ref new; do
    exe="$REF"; [ "$side" = new ] && exe="$NEW"
    echo "############ ncu $side ############"
    CFG_DIR="$STG/configs" KERNEL_RE="step_momentum" LAUNCH_SKIP=12 LAUNCH_COUNT=6 \
        bash "$STG/run_ncu.sh" "$exe" "$RES/ncu_$side"
    python3 "$STG/collect_ncu.py" "$RES/ncu_$side" > "$RES/ncu_$side.md" 2>&1 \
        && { echo; cat "$RES/ncu_$side.md"; } || echo "=== collector failed for $side; CSVs intact ==="
done

# (2) the step, both columns, new first.
export RANKS_SMALL="${RANKS_SMALL:-4,8}" RANKS_BIG="${RANKS_BIG:-4,8}"
CFGS="${CONFIGS:-base_jacobi rect_jacobi refined_yp82_rect_jacobi}"
echo "############ NEW (guard) ############"
CFG_DIR="$STG/configs" CONFIGS="$CFGS" NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$NEW" "$RES/new"
echo "############ REF (no guard) ############"
CFG_DIR="$STG/configs" CONFIGS="$CFGS" NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$REF" "$RES/ref"

python3 "$STG/collect_scaling.py" "$RES/ref" "$RES/new" > "$RES/scaling_predictor.md" 2>&1 \
    && { echo; cat "$RES/scaling_predictor.md"; } || echo "=== collector failed ==="

echo "=== s/step and momentum bucket, every run ==="
for d in "$RES/ref"/*/; do
    n=$(basename "$d"); o="$RES/new/$n"
    [ -f "$d/run.log" ] && [ -f "$o/run.log" ] || continue
    a=$(awk '/^timing:/ {v=$NF} END {print v}' "$d/run.log")
    b=$(awk '/^timing:/ {v=$NF} END {print v}' "$o/run.log")
    ma=$(grep -E "^step_timing: momentum " "$d/run.log" | awk '{print $(NF-2)}')
    mb=$(grep -E "^step_timing: momentum " "$o/run.log" | awk '{print $(NF-2)}')
    python3 -c "
a=$a; b=$b
print(f'  {\"$n\":32s} step {1e3*a:8.3f} -> {1e3*b:8.3f} ms ({100*(b-a)/a:+6.2f}%)   momentum ${ma:-?} -> ${mb:-?}')"
done
echo "=== predictor job finished $(date '+%F %T') ==="
