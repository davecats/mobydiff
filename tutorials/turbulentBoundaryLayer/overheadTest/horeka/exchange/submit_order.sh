#!/bin/bash
#SBATCH --job-name=moby_order
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
# The minimum-surface block order (PREREGISTERED_order.md): since step 7 an
# nb-less case is one block per rank in Morton order with z on top, which put
# base_jacobi's largest faces on the node boundary at 8 ranks (+70 % on the
# step, results_horeka_2026-09-30.md section 3). The case builder now chooses
# the key's bit order from the geometry and records it in the case file.
#
#   ref  REF_DIR   the head before the change (5bdc5eb): legacy order
#   new  CODE_DIR  the order commit
#
# Each column prepares its own case file in every run directory
# (run_matrix.sh, tools/prepare_if_missing.sh beside the binary), so each runs
# its own order from the SAME shipped config, with the DEFAULT placement.
#
# Two measurements, one allocation. (1) The field gate: 20 steps at 8 ranks,
# both binaries, final snapshot through h5maxdiff -- which matches rows on
# (origin, level), since the two columns write their blocks in different row
# orders. Production flags: no expression moves, only ownership. (2) The step
# at 4 and 8 ranks, both columns, new first.
#
# 16 ranks: override RANKS_SMALL=16 RANKS_BIG=16 GATE=0 with --nodes=4
# --partition=accelerated into the same RESDIR (run_matrix.sh skips finished
# runs).
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka"
STG="$RUN_DIR/order_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC/run_matrix.sh" "$SRC/mechanism/collect_scaling.py" "$STG/"
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

# ALWAYS rebuilt from the NEW tree: the row matching is new, and an h5maxdiff
# left over from an earlier job compares rows by position (it would report the
# two orders as different fields). mpicc: this HDF5 is the parallel build.
H5MAXDIFF="$STG/h5maxdiff"
mpicc -O2 -o "$H5MAXDIFF" "$CODE_DIR/tools/h5maxdiff.c" \
    -I"$HDF5_ROOT/include" -L"$HDF5_ROOT/lib" -lhdf5 -Wl,-rpath,"$HDF5_ROOT/lib" \
    || { echo "ERROR: could not build h5maxdiff" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_order}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  (minimum-surface order)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)  (legacy order)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps : ${NSTEPS:-200}"
    echo "layout : --map-by numa --bind-to core (the default placement of every matrix)"
} | tee "$RES/provenance.txt"

# (1) The field gate. A run directory per side, the shipped config with a
# snapshot at the last step; the snapshots (4.4 GB each) are deleted after the
# comparison, the comparison's output is kept.
if [ "${GATE:-1}" = 1 ]; then
    NG="${NSTEPS_GATE:-20}"
    for cfg in ${GATE_CONFIGS:-base_jacobi rect_jacobi}; do
        for side in ref new; do
            exe="$REF"; [ "$side" = new ] && exe="$NEW"
            run="$RES/gate_${cfg}_${side}"
            [ -f "$run/run.log" ] && continue
            rm -rf "$run"; mkdir -p "$run"
            sed -e "s/^nsteps *=.*/nsteps = $NG/" \
                -e "s/^field_interval *=.*/field_interval = $NG/" \
                -e "s/^runtime_interval *=.*/runtime_interval = $NG/" \
                "$STG/configs/$cfg.ini" > "$run/config.ini"
            echo "=== gate $cfg $side ($(date '+%F %T'))"
            pim="$(dirname "$exe")/../tools/prepare_if_missing.sh"
            ( cd "$run" && "$pim" 8 "$exe" config.ini prepare.log ) \
                || { echo "    PREPARE FAILED -- see $run/prepare.log"; continue; }
            ( cd "$run" && mpirun -n 8 --bind-to core --map-by numa "$exe" config.ini > run.log 2>&1 ) \
                || { echo "    FAILED -- see $run/run.log"; mv "$run/run.log" "$run/run.FAILED.log"; }
        done
        a=$(ls "$RES/gate_${cfg}_ref"/overhead_*.h5 2>/dev/null | grep -v case | tail -1)
        b=$(ls "$RES/gate_${cfg}_new"/overhead_*.h5 2>/dev/null | grep -v case | tail -1)
        if [ -n "$a" ] && [ -n "$b" ]; then
            echo "--- h5maxdiff $cfg ($a vs $b)"
            "$H5MAXDIFF" "$a" "$b" | tee "$RES/gate_${cfg}.txt"
            rm -f "$RES/gate_${cfg}_ref"/overhead_[0-9]*.h5 "$RES/gate_${cfg}_new"/overhead_[0-9]*.h5
        else
            echo "--- gate $cfg: NO SNAPSHOT PAIR (a='$a' b='$b')" | tee "$RES/gate_${cfg}.txt"
        fi
    done
fi

# (2) The step, both columns, new first.
export RANKS_SMALL="${RANKS_SMALL:-4,8}" RANKS_BIG="${RANKS_BIG:-4,8}"
CFGS="${CONFIGS:-base_jacobi rect_jacobi refined_yp82_rect_jacobi}"
echo "############ NEW (minimum-surface order) ############"
CFG_DIR="$STG/configs" CONFIGS="$CFGS" NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$NEW" "$RES/new"
echo "############ REF (legacy order) ############"
CFG_DIR="$STG/configs" CONFIGS="$CFGS" NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$REF" "$RES/ref"

python3 "$STG/collect_scaling.py" "$RES/ref" "$RES/new" > "$RES/scaling_order.md" 2>&1 \
    && { echo; cat "$RES/scaling_order.md"; } || echo "=== collector failed ==="

echo "=== step, mpi_wait and the partition line, every run ==="
for d in "$RES/ref"/*/; do
    n=$(basename "$d"); o="$RES/new/$n"
    [ -f "$d/run.log" ] && [ -f "$o/run.log" ] || continue
    a=$(awk '/^timing:/ {v=$NF} END {print v}' "$d/run.log")
    b=$(awk '/^timing:/ {v=$NF} END {print v}' "$o/run.log")
    wa=$(grep -E "^exch_timing: mpi_wait " "$d/run.log" | sed 's/.*seconds_per_step *\([0-9.Ee+-]*\).*/\1/')
    wb=$(grep -E "^exch_timing: mpi_wait " "$o/run.log" | sed 's/.*seconds_per_step *\([0-9.Ee+-]*\).*/\1/')
    python3 -c "
a=$a; b=$b; wa=float('${wa:-nan}'); wb=float('${wb:-nan}')
print(f'  {\"$n\":32s} step {1e3*a:8.3f} -> {1e3*b:8.3f} ms ({100*(b-a)/a:+6.2f}%)   mpi_wait {1e3*wa:7.3f} -> {1e3*wb:7.3f} ms')"
    grep -h "partition: face cells" "$o/run.log" | sed 's/^/      new:/'
done
echo "=== order job finished $(date '+%F %T') ==="
