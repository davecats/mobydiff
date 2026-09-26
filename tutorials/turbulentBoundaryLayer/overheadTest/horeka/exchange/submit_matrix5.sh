#!/bin/bash
#SBATCH --job-name=moby_matrix5
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
# The campaign matrix on the CONSOLIDATED head: everything that landed after
# the last published matrix (job 5159149, results_horeka_2026-09-25.md) --
# the boundaryLayer and scalar merges, and the jacobi-interface feature port.
#
# WHY IT HAD TO BE RE-RUN, and it is not only the usual moved denominator.
# Every cross-config ratio the campaign exists to report -- block tax, strong
# scaling, the 2:1 machinery's coarse-cell-equivalents, red-black against
# Jacobi -- is measured BETWEEN configs at one head, so all of them are stale
# the moment the step time moves. Two things moved it here: the port hardwired
# skew-symmetric convection, which the branch measured at ~+3.4 % s/step, and
# the scalar merge put new (dormant) code on the step path.
#
#   ref  8fa0fc2   the head the published campaign measured
#   new  HEAD      the consolidated head
#
# THE REF COLUMN NEEDS `[flow] convection = skew` INJECTED, and this is
# load-bearing. The shipped configs USED to carry that key; the S3 lockdown made
# it an error, so the key was stripped from them -- which means a pre-lockdown
# binary handed the current configs would silently run DIVERGENCE-form
# convection and the comparison would measure the wrong thing. (config.f90 has
# no `case default`, so the reverse trap is just as quiet: the old binary would
# not complain about a missing key either.) The ref column therefore runs from a
# staged config copy with the key appended, which is the configuration 8fa0fc2
# was published with.
#
# NEW FIRST: if the wall clock bites, the column that describes today's solver
# survives and only the continuity comparison is lost.
#
# SPLIT ACROSS PARTITIONS, because `accelerated` was backlogged to 09-28/09-29
# when this was submitted while `dev_accelerated` (2 nodes, 1 h) schedules
# same-day. 2 nodes = 8 GPUs covers rank counts 1/2/4/8, which is where every
# per-cell effect lives; the 16-rank column needs 4 nodes and is a separate
# submission with RANKS_SMALL=16 RANKS_BIG=16 into the same results directory,
# which works because run_matrix.sh skips any run whose run.log already exists.
# Override RANKS_SMALL / RANKS_BIG to choose.
#
# RESUMABLE: run_matrix.sh skips any run whose run.log exists.
set -uo pipefail
RANKS_SMALL="${RANKS_SMALL:-1,2,4,8}"
RANKS_BIG="${RANKS_BIG:-4,8}"
export RANKS_SMALL RANKS_BIG

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
LBL_REF="${LBL_REF:-ref}"; LBL_NEW="${LBL_NEW:-new}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka"

# Stage the driver AND the configs it resolves relative to itself, then run the
# copy -- bash reads a script by byte offset and an edit mid-run splices two
# versions together (job 5139461).
STG="$RUN_DIR/matrix5_staged"
rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC/run_matrix.sh" "$SRC/mechanism/collect_scaling.py" "$STG/"
cp -r "$SRC/configs" "$STG/configs"
cp -r "$SRC/configs" "$STG/configs_skew"
for f in "$STG"/configs_skew/*.ini; do printf '\n[flow]\nconvection = skew\n' >> "$f"; done

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

# A CMake cache records the ABSOLUTE source path, so a worktree that was renamed
# carries a cache pointing at its old name and cmake refuses to configure --
# which killed job 5145816 after a day in the queue. Check and wipe.
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

RES="${RESDIR:-$RUN_DIR/results_matrix5}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  ($LBL_NEW, configs as shipped)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)  ($LBL_REF, configs + convection = skew)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps : ${NSTEPS:-200}"
    echo "layout : --map-by numa --bind-to core (as every earlier matrix)"
} | tee "$RES/provenance.txt"

echo "############ NEW ($LBL_NEW) -- the consolidated head, all 6 configs ############"
CFG_DIR="$STG/configs" NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$NEW" "$RES/new"
echo "############ REF ($LBL_REF) -- 8fa0fc2, configs WITH convection = skew ############"
CFG_DIR="$STG/configs_skew" NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$REF" "$RES/ref"

python3 "$STG/collect_scaling.py" "$RES/ref" "$RES/new" > "$RES/scaling_consolidation.md" 2>&1 \
    && { echo; cat "$RES/scaling_consolidation.md"; } || echo "=== collector failed ==="

echo "=== s/step, every run ==="
for col in new ref; do
    for d in "$RES/$col"/*/; do
        [ -f "$d/run.log" ] || continue
        printf "  %-8s %-46s %s\n" "$col" "$(basename "$d")" \
            "$(grep -oE 'seconds_per_step +[0-9.E+-]+' "$d/run.log" | tail -1 | awk '{print $2}')"
    done
done
echo "=== matrix5 finished $(date '+%F %T') ==="
