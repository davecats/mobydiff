#!/bin/bash
#SBATCH --job-name=moby_matrix4
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=04:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# --no-requeue is deliberate. On 2026-09-19 job 5150030 hit
# "user env retrieval failed" on the node it was allocated, and Slurm requeued
# it AND HELD it at priority 0 -- where it sat, invisibly, for two days, because
# a held job looks exactly like a queued one in squeue. run_matrix.sh is
# resumable (it skips any run whose run.log exists), so a visible failure and a
# manual resubmit is strictly better than a silent hold.
#
# The rdenom work and the rough-wall benchmark, across the full rank sweep.
#
# results_horeka_2026-09-23.md section 4 named two gaps: that matrix predates
# the rdenom narrowing, and it has no case with a BODY at all. This closes both.
#
#   ref  b9414bd   the published campaign point. Runs the 5 body-free configs
#                  ONLY -- see CFG_REF below, this is load-bearing.
#   mid  809759e   HEAD with the rdenom narrowing alone disabled
#                  (branch bench/rdenom-always). Runs all 6.
#   new  8fa0fc2   HEAD. Runs all 6.
#
# Pairings: mid->new isolates the narrowing (on body-free AND body cases);
# ref->mid is a CONTROL for the roughness feature being dormant plus the cost of
# the rdenomBlocks indirection that mid carries and never benefits from.
#
# WHY ref IS GATED. config.f90 has no `case default`: an unknown key in a known
# section is SILENTLY IGNORED. b9414bd therefore runs rough_jacobi with
# `wall_shape = eggcarton` discarded, falling back to the hardcoded 2D wavy
# wall, and reports success. It must never be handed that config.
#
# Readings pre-registered in PREREGISTERED_matrix4.md.
#
# NEW then MID then REF: if the wall clock bites, mid->new (the increment under
# test, on every config) survives and only the continuity column is lost.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"; MID_DIR="${MID_DIR:?}"
# The body-free five. b9414bd predates [ibm] wall_shape and must get only these.
CFG_REF="base_jacobi rect_jacobi refined_yp82_rect_jacobi refined_yp82_rect_redblack refined_big_rect_jacobi"
# What each column IS varies by run; the collector writes scaling_<pair>.md
# named from these, so a table can never be quoted against the wrong question.
LBL_REF="${LBL_REF:-ref}"; LBL_MID="${LBL_MID:-mid}"; LBL_NEW="${LBL_NEW:-new}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka"

# Stage the driver AND the configs it resolves relative to itself, then run the
# copy -- bash reads a script by byte offset and an edit mid-run splices two
# versions together (job 5139461).
STG="$RUN_DIR/matrix4_staged"
rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC/run_matrix.sh" "$SRC/mechanism/collect_scaling.py" "$STG/"
cp -r "$SRC/configs" "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

# A CMake cache records the ABSOLUTE source path, so a worktree that was
# renamed (git worktree move) carries a cache pointing at its old name and
# cmake refuses to configure -- which killed job 5145816 after it had waited a
# day in the queue. Check and wipe rather than trust.
for d in "$REF_DIR" "$MID_DIR" "$CODE_DIR"; do
    cache="$d/build_gpu/CMakeCache.txt"
    if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$d" "$cache"; then
        echo "=== $d: build_gpu cache was configured for a different source dir -- discarding"
        grep -m1 "^CMAKE_HOME_DIRECTORY:" "$cache" | sed 's/^/    /'
        rm -rf "$d/build_gpu"
    fi
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
NEW="$CODE_DIR/build_gpu/moby_solve"; REF="$REF_DIR/build_gpu/moby_solve"
MID="$MID_DIR/build_gpu/moby_solve"
[ -x "$NEW" ] && [ -x "$REF" ] && [ -x "$MID" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_matrix4}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  ($LBL_NEW)"
    echo "mid    : $(git -C "$MID_DIR" rev-parse HEAD)  ($LBL_MID)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)  ($LBL_REF)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps : ${NSTEPS:-200}"
    echo "layout : --map-by numa --bind-to core (as every earlier matrix)"
} | tee "$RES/provenance.txt"

echo "############ NEW ($LBL_NEW) -- all 6 configs ############"
NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$NEW" "$RES/new"
echo "############ MID ($LBL_MID) -- all 6 configs ############"
NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$MID" "$RES/mid"
echo "############ REF ($LBL_REF) -- body-free five ONLY ############"
CONFIGS="$CFG_REF" NSTEPS="${NSTEPS:-200}" bash "$STG/run_matrix.sh" "$REF" "$RES/ref"

# Three pairings: the total, and each step of the way on its own.
for pair in "mid new rdenom" "ref mid roughness_dormant_control" "ref new total"; do
    set -- $pair
    python3 "$STG/collect_scaling.py" "$RES/$1" "$RES/$2" > "$RES/scaling_$3.md" 2>&1 \
        && { echo; echo "===== $3 ($1 -> $2)"; cat "$RES/scaling_$3.md"; } \
        || echo "=== collector failed for $3; raw logs intact ==="
done
# The roughness cost, like for like: rough_jacobi against its body-free twin.
echo "=== Q4: what the roughness costs (same grid, blocks and flow as rect) ==="
for col in new mid; do for n in 1 2 4 8 16; do
  a=$(grep -hE "^timing: nsteps" "$RES/$col/rect_jacobi_n$n/run.log" 2>/dev/null | awk '{print $NF}')
  b=$(grep -hE "^timing: nsteps" "$RES/$col/rough_jacobi_n$n/run.log" 2>/dev/null | awk '{print $NF}')
  [ -n "$a" ] && [ -n "$b" ] && python3 -c "
a=$a; b=$b
print(f'  $col  n=$n  rect {1e3*a:8.3f}  rough {1e3*b:8.3f} ms  ->  {100*(b-a)/a:+6.2f}%')"
done; done
grep -h "rdenom recomputed" "$RES/new"/rough_jacobi_n*/run.log 2>/dev/null | sort -u | sed 's/^/  /'
echo "=== matrix4 job finished ==="
