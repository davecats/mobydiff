#!/bin/bash
#SBATCH --job-name=moby_rough_tbl
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
# The elected turbulentBoundaryLayer case (CaNS grid, CONS convection) with the
# PREDICTED outlet face, smooth and rough-walled, restarted from the smooth
# developed field of the 2026-10-02 campaign (cans_match_run/cons_p2_675000.h5,
# a pre-step-7 global-3D file, readable on any rank count), 4 ranks per side,
# ONE NODE PER SIDE in the same allocation, statistics on.
#
#     export CODE_DIR=<worktree> RUN_DIR=<run dir> HDF5_ROOT=$HOME/hdf5 \
#            [NSTEPS=12500] [END_STEP=775000] [SIDES="smooth rough"]
#     cd $RUN_DIR && sbatch --export=ALL submit_rough_tbl.sh
#
#   smooth   tutorials/turbulentBoundaryLayer/production_stats.ini
#   rough    tutorials/turbulentBoundaryLayer/production_stats_rough.ini
#
# RESUMABLE AND SELF-CHAINING, in chunks of NSTEPS: each side restarts from its
# own latest snapshot (<side>_<step>.h5; the first chunk from
# $RUN_DIR/start_field.h5), its statistics file continues (raw sums + counts:
# a window is the difference of two chunk-end copies, kept as
# stats_<side>_<step>.h5), and a chunk writes ONE snapshot, at its end. A
# chunk that hits the wall clock leaves no snapshot and is simply run again.
# The successor job is submitted at the START of each job (afterany on this
# one), so a wall-clock kill still chains; a job whose every side has reached
# END_STEP, or whose previous two chunks of one side failed, exits at once
# without a successor. NO_CHAIN=1 runs one chunk. Keep NSTEPS x
# seconds-per-step under ~45 minutes (0.18 s/step on 4 A100 -> 12500 steps).
# RUN_DIR holds start_field.h5 (or a symlink); the inis are generated here
# from the tutorial's.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
NSTEPS="${NSTEPS:-12500}"; END_STEP="${END_STEP:-775000}"; SIDES="${SIDES:-smooth rough}"
TUT="$CODE_DIR/tutorials/turbulentBoundaryLayer"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

cache="$CODE_DIR/build_gpu/CMakeCache.txt"
if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$CODE_DIR" "$cache"; then
    rm -rf "$CODE_DIR/build_gpu"
fi
if [ ! -x "$CODE_DIR/build_gpu/moby_solve" ]; then
    echo "=== building $CODE_DIR ($(git -C "$CODE_DIR" rev-parse --short HEAD))"
    ( cd "$CODE_DIR" && ./compile.sh gpu ) > "$RUN_DIR/build.log" 2>&1 || { echo "BUILD FAILED"; exit 1; }
fi

cd "$RUN_DIR"
[ -f start_field.h5 ] || { echo "no start_field.h5 in $RUN_DIR"; exit 1; }
mapfile -t NODES < <(scontrol show hostnames "$SLURM_JOB_NODELIST")
{
    echo "job    : ${SLURM_JOB_ID:-none}   $(date '+%F %T %Z')"
    echo "nodes  : ${NODES[*]}"
    echo "code   : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "chunk  : $NSTEPS steps, end $END_STEP, sides: $SIDES"
} | tee -a provenance.txt

latest_of() { ls -1 ${1}_[0-9]*.h5 2>/dev/null | sed 's/.*_\([0-9]*\)\.h5/\1 &/' | sort -n | tail -1 | cut -d' ' -f2; }
step_of() { echo "$1" | sed 's/.*_\([0-9]*\)\.h5/\1/'; }

# Chain decision, BEFORE running: done, or broken, or queue the successor.
todo=0
for side in $SIDES; do
    new=$(latest_of "$side"); step=$(step_of "${new:-x_0.h5}")
    [ "$step" -lt "$END_STEP" ] && todo=1
    if [ "$(grep "^exit .* $side " status.txt 2>/dev/null | tail -2 | grep -c '^exit [^0]')" -ge 2 ]; then
        echo "$side: the last two chunks failed -- stopping the chain"; exit 1
    fi
done
[ "$todo" = 1 ] || { echo "every side at END_STEP $END_STEP: nothing to do"; exit 0; }
if [ -z "${NO_CHAIN:-}" ]; then
    sbatch --dependency=afterany:$SLURM_JOB_ID --export=ALL \
        "$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange/submit_rough_tbl.sh" \
        | tee -a provenance.txt
fi

n=0
for side in $SIDES; do
    src="$TUT/production_stats.ini"; [ "$side" = rough ] && src="$TUT/production_stats_rough.ini"
    node="${NODES[$((n % ${#NODES[@]}))]}"; n=$((n + 1))
    latest=$(latest_of "$side"); restart="${latest:-start_field.h5}"
    [ "$(step_of "${latest:-x_0.h5}")" -ge "$END_STEP" ] && { echo "=== $side at END_STEP"; continue; }
    ini="${side}.ini"
    sed -e "s/^nsteps = .*/nsteps = $NSTEPS/" \
        -e "s/^field_interval = .*/field_interval = 0/" \
        -e "s/^field_prefix = .*/field_prefix = ${side}/" \
        -e "s/^runtime_file = .*/runtime_file = ${side}_runtime.txt/" \
        -e "s/^stats_file = .*/stats_file = ${side}_stats.h5/" \
        -e "s/^stats_write_interval = .*/stats_write_interval = $NSTEPS/" \
        -e "/^\[restart\]/,/^$/ s|^file = .*|file = $restart|" "$src" > "$ini"
    MPIRUN="mpirun -n 4 --host $node:4 --bind-to core --map-by numa"
    if [ ! -f "${side}.case.h5" ]; then
        $MPIRUN "$CODE_DIR/build_gpu/moby_prepare" "$ini" > "prep_$side.log" 2>&1 \
            || { echo "PREPARE FAILED $side"; tail -5 "prep_$side.log"; exit 1; }
    fi
    log="${side}_from_$(basename "$restart" .h5).log"
    echo "=== $side on $node: $NSTEPS steps from $restart"
    ( echo "start $side $restart $(date '+%F %T')" >> status.txt
      $MPIRUN "$CODE_DIR/build_gpu/moby_solve" "$ini" > "$log" 2>&1
      rc=$?
      echo "exit $rc $side $restart $(date '+%F %T')" >> status.txt
      new=$(latest_of "$side")
      [ -n "$new" ] && [ "$new" != "$latest" ] && cp "${side}_stats.h5" "stats_${side}_$(step_of "$new").h5"
    ) &
done
wait
tail -4 status.txt
for side in $SIDES; do tail -1 ${side}_runtime.txt 2>/dev/null; done

