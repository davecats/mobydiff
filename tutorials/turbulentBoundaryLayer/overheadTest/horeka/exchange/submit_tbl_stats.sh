#!/bin/bash
#SBATCH --job-name=tbl_stats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=36:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# Phase 2 of the turbulent boundary-layer tutorial (statistics, the committed
# production_stats.ini: 500000 steps = 10000 time units) with the present
# code, to regenerate assets/mobydiff/.../data.nc after the outlet change of
# 2026-10-01 (README "The outlet zone"). About 27 h on 4 A100 (0.194 s/step).
#
#     export CODE_DIR=<worktree> RUN_DIR=<run dir> HDF5_ROOT=$HOME/hdf5
#     cd $RUN_DIR && sbatch --export=ALL submit_tbl_stats.sh
#
# RUN_DIR holds production_stats.ini (its [restart] file is rewritten here)
# and start_field.h5, a field settled with the present code (a symlink to
# outlet_tbl_run/tbl_new_1087500.h5: 750 time units with the predicted outlet
# from the tutorial's developed field).
# RESUMABLE: submit again and it continues from the latest production_p2_*.h5
# with the remaining steps; the statistics file continues by itself.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
TOTAL="${TOTAL:-500000}"          # steps of the whole phase
START_STEP="${START_STEP:-1087500}"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

cache="$CODE_DIR/build_gpu/CMakeCache.txt"
if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$CODE_DIR" "$cache"; then
    rm -rf "$CODE_DIR/build_gpu"
fi
echo "=== building $CODE_DIR ($(git -C "$CODE_DIR" rev-parse --short HEAD))"
( cd "$CODE_DIR" && ./compile.sh gpu ) > "$RUN_DIR/build.log" 2>&1 || { echo "BUILD FAILED"; exit 1; }

cd "$RUN_DIR"
latest=$(ls -1 production_p2_[0-9]*.h5 2>/dev/null | sort -t_ -k3 -n | tail -1)
restart="${latest:-start_field.h5}"
done_steps=0
if [ -n "$latest" ]; then
    step=${latest#production_p2_}; step=${step%.h5}
    done_steps=$((step - START_STEP))
fi
left=$((TOTAL - done_steps))
[ "$left" -gt 0 ] || { echo "phase complete: $done_steps steps"; exit 0; }
sed -i -e "s/^nsteps = .*/nsteps = $left/" \
       -e "/^\[restart\]/,/^$/ s|^file = .*|file = $restart|" production_stats.ini
{
    echo "job    : ${SLURM_JOB_ID:-none}   $(date '+%F %T %Z')   node ${SLURM_JOB_NODELIST:-?}"
    echo "code   : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "from   : $restart, $left steps left of $TOTAL"
} | tee -a provenance.txt

MPIRUN="mpirun -n 4 --bind-to core --map-by numa"
if [ ! -f production_p2.case.h5 ]; then
    $MPIRUN "$CODE_DIR/build_gpu/moby_prepare" production_stats.ini > prepare.log 2>&1 \
        || { echo "PREPARE FAILED"; tail -5 prepare.log; exit 1; }
fi
echo "start $restart $(date '+%F %T')" >> status.txt
$MPIRUN "$CODE_DIR/build_gpu/moby_solve" production_stats.ini > "production_stats_from_$(basename "$restart" .h5).log" 2>&1
rc=$?
echo "exit $rc $restart $(date '+%F %T')" >> status.txt
tail -3 status.txt
tail -2 production_stats_runtime.txt 2>/dev/null   # absent on a run shorter than runtime_interval
exit $rc
