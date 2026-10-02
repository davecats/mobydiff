#!/bin/bash
#SBATCH --job-name=moby_outlet_tbl
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
# The turbulent boundary-layer tutorial across the outlet change
# (docs/next_session_outlet.md, follow-up of 2026-10-02): the production case
# (4096 x 224 x 192, red-black 6) restarted from the same developed field
# with the two binaries, 4 ranks each, ONE NODE PER SIDE in the same
# allocation, statistics on.
#
#     export CODE_DIR=<worktree at 949148f> REF_DIR=<worktree at 154e48f> \
#            RUN_DIR=<run dir> HDF5_ROOT=$HOME/hdf5 [NSTEPS=12500] [SIDES="ref new"]
#     cd $RUN_DIR && sbatch --export=ALL submit_outlet_tbl.sh
#
#   ref  REF_DIR   the outlet face reset by apply_bc
#   new  CODE_DIR  the outlet face predicted
#
# RESUMABLE, in chunks of NSTEPS: each side restarts from its own latest
# snapshot (tbl_<side>_<step>.h5; the first chunk from
# reftrip_restart_global.h5, a global-3D file readable on any rank count),
# its statistics file continues, and a chunk writes ONE snapshot, at its end.
# A chunk that hits the wall clock leaves no snapshot and is simply run again:
# keep NSTEPS x seconds-per-step under ~45 minutes. SIDES=smoke runs 100 steps
# of the new binary on the first node, for the step time.
# RUN_DIR holds reftrip_restart_global.h5 and tbl.ini (the production ini with
# the statistics keys; prefix, restart and step count are set here).
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; REF_DIR="${REF_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
NSTEPS="${NSTEPS:-12500}"; SIDES="${SIDES:-ref new}"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

for d in "$REF_DIR" "$CODE_DIR"; do
    cache="$d/build_gpu/CMakeCache.txt"
    if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$d" "$cache"; then
        rm -rf "$d/build_gpu"
    fi
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) > "$RUN_DIR/build_$(basename "$d").log" 2>&1 || { echo "BUILD FAILED $d"; exit 1; }
done

cd "$RUN_DIR"
mapfile -t NODES < <(scontrol show hostnames "$SLURM_JOB_NODELIST")
{
    echo "job    : ${SLURM_JOB_ID:-none}   $(date '+%F %T %Z')"
    echo "nodes  : ${NODES[*]}"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "chunk  : $NSTEPS steps, sides: $SIDES"
} | tee -a provenance.txt

n=0
for side in $SIDES; do
    dir="$CODE_DIR"; steps="$NSTEPS"
    [ "$side" = ref ] && dir="$REF_DIR"
    [ "$side" = smoke ] && steps=100
    node="${NODES[$((n % ${#NODES[@]}))]}"; n=$((n + 1))
    # The latest snapshot of this side: tbl_<side>_<step>.h5, largest step.
    latest=$(ls -1 tbl_${side}_[0-9]*.h5 2>/dev/null | sort -t_ -k3 -n | tail -1)
    restart="${latest:-reftrip_restart_global.h5}"
    ini="tbl_${side}.ini"
    sed -e "s/^nsteps = .*/nsteps = $steps/" \
        -e "s/^field_interval = .*/field_interval = 0/" \
        -e "s/^field_prefix = .*/field_prefix = tbl_${side}/" \
        -e "s/^runtime_file = .*/runtime_file = tbl_${side}_runtime.txt/" \
        -e "s/^stats_file = .*/stats_file = tbl_${side}_stats.h5/" \
        -e "s/^stats_write_interval = .*/stats_write_interval = $steps/" \
        -e "/^\[restart\]/,/^$/ s|^file = .*|file = $restart|" tbl.ini > "$ini"
    MPIRUN="mpirun -n 4 --host $node:4 --bind-to core --map-by numa"
    if [ ! -f "tbl_${side}.case.h5" ]; then
        $MPIRUN "$dir/build_gpu/moby_prepare" "$ini" > "prep_$side.log" 2>&1 \
            || { echo "PREPARE FAILED $side"; tail -5 "prep_$side.log"; exit 1; }
    fi
    log="tbl_${side}_from_$(basename "$restart" .h5).log"
    echo "=== $side on $node: $steps steps from $restart"
    ( echo "start $side $restart $(date '+%F %T')" >> status.txt
      $MPIRUN "$dir/build_gpu/moby_solve" "$ini" > "$log" 2>&1
      echo "exit $? $side $restart $(date '+%F %T')" >> status.txt ) &
done
wait
tail -4 status.txt
for side in $SIDES; do tail -2 tbl_${side}_runtime.txt 2>/dev/null; done
