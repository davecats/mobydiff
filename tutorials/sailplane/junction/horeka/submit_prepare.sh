#!/bin/bash
#SBATCH --job-name=moby_prep
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=480000
#SBATCH --time=06:00:00
#SBATCH --partition=cpuonly
#SBATCH --account=hk-project-p0024333
#SBATCH --no-requeue
#
# moby_prepare of a large case on HoreKa CPU nodes (the c/2048 calibration
# case, 2026-10-29). Builds build_cpu in the worktree if missing or from
# another commit. With an
# explicit [blocks] nb the case file does not depend on the rank count.
#
# RANKS PER NODE: the per-level touch/buried masks of refine_body are
# dense rasters REPLICATED on every rank (classify_refinement_masks; the
# ranks split the scan and merge by allreduce), every level column sized
# to the LARGEST window. Windowed to the STL bounding box the level-7
# window of the c/2048 case was 1e9 blocks: 2 masks x 8 columns x 4 B =
# 64 GB per rank, OOM-killed at 64 AND at 32 ranks per 480 GB node (jobs
# 5217147 and 5217871, 2026-10-29/30). Since 2026-10-30 the levels above
# refine_body_levels are windowed to the refine_body_box hull: the
# prepare prints "refine_body masks: level l window ..." and the GB per
# rank before allocating -- READ IT in the log; the c/2048 case is 72 M
# blocks at level 7, 4.3 GB per rank. 32 ranks per node leave the margin
# for the coefficients and the wall distance; add nodes for speed. The
# time goes to the per-leaf coefficients and wall distance (0.47 rank-s
# per leaf measured at c/512: 3-4 M leaves = 400-500 rank-h, 2-3 h on
# 6 x 32 ranks).
#     export CODE_DIR=$W/moby-junction RUN_DIR=$W/calib_run HDF5_ROOT=$HOME/hdf5 PREP_INI=calib_a1.5_y3.5_l7.ini
#     cd $RUN_DIR && sbatch --export=ALL submit_prepare.sh        # --nodes=N --ntasks-per-node=32 to resize
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; PREP_INI="${PREP_INI:?}"
module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
# Rebuild when the worktree moved to another commit (the binary records
# the commit it was built from).
head=$(git -C "$CODE_DIR" rev-parse HEAD)
if [ ! -x "$CODE_DIR/build_cpu/moby_prepare" ] || [ "$(cat "$CODE_DIR/build_cpu/.commit" 2>/dev/null)" != "$head" ]; then
    echo "=== building the CPU binaries in $CODE_DIR (${head:0:7})"
    ( cd "$CODE_DIR" && ./compile.sh cpu ) > "$RUN_DIR/build_cpu.log" 2>&1 || { echo "CPU BUILD FAILED"; exit 1; }
    echo "$head" > "$CODE_DIR/build_cpu/.commit"
fi
cd "$RUN_DIR"
# The STL is git-ignored (14 MB per pitched copy): ship it into RUN_DIR by
# hand. Stop here rather than let moby_prepare stop on every rank after the
# queue wait (the first run of the l7 case lost its slot exactly that way).
stl=$(sed -n 's/^stl_file *= *"\{0,1\}\([^"]*\)"\{0,1\}.*/\1/p' "$PREP_INI" | head -1)
if [ -n "$stl" ] && [ ! -f "$stl" ]; then echo "no STL $stl in $RUN_DIR (ship it: ssh horeka \"cat > $RUN_DIR/$stl\" < calib/$stl)"; exit 1; fi
log="${PREP_INI%.ini}.prep.log"
echo "job ${SLURM_JOB_ID:-none} $(date '+%F %T') prepare $PREP_INI on $SLURM_NNODES nodes x $SLURM_NTASKS_PER_NODE ranks, code $(git -C "$CODE_DIR" rev-parse --short HEAD)" | tee -a provenance.txt
/usr/bin/time -f "prepare wall %e s, max RSS %M kB (rank 0 of the time wrapper)" \
    mpirun -n $SLURM_NTASKS "$CODE_DIR/build_cpu/moby_prepare" "$PREP_INI" > "$log" 2>&1
rc=$?
grep -h "leaves\|prepare wall" "$log" | tail -3
echo "exit $rc prepare $PREP_INI $(date '+%F %T')" | tee -a provenance.txt
