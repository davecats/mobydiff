#!/bin/bash
#SBATCH --job-name=moby_prep
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=480000
#SBATCH --time=04:00:00
#SBATCH --partition=cpuonly
#SBATCH --account=hk-project-p0024333
#SBATCH --no-requeue
#
# moby_prepare of a large case on HoreKa CPU nodes (the c/2048 calibration
# case, 2026-10-29: about 1.5 M leaves, 450 GB resident -- the ISTM hosts
# hold 251 GB). Builds build_cpu in the worktree if missing. With an
# explicit [blocks] nb the case file does not depend on the rank count.
#     export CODE_DIR=$W/moby-junction RUN_DIR=$W/calib_run HDF5_ROOT=$HOME/hdf5 PREP_INI=calib_a1.5_y3.5_l7.ini
#     cd $RUN_DIR && sbatch --export=ALL submit_prepare.sh
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; PREP_INI="${PREP_INI:?}"
module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
if [ ! -x "$CODE_DIR/build_cpu/moby_prepare" ]; then
    echo "=== building the CPU binaries in $CODE_DIR ($(git -C "$CODE_DIR" rev-parse --short HEAD))"
    ( cd "$CODE_DIR" && ./compile.sh cpu ) > "$RUN_DIR/build_cpu.log" 2>&1 || { echo "CPU BUILD FAILED"; exit 1; }
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
