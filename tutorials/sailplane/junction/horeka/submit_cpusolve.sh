#!/bin/bash
#SBATCH --job-name=moby_cpusolve
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=480000
#SBATCH --time=01:00:00
#SBATCH --partition=cpuonly
#SBATCH --account=hk-project-p0024333
#SBATCH --no-requeue
#
# ONE short run of an ini on the CPU build (the reference path), on a
# cpuonly node: the first-step blow-up of the c/2048 calibration box
# reproduced one level lower (2026-11-01) -- host and device memory
# differ in what an unwritten halo or an unset coefficient holds, so a
# CPU run that does NOT blow up where the GPU run does points at
# uninitialised data, not at the numerics. The binary comes from the
# prepare job's build (build_cpu, rebuilt here when the commit moved).
#     export CODE_DIR=$W/moby-junction RUN_DIR=$W/calib_run RUN_INI=.box6cpu.ini
#     cd $RUN_DIR && sbatch --export=ALL submit_cpusolve.sh
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; RUN_INI="${RUN_INI:?}"
module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
head=$(git -C "$CODE_DIR" rev-parse HEAD)
if [ ! -x "$CODE_DIR/build_cpu/moby_solve" ] || [ "$(cat "$CODE_DIR/build_cpu/.commit" 2>/dev/null)" != "$head" ]; then
    echo "=== building the CPU binaries in $CODE_DIR (${head:0:7})"
    ( cd "$CODE_DIR" && ./compile.sh cpu ) > "$RUN_DIR/build_cpu.log" 2>&1 || { echo "CPU BUILD FAILED"; exit 1; }
    echo "$head" > "$CODE_DIR/build_cpu/.commit"
fi
cd "$RUN_DIR"
log="${RUN_INI#.}"; log="${log%.ini}.log"
echo "job ${SLURM_JOB_ID:-none} $(date '+%F %T') $RUN_INI -> $log (CPU, $SLURM_NTASKS ranks), code $(git -C "$CODE_DIR" rev-parse --short HEAD)" | tee -a provenance.txt
mpirun -n "${SLURM_NTASKS:-64}" --bind-to core "$CODE_DIR/build_cpu/moby_solve" "$RUN_INI" > "$log" 2>&1
echo "exit $? $RUN_INI $(date '+%F %T')" | tee -a provenance.txt
