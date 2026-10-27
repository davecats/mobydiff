#!/bin/bash
#SBATCH --job-name=moby_diag
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:30:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-p0024333
#SBATCH --no-requeue
#
# ONE short run of an ini on 4 A100 ranks with the already-built worktree:
# the diagnostics of the calibration box (2026-10-27, the cold-start
# divergence). No chunking, no restart logic, the ini says everything.
#     export CODE_DIR=$W/moby-junction RUN_DIR=$W/calib_run RUN_INI=.diag_base.ini
#     cd $RUN_DIR && sbatch --export=ALL submit_diag.sh
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; RUN_INI="${RUN_INI:?}"
module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1
cd "$RUN_DIR"
log="${RUN_INI#.}"; log="${log%.ini}.log"
echo "job ${SLURM_JOB_ID:-none} $(date '+%F %T') $RUN_INI -> $log, code $(git -C "$CODE_DIR" rev-parse --short HEAD)" | tee -a provenance.txt
mpirun -n 4 --bind-to core --map-by numa "$CODE_DIR/build_gpu/moby_solve" "$RUN_INI" > "$log" 2>&1
echo "exit $? $RUN_INI $(date '+%F %T')" | tee -a provenance.txt
