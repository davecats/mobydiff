#!/bin/bash
#SBATCH --job-name=moby_nsys
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=01:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# Timeline probe. Runs alongside (not inside) the timing job: nothing here is a
# timing comparison, so it does not need to share that allocation -- the
# measurement is the STRUCTURE of the timeline, which is a property of the
# configuration, not of the node set.
#
# Why this exists at all: three mechanisms have been proposed and refuted from
# aggregate mpi_wait, which can only report how long a call took. It cannot say
# whether the Waitall is blocked on the wire or on an unfinished CUDA stream,
# whether all 39 rounds are slow or one is catastrophic, or which rank posts
# late. A trace can, and Score-P/Scalasca -- the tool built for exactly this --
# is unavailable here (it needs compiler/gnu + mpi/openmpi, and the solver needs
# nvfortran for OpenMP target offload).
set -uo pipefail

CODE_DIR="${CODE_DIR:-$HOME/moby-2to1-code}"
RUN_DIR="${RUN_DIR:-$SLURM_SUBMIT_DIR}"
EXE="$CODE_DIR/build_gpu/moby_solve"

module purge
module load toolkit/nvidia-hpc-sdk/25.3

export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n
export OMP_NUM_THREADS=1

cd "$CODE_DIR" || exit 1
echo "commit under test: $(git rev-parse HEAD)"
./compile.sh gpu || exit 1
[ -x "$EXE" ] || { echo "ERROR: no $EXE" >&2; exit 1; }

cd "$RUN_DIR" || exit 1
RES="$RUN_DIR/results_nsys"
mkdir -p "$RES"
{
    echo "job          : ${SLURM_JOB_ID:-none}"
    echo "date         : $(date '+%F %T %Z')"
    echo "nodes        : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "commit       : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "dirty        : $(git -C "$CODE_DIR" status --porcelain --untracked-files=no | wc -l)"
    echo "gpu          : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps       : ${NSTEPS:-40}"
} | tee "$RES/provenance.txt"

bash "$RUN_DIR/mechanism/run_nsys.sh" "$EXE" "$RES"

python3 "$RUN_DIR/mechanism/collect_nsys.py" "$RES" > "$RES/nsys.md" 2>&1 \
    && { echo; cat "$RES/nsys.md"; } \
    || echo "=== collect_nsys failed; CSVs and .nsys-rep intact ==="

echo "=== nsys job finished ==="
