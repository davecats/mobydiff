#!/bin/bash
#SBATCH --job-name=moby_mechanism
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=02:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# Settle the mechanism behind the blocked exchange's 10x mpi_wait jump.
#
# 4 nodes are allocated so that a FIXED rank count can be spread over 1, 2 or 4
# of them. Most runs use a fraction of the allocation on purpose: holding the
# rank count fixed while varying only the node span is the whole point, and
# doing it inside one allocation keeps the machine state fixed.
#
# RESUMABLE: run_placement.sh skips any run whose run.log exists.
set -uo pipefail

CODE_DIR="${CODE_DIR:-$HOME/moby-2to1-code}"
RUN_DIR="${RUN_DIR:-$SLURM_SUBMIT_DIR}"
EXE="$CODE_DIR/build_gpu/moby_solve"

# --- environment -----------------------------------------------------------
# module purge FIRST: HoreKa's default Intel toolchain exports
# FFLAGS="-O2 -xCORE-AVX2", which CMake picks up as CMAKE_Fortran_FLAGS and
# nvfortran rejects. Purging clears it coherently, along with the PATH and
# library entries that would interfere just as easily.
module purge
module load toolkit/nvidia-hpc-sdk/25.3

export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
if [ ! -e "$HDF5_ROOT/lib/libhdf5.so" ] && [ ! -e "$HDF5_ROOT/lib/libhdf5.a" ]; then
    bash "$RUN_DIR/build_hdf5.sh" "$HDF5_ROOT" || exit 1
fi
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n
export OMP_NUM_THREADS=1

# --- build (reused from the main campaign if already present) --------------
# A build dir with no executable is a failed configure whose CMakeCache has
# already captured the environment's flags; reusing it re-applies them.
if [ ! -x "$EXE" ] && [ -d "$CODE_DIR/build_gpu" ]; then
    echo "=== discarding incomplete build dir $CODE_DIR/build_gpu ==="
    rm -rf "$CODE_DIR/build_gpu"
fi
if [ ! -x "$EXE" ]; then
    cd "$CODE_DIR" || exit 1
    echo "commit under test: $(git rev-parse HEAD)"
    ./compile.sh gpu || exit 1
fi
[ -x "$EXE" ] || { echo "ERROR: no $EXE" >&2; exit 1; }

cd "$RUN_DIR" || exit 1
RESULTS="$RUN_DIR/results_mechanism"
mkdir -p "$RESULTS"
{
    echo "job          : ${SLURM_JOB_ID:-none}"
    echo "date         : $(date '+%F %T %Z')"
    echo "nodes        : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "commit       : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "dirty        : $(git -C "$CODE_DIR" status --porcelain --untracked-files=no | wc -l)"
    echo "gpu          : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps       : ${NSTEPS:-200}"
} | tee "$RESULTS/provenance.txt"

bash "$RUN_DIR/mechanism/run_placement.sh" "$EXE" "$RESULTS"

python3 "$RUN_DIR/mechanism/collect_mechanism.py" "$RESULTS" > "$RESULTS/mechanism.md" 2>&1 \
    && { echo; cat "$RESULTS/mechanism.md"; } \
    || echo "=== collector failed; raw logs intact ==="

echo "=== mechanism job finished ==="
