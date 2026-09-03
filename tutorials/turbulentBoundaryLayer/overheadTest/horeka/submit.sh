#!/bin/bash
#SBATCH --job-name=moby_2to1_prof
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=06:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# Timing + profiling matrix for the 2:1 block-refinement machinery.
#
# 4 HoreKa Green nodes = 16x A100-40. The matrix runs each config at 1/2/4
# (intra-node, NVLink) and 8/16 (inter-node) ranks; a small run inside a large
# allocation is deliberate -- one allocation gives a rank sweep with the machine
# state held fixed, which is the only way the sweep is internally comparable.
#
# RESUMABLE: run_matrix.sh skips any run whose run.log exists, so if the job
# hits the wall clock just `sbatch submit.sh` again.
set -uo pipefail

CODE_DIR="${CODE_DIR:-$HOME/moby-2to1-code}"
RUN_DIR="${RUN_DIR:-$SLURM_SUBMIT_DIR}"
EXE="$CODE_DIR/build_gpu/moby_solve"

# --- environment -----------------------------------------------------------
module purge
module load toolkit/nvidia-hpc-sdk/25.3
# cmake (3.26.5) is on PATH by default on HoreKa -- no module needed.

# Parallel (MPI-IO) HDF5, built from source with the same nvhpc toolchain:
# HoreKa ships only a SERIAL HDF5 module and mobydiff's field I/O is collective.
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
if [ ! -e "$HDF5_ROOT/lib/libhdf5.so" ] && [ ! -e "$HDF5_ROOT/lib/libhdf5.a" ]; then
    echo "=== parallel HDF5 not found at $HDF5_ROOT -- building it ==="
    bash "$RUN_DIR/build_hdf5.sh" "$HDF5_ROOT" || exit 1
fi
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"

# CUDA-aware MPI. A100 nodes have NVLink so cuda_ipc/P2P work; UCX_MEMTYPE_CACHE=n
# is the cheap safeguard against the UCX memtype-cache segfault. Compute is on the
# GPU, so one OMP thread per rank.
export UCX_MEMTYPE_CACHE=n
export OMP_NUM_THREADS=1

# --- build once, on the compute node (A100 cc80 auto-detected) -------------
if [ ! -x "$EXE" ]; then
    echo "=== building the GPU solver in $CODE_DIR ==="
    cd "$CODE_DIR" || exit 1
    echo "commit under test: $(git rev-parse HEAD)"
    git --no-pager log --oneline -1
    ./compile.sh gpu || exit 1
fi
[ -x "$EXE" ] || { echo "ERROR: build did not produce $EXE" >&2; exit 1; }

# --- provenance, recorded next to the results ------------------------------
cd "$RUN_DIR" || exit 1
RESULTS="$RUN_DIR/results"
mkdir -p "$RESULTS"
{
    echo "job          : ${SLURM_JOB_ID:-none} on $(hostname)"
    echo "date         : $(date '+%F %T %Z')"
    echo "nodes        : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "commit       : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "commit_subj  : $(git -C "$CODE_DIR" log --oneline -1)"
    echo "dirty        : $(git -C "$CODE_DIR" status --porcelain --untracked-files=no | wc -l) tracked changes"
    echo "compiler     : $(nvfortran --version 2>/dev/null | sed -n 2p)"
    echo "hdf5_root    : $HDF5_ROOT"
    echo "gpu          : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RESULTS/provenance.txt"

# --- the matrix ------------------------------------------------------------
bash "$RUN_DIR/run_matrix.sh" "$EXE" "$RESULTS"

# --- collect ---------------------------------------------------------------
# Best-effort: needs a python3, which HoreKa has on PATH. The collector uses
# only the standard library precisely so it cannot fail for a missing module.
python3 "$RUN_DIR/collect_profile.py" "$RESULTS" > "$RESULTS/summary.md" 2>&1 \
    && echo "=== summary written to $RESULTS/summary.md ===" \
    || echo "=== collector failed; the raw run.log files are intact ==="

echo "=== job finished (or hit the wall clock) ==="
