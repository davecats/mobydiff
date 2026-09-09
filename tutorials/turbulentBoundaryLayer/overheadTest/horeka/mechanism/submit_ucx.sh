#!/bin/bash
#SBATCH --job-name=moby_ucx
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:30:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# Environment-only sweep over UCX's handling of GPU memory, on the exact
# configuration and placement of the 738 us anomaly. No solver code, no rebuild.
# See run_ucx.sh for what the traces showed and why these knobs.
set -uo pipefail

CODE_DIR="${CODE_DIR:-$HOME/moby-2to1-code}"
RUN_DIR="${RUN_DIR:-$SLURM_SUBMIT_DIR}"
EXE="$CODE_DIR/build_gpu/moby_solve"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
# deliberately NOT setting UCX_MEMTYPE_CACHE here: it is one of the variables
# under test, and run_ucx.sh sets it per variant.

[ -x "$EXE" ] || { echo "ERROR: no $EXE -- build it in a timing job first" >&2; exit 1; }

cd "$RUN_DIR" || exit 1
RES="${RESDIR:-$RUN_DIR/results_ucx}"
mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "commit : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "dirty  : $(git -C "$CODE_DIR" status --porcelain --untracked-files=no | wc -l)"
} | tee "$RES/provenance.txt"

echo "=== HCA / GPU topology on this node ==="
{ ibstat -l; echo "--- nvidia-smi topo ---"; nvidia-smi topo -m; } 2>&1 | tee "$RES/topology.txt"

bash "$RUN_DIR/mechanism/run_ucx.sh" "$EXE" "$RES"
echo "=== ucx job finished ==="
