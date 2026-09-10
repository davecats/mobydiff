#!/bin/bash
#SBATCH --job-name=moby_mapgate
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=01:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The 7-case suite for select_target_device, against the pre-change binary that
# the campaign left in build_gpu/moby_solve.ref (verified: no MOBY_GPU_ORDER
# string, so it predates 059248e; it does carry be48d44's exchange diagnostics).
# See run_mapgate.sh for why nofma is deliberately not used.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
# Run a STAGED copy -- see submit_exchange.sh for why.
EXCH="$RUN_DIR/mapgate_staged"
rm -rf "$EXCH"; mkdir -p "$EXCH"
cp "$SRC"/run_mapgate.sh "$EXCH/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

cd "$CODE_DIR" || exit 1
./compile.sh gpu || exit 1
EXE="$CODE_DIR/build_gpu/moby_solve"
REF="$CODE_DIR/build_gpu/moby_solve.ref"
[ -x "$EXE" ] && [ -x "$REF" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_mapgate}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "commit : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "dirty  : $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "ref    : build_gpu/moby_solve.ref (pre-059248e, production flags)"
    echo "nsteps : ${NSTEPS:-200}"
} | tee "$RES/provenance.txt"

MOBY_ROOT="$CODE_DIR" bash "$EXCH/run_mapgate.sh" "$EXE" "$REF" "$RES"
echo "=== summary ==="
for f in "$RES"/*.txt; do
    [ "$(basename "$f")" = provenance.txt ] && continue
    echo "--- $(basename "$f" .txt)"; cat "$f"
done
echo "=== map gate job finished ==="
