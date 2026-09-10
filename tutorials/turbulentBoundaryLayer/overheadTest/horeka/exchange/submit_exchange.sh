#!/bin/bash
#SBATCH --job-name=moby_exch
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
# Phases 1 and 2 of docs/next_session_2to1_performance.md in ONE allocation,
# plus the bit-exactness gate that the op-split diagnostic owes.
#
# TWO binaries are built here rather than reusing whatever build_gpu/ holds:
# CODE_DIR is the working tree (op split + A0 probe), REF_DIR is a git worktree
# detached at the pre-change commit. Same node, same module, same CMake
# configuration -- so the only difference between them is the diagnostic.
# (The earlier campaign was bitten twice by a stale shared build_gpu/.)
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"

# STAGE the driver, and run the COPY. bash reads a script incrementally by byte
# offset, so editing the file a running job is executing makes it resume at the
# wrong place -- in job 5139461 that silently re-launched the last A0 case
# without its probe. The campaign learned this once already (c1d5910); this is
# the fix rather than the discipline.
EXCH="$RUN_DIR/exchange_staged"
rm -rf "$EXCH"; mkdir -p "$EXCH"
cp "$SRC"/run_exchange.sh "$SRC"/run_mapgate.sh "$SRC"/collect_exchange.py "$EXCH/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

# Build on the compute node: the login node has no GPU, and nvfortran then emits
# every CC target instead of cc80 (and fails to link the device code at all).
for d in "$REF_DIR" "$CODE_DIR"; do
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
EXE="$CODE_DIR/build_gpu/moby_solve"
REF="$REF_DIR/build_gpu/moby_solve"
[ -x "$EXE" ] && [ -x "$REF" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_exchange}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "commit : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "dirty  : $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

REF="$REF" H5MAXDIFF="$CODE_DIR/tools/h5maxdiff" \
    CFG_DIR="$SRC/../configs" \
    bash "$EXCH/run_exchange.sh" "$EXE" "$RES"

python3 "$EXCH/collect_exchange.py" "$RES" > "$RES/exchange.md" 2>&1 \
    && { echo; cat "$RES/exchange.md"; } || echo "=== collector failed; raw logs intact ==="
echo "=== exchange job finished ==="
