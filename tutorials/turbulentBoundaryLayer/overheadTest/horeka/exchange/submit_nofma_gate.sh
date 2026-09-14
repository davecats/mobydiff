#!/bin/bash
#SBATCH --job-name=moby_nofma
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:40:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# THE bit-exactness gate of CLAUDE.md's "Verification": both sides built with
# -Mnofma -gpu=nofma, so the compiler cannot fuse a*b+c differently on either
# side of a refactor, and a surviving difference is a real one.
#
# Use this, not submit_apply.sh's production-flag gate, whenever the change MOVES
# an expression -- a hoisted metric, a precomputed table, a reordered sum. The
# production-flag comparison is the tighter test only when the arithmetic source
# is byte-identical between the two binaries.
#
# Timings must NOT be quoted from these binaries: -Mnofma is not how production
# is built.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/nofma_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC"/run_exchange.sh "$SRC"/run_mapgate.sh "$SRC"/collect_exchange.py "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

for d in "$REF_DIR" "$CODE_DIR"; do
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD)) -- gpu_nofma"
    ( cd "$d" && ./compile.sh gpu_nofma ) || exit 1
done
NEW="$CODE_DIR/build_gpu_nofma/moby_solve"; REF="$REF_DIR/build_gpu_nofma/moby_solve"
[ -x "$NEW" ] && [ -x "$REF" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_nofma}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  dirty $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
    echo "flags  : -Mnofma -gpu=nofma, BOTH sides"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

echo "############ Pass G -- production cases ############"
REF="$REF" H5MAXDIFF="$CODE_DIR/tools/h5maxdiff" CFG_DIR="$SRC/../configs" \
    PASSES="G" NSTEPS_GATE=20 bash "$STG/run_exchange.sh" "$NEW" "$RES/gate"

echo "############ the 7-case suite ############"
MOBY_ROOT="$CODE_DIR" NSTEPS=200 \
    bash "$STG/run_mapgate.sh" "$NEW" "$REF" "$RES/suite"

echo "=== verdict ==="
grep -h -E "^(OK|FAIL)" "$RES"/gate/*.txt "$RES"/suite/*.txt 2>/dev/null | sort | uniq -c
echo "=== nofma gate job finished ==="
