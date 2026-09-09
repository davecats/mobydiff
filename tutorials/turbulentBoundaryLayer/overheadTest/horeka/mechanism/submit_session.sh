#!/bin/bash
#SBATCH --job-name=moby_balance
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
# Phase 1 + Phase 2 of HANDOUT_cluster_session.md in ONE allocation.
#
# Deliberate deviation from the handout's "resume and add three runs": the whole
# 26-spec matrix is re-run into a FRESH directory with the be48d44 binary. The
# 23 existing runs in results_mechanism/ came from two different jobs, two node
# sets and the pre-diagnostics binary; re-running them costs ~22 min of an
# otherwise idle 2 h allocation and removes the cross-job caveat from EVERY
# comparison at once, not just from rect_jacobi:8:2. It also delivers the free
# `exchange balance:` line for all 26 rather than for the four of Phase 2(a).
#
#   pass A -> results_balance : normal timing + the balance line (quotable)
#   pass B -> results_barrier : exchange_barrier=true, skew/transfer split
#                               (SERIALISING -- attribution only, never quote
#                                its step times)
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

# ALWAYS rebuild: the point of this session is the be48d44 diagnostics, and the
# build_gpu/ left by the earlier jobs is the a11e355 binary. Building on the
# compute node keeps the -mp=gpu target selection identical to the earlier
# campaign's (the login node has no GPU).
cd "$CODE_DIR" || exit 1
echo "commit under test: $(git rev-parse HEAD)"
./compile.sh gpu || exit 1
[ -x "$EXE" ] || { echo "ERROR: no $EXE" >&2; exit 1; }

cd "$RUN_DIR" || exit 1
BAL="$RUN_DIR/results_balance"
BAR="$RUN_DIR/results_barrier"
mkdir -p "$BAL" "$BAR"
{
    echo "job          : ${SLURM_JOB_ID:-none}"
    echo "date         : $(date '+%F %T %Z')"
    echo "nodes        : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "commit       : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "dirty        : $(git -C "$CODE_DIR" status --porcelain --untracked-files=no | wc -l)"
    echo "gpu          : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "nsteps       : ${NSTEPS:-200}"
} | tee "$BAL/provenance.txt"
cp "$BAL/provenance.txt" "$BAR/provenance.txt"

echo "############ PASS A -- timing + exchange balance ############"
bash "$RUN_DIR/mechanism/run_placement.sh" "$EXE" "$BAL"

echo "############ PASS B -- barrier: skew vs transfer ############"
BARRIER=1 SPECS="rect_jacobi:8:2 rect_jacobi:8:4 base_jacobi:8:2 nb16_jacobi:8:2 blk8_jacobi:8:2 rect_jacobi:4:1" \
    bash "$RUN_DIR/mechanism/run_placement.sh" "$EXE" "$BAR"

python3 "$RUN_DIR/mechanism/collect_mechanism.py" "$BAL" > "$BAL/mechanism.md" 2>&1 \
    || echo "=== collect_mechanism failed; raw logs intact ==="
python3 "$RUN_DIR/mechanism/collect_balance.py" "$BAL" "$BAR" > "$RUN_DIR/balance.md" 2>&1 \
    && { echo; cat "$RUN_DIR/balance.md"; } \
    || echo "=== collect_balance failed; raw logs intact ==="

echo "=== balance job finished ==="
