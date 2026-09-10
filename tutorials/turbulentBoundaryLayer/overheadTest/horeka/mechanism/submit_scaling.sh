#!/bin/bash
#SBATCH --job-name=moby_scaling
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
# Re-measure the 2026-09-07 campaign at the corrected rank-to-GPU mapping.
#
# Every scaling number, block tax and strong-scaling efficiency in
# results_horeka_2026-09-07.md was taken with `device = local_rank mod ndev`,
# which pairs a node's LAST local rank with the next node's FIRST across every
# node boundary -- the mismatched-affinity case now known to cost 20-27 % of the
# step. Those tables therefore measure the mapping as much as they measure
# blocking or refinement.
#
# The SAME 23-run matrix (run_matrix.sh, unchanged, --map-by numa --bind-to core,
# 400 steps) is run twice in ONE allocation: once with the pre-change binary and
# once with the new one. Same nodes, same machine state, so the difference is the
# mapping and nothing else -- which is what the campaign's own rule demands and
# what comparing against the Sep-7 numbers (different nodes, two days apart, ~6 %
# node variation) could not give.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
EXE="$CODE_DIR/build_gpu/moby_solve"; REF="$CODE_DIR/build_gpu/moby_solve.ref"
HOREKA="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka"

module purge; module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

[ -x "$EXE" ] || { echo "ERROR: no $EXE" >&2; exit 1; }
[ -x "$REF" ] || { echo "ERROR: no reference binary $REF" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_scaling}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "commit : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "dirty  : $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "nsteps : ${NSTEPS:-400}"
    echo "layout : --map-by numa --bind-to core (as the 2026-09-07 campaign)"
} | tee "$RES/provenance.txt"

# NEW first: if the wall clock bites, the corrected numbers are the deliverable
# and the reference column is the nice-to-have.
echo "############ NEW binary (topology-aware mapping) ############"
NSTEPS="${NSTEPS:-400}" bash "$HOREKA/run_matrix.sh" "$EXE" "$RES/new"
echo "############ REFERENCE binary (local_rank mod ndev) ############"
NSTEPS="${NSTEPS:-400}" bash "$HOREKA/run_matrix.sh" "$REF" "$RES/ref"

python3 "$RUN_DIR/mechanism/collect_scaling.py" "$RES/ref" "$RES/new" > "$RES/scaling.md" 2>&1 \
    && { echo; cat "$RES/scaling.md"; } || echo "=== collector failed; raw logs intact ==="
echo "=== scaling job finished ==="
