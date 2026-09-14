#!/bin/bash
#SBATCH --job-name=moby_apply8
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:40:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The 8-rank (and a same-allocation 4-rank) point, on dev_accelerated: `accelerated`
# is three days out for the 4-node job, and the dev partition schedules same-day.
# Pass A only -- untraced brackets and step times, before and after, in ONE
# allocation (the only way the two sides are comparable).
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/apply8_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC"/run_timeline.sh "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1
for d in "$REF_DIR" "$CODE_DIR"; do ( cd "$d" && ./compile.sh gpu ) || exit 1; done
NEW="$CODE_DIR/build_gpu/moby_solve"; REF="$REF_DIR/build_gpu/moby_solve"

RES="${RESDIR:-$RUN_DIR/results_apply8}"; mkdir -p "$RES"
{ echo "job    : ${SLURM_JOB_ID:-none}"; echo "date   : $(date '+%F %T %Z')"
  echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
  echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)"
  echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
} | tee "$RES/provenance.txt"

SPECS="refined_yp82_rect_jacobi:8:2 rect_jacobi:8:2 refined_yp82_rect_jacobi:4:1 rect_jacobi:4:1"
for side in before after; do
    exe="$REF"; [ "$side" = after ] && exe="$NEW"
    echo "############ $side ############"
    CFG_DIR="$SRC/../configs" TL_PASSES="A" SPECS="$SPECS" NSTEPS="${NSTEPS:-100}" \
        bash "$STG/run_timeline.sh" "$exe" "$RES/$side"
done

echo "=== summary: s/step and the projection buckets ==="
for side in before after; do
  for d in "$RES/$side"/host_*; do
    [ -d "$d" ] || continue
    printf "%-7s %-34s " "$side" "$(basename "$d" | sed 's/^host_//')"
    grep -E "^timing: nsteps" "$d/run.log" | awk '{printf "s/step %s  ", $NF}'
    grep -E "proj_timing: (sweep|apply) " "$d/run.log" \
      | awk '{printf "%s=%s ", $2, $(NF-2)}'
    echo
  done
done
echo "=== apply8 job finished ==="
