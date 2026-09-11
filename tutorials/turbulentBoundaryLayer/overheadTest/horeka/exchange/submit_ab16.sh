#!/bin/bash
#SBATCH --job-name=moby_ab16
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:30:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The one number the A/B did not measure: what `map(to: c)` is worth at 8 and 16
# ranks, where the launch bill lands on a step three to five times smaller. Pass A
# only -- untraced brackets and step times, before and after, in one allocation.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
EXCH="$RUN_DIR/ab16_staged"; rm -rf "$EXCH"; mkdir -p "$EXCH"
cp "$SRC"/run_timeline.sh "$EXCH/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

for d in "$REF_DIR" "$CODE_DIR"; do
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
NEW="$CODE_DIR/build_gpu/moby_solve"; REF="$REF_DIR/build_gpu/moby_solve"

RES="${RESDIR:-$RUN_DIR/results_ab16}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)"
} | tee "$RES/provenance.txt"

SPECS="refined_yp82_rect_jacobi:16:4 rect_jacobi:16:4 refined_yp82_rect_jacobi:8:2 rect_jacobi:8:2"
for side in before after; do
    exe="$REF"; [ "$side" = after ] && exe="$NEW"
    echo "############ $side ############"
    CFG_DIR="$SRC/../configs" TL_PASSES="A" SPECS="$SPECS" NSTEPS="${NSTEPS:-100}" \
        bash "$EXCH/run_timeline.sh" "$exe" "$RES/$side"
done

echo "=== summary: s/step and per-round buckets ==="
for side in before after; do
  for d in "$RES/$side"/host_*; do
    printf "%-7s %-34s " "$side" "$(basename "$d" | sed 's/^host_//')"
    grep -E "^timing: nsteps" "$d/run.log" | awk '{printf "s/step %s  ", $NF}'
    grep -E "exch_timing: (pack|unpack|local_copy|copy_cross) " "$d/run.log" \
      | awk '{printf "%s=%.1f ", $2, $(NF-2)*1e6/39}'
    echo
  done
done
echo "=== ab16 job finished ==="
