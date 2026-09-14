#!/bin/bash
#SBATCH --job-name=moby_divgate
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:30:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# Pass G alone, re-run: job 5145805's Pass G produced no comparison because
# h5maxdiff is a BUILD PRODUCT and the pinned worktree did not have it -- and
# run_exchange.sh deletes the snapshots whether or not the comparison ran, so
# the gate was lost rather than merely unreported. submit_divhalo.sh now builds
# it; this job re-does the pass on the same two production cases.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/divgate_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC"/run_exchange.sh "$SRC"/collect_exchange.py "$STG/"

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
H5MAXDIFF="$CODE_DIR/tools/h5maxdiff"
[ -x "$H5MAXDIFF" ] || gcc -O2 -o "$H5MAXDIFF" "$CODE_DIR/tools/h5maxdiff.c" \
    -I"$HDF5_ROOT/include" -L"$HDF5_ROOT/lib" -lhdf5 -Wl,-rpath,"$HDF5_ROOT/lib" || exit 1
[ -x "$NEW" ] && [ -x "$REF" ] || { echo "ERROR: missing binary" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_divhalo/gate2}"; mkdir -p "$RES"
echo "job $SLURM_JOB_ID  new $(git -C "$CODE_DIR" rev-parse --short HEAD)  ref $(git -C "$REF_DIR" rev-parse --short HEAD)" \
    | tee "$RES/provenance.txt"

REF="$REF" H5MAXDIFF="$H5MAXDIFF" CFG_DIR="$SRC/../configs" \
    PASSES="G" NSTEPS_GATE=20 bash "$STG/run_exchange.sh" "$NEW" "$RES"
echo "=== divhalo Pass G re-run finished ==="
