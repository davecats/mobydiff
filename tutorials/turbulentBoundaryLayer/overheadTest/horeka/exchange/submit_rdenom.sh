#!/bin/bash
#SBATCH --job-name=moby_rdenom
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:55:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# rdenom formed ONCE on a body-free rank instead of every substage.
#
# TWO columns; the third slot is unused and set to the same as `prefix`.
#   base   b9414bd  -- the index list, i.e. the current published state
#   list   HEAD     -- + rdenom static
#
# Rank counts 1, 4, 8, and RED-BLACK included as a CONTROL: redblack_projection
# never calls compute_rdenom, so its `setup` must not move.
#
# Readings pre-registered in PREREGISTERED_rdenom_static.md.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"; BASE_DIR="${BASE_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/rdenom_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC"/run_timeline.sh "$SRC"/run_exchange.sh "$SRC"/run_mapgate.sh \
   "$SRC"/collect_exchange.py "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

# A CMake cache records the ABSOLUTE source path; a renamed worktree carries a
# stale one and cmake refuses (job 5145816).
for d in "$BASE_DIR" "$REF_DIR" "$CODE_DIR"; do
    cache="$d/build_gpu/CMakeCache.txt"
    if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$d" "$cache"; then
        echo "=== $d: stale build_gpu cache -- discarding"; rm -rf "$d/build_gpu"
    fi
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
LIST="$CODE_DIR/build_gpu/moby_solve"
PREFIX="$REF_DIR/build_gpu/moby_solve"
BASE="$BASE_DIR/build_gpu/moby_solve"
for x in "$LIST" "$PREFIX" "$BASE"; do [ -x "$x" ] || { echo "ERROR: missing $x" >&2; exit 1; }; done
# h5maxdiff is a build product; a fresh worktree lacks it and Pass G deletes its
# snapshots regardless, so a missing comparator LOSES the gate (job 5145805).
# mpicc, NOT gcc: this HDF5 is the PARALLEL build, so H5public.h includes
# mpi.h and `module purge` has removed it from the default include path.
H5MAXDIFF="$CODE_DIR/tools/h5maxdiff"
[ -x "$H5MAXDIFF" ] || mpicc -O2 -o "$H5MAXDIFF" "$CODE_DIR/tools/h5maxdiff.c" \
    -I"$HDF5_ROOT/include" -L"$HDF5_ROOT/lib" -lhdf5 -Wl,-rpath,"$HDF5_ROOT/lib" || exit 1

RES="${RESDIR:-$RUN_DIR/results_rdenom}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  (rdenom static)"
    echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)  (index list, b9414bd)"
    echo "(two columns only)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

SPECS="rect_jacobi:8:2 refined_yp82_rect_jacobi:8:2 refined_yp82_rect_redblack:8:2 \
rect_jacobi:4:1 refined_yp82_rect_jacobi:4:1 refined_yp82_rect_redblack:4:1 \
rect_jacobi:1:1 refined_yp82_rect_jacobi:1:1 refined_yp82_rect_redblack:1:1"
for side in list prefix; do
    case $side in list) exe="$LIST";; prefix) exe="$PREFIX";; base) exe="$BASE";; esac
    echo "############ timeline: $side ############"
    CFG_DIR="$SRC/../configs" TL_PASSES="A" SPECS="$SPECS" NSTEPS="${NSTEPS:-100}" \
        bash "$STG/run_timeline.sh" "$exe" "$RES/$side"
done

echo "############ Pass G vs the PREFIX form -- must be max_abs 0 ############"
REF="$PREFIX" H5MAXDIFF="$H5MAXDIFF" CFG_DIR="$SRC/../configs" \
    PASSES="G" NSTEPS_GATE=20 bash "$STG/run_exchange.sh" "$LIST" "$RES/gate_vs_prefix"

echo "############ the 7-case suite vs the PREFIX form ############"
MOBY_ROOT="$CODE_DIR" H5MAXDIFF="$H5MAXDIFF" NSTEPS=200 \
    bash "$STG/run_mapgate.sh" "$LIST" "$PREFIX" "$RES/suite_vs_prefix"

echo "=== divlist job finished ==="
