#!/bin/bash
#SBATCH --job-name=moby_divlist
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
# The divergence subset as an INDEX LIST instead of a prefix.
#
# THREE columns, because the two claims need different references:
#   base   3c2903a  -- before any divergence work. Red-black must return to THIS.
#   prefix 95312d7  -- the prefix form. Jacobi must beat THIS.
#   list   HEAD     -- the index list.
#
# Rank counts 1, 4, 8 and the RED-BLACK config included, because that is where
# the prefix's cost showed: 1 rank has no message saving to hide it and
# red-black never runs the reduced round at all
# (results_horeka_2026-09-17.md section 4).
#
# Readings pre-registered in PREREGISTERED_divlist.md.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"; BASE_DIR="${BASE_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/divlist_staged"; rm -rf "$STG"; mkdir -p "$STG"
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
H5MAXDIFF="$CODE_DIR/tools/h5maxdiff"
[ -x "$H5MAXDIFF" ] || gcc -O2 -o "$H5MAXDIFF" "$CODE_DIR/tools/h5maxdiff.c" \
    -I"$HDF5_ROOT/include" -L"$HDF5_ROOT/lib" -lhdf5 -Wl,-rpath,"$HDF5_ROOT/lib" || exit 1

RES="${RESDIR:-$RUN_DIR/results_divlist}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "list   : $(git -C "$CODE_DIR" rev-parse HEAD)  (index list)"
    echo "prefix : $(git -C "$REF_DIR" rev-parse HEAD)  (prefix form, 95312d7)"
    echo "base   : $(git -C "$BASE_DIR" rev-parse HEAD)  (before any divergence work)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

SPECS="rect_jacobi:8:2 refined_yp82_rect_jacobi:8:2 refined_yp82_rect_redblack:8:2 \
rect_jacobi:4:1 refined_yp82_rect_jacobi:4:1 refined_yp82_rect_redblack:4:1 \
rect_jacobi:1:1 refined_yp82_rect_jacobi:1:1 refined_yp82_rect_redblack:1:1"
for side in list prefix base; do
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
