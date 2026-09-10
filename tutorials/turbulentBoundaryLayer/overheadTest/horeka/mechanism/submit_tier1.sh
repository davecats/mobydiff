#!/bin/bash
#SBATCH --job-name=moby_tier1
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:45:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# Gate + measure the topology-aware device selection (comm.f90
# select_target_device). Three questions, in order:
#
#   1. Is it bit-exact?  Same case, same placement, reference binary (built
#      before the change, kept as build_gpu/moby_solve.ref) against the new one.
#      Different ranks land on different physical GPUs, but the cards are
#      identical A100s and the arithmetic does not move, so the fields must
#      match EXACTLY. Compared with tools/h5maxdiff (HoreKa has no h5diff/h5py).
#   2. Does the automatic map reproduce the hand-tuned optimum?  The new binary
#      with NO environment set should land on the map that measured 0,2,3,1.
#   3. Is a single-node run unchanged?  With no off-node peers every rank has
#      degree 0, the order collapses to local-rank order, and the mapping is the
#      old round-robin -- so single-node users see no change at all.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
EXE="$CODE_DIR/build_gpu/moby_solve"
REF="$CODE_DIR/build_gpu/moby_solve.ref"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n
export OMP_NUM_THREADS=1

[ -x "$REF" ] || { echo "ERROR: no reference binary $REF" >&2; exit 1; }
cd "$CODE_DIR" || exit 1
echo "commit under test: $(git rev-parse HEAD)  dirty=$(git status --porcelain -uno | wc -l)"
./compile.sh gpu || exit 1
[ -x "$EXE" ] || { echo "ERROR: build produced no $EXE" >&2; exit 1; }

RES="${RESDIR:-$RUN_DIR/results_tier1}"; mkdir -p "$RES"
CFG="$RUN_DIR/configs"
cd "$RES" || exit 1
{ echo "job    : ${SLURM_JOB_ID:-none}"; echo "date   : $(date '+%F %T %Z')"
  echo "nodes  : ${SLURM_JOB_NUM_NODES:-?} (${SLURM_JOB_NODELIST:-?})"
  echo "commit : $(git -C "$CODE_DIR" rev-parse HEAD)"; } | tee provenance.txt

mkcfg() {  # mkcfg <name> <config> <nsteps>
    mkdir -p "$1"; cp "$CFG/$2.ini" "$1/config.ini"
    sed -i 's/^profile *=.*/profile = true/' "$1/config.ini"
    sed -i -e "s/^nsteps *=.*/nsteps = $3/" -e "s/^runtime_interval *=.*/runtime_interval = $3/" \
           -e "s/^field_interval *=.*/field_interval = $3/" "$1/config.ini"
}
run() {    # run <dir> <exe> <ranks> <nodes> [env...]
    local d="$1" x="$2" n="$3" N="$4"; shift 4
    ( cd "$d" && env "$@" mpirun -n "$n" --map-by "ppr:$((n/N)):node" --bind-to core \
        --display-map "$x" config.ini > run.log 2>&1 < /dev/null )
    printf "  %-22s " "$(basename "$d")"
    grep -m1 "gpu binding" "$d/run.log" | sed 's/^ *//' || echo "(no binding line)"
    grep -m1 "^timing: nsteps" "$d/run.log" | awk '{printf "      s/step %s", $NF}'
    grep -m1 "^exch_timing: mpi_wait" "$d/run.log" | awk '{printf "   wait/round %.1f us\n", $8/$4*1e6}'
}

echo "=== 1. BIT-EXACTNESS: reference vs new, rect_jacobi 8x2, 50 steps ==="
mkcfg exact_ref rect_jacobi 50; mkcfg exact_new rect_jacobi 50
run exact_ref "$REF" 8 2
run exact_new "$EXE" 8 2
"$CODE_DIR/tools/h5maxdiff" exact_ref/overhead_50.h5 exact_new/overhead_50.h5

echo "=== 2. PERFORMANCE, reference vs new, back to back in one allocation ==="
for spec in "r8_rect rect_jacobi 8 2" "r16_rect rect_jacobi 16 4" \
            "r16_refined refined_yp82_rect_jacobi 16 4" "r16_base base_jacobi 16 4"; do
    set -- $spec
    mkcfg "${1}_ref" "$2" 200; run "${1}_ref" "$REF" "$3" "$4"
    mkcfg "${1}_new" "$2" 200; run "${1}_new" "$EXE" "$3" "$4"
done

echo "=== 3. SINGLE NODE must be unchanged (identity mapping, bit-exact) ==="
mkcfg single_ref rect_jacobi 50; mkcfg single_new rect_jacobi 50
run single_ref "$REF" 4 1
run single_new "$EXE" 4 1
"$CODE_DIR/tools/h5maxdiff" single_ref/overhead_50.h5 single_new/overhead_50.h5

echo "=== 4. MOBY_GPU_ORDER override takes effect ==="
mkcfg override rect_jacobi 200
run override "$EXE" 8 2 MOBY_GPU_ORDER=3,2,1,0

echo "=== tier1 job finished ==="
