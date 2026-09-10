#!/bin/bash
#SBATCH --job-name=moby_class
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:30:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#
# TEST THE CORRECTED MODEL. The confound matrix showed the slow cases are the
# ones where the two ends of a cross-node link sit in DIFFERENT NUMA/NIC
# affinity classes (dev3<->dev0, dev0<->dev2), while BOTH-good (dev1<->dev0,
# dev0<->dev0) and BOTH-BAD (dev3<->dev3) are all fast. So the rule is not
# "both ends need a local HCA" -- it is "both ends must MATCH".
#
# Predictions, written before the runs:
#   ref  GPUS=1,0,3,2  -> A:dev2 B:dev1  MIXED  -> SLOW
#   ref  GPUS=2,0,1,3  -> A:dev3 B:dev2  both NUMA1 -> FAST
#   new  MOBY_GPU_ORDER=k,... -> at 8x2 BOTH cross-node ranks take pref[0]=k,
#        so every k is matched -> all four FAST. If some k is slow, "matched" is
#        not the whole rule either.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
EXE="$CODE_DIR/build_gpu/moby_solve"; REF="$CODE_DIR/build_gpu/moby_solve.ref"
module purge; module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1
RES="${RESDIR:-$RUN_DIR/results_class}"; mkdir -p "$RES"; cd "$RES" || exit 1
cat > pin.sh <<'PIN'
#!/bin/bash
if [ -n "${GPUS:-}" ]; then
    IFS=',' read -r -a d <<< "$GPUS"
    lr="${OMPI_COMM_WORLD_LOCAL_RANK:-${SLURM_LOCALID:-0}}"
    export CUDA_VISIBLE_DEVICES="${d[$(( lr % ${#d[@]} ))]}"
fi
exec "$@"
PIN
chmod +x pin.sh

go() {  # go <name> <exe> <env assignments...>
    local d="$1" x="$2"; shift 2
    mkdir -p "$d"; cp "$RUN_DIR/configs/rect_jacobi.ini" "$d/config.ini"
    sed -i -e 's/^profile *=.*/profile = true/' -e 's/^nsteps *=.*/nsteps = 200/' \
           -e 's/^runtime_interval *=.*/runtime_interval = 200/' "$d/config.ini"
    ( cd "$d" && env "$@" mpirun -n 8 --map-by ppr:4:node --bind-to core \
        -x GPUS -x MOBY_GPU_ORDER "$RES/pin.sh" "$x" config.ini > run.log 2>&1 < /dev/null )
    printf "%-26s " "$d"
    grep -m1 "^timing: nsteps" "$d/run.log" | awk '{printf "%s  ", $NF}'
    grep -m1 "^exch_timing: mpi_wait" "$d/run.log" | awk '{printf "wait %7.1f us  ", $8/$4*1e6}'
    grep -m1 "gpu binding" "$d/run.log" | sed 's/.*node)://;s/MOBY.*//' | tr -d '\n'
    echo
    rm -f "$d"/overhead_*.h5
}

echo "--- pinned pairs (ref binary): does MIXED class predict slow? ---"
go pin_mixed_2v1    "$REF" GPUS=1,0,3,2
go pin_bothbad_3v2  "$REF" GPUS=2,0,1,3
go pin_bothgood_1v0 "$REF" GPUS=0,2,3,1
echo "--- new binary, first preference swept (both ends always take pref[0]) ---"
for k in 0 1 2 3; do go "auto_pref$k" "$EXE" MOBY_GPU_ORDER=$k; done
echo "=== class job finished ==="
