#!/bin/bash
#SBATCH --job-name=moby_confound
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:30:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#
# SEPARATE THE CONFOUND. The 2026-09-10 sweep compared a run with 4 devices
# visible and no pinning ("default", 745.7 us) against runs that set
# CUDA_VISIBLE_DEVICES per rank ("gpumap", 70.2 us). That changed the PERMUTATION
# and the DEVICE VISIBILITY together, and attributed the whole 10.6x to the
# permutation. Then MOBY_GPU_ORDER=3,2,1,0 -- which puts BOTH cross-node ranks on
# the supposedly bad cards -- came back fast, which the permutation story cannot
# explain.
#
# 2 binaries x 3 pinnings, one allocation, identical case and placement:
#   none      4 devices visible, the solver picks
#   identity  1 visible, local rank k -> physical k  (the old mapping, pinned)
#   boundary  1 visible, permutation 0,2,3,1
# If identity-pinned is FAST, the permutation was never the cause and visibility
# is. If it is SLOW, the permutation is real and the override result needs its
# own explanation.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
EXE="$CODE_DIR/build_gpu/moby_solve"; REF="$CODE_DIR/build_gpu/moby_solve.ref"
module purge; module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1
RES="${RESDIR:-$RUN_DIR/results_confound}"; mkdir -p "$RES"; cd "$RES" || exit 1
{ echo "job   : ${SLURM_JOB_ID:-none}"; echo "nodes : ${SLURM_JOB_NODELIST:-?}"
  echo "commit: $(git -C "$CODE_DIR" rev-parse HEAD)"; } | tee provenance.txt

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

printf "%-22s %-16s %-14s %s\n" run binary pinning "s/step / wait per round"
for b in ref new; do
  x="$REF"; [ "$b" = new ] && x="$EXE"
  for p in none identity boundary; do
    d="${b}_${p}"; mkdir -p "$d"; cp "$RUN_DIR/configs/rect_jacobi.ini" "$d/config.ini"
    sed -i -e 's/^profile *=.*/profile = true/' -e 's/^nsteps *=.*/nsteps = 200/' \
           -e 's/^runtime_interval *=.*/runtime_interval = 200/' "$d/config.ini"
    g=""; [ "$p" = identity ] && g="0,1,2,3"; [ "$p" = boundary ] && g="0,2,3,1"
    ( cd "$d" && env ${g:+GPUS=$g} mpirun -n 8 --map-by ppr:4:node --bind-to core \
        -x GPUS --display-map "$RES/pin.sh" "$x" config.ini > run.log 2>&1 < /dev/null )
    printf "%-22s %-16s %-14s " "$d" "$b" "${g:-4 visible}"
    grep -m1 "^timing: nsteps" "$d/run.log" | awk '{printf "%s  ", $NF}'
    grep -m1 "^exch_timing: mpi_wait" "$d/run.log" | awk '{printf "wait %.1f us  ", $8/$4*1e6}'
    grep -m1 "devices available" "$d/run.log" | awk '{printf "(ndev %s)", $NF}'
    grep -m1 "gpu binding" "$d/run.log" | sed 's/.*node)://' | tr -d '\n'
    echo
    rm -f "$d"/overhead_*.h5
  done
done
echo "=== confound job finished ==="
