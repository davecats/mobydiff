#!/bin/bash
#SBATCH --job-name=moby_rdblk
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:40:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# rdenom recomputed only on the blocks that hold a body.
#
# The overheadTest configs are all BODY-FREE, so they cannot measure this: the
# cases here are the two body geometries that bracket production --
#   les_ibm    plane walls spanning the domain  (256/640 blocks on CPU)
#   sailplane  a compact body in a large domain (48/4500 at nb=10)
# ref = b9414bd (before any rdenom work), new = HEAD.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
module purge; module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

for d in "$REF_DIR" "$CODE_DIR"; do
    cache="$d/build_gpu/CMakeCache.txt"
    if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$d" "$cache"; then
        echo "=== $d: stale build_gpu cache -- discarding"; rm -rf "$d/build_gpu"; fi
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD))"
    ( cd "$d" && ./compile.sh gpu ) || exit 1
done
NEW="$CODE_DIR/build_gpu/moby_solve"; REF="$REF_DIR/build_gpu/moby_solve"
H5MAXDIFF="$CODE_DIR/tools/h5maxdiff"
[ -x "$H5MAXDIFF" ] || mpicc -O2 -o "$H5MAXDIFF" "$CODE_DIR/tools/h5maxdiff.c" \
    -I"$HDF5_ROOT/include" -L"$HDF5_ROOT/lib" -lhdf5 -Wl,-rpath,"$HDF5_ROOT/lib" || exit 1

RES="${RESDIR:-$RUN_DIR/results_rdblk}"; mkdir -p "$RES"
echo "job ${SLURM_JOB_ID:-none}  new $(git -C "$CODE_DIR" rev-parse --short HEAD)  ref $(git -C "$REF_DIR" rev-parse --short HEAD)" | tee "$RES/provenance.txt"

# case : ini (relative to the repo) : ranks : extra ini lines
run_case() {
    local name="$1" ini="$2" ranks="$3" extra="$4" side="$5" exe="$6"
    local run="$RES/${name}_${side}"
    rm -rf "$run"; cp -a "$CODE_DIR/$(dirname "$ini")" "$run"
    cp "$CODE_DIR/$ini" "$run/config.ini"
    sed -i -e "s/^nsteps *=.*/nsteps = ${NSTEPS:-50}/" -e "s/^t_final *=.*/t_final = 0.0/" "$run/config.ini"
    printf '%b' "$extra" >> "$run/config.ini"
    ( cd "$run" && mpirun -n "$ranks" --bind-to core --map-by numa "$exe" config.ini > run.log 2>&1 ) \
        || echo "    $name $side FAILED"
}
PROF='\n[output]\nprofile = true\n'
for side in ref new; do
    exe="$REF"; [ "$side" = new ] && exe="$NEW"
    run_case les_ibm  validation/channel_interface/les_ibm/channel_ibm.ini 1 "$PROF" "$side" "$exe"
    run_case sailplane tutorials/sailplane/input.ini 4 "\n[blocks]\nnb = 10\n$PROF" "$side" "$exe"
done

echo "=== body-block fractions, buckets and bit-exactness ==="
for name in les_ibm sailplane; do
    echo "--- $name"
    grep -h "rdenom recomputed" "$RES/${name}_new/run.log" | sed 's/^/    /'
    for side in ref new; do
        printf "    %-4s " "$side"
        grep -E "^timing: nsteps" "$RES/${name}_${side}/run.log" | awk '{printf "s/step %s  ", $NF}'
        grep -E "^proj_timing: setup " "$RES/${name}_${side}/run.log" | awk '{printf "setup/step %s", $(NF-2)}'
        echo
    done
    a=$(ls "$RES/${name}_ref"/*_[0-9]*.h5 2>/dev/null | tail -1)
    b=$(ls "$RES/${name}_new"/*_[0-9]*.h5 2>/dev/null | tail -1)
    [ -n "$a" ] && [ -n "$b" ] && "$H5MAXDIFF" "$a" "$b" | tail -3 | sed 's/^/    /'
done
echo "=== rdblk job finished ==="
