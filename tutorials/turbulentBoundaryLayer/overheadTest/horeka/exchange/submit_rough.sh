#!/bin/bash
#SBATCH --job-name=moby_rough
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:50:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The rdenom narrowing on a PRODUCTION-SHAPED BODY case: configs/rough_jacobi,
# rect_jacobi's grid/blocks/flow plus a 3D egg-carton wall (25% body blocks).
#
# ref is bench/rdenom-always -- HEAD with the narrowing alone disabled. The
# pre-rdenom commits cannot be the reference here: they have no [ibm] wall_shape
# key and cannot run the case at all.
#
# rect_jacobi runs alongside as the BODY-FREE control: it must be unchanged
# between the two binaries except for the narrowing it already had.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/rough_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC"/run_timeline.sh "$SRC"/run_mapgate.sh "$STG/"

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

RES="${RESDIR:-$RUN_DIR/results_rough}"; mkdir -p "$RES"
{ echo "job    : ${SLURM_JOB_ID:-none}"
  echo "date   : $(date '+%F %T %Z')"
  echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  (rdenom narrowed)"
  echo "ref    : $(git -C "$REF_DIR" rev-parse HEAD)  (bench/rdenom-always)"
} | tee "$RES/provenance.txt"

SPECS="rough_jacobi:8:2 rect_jacobi:8:2 rough_jacobi:4:1 rect_jacobi:4:1"
for side in ref new; do
    exe="$REF"; [ "$side" = new ] && exe="$NEW"
    echo "############ timeline: $side ############"
    CFG_DIR="$SRC/../configs" TL_PASSES="A" SPECS="$SPECS" NSTEPS="${NSTEPS:-200}" \
        bash "$STG/run_timeline.sh" "$exe" "$RES/$side"
done

echo "=== body-block fraction, buckets, and the body-free control ==="
grep -h "rdenom recomputed" "$RES/new"/host_rough_jacobi_r4N1/run.log | sed 's/^/    /'
for d in "$RES"/ref/host_*; do
    n=$(basename "$d"); o="$RES/new/$n"
    printf "%-38s " "${n#host_}"
    a=$(grep -E "^timing: nsteps" "$d/run.log" | awk '{print $NF}')
    b=$(grep -E "^timing: nsteps" "$o/run.log" | awk '{print $NF}')
    sa=$(grep -E "^proj_timing: setup " "$d/run.log" | awk '{print $(NF-2)}')
    sb=$(grep -E "^proj_timing: setup " "$o/run.log" | awk '{print $(NF-2)}')
    python3 -c "
a=$a; b=$b; sa=$sa; sb=$sb
print(f'step {1e3*a:8.3f} -> {1e3*b:8.3f} ms ({100*(b-a)/a:+6.2f}%)   setup {1e3*sa:7.3f} -> {1e3*sb:7.3f} ({100*(sb-sa)/sa:+7.2f}%)')"
done
echo "=== the 7-case suite (must be untouched by all of this) ==="
MOBY_ROOT="$CODE_DIR" H5MAXDIFF="$H5MAXDIFF" NSTEPS=200 \
    bash "$STG/run_mapgate.sh" "$NEW" "$REF" "$RES/suite"
echo "=== rough job finished ==="
