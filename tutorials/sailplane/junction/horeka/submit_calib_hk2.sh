#!/bin/bash
#SBATCH --job-name=moby_calib
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:4
#SBATCH --time=2-00:00:00
#SBATCH --partition=gpu-h100
#SBATCH --account=hk-project-p0024333
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# HoreKa 2 twin of submit_calib.sh (2026-11-12): the same resumable chunk
# chain -- restart from the latest <PREFIX>_<step>.h5, NSTEPS steps, done
# when a chunk ends on "main loop ended" before a full chunk, successor
# submitted at job start (afterany) -- on the `gpu-h100` partition (Teal
# nodes: 4 H100-94 GB, 64 cores; 16 cores per rank), the EasyBuild NVHPC
# 25.9 module (nvfortran, HPC-X mpirun) and its parallel HDF5 module, built
# for cc90 (build_hk2.sh builds on a CPU node ahead of the first chunk, so
# a queued GPU job never dies on a compile error). The c/2048 calibration
# box continued from the legacy cluster's t 18.76 snapshot:
#
#     W2=$(ws_find junction)
#     export CODE_DIR=$PROJECT/moby-junction RUN_DIR=$W2/calib_run \
#            STAGE_INI=calib_a1.5_y3.5_l7.ini PREFIX=calib_a1.5_y3.5_l7 NSTEPS=34000
#     cd $RUN_DIR && sbatch --nodes=6 --ntasks-per-node=4 --export=ALL submit_calib_hk2.sh
#
# (33 750 steps to t 20 fit in one 2-day chunk: 27 h at the A100 rate of
# 2.9 s/step, less on H100; the successor finds the run done and exits.)
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; STAGE_INI="${STAGE_INI:?}"; PREFIX="${PREFIX:?}"
NSTEPS="${NSTEPS:-20000}"

module purge
module load NVHPC/25.9-CUDA-12.9.1 HDF5/1.14.6-NVHPC-25.9-CUDA-12.9.1   # nvfortran + HPC-X (CUDA-aware Open MPI 4.1) + parallel HDF5 built with them
export HDF5_ROOT="${HDF5_ROOT:-$EBROOTHDF5}"
export OPENMP_OFFLOAD_FLAGS="${OPENMP_OFFLOAD_FLAGS:--mp=gpu -gpu=cc90}"   # H100 / H200 (the build job runs on a CPU node)
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

cache="$CODE_DIR/build_gpu/CMakeCache.txt"
if [ -f "$cache" ] && ! grep -qx "CMAKE_HOME_DIRECTORY:INTERNAL=$CODE_DIR" "$cache"; then
    rm -rf "$CODE_DIR/build_gpu"
fi
if [ ! -x "$CODE_DIR/build_gpu/moby_solve" ]; then
    echo "=== building $CODE_DIR ($(git -C "$CODE_DIR" rev-parse --short HEAD))"
    ( cd "$CODE_DIR" && ./compile.sh gpu ) > "$RUN_DIR/build.log" 2>&1 || { echo "BUILD FAILED"; exit 1; }
fi

cd "$RUN_DIR"
[ -f "$STAGE_INI" ] || { echo "no $STAGE_INI in $RUN_DIR"; exit 1; }
{
    echo "job    : ${SLURM_JOB_ID:-none}   $(date '+%F %T %Z')"
    echo "node   : $(hostname)"
    echo "code   : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "run    : $STAGE_INI, prefix $PREFIX, chunk $NSTEPS steps"
} | tee -a provenance.txt

latest_of() { ls -1 ${1}_[0-9]*.h5 2>/dev/null | sed 's/.*_\([0-9]*\)\.h5/\1 &/' | sort -n | tail -1 | cut -d' ' -f2; }
done_already() {  # the run reached t_final: the newest chunk ended ("main loop ended")
    # after FEWER than NSTEPS steps. The generic case prints no time rows, so
    # the step count comes from the snapshot names: the newest snapshot's
    # step minus the step the last chunk restarted from (status.txt holds
    # "exit <rc> <prefix> <restart file>"; a cold start is step 0).
    local log; log=$(ls -t ${PREFIX}_from_*.log 2>/dev/null | head -1); [ -n "$log" ] || return 1
    grep -q "main loop ended" "$log" || return 1
    local from; from=$(grep "^exit 0 $PREFIX " status.txt 2>/dev/null | tail -1 | awk '{print $4}')
    [ -n "$from" ] || return 1
    local s0=0; [ "$from" != "cold-start" ] && s0=$(echo "$from" | sed 's/.*_\([0-9]*\)\.h5/\1/')
    local s1; s1=$(latest_of "$PREFIX" | sed 's/.*_\([0-9]*\)\.h5/\1/'); [ -n "$s1" ] || return 1
    [ $((s1 - s0)) -lt "$NSTEPS" ] && [ $((s1 - s0)) -ge 0 ]
}
if done_already; then echo "$STAGE_INI is done: nothing to do"; exit 0; fi
if [ "$(grep "^exit [^0].* $PREFIX " status.txt 2>/dev/null | tail -2 | wc -l)" -ge 2 ] \
   && [ "$(grep " $PREFIX " status.txt | grep '^exit' | tail -2 | grep -c '^exit [^0]')" -ge 2 ]; then
    echo "the last two chunks of $PREFIX failed -- stopping the chain"; exit 1
fi

SCRIPT="$RUN_DIR/submit_calib.sh"
submit_successor() {
    local k=0 out
    while [ $k -lt 25 ]; do
        if out=$(sbatch --dependency=afterany:$SLURM_JOB_ID --job-name="$SLURM_JOB_NAME" \
                     --partition="${SLURM_JOB_PARTITION:-gpu-h100}" \
                     --nodes="${SLURM_NNODES:-1}" --ntasks-per-node="${SLURM_NTASKS_PER_NODE:-4}" \
                     --time="$(squeue -h -j "$SLURM_JOB_ID" -o %l)" --export=ALL "$SCRIPT" 2>&1); then
            echo "successor of $SLURM_JOB_ID: $out (attempt $((k + 1)))" | tee -a provenance.txt; return 0
        fi
        k=$((k + 1)); sleep 120
    done
    echo "CHAIN BROKEN after $k attempts ($out) -- resubmit by hand: cd $RUN_DIR && sbatch --export=ALL $SCRIPT" | tee -a provenance.txt
}
[ -n "${NO_CHAIN:-}" ] || submit_successor &

latest=$(latest_of "$PREFIX")
chunk="chunk_$PREFIX.ini"
sed -e "s/^nsteps = .*/nsteps = $NSTEPS/" "$STAGE_INI" > "$chunk"
if [ -n "$latest" ]; then
    if grep -q '^\[restart\]' "$chunk"; then sed -i "/^\[restart\]/,/^$/ s|^file = .*|file = $latest|" "$chunk"
    else printf '\n[restart]\nfile = %s\n' "$latest" >> "$chunk"; fi
fi
restart=${latest:-cold-start}
log="${PREFIX}_from_$(basename "$restart" .h5).log"
echo "=== $NSTEPS steps from $restart"
echo "start $PREFIX $restart $(date '+%F %T')" >> status.txt
mpirun -n "${SLURM_NTASKS:-4}" --bind-to core --map-by numa "$CODE_DIR/build_gpu/moby_solve" "$chunk" > "$log" 2>&1
rc=$?
grep -q "main loop ended" "$log" || rc=$(( rc == 0 ? 99 : rc ))
echo "exit $rc $PREFIX $restart $(date '+%F %T')" >> status.txt
grep -E "seconds_per_step" "$log" | tail -1
tail -3 status.txt
