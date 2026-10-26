#!/bin/bash
#SBATCH --job-name=moby_calib
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=05:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-p0024333
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# ONE calibration-box run of the sailplane junction study (docs/
# next_session_sailplane_junction.md section 4.2) on 4 A100 ranks, as a
# resumable chunk of NSTEPS steps (one `accelerated` job of 5 h per case:
# the c/512 band makes 190 M cells, about 0.7 s/step on 4 A100 by the NACA
# scaling, 15 000 steps to t_final; `sbatch --partition=dev_accelerated
# --time=01:00:00` for a short test): tutorials/naca/rans/horeka/
# submit_naca_stage.sh with the airfoil-specific steady-state test replaced
# by "the run reached the ini's t_final" (the generic case has no runtime
# force output), the project account of decision D5, and no second account
# to alternate with. The ini (field_prefix = PREFIX) and the case file it
# names are prepared on the workstation and shipped to RUN_DIR:
#
#     export CODE_DIR=$W/moby-junction RUN_DIR=$W/calib_run HDF5_ROOT=$HOME/hdf5 \
#            STAGE_INI=calib_a1.5_y3.5.ini PREFIX=calib_a1.5_y3.5 NSTEPS=20000 NO_CHAIN=1
#     cd $RUN_DIR && sbatch --export=ALL submit_calib.sh
#
# A chunk restarts from the latest <PREFIX>_<step>.h5 (the first cold-starts),
# runs NSTEPS steps keeping the ini's field_interval (the read-out wants the
# last few snapshots, not only the final one); the run is DONE when
# the newest chunk log says "main loop ended" before a full chunk (t_final
# trims the last step). Without NO_CHAIN the successor is submitted at the
# start of each job (afterany, retried against the dev submit cap); two
# failed chunks in a row stop the chain.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; STAGE_INI="${STAGE_INI:?}"; PREFIX="${PREFIX:?}"
NSTEPS="${NSTEPS:-20000}"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
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
done_already() {  # the run reached t_final: the newest chunk log ended before a full chunk
    local log; log=$(ls -t ${PREFIX}_from_*.log 2>/dev/null | head -1); [ -n "$log" ] || return 1
    grep -q "main loop ended" "$log" || return 1
    local tf; tf=$(grep '^t_final' "$STAGE_INI" | awk '{print $3}')
    local tl; tl=$(grep -E '^ *[0-9]+ +[0-9.E+-]+ +[0-9.E+-]+ +[0-9.E+-]+ +[0-9.E+-]+ *$' "$log" | tail -1 | awk '{print $2}')
    [ -n "$tl" ] && awk -v a="$tl" -v b="$tf" 'BEGIN{exit !(a >= b - 1e-6)}' && return 0
    return 1
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
                     --partition="${SLURM_JOB_PARTITION:-dev_accelerated}" --export=ALL "$SCRIPT" 2>&1); then
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
mpirun -n 4 --bind-to core --map-by numa "$CODE_DIR/build_gpu/moby_solve" "$chunk" > "$log" 2>&1
rc=$?
grep -q "main loop ended" "$log" || rc=$(( rc == 0 ? 99 : rc ))
echo "exit $rc $PREFIX $restart $(date '+%F %T')" >> status.txt
grep -E "seconds_per_step" "$log" | tail -1
tail -3 status.txt
