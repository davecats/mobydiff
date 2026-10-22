#!/bin/bash
#SBATCH --job-name=moby_naca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=01:00:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# ONE STAGE of the NACA 0012 from-scratch protocol (tutorials/naca/rans/run_case.sh
# scratch: L10 transient -> L11 -> the nose band) as a RESUMABLE, SELF-CHAINING
# chunk chain on 4 A100 ranks, for the 2026-10-21 re-validation (decision 2:
# re-converge from scratch on the nose grid). The stage ini is composed on the
# workstation exactly as run_case.sh composes it (the steady_tol / steady_samples /
# force_sample_interval keys of the stage, t_final as the safety net) and shipped
# with the case file it names; the interpolation between stages
# (interp_restart.py, h5py) runs on the workstation, so each stage is its own chain:
#
#     export CODE_DIR=<worktree> RUN_DIR=<run dir> HDF5_ROOT=$HOME/hdf5 \
#            STAGE_INI=stage1.ini PREFIX=.l10_stage [NSTEPS=15000]
#     cd $RUN_DIR && sbatch --export=ALL submit_naca_stage.sh
#
# A chunk restarts from the latest <PREFIX>_<step>.h5 (the first chunk of stage 1
# cold-starts; stages 2/3 start from the shipped interpolated restart the ini
# names), runs NSTEPS steps and writes one snapshot at its end. The stage is DONE
# when a chunk's log says "airfoil: steady state reached" (the solver's own
# stop, which also writes the field) or the run reached the ini's t_final
# (the log's "main loop ended" without a full chunk). NOTE the steady window
# (steady_samples consecutive samples) restarts at every chunk boundary, so a
# chunk boundary can delay the stop by up to one window. The successor is
# submitted at the START of each job (afterany) by a retry loop that alternates
# the two accounts; two failed chunks in a row stop the chain. NO_CHAIN=1 runs
# one chunk.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; STAGE_INI="${STAGE_INI:?}"; PREFIX="${PREFIX:?}"
NSTEPS="${NSTEPS:-15000}"

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
    echo "stage  : $STAGE_INI, prefix $PREFIX, chunk $NSTEPS steps"
} | tee -a provenance.txt

latest_of() { ls -1 ${1}_[0-9]*.h5 2>/dev/null | sed 's/.*_\([0-9]*\)\.h5/\1 &/' | sort -n | tail -1 | cut -d' ' -f2; }
done_already() {  # the stage ended: the steady stop or t_final, in the newest chunk log
    local log; log=$(ls -t ${PREFIX}_from_*.log 2>/dev/null | head -1); [ -n "$log" ] || return 1
    grep -q "steady state reached" "$log" && return 0
    # t_final reached: the loop ended before the chunk's nsteps
    grep -q "main loop ended" "$log" && [ "$(grep -c '^ *[0-9]\+ \+[0-9.E+-]\+ \+' "$log")" -lt 2 ] && return 0
    local tf; tf=$(grep '^t_final' "$STAGE_INI" | awk '{print $3}')
    local tl; tl=$(grep -E '^ *[0-9]+ +[0-9.E+-]+ +[0-9.E+-]+ +[0-9.E+-]+ +[0-9.E+-]+ *$' "$log" | tail -1 | awk '{print $2}')
    [ -n "$tl" ] && awk -v a="$tl" -v b="$tf" 'BEGIN{exit !(a >= b - 1e-6)}' && return 0
    return 1
}
if done_already; then echo "stage $STAGE_INI is done: nothing to do"; exit 0; fi
if [ "$(grep '^exit' status.txt 2>/dev/null | tail -2 | grep -c '^exit [^0]')" -ge 2 ]; then
    echo "the last two chunks failed -- stopping the chain"; exit 1
fi

SCRIPT="$RUN_DIR/submit_naca_stage.sh"
submit_successor() {
    local k=0 out acct
    while [ $k -lt 25 ]; do
        acct=hk-project-exasim; [ $((k % 2)) -eq 1 ] && acct=hk-project-p0024333
        # the successor inherits the partition this job ran on (dev_accelerated 1 h
        # chunks, or accelerated 5 h chunks when the dev submit cap is in the way)
        local part="${SLURM_JOB_PARTITION:-dev_accelerated}" tlim=01:00:00
        [ "$part" = accelerated ] && tlim=05:00:00
        if out=$(sbatch --dependency=afterany:$SLURM_JOB_ID --job-name="$SLURM_JOB_NAME" \
                     --partition="$part" --time="$tlim" --account="$acct" --export=ALL "$SCRIPT" 2>&1); then
            echo "successor of $SLURM_JOB_ID: $out (attempt $((k + 1)), account $acct)" | tee -a provenance.txt; return 0
        fi
        k=$((k + 1)); sleep 120
    done
    echo "CHAIN BROKEN after $k attempts ($out) -- resubmit by hand: cd $RUN_DIR && sbatch --export=ALL $SCRIPT" | tee -a provenance.txt
}
[ -n "${NO_CHAIN:-}" ] || submit_successor &

latest=$(latest_of "$PREFIX")
sed -e "s/^nsteps = .*/nsteps = $NSTEPS/" -e "s/^field_interval = .*/field_interval = 0/" "$STAGE_INI" > chunk.ini
if [ -n "$latest" ]; then
    if grep -q '^\[restart\]' chunk.ini; then sed -i "/^\[restart\]/,/^$/ s|^file = .*|file = $latest|" chunk.ini
    else printf '\n[restart]\nfile = %s\n' "$latest" >> chunk.ini; fi
fi
restart=${latest:-cold-start}
log="${PREFIX}_from_$(basename "$restart" .h5).log"
echo "=== $NSTEPS steps from $restart"
echo "start $PREFIX $restart $(date '+%F %T')" >> status.txt
mpirun -n 4 --bind-to core --map-by numa "$CODE_DIR/build_gpu/moby_solve" chunk.ini > "$log" 2>&1
rc=$?
grep -q "main loop ended" "$log" || rc=$(( rc == 0 ? 99 : rc ))
echo "exit $rc $PREFIX $restart $(date '+%F %T')" >> status.txt
grep -E "seconds_per_step|steady state reached" "$log" | tail -2
tail -3 status.txt
