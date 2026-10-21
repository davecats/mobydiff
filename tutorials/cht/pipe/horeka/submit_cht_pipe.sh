#!/bin/bash
#SBATCH --job-name=moby_cht_pipe
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=05:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The conjugate pipe (tutorials/cht/pipe), leg D of run_pipe.sh -- the
# statistics measurement -- re-run on HoreKa for the 2026-10-21 re-validation
# from the campaign's settled state (p_pr_settle2_30000.h5 with its pressure
# zeroed: that snapshot was converged at niter 6), 4 A100 ranks, in RESUMABLE,
# SELF-CHAINING chunks of NSTEPS steps. The projection is the one the
# istmcorax timing elected (decision 4 of 2026-10-21): the chunk ini is
# generated on the workstation (run_pipe.sh's `thermal` composition, the
# sink scale of measure_generation.py applied, the projection keys set) and
# shipped as chunk_template.ini with @NSTEPS@ and @RESTART@ left open.
#
#     export CODE_DIR=<worktree> RUN_DIR=<run dir> HDF5_ROOT=$HOME/hdf5 [NSTEPS=36000]
#     cd $RUN_DIR && sbatch --export=ALL submit_cht_pipe.sh
#
# RUN_DIR holds start_field.h5 (step 30000), pipe_prod.h5 (the case file) and
# chunk_template.ini. Steps continue across restarts: the leg runs from step
# 30000 to END (198000 = the 168 000-step window of the README); each chunk
# restarts from the latest p_pr_stat_<step>.h5 and writes snapshots every
# FIELD steps (12 000, as the campaign did, for pipe_stats.py) plus one at
# its end; the scalar statistics (p_pr_stats.h5) and the heat file continue
# across chunks as long as they stay. The successor is submitted at the START
# of each job (afterany), by a retry loop; a job that finds END reached, or two
# failed chunks in a row, exits without a successor. NO_CHAIN=1 runs one chunk.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
NSTEPS="${NSTEPS:-36000}"; FIELD="${FIELD:-12000}"
START=30000; END=198000; PREFIX=p_pr_stat

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
for f in start_field.h5 pipe_prod.h5 chunk_template.ini; do
    [ -f "$f" ] || { echo "no $f in $RUN_DIR"; exit 1; }
done
{
    echo "job    : ${SLURM_JOB_ID:-none}   $(date '+%F %T %Z')"
    echo "node   : $(hostname)"
    echo "code   : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "chunk  : $NSTEPS steps, snapshots every $FIELD; $PREFIX $START -> $END"
    grep -A4 '^\[pressure\]' chunk_template.ini | tr '\n' ' '; echo
} | tee -a provenance.txt

latest_of() { ls -1 ${1}_[0-9]*.h5 2>/dev/null | sed 's/.*_\([0-9]*\)\.h5/\1 &/' | sort -n | tail -1 | cut -d' ' -f2; }
step_of() { echo "$1" | sed 's/.*_\([0-9]*\)\.h5/\1/'; }

latest=$(latest_of $PREFIX); restart=${latest:-start_field.h5}
step=$START; [ -n "$latest" ] && step=$(step_of "$latest")
if [ "$step" -ge "$END" ]; then echo "$PREFIX at END $END: nothing to do"; exit 0; fi
if [ "$(grep '^exit' status.txt 2>/dev/null | tail -2 | grep -c '^exit [^0]')" -ge 2 ]; then
    echo "the last two chunks failed -- stopping the chain"; exit 1
fi
n=$(( END - step )); [ "$n" -gt "$NSTEPS" ] && n=$NSTEPS

SCRIPT="$RUN_DIR/submit_cht_pipe.sh"
submit_successor() {
    # the project's GPU hours may run out mid-chain: odd attempts use the
    # script's account, even ones the fallback the user named on 2026-10-21
    local k=0 out acct
    while [ $k -lt 25 ]; do
        acct=hk-project-exasim; [ $((k % 2)) -eq 1 ] && acct=hk-project-p0024333
        if out=$(sbatch --dependency=afterany:$SLURM_JOB_ID --job-name="$SLURM_JOB_NAME" \
                     --account="$acct" --export=ALL "$SCRIPT" 2>&1); then
            echo "successor of $SLURM_JOB_ID: $out (attempt $((k + 1)), account $acct)" | tee -a provenance.txt; return 0
        fi
        k=$((k + 1)); sleep 120
    done
    echo "CHAIN BROKEN after $k attempts ($out) -- resubmit by hand: cd $RUN_DIR && sbatch --export=ALL $SCRIPT" | tee -a provenance.txt
}
[ -n "${NO_CHAIN:-}" ] || submit_successor &

sed -e "s|@NSTEPS@|$n|" -e "s|@RESTART@|$restart|" -e "s|^field_interval = .*|field_interval = $FIELD|" \
    chunk_template.ini > chunk.ini
if grep -v '^;' chunk.ini | grep -q '@'; then echo "unsubstituted placeholder"; exit 1; fi

log="${PREFIX}_from_$(basename "$restart" .h5).log"
echo "=== $n steps from $restart (step $step -> $((step + n)) of $END)"
echo "start $PREFIX $restart $(date '+%F %T')" >> status.txt
mpirun -n 4 --bind-to core --map-by numa "$CODE_DIR/build_gpu/moby_solve" chunk.ini > "$log" 2>&1
rc=$?
grep -q "main loop ended" "$log" || rc=$(( rc == 0 ? 99 : rc ))   # the solver's own word, not mpirun's
echo "exit $rc $PREFIX $restart $(date '+%F %T')" >> status.txt
grep -E "seconds_per_step" "$log" | tail -1
tail -3 status.txt
