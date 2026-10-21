#!/bin/bash
#SBATCH --job-name=moby_cht_channel
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
# The Flageul-matched conjugate channel (tutorials/cht/channel) re-run on
# HoreKa for the 2026-10-21 re-validation: the `settle` and `clean` legs of
# run_flageul.sh, from the campaign's settled state (Cdev_65000.h5 with its
# pressure zeroed, since that snapshot was converged at niter 6), on 4 A100
# ranks, in RESUMABLE, SELF-CHAINING chunks of NSTEPS steps.
#
#     export CODE_DIR=<worktree> RUN_DIR=<run dir> HDF5_ROOT=$HOME/hdf5 [NSTEPS=20000]
#     cd $RUN_DIR && sbatch --export=ALL submit_cht_channel.sh
#
# RUN_DIR holds start_field.h5 (step 65000), cht149_flageul.h5 (the case
# file, prepared on the workstation) and ynodes_f149.dat; the chunk inis are
# generated from the tutorial's cht149_flageul.ini. Step numbers continue
# across restarts, so the legs are ABSOLUTE step targets:
#
#   leg 1  Cdev2   settle, statistics OFF        65000 -> 130000
#   leg 2  Cstat   clean, statistics FROM ZERO  130000 -> 780000
#          (scalar stats flageul_clean.h5, channel stats Cstat_vel.h5; the
#          accumulators continue across chunks as long as the files stay)
#
# A chunk restarts from the latest snapshot of the current leg (the first
# clean chunk from leg 1's final field), runs min(NSTEPS, remaining) steps
# and writes ONE snapshot at its end; a chunk that hits the wall clock leaves
# none and is simply run again. The successor is submitted at the START of
# each job (afterany), by a retry loop because the per-user submit limit
# counts the user's other jobs. A job that finds leg 2 complete, or two
# failed chunks in a row, exits without a successor. NO_CHAIN=1 runs one chunk.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
NSTEPS="${NSTEPS:-20000}"
END1=130000; END2=780000
TUT="$CODE_DIR/tutorials/cht/channel"

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
for f in start_field.h5 cht149_flageul.h5 ynodes_f149.dat; do
    [ -f "$f" ] || { echo "no $f in $RUN_DIR"; exit 1; }
done
{
    echo "job    : ${SLURM_JOB_ID:-none}   $(date '+%F %T %Z')"
    echo "node   : $(hostname)"
    echo "code   : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
    echo "chunk  : $NSTEPS steps; legs Cdev2 -> $END1, Cstat -> $END2"
} | tee -a provenance.txt

latest_of() { ls -1 ${1}_[0-9]*.h5 2>/dev/null | sed 's/.*_\([0-9]*\)\.h5/\1 &/' | sort -n | tail -1 | cut -d' ' -f2; }
step_of() { echo "$1" | sed 's/.*_\([0-9]*\)\.h5/\1/'; }

# Which leg, and from which snapshot.
s2=$(latest_of Cstat); s1=$(latest_of Cdev2)
if [ -n "$s2" ] && [ "$(step_of "$s2")" -ge "$END2" ]; then
    echo "leg 2 at END2 $END2: nothing to do"; exit 0
fi
if [ "$(grep '^exit' status.txt 2>/dev/null | tail -2 | grep -c '^exit [^0]')" -ge 2 ]; then
    echo "the last two chunks failed -- stopping the chain"; exit 1
fi
if [ -n "$s2" ]; then
    leg=2; restart=$s2
elif [ -n "$s1" ] && [ "$(step_of "$s1")" -ge "$END1" ]; then
    leg=2; restart=$s1
    rm -f flageul_clean.h5 Cstat_vel.h5          # statistics FROM ZERO
else
    leg=1; restart=${s1:-start_field.h5}
fi
step=$(step_of "${restart/start_field.h5/x_65000.h5}")
if [ $leg = 1 ]; then
    prefix=Cdev2; end=$END1; sample=0; flush=0; statsfile=x
else
    prefix=Cstat; end=$END2; sample=25; flush=2000; statsfile=flageul_clean.h5
fi
n=$(( end - step )); [ "$n" -gt "$NSTEPS" ] && n=$NSTEPS

SCRIPT="$RUN_DIR/submit_cht_channel.sh"
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

# The chunk ini: run_flageul.sh's `ini` substitution plus the rerun settings
# (cflmax 1.0; niter 12 and sor 0.8 are explicit in the template, as a restart
# ini must be -- the start field carries niter 6).
sed -e "s|@CASE@|cht149_flageul.h5|" -e "s|@PREFIX@|$prefix|g" -e "s|@NSTEPS@|$n|" \
    -e "s|@WRITE@|0|" -e "s|@RESTART@|$restart|" -e "s|@SAMPLE@|$sample|g" \
    -e "s|@STATS@|$flush|g" -e "s|@STATSFILE@|$statsfile|g" \
    -e "s|^cflmax = 0.8|cflmax = 1.0|" "$TUT/cht149_flageul.ini" > "chunk_$prefix.ini"
if grep -v '^;' "chunk_$prefix.ini" | grep -q '@'; then echo "unsubstituted placeholder"; exit 1; fi

log="${prefix}_from_$(basename "$restart" .h5).log"
echo "=== leg $leg ($prefix): $n steps from $restart (step $step -> $((step + n)) of $end)"
echo "start $prefix $restart $(date '+%F %T')" >> status.txt
mpirun -n 4 --bind-to core --map-by numa "$CODE_DIR/build_gpu/moby_solve" "chunk_$prefix.ini" > "$log" 2>&1
rc=$?
grep -q "main loop ended" "$log" || rc=$(( rc == 0 ? 99 : rc ))   # the solver's own word, not mpirun's
echo "exit $rc $prefix $restart $(date '+%F %T')" >> status.txt
grep -E "seconds_per_step" "$log" | tail -1
tail -3 status.txt
