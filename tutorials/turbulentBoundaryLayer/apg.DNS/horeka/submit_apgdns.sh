#!/bin/bash
#SBATCH --job-name=apgdns
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00
#SBATCH --partition=accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The apg.DNS run (../PLAN.md) on legacy HoreKa: 8 A100 ranks (2 nodes), in
# chunks of NSTEPS steps run back to back INSIDE the job while its wall clock
# allows another one, and self-chained across jobs until DEADLINE (the legacy
# job stop, 20 Nov 2026 noon).
#
#     export CODE_DIR=<clone at the run's commit> RUN_DIR=<run dir> HDF5_ROOT=$HOME/hdf5 \
#            [NSTEPS=12500] [DEADLINE="2026-11-20 11:45"] [NO_CHAIN=1] [MAX_CHUNKS=n] [CHUNK_SECS=3600]
#     cd $RUN_DIR && sbatch --export=ALL submit_apgdns.sh
#
# RUN_DIR holds apg_dns.ini, ptop.txt and start_field.h5 (assets/tile_restart.py).
# The case file apg.case.h5 is prepared on the first chunk (8 ranks).
#
# A chunk restarts from the snapshot of the last chunk that ENDED CLEANLY (an
# "ok <step>" line in status.txt; a chunk killed while writing leaves a file
# that is never read), runs NSTEPS steps, writes one snapshot and the
# cumulative statistics at its end. Before it starts, apg_stats.h5 is reset to
# the copy taken at that same step (stats_apg_<step>.h5), so a killed chunk
# never leaves samples in the accumulators that its restart field does not
# have. Windows are differences of these copies. The two newest snapshots are
# kept (10.7 GB each).
#
# Steering between chunks: the chunk ini is regenerated from apg_dns.ini
# every chunk, so an edited apg_dns.ini or a replaced ptop.txt (write it
# under another name, then mv) takes effect at the next chunk; a file STOP in
# RUN_DIR ends the job after the running chunk, without a successor.
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
NSTEPS="${NSTEPS:-12500}"; DEADLINE="${DEADLINE:-2026-11-20 11:45}"
MAX_CHUNKS="${MAX_CHUNKS:-1000}"; CHUNK_SECS="${CHUNK_SECS:-3600}"

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
for f in apg_dns.ini ptop.txt start_field.h5; do [ -f "$f" ] || { echo "no $f in $RUN_DIR"; exit 1; }; done
[ -f STOP ] && { echo "STOP present: nothing to do"; exit 0; }
{
    echo "job    : ${SLURM_JOB_ID:-none}   $(date '+%F %T %Z')   limit $(squeue -h -j "${SLURM_JOB_ID:-0}" -o %l 2>/dev/null)"
    echo "nodes  : ${SLURM_JOB_NODELIST:-$(hostname)}"
    echo "code   : $(git -C "$CODE_DIR" rev-parse HEAD)"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1) x ${SLURM_NTASKS:-?}"
    echo "chunk  : $NSTEPS steps, deadline $DEADLINE"
} | tee -a provenance.txt

if [ "$(grep '^exit' status.txt 2>/dev/null | tail -2 | grep -c '^exit [^0]')" -ge 2 ]; then
    echo "the last two chunks failed -- stopping the chain"; exit 1
fi

secs() { date -d "$1" +%s; }
to_secs() {  # a Slurm duration ([d-]hh:mm:ss, mm:ss) in seconds
    local s=$1 d=0 a b c
    [[ $s == *-* ]] && { d=${s%%-*}; s=${s#*-}; }
    IFS=: read -r a b c <<< "$s"
    if [ -n "$c" ]; then echo $(( d * 86400 + 3600 * 10#$a + 60 * 10#$b + 10#$c ))
    elif [ -n "$b" ]; then echo $(( d * 86400 + 60 * 10#$a + 10#$b ))
    else echo $(( d * 86400 + 10#$a )); fi
}
# The successor waits for this job (afterany) and gets the time left between
# this job's END (at the latest) and DEADLINE, capped at this job's limit; no
# successor when that is under 90 minutes.
SCRIPT="$RUN_DIR/submit_apgdns.sh"
submit_successor() {
    local limit_s left_s end_s t k=0 out
    limit_s=$(to_secs "$(squeue -h -j "$SLURM_JOB_ID" -o %l)")
    end_s=$(( $(date +%s) + limit_s ))
    left_s=$(( $(secs "$DEADLINE") - end_s ))
    [ "$left_s" -gt "$limit_s" ] && left_s=$limit_s
    if [ "$left_s" -lt 5400 ]; then echo "no successor: $((left_s / 60)) min left before $DEADLINE" | tee -a provenance.txt; return 0; fi
    t=$(printf '%d:%02d:00' $((left_s / 3600)) $(((left_s % 3600) / 60)))
    while [ $k -lt 25 ]; do
        if out=$(sbatch --dependency=afterany:$SLURM_JOB_ID --job-name="$SLURM_JOB_NAME" \
                     --partition="$SLURM_JOB_PARTITION" --nodes="$SLURM_NNODES" \
                     --time="$t" --export=ALL "$SCRIPT" 2>&1); then
            echo "successor of $SLURM_JOB_ID: $out, limit $t (attempt $((k + 1)))" | tee -a provenance.txt; return 0
        fi
        k=$((k + 1)); sleep 120
    done
    echo "CHAIN BROKEN after $k attempts ($out) -- resubmit by hand: cd $RUN_DIR && sbatch --export=ALL $SCRIPT" | tee -a provenance.txt
}
[ -n "${NO_CHAIN:-}" ] || submit_successor &

MPIRUN="mpirun -n ${SLURM_NTASKS:-8} --bind-to core --map-by numa"
if [ ! -f apg.case.h5 ]; then
    $MPIRUN "$CODE_DIR/build_gpu/moby_prepare" apg_dns.ini > prep.log 2>&1 \
        || { echo "PREPARE FAILED"; tail -5 prep.log; exit 1; }
fi

job_left() { to_secs "$(squeue -h -j "$SLURM_JOB_ID" -o %L)"; }   # seconds left in this job
need=$CHUNK_SECS   # a chunk's wall time with margin; then 1.25 x the last chunk + 5 min
chunks=0
while [ $chunks -lt "$MAX_CHUNKS" ]; do
    [ -f STOP ] && { echo "STOP found"; break; }
    left=$(job_left)
    [ "$left" -gt "$need" ] || { echo "$left s left in the job, a chunk needs $need: ending"; break; }
    last=$(grep '^ok ' status.txt 2>/dev/null | tail -1 | awk '{print $2}')
    if [ -n "$last" ]; then
        restart="apg_${last}.h5"
        cp "stats_apg_${last}.h5" apg_stats.h5 || { echo "no stats copy for step $last"; exit 1; }
    else
        restart=start_field.h5
        rm -f apg_stats.h5
    fi
    sed -e "s/^nsteps = .*/nsteps = $NSTEPS/" \
        -e "s/^stats_write_interval = .*/stats_write_interval = $NSTEPS/" \
        -e "/^\[restart\]/,/^$/ s|^file = .*|file = $restart|" apg_dns.ini > chunk.ini
    log="apg_from_$(basename "$restart" .h5).log"
    echo "start $restart $(date '+%F %T') job ${SLURM_JOB_ID:-none} ptop $(md5sum < ptop.txt | cut -c1-8)" >> status.txt
    t0=$(date +%s)
    $MPIRUN "$CODE_DIR/build_gpu/moby_solve" chunk.ini > "$log" 2>&1
    rc=$?
    grep -q "main loop ended" "$log" || rc=$(( rc == 0 ? 99 : rc ))
    echo "exit $rc $restart $(date '+%F %T')" >> status.txt
    [ $rc -eq 0 ] || { echo "chunk failed (rc $rc), see $log"; tail -5 "$log"; exit $rc; }
    step=$(grep -o 'field filename: apg_[0-9]*\.h5' "$log" | tail -1 | sed 's/.*apg_\([0-9]*\)\.h5/\1/')
    [ -n "$step" ] && [ -f "apg_${step}.h5" ] || { echo "no final snapshot in $log"; exit 1; }
    cp apg_stats.h5 "stats_apg_${step}.h5"
    echo "ok $step $(date '+%F %T')" >> status.txt
    need=$(( ($(date +%s) - t0) * 5 / 4 + 300 ))
    chunks=$((chunks + 1))
    # keep the two newest snapshots of finished chunks
    for f in $(ls -1 apg_[0-9]*.h5 | sed 's/apg_\([0-9]*\)\.h5/\1 &/' | sort -n | head -n -2 | cut -d' ' -f2); do rm -f "$f"; done
    grep seconds_per_step "$log" | tail -1
    tail -1 apg_runtime.txt
done
tail -4 status.txt 2>/dev/null || true
