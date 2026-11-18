#!/bin/bash
#SBATCH --job-name=sd7003
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:4
#SBATCH --time=01:00:00
#SBATCH --partition=dev-gpu-h100
#SBATCH --account=hk-project-p0024333
#SBATCH --no-requeue
#
# Step 3b, the SD7003 chain (tutorials/sailplane/junction/sd7003) on HoreKa 2
# in one-hour jobs of the dev partition (the regular gpu-h100 queue estimated
# a week on 2026-11-16, the dev node an hour). Every job runs ONE chunk of
# the first unfinished stage -- seed (to step 7500), dev (30 000), then the
# branches c12, c24, rb (55 000 each, from sd_30000.h5) -- from that stage's
# newest snapshot, reduces a branch's snapshots to span-averaged moments
# (sd_reduce.py) and deletes all but the newest, and submits its successor
# at the start (afterany). A stage is finished when its end snapshot exists.
#
#     export CODE_DIR=$PROJECT/moby-junction RUN_DIR=$(ws_find junction)/sd7003_run
#     cd $RUN_DIR && sbatch --export=ALL $CODE_DIR/tutorials/sailplane/junction/horeka/chain_sd7003_hk2.sh
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
SECONDS_PER_STEP_C12="${SECONDS_PER_STEP_C12:?measured step time of the 4-rank Chebyshev 12 run}"
BUDGET="${BUDGET:-2700}"                      # seconds of stepping per job (of 3600)
module purge
module load NVHPC/25.9-CUDA-12.9.1 HDF5/1.14.6-NVHPC-25.9-CUDA-12.9.1
export LD_LIBRARY_PATH="$EBROOTHDF5/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1
BIN="$CODE_DIR/build_gpu/moby_solve"
SELF="$CODE_DIR/tutorials/sailplane/junction/horeka/chain_sd7003_hk2.sh"
cd "$RUN_DIR"
newest() { ls -1 ${1}_[0-9]*.h5 2>/dev/null | sed 's/.*_\([0-9]*\)\.h5/\1 &/' | sort -n | tail -1 | cut -d' ' -f2; }
step_of() { echo "$1" | sed 's/.*_\([0-9]*\)\.h5/\1/'; }

# the first unfinished stage: name, prefix, end step, cost factor per step (x 10)
stage=""
for s in "seed sd 7500 10" "dev sd 30000 10" "c12 sd_c12 55000 10" "c24 sd_c24 55000 ${FACTOR_C24:-19}" "rb sd_rb 55000 ${FACTOR_RB:-19}" "c48 sd_c48 55000 ${FACTOR_C48:-35}" "rb24 sd_rb24 55000 ${FACTOR_RB24:-35}"; do
    set -- $s
    [ -f "$2_$3.h5" ] && continue
    stage=$1; prefix=$2; end=$3; factor=$4; break
done
[ -n "$stage" ] || { echo "all stages done"; exit 0; }
if [ "$(grep -c '^exit' status.txt 2>/dev/null)" -ge 2 ] && [ "$(grep '^exit' status.txt | tail -2 | grep -c '^exit [^0]')" -ge 2 ]; then
    echo "the last two chunks failed -- stopping the chain"; exit 1
fi
if [ -z "${NO_CHAIN:-}" ]; then
    ( k=0; until out=$(sbatch --dependency=afterany:$SLURM_JOB_ID --export=ALL -o sd7003-%j.out "$SELF" 2>&1); do
          k=$((k + 1)); [ $k -ge 25 ] && { echo "CHAIN BROKEN: $out" >> status.txt; exit; }; sleep 120; done
      echo "successor of $SLURM_JOB_ID: $out" >> status.txt ) &
fi

latest=$(newest "$prefix")
case $stage in seed) : ;; dev) [ "$(step_of "${latest:-x_0.h5}")" -lt 7500 ] && latest="" ; [ -z "$latest" ] && latest=sd_7500.h5 ;;
               *) [ -z "$latest" ] && latest=sd_30000.h5 ;; esac
from=0; [ -n "$latest" ] && from=$(step_of "$latest")
chunk=$(awk -v b="$BUDGET" -v s="$SECONDS_PER_STEP_C12" -v f="$factor" 'BEGIN { printf "%d", int(b / (s * f / 10) / 250) * 250 }')
[ $((from + chunk)) -gt "$end" ] && chunk=$((end - from))
interval=$chunk; case $stage in c12|c24|rb|c48|rb24) interval=250 ;; esac
ini=chunk_$stage.ini
sed -e "s/^nsteps = .*/nsteps = $chunk/" -e "s/^field_interval = .*/field_interval = $interval/" \
    -e '/^\[restart\]/,$d' "sd_$stage.ini" > "$ini"
[ -n "$latest" ] && printf '[restart]\nfile = %s\n' "$latest" >> "$ini"
log="sd_${stage}_from_$from.log"
echo "start $stage from $from, $chunk steps, job ${SLURM_JOB_ID:-none} $(date '+%F %T') $(hostname) code $(git -C "$CODE_DIR" rev-parse --short HEAD)" >> status.txt
mpirun -n "${SLURM_NTASKS:-4}" --bind-to none "$BIN" "$ini" > "$log" 2>&1; rc=$?
grep -q "main loop ended" "$log" || rc=$(( rc == 0 ? 99 : rc ))
echo "exit $rc $stage to $(step_of "$(newest "$prefix")") $(grep -o 'seconds_per_step .*' "$log" | head -1) $(date '+%F %T')" >> status.txt
[ $rc -eq 0 ] || exit 1

case $stage in c12|c24|rb|c48|rb24)
    ( module purge; module load h5py/3.12.1-foss-2024a
      files=$(ls -1 ${prefix}_[0-9]*.h5 | sed 's/.*_\([0-9]*\)\.h5/\1 &/' | sort -n | cut -d' ' -f2)
      todo=""; for f in $files; do [ -s "${f%.h5}.red.npz" ] || todo="$todo $f"; done
      python3 sd_reduce.py $todo >> reduce_$stage.log 2>&1 || { echo "reduce FAILED $stage" >> status.txt; exit 1; }
      for f in $(echo "$files" | head -n -1); do [ -s "${f%.h5}.red.npz" ] && rm -f "$f"; done ) ;;
esac
tail -4 status.txt
