#!/bin/bash
# The timing/profiling matrix, run INSIDE the SLURM job by submit.sh.
#
#   run_matrix.sh <exe> <results_dir>
#
# Every run is profiled ([output] profile = true): the profiler only reads
# clocks, so the fields are bit-identical to an unprofiled run and there is no
# reason not to. Each run gets its own directory holding the exact config it
# used, runtime.txt and run.log, so a number can always be traced back to the
# input that produced it.
#
# Runs are SKIPPED if their run.log already exists, which makes the whole matrix
# resumable: if the job hits the wall clock, `sbatch submit.sh` again and it
# picks up where it stopped.
set -uo pipefail

EXE="${1:?usage: run_matrix.sh <exe> <results_dir>}"
RES="${2:?usage: run_matrix.sh <exe> <results_dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="$HERE/configs"

# Rank counts to sweep. On HoreKa accelerated (4x A100 per node): 1/2/4 are
# intra-node (NVLink P2P), 8 and 16 cross the interconnect. The 4 -> 8 step is
# the one that matters -- every local measurement to date was intra-node, so the
# inter-node behaviour of the 2:1 exchange is entirely unmeasured.
RANKS_SMALL="${RANKS_SMALL:-1,2,4,8,16}"
# The production-size case holds ~39 GB of field state: 4 ranks is the smallest
# safe allocation on A100-40 (~9.7 GB/GPU).
RANKS_BIG="${RANKS_BIG:-4,8,16}"

# config:rank-list, rank lists COMMA-separated -- this table is iterated with
# word splitting, so a space inside an entry would split it into bogus configs.  base_jacobi is the unblocked reference that turns a
# s/step into a block tax; rect_jacobi is BOTH the single-level baseline and the
# exact twin of refined_big_rect_jacobi.
# CONFIGS restricts the matrix to a subset, e.g.
#   CONFIGS="refined_big_rect_jacobi rect_jacobi" RANKS_BIG=8 run_matrix.sh ...
# Useful for a targeted re-measurement without re-running the whole sweep.
CONFIGS="${CONFIGS:-}"

# Extra mpirun flags. Add --report-bindings on the first multi-node run to
# confirm each rank landed on the socket nearest its GPU; a silent fallback to
# one GPU per node reads as catastrophic scaling rather than as an error.
MPIRUN_EXTRA="${MPIRUN_EXTRA:-}"

MATRIX="
base_jacobi:$RANKS_SMALL
rect_jacobi:$RANKS_SMALL
refined_yp82_rect_jacobi:$RANKS_SMALL
refined_yp82_rect_redblack:$RANKS_SMALL
refined_big_rect_jacobi:$RANKS_BIG
"

command -v mpirun >/dev/null || {
    echo "ERROR: mpirun not on PATH -- load the toolchain module first" >&2
    echo "       (HoreKa: module load toolkit/nvidia-hpc-sdk/25.3)" >&2
    exit 1
}
[ -x "$EXE" ] || { echo "ERROR: solver not executable: $EXE" >&2; exit 1; }

mkdir -p "$RES"
echo "=== matrix start $(date '+%F %T')  exe=$EXE"

for entry in $MATRIX; do
    cfg="${entry%%:*}"
    ranks="${entry#*:}"
    ranks="${ranks//,/ }"
    [ -f "$CFG/$cfg.ini" ] || { echo "MISSING CONFIG: $cfg.ini" >&2; continue; }
    if [ -n "$CONFIGS" ] && ! printf '%s\n' $CONFIGS | grep -qx "$cfg"; then
        continue
    fi

    for n in $ranks; do
        run="$RES/${cfg}_n${n}"
        if [ -f "$run/run.log" ]; then
            echo "--- skip ${cfg} n=${n} (already done)"
            continue
        fi
        rm -rf "$run"; mkdir -p "$run"

        # Copy the config next to its output and switch profiling on, so the
        # run directory is self-describing.
        cp "$CFG/$cfg.ini" "$run/config.ini"
        if grep -qE '^[[:space:]]*profile[[:space:]]*=' "$run/config.ini"; then
            sed -i 's/^profile *=.*/profile = true/' "$run/config.ini"
        else
            printf '\n[output]\nprofile = true\n' >> "$run/config.ini"
        fi
        [ -n "${NSTEPS:-}" ] && sed -i \
            -e "s/^nsteps *=.*/nsteps = $NSTEPS/" \
            -e "s/^runtime_interval *=.*/runtime_interval = $(( NSTEPS/4 > 0 ? NSTEPS/4 : 1 ))/" \
            "$run/config.ini"

        echo "=== ${cfg}  n=${n}  ($(date '+%F %T'))"
        # --map-by numa --bind-to core keeps each rank on the socket nearest its
        # GPU; mobydiff binds rank -> device as local_rank mod n_devices.
        ( cd "$run" && mpirun -n "$n" --bind-to core --map-by numa $MPIRUN_EXTRA \
              "$EXE" config.ini > run.log 2>&1 )
        rc=$?
        if [ $rc -ne 0 ]; then
            echo "    FAILED (exit $rc) -- see $run/run.log"
            # Keep the log but drop nothing else; a failed run must NOT look
            # like a completed one to the resume check above.
            mv "$run/run.log" "$run/run.FAILED.log"
        else
            tail -2 "$run/runtime.txt" 2>/dev/null | sed 's/^/    /'
        fi
        # The unconditional end-of-run snapshot is written outside the timed
        # region but is GB-sized; the matrix only needs the timings.
        rm -f "$run"/overhead_*.h5
    done
done

echo "=== matrix done $(date '+%F %T')"
