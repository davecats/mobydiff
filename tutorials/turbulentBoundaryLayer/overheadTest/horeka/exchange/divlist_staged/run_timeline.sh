#!/bin/bash
# Per-kernel timeline: WHAT are the 56-114 us of per-launch fixed cost?
#
#     run_timeline.sh <exe> <results_dir>
#
# results_horeka_exchange_2026-09-10.md measured a fixed cost per exchange kernel
# LAUNCH -- same-level copy 56-67 us, cross-level copy 85-87, pack ~92-114,
# unpack ~83-100 -- against jacobi_compute_phi's 24-26, and refused to name what
# it consists of. This probe names it.
#
# TWO PASSES over the same two configs, in one allocation:
#
#   A (untraced)  the exch_timing/proj_timing brackets, i.e. the HOST-side cost
#                 per launch that the report quotes. Nothing here is new; it is
#                 re-measured on these nodes so pass B compares like with like.
#   B (nsys)      per-kernel DEVICE duration, launch->start latency, the CUDA API
#                 call mix and the H2D memcpy count.
#
# WHY THE COMPARISON IS SOUND DESPITE TRACING OVERHEAD. nsys inflates host-side
# timings, and unevenly (results_horeka_2026-09-10.md section 5: "READ SHAPE, NOT
# MAGNITUDE"). It does NOT inflate the GPU execution time of a kernel. So the
# load-bearing number here is
#       (untraced host bracket)  minus  (traced device duration),
# where the first term comes from pass A and the second is tracing-independent:
# whatever that difference is, it is host-side work that is not the kernel. The
# API-call durations from pass B are used only to say WHICH calls, never how
# long they take.
set -uo pipefail

EXE="${1:?usage: run_timeline.sh <exe> <results_dir>}"
RES="${2:?usage: run_timeline.sh <exe> <results_dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="${CFG_DIR:-$HERE/../configs}"

NSYS="${NSYS:-/software/all/toolkit/nvidia_hpc_sdk/25.3/Linux_x86_64/25.3/profilers/Nsight_Systems/bin/nsys}"
NSTEPS="${NSTEPS:-30}"

# 4 ranks on ONE node. The fixed cost is a per-LAUNCH constant -- the fits put it
# within a few percent at 2, 4, 8 and 16 ranks -- so the smallest placement that
# still has peers (and therefore pack/unpack) is the cleanest place to dissect
# it, with no cross-node traffic in the timeline at all.
SPECS="${SPECS:-refined_yp82_rect_jacobi:4:1 rect_jacobi:4:1}"
# "A B" by default; "A" alone skips nsys, which is what a plain
# before/after step-time comparison at another rank count needs.
TL_PASSES="${TL_PASSES:-A B}"

command -v mpirun >/dev/null || { echo "ERROR: mpirun not on PATH" >&2; exit 1; }
[ -x "$EXE" ]  || { echo "ERROR: solver not executable: $EXE" >&2; exit 1; }
[ -x "$NSYS" ] || { echo "ERROR: nsys not found at $NSYS" >&2; exit 1; }
"$NSYS" --version
mkdir -p "$RES"

stage() {   # stage <cfg> <run_dir>
    local cfg="$1" run="$2"
    rm -rf "$run"; mkdir -p "$run"
    cp "$CFG/$cfg.ini" "$run/config.ini" || return 1
    if grep -qE '^[[:space:]]*profile[[:space:]]*=' "$run/config.ini"; then
        sed -i 's/^profile *=.*/profile = true/' "$run/config.ini"
    else
        printf '\n[output]\nprofile = true\n' >> "$run/config.ini"
    fi
    sed -i -e "s/^nsteps *=.*/nsteps = $NSTEPS/" \
           -e "s/^runtime_interval *=.*/runtime_interval = $NSTEPS/" "$run/config.ini"
}

echo "=== timeline start $(date '+%F %T')  nsteps=$NSTEPS"
for spec in $SPECS; do
    cfg="${spec%%:*}"; rest="${spec#*:}"; ranks="${rest%%:*}"; nodes="${rest##*:}"
    per_node=$(( ranks / nodes ))

    ############ pass A -- untraced brackets ############
    run="$RES/host_${cfg}_r${ranks}N${nodes}"
    if [ ! -f "$run/run.log" ]; then
        stage "$cfg" "$run" || { echo "MISSING CONFIG $cfg"; continue; }
        echo "=== host $cfg ${ranks}x${nodes} ($(date '+%F %T'))"
        ( cd "$run" && mpirun -n "$ranks" --map-by "ppr:${per_node}:node" --bind-to core \
            --display-map "$EXE" config.ini > run.log 2>&1 )
        grep -E "exch_timing|proj_timing" "$run/run.log" | sed 's/^/    /' | head -4
        rm -f "$run"/overhead_*.h5
    fi

    ############ pass B -- nsys ############
    case " $TL_PASSES " in *" B "*) ;; *) continue;; esac
    run="$RES/nsys_${cfg}_r${ranks}N${nodes}"
    [ -f "$run/run.log" ] && { echo "--- skip nsys $spec (done)"; continue; }
    stage "$cfg" "$run" || continue
    # Paths baked in at write time: OpenMPI forwards only OMPI_* to remote nodes.
    # exit 0 because the solver segfaults during process EXIT under nsys, after
    # the main loop and every timing line; without this mpirun kills the siblings
    # mid-report and only some ranks' traces survive.
    cat > "$run/nsys_wrap.sh" <<WRAP
#!/bin/bash
"$NSYS" profile --trace=cuda --sample=none --cpuctxsw=none \
    --force-overwrite=true -o "rep_\${OMPI_COMM_WORLD_RANK}" "$EXE" config.ini
exit 0
WRAP
    chmod +x "$run/nsys_wrap.sh"
    echo "=== nsys $cfg ${ranks}x${nodes} ($(date '+%F %T'))"
    ( cd "$run" && mpirun -n "$ranks" --map-by "ppr:${per_node}:node" --bind-to core \
        --display-map ./nsys_wrap.sh > run.log 2>&1 )
    rc=$?
    grep -q "main loop ended" "$run/run.log" || {
        echo "    FAILED (exit $rc, no completed main loop)"
        mv "$run/run.log" "$run/run.FAILED.log"; continue; }
    [ $rc -ne 0 ] && echo "    (mpirun exit $rc -- teardown only)"

    # Rank 0 only: its exch_timing is what every report in this campaign quotes,
    # and a full stats pass per rank costs minutes for no extra information.
    rep="$run/rep_0.nsys-rep"
    if [ -e "$rep" ]; then
        # ONE invocation: each `nsys stats` call re-exports the report to sqlite,
        # so asking for the reports separately pays that export once each.
        "$NSYS" stats --format csv -o "$run/rep_0" \
            --report cuda_gpu_kern_sum --report cuda_kern_exec_sum \
            --report cuda_api_sum --report cuda_gpu_mem_size_sum \
            --report cuda_gpu_mem_time_sum --report cuda_gpu_trace \
            "$rep" > "$run/stats.log" 2>&1 \
            || echo "    stats failed -- see $run/stats.log"
        echo "    csv: $(ls "$run"/rep_0_*.csv 2>/dev/null | wc -l) reports"
    else
        echo "    *** no rep_0.nsys-rep"
    fi
    rm -f "$run"/overhead_*.h5
done
echo "=== timeline done $(date '+%F %T')"
