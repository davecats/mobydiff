#!/bin/bash
# Timeline probe: what is a rank ACTUALLY waiting for inside MPI_Waitall?
#
#     run_nsys.sh <exe> <results_dir>
#
# Same spec grammar as run_placement.sh (config:ranks:nodes, same --map-by, same
# placement check), but each rank runs under Nsight Systems with MPI and CUDA
# tracing, so every Isend/Irecv/Waitall lands on the same timeline as the GPU
# kernels. That is the one thing the aggregate exch_timing buckets cannot show:
# a Waitall that overlaps CUDA activity is waiting for the DEVICE, not the wire,
# and no partitioning or transport change would touch it.
#
# No solver code is involved -- this profiles the same binary the timing passes
# run, so it needs no bit-exactness gate.
set -uo pipefail

EXE="${1:?usage: run_nsys.sh <exe> <results_dir>}"
RES="${2:?usage: run_nsys.sh <exe> <results_dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="$HERE/../configs"

# nsys is in the nvhpc tree but NOT on the PATH the modulefile sets.
NSYS="${NSYS:-/software/all/toolkit/nvidia_hpc_sdk/25.3/Linux_x86_64/25.3/profilers/Nsight_Systems/bin/nsys}"

# 40 steps: enough that the last 20 are past every start-up transient, few
# enough that a trace stays small. Step TIMES from a traced run are not quoted
# anywhere -- the timeline structure is the measurement.
NSTEPS="${NSTEPS:-40}"
MPIRUN_EXTRA="${MPIRUN_EXTRA:-}"

# NIC=1 adds --nic-metrics=true (HCA counters). OFF by default and deliberately
# a SECOND pass, not part of the first: it is a system-scope collector that can
# fail on permissions and take the profile session with it, and it is only worth
# having once the first pass has said the wait is genuinely MPI rather than
# unfinished device work. The nodes carry bonded mlx5 HCAs (ibstat -l), so this
# is what would show 4 ranks/node saturating a link that 2 ranks/node does not.
NIC="${NIC:-0}"
NIC_FLAG=""
[ "$NIC" = 1 ] && NIC_FLAG="--nic-metrics=true"

# The four cells the timeline has to explain, all at 8 ranks except the last:
#   rect 8x2  the reproducible 752 us anomaly (4 ranks/node, 1 crossing)
#   base 8x2  the 101 us control at the SAME placement
#   rect 8x4  the same blocked config, FAST (182 us) spread 2 ranks/node
#   rect 4x1  the mild regime, no node crossing at all
SPECS="${SPECS:-
rect_jacobi:8:2
base_jacobi:8:2
rect_jacobi:8:4
rect_jacobi:4:1
}"

command -v mpirun >/dev/null || { echo "ERROR: mpirun not on PATH" >&2; exit 1; }
[ -x "$EXE" ]  || { echo "ERROR: solver not executable: $EXE" >&2; exit 1; }
[ -x "$NSYS" ] || { echo "ERROR: nsys not found at $NSYS" >&2; exit 1; }
"$NSYS" --version
mkdir -p "$RES"
echo "=== nsys matrix start $(date '+%F %T')  nsteps=$NSTEPS"

for spec in $SPECS; do
    cfg="${spec%%:*}"; rest="${spec#*:}"
    ranks="${rest%%:*}"; nodes="${rest##*:}"
    [ -f "$CFG/$cfg.ini" ] || { echo "MISSING CONFIG: $cfg.ini" >&2; continue; }
    per_node=$(( ranks / nodes ))
    if [ $(( per_node * nodes )) -ne "$ranks" ] || [ "$per_node" -lt 1 ] || [ "$per_node" -gt 4 ]; then
        echo "SKIP $spec: $ranks ranks do not lay out over $nodes nodes" >&2; continue
    fi

    run="$RES/${cfg}_r${ranks}_N${nodes}"
    [ -f "$run/run.log" ] && { echo "--- skip $cfg ${ranks}x${nodes} (done)"; continue; }
    rm -rf "$run"; mkdir -p "$run"

    cp "$CFG/$cfg.ini" "$run/config.ini"
    if grep -qE '^[[:space:]]*profile[[:space:]]*=' "$run/config.ini"; then
        sed -i 's/^profile *=.*/profile = true/' "$run/config.ini"
    else
        printf '\n[output]\nprofile = true\n' >> "$run/config.ini"
    fi
    sed -i -e "s/^nsteps *=.*/nsteps = $NSTEPS/" \
           -e "s/^runtime_interval *=.*/runtime_interval = $NSTEPS/" "$run/config.ini"

    echo "=== $cfg  ${ranks} ranks on ${nodes} node(s) = ${per_node}/node  ($(date '+%F %T'))"
    # One report per rank. --sample=none/--cpuctxsw=none keep the trace to the
    # two event streams the question needs; sampling would add noise and bulk
    # without answering anything here.
    ( cd "$run" && mpirun -n "$ranks" --map-by "ppr:${per_node}:node" --bind-to core \
        --display-map $MPIRUN_EXTRA \
        "$NSYS" profile --trace=mpi,cuda --mpi-impl=openmpi \
                --sample=none --cpuctxsw=none --force-overwrite=true $NIC_FLAG \
                -o "rep_%q{OMPI_COMM_WORLD_RANK}" \
                "$EXE" config.ini > run.log 2>&1 )
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "    FAILED (exit $rc)"; mv "$run/run.log" "$run/run.FAILED.log"; continue
    fi
    awk '/Data for node:/ {print $4}' "$run/run.log" | sort -u > "$run/hosts.txt"
    got=$(wc -l < "$run/hosts.txt")
    echo "    hosts ($got): $(tr '\n' ' ' < "$run/hosts.txt")"
    [ "$got" -ne "$nodes" ] && {
        echo "    *** PLACEMENT MISMATCH: asked $nodes, got $got -- VOID"
        echo "PLACEMENT MISMATCH: requested $nodes nodes, used $got" > "$run/VOID"; }

    # Reduce to CSV inside the job: the .nsys-rep files stay for a GUI, but the
    # analysis is headless and must not depend on one.
    for rep in "$run"/rep_*.nsys-rep; do
        [ -e "$rep" ] || continue
        b="${rep%.nsys-rep}"
        # BOTH reports in one call: each invocation exports the report to sqlite
        # first, so two calls would pay that export twice for the same trace.
        # Writes ${b}_mpi_event_trace.csv and ${b}_cuda_gpu_trace.csv.
        "$NSYS" stats --report mpi_event_trace,cuda_gpu_trace --format csv \
                -o "$b" "$rep" > "${b}_stats.log" 2>&1 \
            || echo "    stats failed for $(basename "$rep") -- see ${b}_stats.log"
    done
    echo "    csv: $(ls "$run"/*mpi_event_trace*.csv 2>/dev/null | wc -l) mpi, $(ls "$run"/*cuda_gpu_trace*.csv 2>/dev/null | wc -l) gpu"
    rm -f "$run"/overhead_*.h5
done
echo "=== nsys matrix done $(date '+%F %T')"
