#!/bin/bash
# Placement matrix: the same job at the same RANK COUNT, laid out on a
# different number of nodes. Rank count, decomposition, peer count and message
# volume are all held FIXED; the only thing that changes is how many of the
# links cross a node boundary.
#
#   run_placement.sh <exe> <results_dir>
#
# Each spec is  config:ranks:nodes  and is launched with
#   mpirun -n <ranks> --map-by ppr:<ranks/nodes>:node
# so OpenMPI puts exactly ranks/nodes tasks on each of the first <nodes> nodes
# of the allocation. --display-map is captured in every log: the placement is
# EVIDENCE, not an assumption, and must be read back before the timings.
set -uo pipefail

EXE="${1:?usage: run_placement.sh <exe> <results_dir>}"
RES="${2:?usage: run_placement.sh <exe> <results_dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="$HERE/../configs"

# 200 steps is ample: the per-round exchange rate converges in a few tens of
# steps, and every number this matrix produces is a RATIO within it.
NSTEPS="${NSTEPS:-200}"
MPIRUN_EXTRA="${MPIRUN_EXTRA:-}"

# BARRIER=1 turns on [output] exchange_barrier for every run in this pass,
# splitting the exchange's arrival SKEW (bucket skew_barrier) from its
# TRANSFER (what is left in mpi_wait). It serialises the exchange, so a
# barrier pass is for ATTRIBUTION ONLY -- its step times mean nothing. Write
# it to its own results directory so the two passes never mix.
BARRIER="${BARRIER:-0}"

# config:ranks:nodes -- comma-separated, no spaces (word splitting).
#
# ROUND 2 (2026-09-09) added nb16_jacobi:8:2, nb16_jacobi:8:4 and
# blk8_jacobi:4:2. Round 1 refuted all three original hypotheses; nb16 at 8
# ranks / 4 per node is the missing direct test, since it carries 7.6x
# rect's copy load per rank on the SAME 2-peer chain at the SAME placement as
# the reproducible 752 us anomaly. blk8:4:2 completes the (copy, placement)
# 2x2. Re-running rect_jacobi:8:2 in the same allocation (delete its result
# directory to force it) removes the cross-job caveat.
#
# TIER 5 (blk8_jacobi) is the second discriminator: node crossing WITHOUT the
# device-local copy load, one block per rank. See its config header.
# TIER 4 (nb16_jacobi) varies the exchange VOLUME several-fold at a fixed
# placement: if the wait is flat against a large change in send points it is a
# latency/serialization limit, if it tracks the bytes it is bandwidth.
#
# TIER 1 is the decisive pair. At 2 ranks each rank has exactly ONE peer, the
# two are symmetric (36800 pts each), and there is no chain to propagate a
# stall. 2x1 vs 2x2 therefore isolates a SINGLE cross-node link with nothing
# else changing anywhere in the problem.
SPECS="${SPECS:-
rect_jacobi:2:1
rect_jacobi:2:2
base_jacobi:2:1
base_jacobi:2:2
rect_jacobi:4:1
rect_jacobi:4:2
rect_jacobi:4:4
base_jacobi:4:1
base_jacobi:4:2
base_jacobi:4:4
rect_jacobi:8:2
rect_jacobi:8:4
base_jacobi:8:2
base_jacobi:8:4
refined_yp82_rect_jacobi:2:1
refined_yp82_rect_jacobi:2:2
refined_yp82_rect_jacobi:4:1
refined_yp82_rect_jacobi:4:2
nb16_jacobi:4:1
nb16_jacobi:4:2
blk8_jacobi:4:1
blk8_jacobi:8:2
blk8_jacobi:8:4
nb16_jacobi:8:2
nb16_jacobi:8:4
blk8_jacobi:4:2
}"

command -v mpirun >/dev/null || { echo "ERROR: mpirun not on PATH" >&2; exit 1; }
[ -x "$EXE" ] || { echo "ERROR: solver not executable: $EXE" >&2; exit 1; }
mkdir -p "$RES"
echo "=== placement matrix start $(date '+%F %T')  nsteps=$NSTEPS"

for spec in $SPECS; do
    cfg="${spec%%:*}"; rest="${spec#*:}"
    ranks="${rest%%:*}"; nodes="${rest##*:}"
    [ -f "$CFG/$cfg.ini" ] || { echo "MISSING CONFIG: $cfg.ini" >&2; continue; }
    per_node=$(( ranks / nodes ))
    if [ $(( per_node * nodes )) -ne "$ranks" ] || [ "$per_node" -lt 1 ]; then
        echo "SKIP $spec: $ranks ranks do not divide over $nodes nodes" >&2; continue
    fi
    if [ "$per_node" -gt 4 ]; then
        echo "SKIP $spec: $per_node ranks/node exceeds the 4 GPUs on a node" >&2; continue
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
    if [ "$BARRIER" = 1 ]; then
        # Anchor to the profile line, which the block above has just guaranteed
        # exists inside [output]. A blind append to the end of the file would
        # land in whatever section happens to be last.
        if grep -qE '^[[:space:]]*exchange_barrier[[:space:]]*=' "$run/config.ini"; then
            sed -i 's/^exchange_barrier *=.*/exchange_barrier = true/' "$run/config.ini"
        else
            sed -i '0,/^profile *=.*/s//&\nexchange_barrier = true/' "$run/config.ini"
        fi
        grep -qE '^exchange_barrier *= *true' "$run/config.ini" || {
            echo "    *** could not set exchange_barrier -- skipping"; continue; }
    fi
    sed -i -e "s/^nsteps *=.*/nsteps = $NSTEPS/" \
           -e "s/^runtime_interval *=.*/runtime_interval = $(( NSTEPS/4 > 0 ? NSTEPS/4 : 1 ))/" \
           "$run/config.ini"

    echo "=== $cfg  ${ranks} ranks on ${nodes} node(s) = ${per_node}/node  ($(date '+%F %T'))"
    ( cd "$run" && mpirun -n "$ranks" --map-by "ppr:${per_node}:node" --bind-to core \
          --display-map $MPIRUN_EXTRA "$EXE" config.ini > run.log 2>&1 )
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "    FAILED (exit $rc)"; mv "$run/run.log" "$run/run.FAILED.log"
    else
        # Placement is EVIDENCE. --display-map prints one "Data for node: <host>"
        # line per node actually used; record it and CHECK it, so a run that
        # OpenMPI laid out differently from the request cannot be read as if it
        # had obeyed.
        awk '/Data for node:/ {print $4}' "$run/run.log" | sort -u > "$run/hosts.txt"
        got=$(wc -l < "$run/hosts.txt")
        echo "    hosts ($got): $(tr '\n' ' ' < "$run/hosts.txt")"
        if [ "$got" -ne "$nodes" ]; then
            echo "    *** PLACEMENT MISMATCH: asked for $nodes node(s), got $got ***"
            echo "    *** this run is VOID -- --map-by ppr:$per_node:node was not honoured"
            echo "PLACEMENT MISMATCH: requested $nodes nodes, used $got" > "$run/VOID"
        fi
        grep -m1 "exchange sizes" "$run/run.log" | sed 's/^/    /'
    fi
    rm -f "$run"/overhead_*.h5
done
echo "=== placement matrix done $(date '+%F %T')"
