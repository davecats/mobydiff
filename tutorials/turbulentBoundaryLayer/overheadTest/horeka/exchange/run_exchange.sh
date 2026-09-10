#!/bin/bash
# The three passes of docs/next_session_2to1_performance.md, run inside one
# SLURM allocation.
#
#   run_exchange.sh <exe> <results_dir>
#
# PASS G (gate)   -- the new binary against the pre-change one, same nodes, same
#                    step count: the op-split counters are init-time integers and
#                    prints, so the fields must be max_abs 0.
# PASS 1 (op split) -- `exchange by op` at 1/4/16 ranks for the refined case and
#                    its single-level twin. The quantity is a property of the leaf
#                    table and the rank split, not of the timing, but the runs are
#                    profiled anyway so the split and the exchange buckets come
#                    from the SAME run.
# PASS 2 (A0)     -- the overlap probe repeated at 8 and 16 ranks. A compute-bound
#                    kernel between the Isend/Irecv posts and the Waitall: if the
#                    transfer progresses during it, mpi_wait collapses.
#                    Baseline and probe differ ONLY by MOBY_A0_PROBE.
set -uo pipefail

EXE="${1:?usage: run_exchange.sh <exe> <results_dir>}"
RES="${2:?usage: run_exchange.sh <exe> <results_dir>}"
REF="${REF:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="$HERE/../configs"
H5MAXDIFF="${H5MAXDIFF:-$HERE/../../../../../tools/h5maxdiff}"

NSTEPS_P1="${NSTEPS_P1:-200}"      # as the re-measured campaign
NSTEPS_P2="${NSTEPS_P2:-100}"      # the probe multiplies the step time
NSTEPS_GATE="${NSTEPS_GATE:-20}"
A0_ITER="${A0_ITER:-100000}"       # ~5 ms/round on an A100; reported, not assumed
# Which passes to run: "G 1 2" by default. The gate and the A0 smoke fit a
# single-node dev allocation, so they are usually run there first and the
# 4-node job then does 1 and 2 into the same results directory.
PASSES="${PASSES:-G 1 2}"
has_pass() { case " $PASSES " in *" $1 "*) return 0;; esac; return 1; }

mkdir -p "$RES"

# Stage a config: its own directory, profiling on, nsteps overridden.
stage() {   # stage <cfg> <run_dir> <nsteps>
    local cfg="$1" run="$2" n="$3"
    rm -rf "$run"; mkdir -p "$run"
    cp "$CFG/$cfg.ini" "$run/config.ini" || return 1
    if grep -qE '^[[:space:]]*profile[[:space:]]*=' "$run/config.ini"; then
        sed -i 's/^profile *=.*/profile = true/' "$run/config.ini"
    else
        printf '\n[output]\nprofile = true\n' >> "$run/config.ini"
    fi
    sed -i -e "s/^nsteps *=.*/nsteps = $n/" \
           -e "s/^runtime_interval *=.*/runtime_interval = $(( n/4 > 0 ? n/4 : 1 ))/" \
           "$run/config.ini"
}

# Launch. Placement is EVIDENCE: --display-map is captured and the node count
# checked back, exactly as run_placement.sh does.
launch() {  # launch <exe> <ranks> <nodes|0 for numa> <run_dir>
    local exe="$1" ranks="$2" nodes="$3" run="$4" rc
    local mapflag="--map-by numa"
    if [ "$nodes" != 0 ]; then
        mapflag="--map-by ppr:$(( ranks / nodes )):node"
    fi
    ( cd "$run" && mpirun -n "$ranks" $mapflag --bind-to core --display-map \
          "$exe" config.ini > run.log 2>&1 )
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "    FAILED (exit $rc)"; mv "$run/run.log" "$run/run.FAILED.log"; return 1
    fi
    awk '/Data for node:/ {print $4}' "$run/run.log" | sort -u > "$run/hosts.txt"
    echo "    hosts: $(tr '\n' ' ' < "$run/hosts.txt")"
    return 0
}

############################  PASS G -- the gate  ############################
if has_pass G && [ -n "$REF" ] && [ -x "$REF" ]; then
    echo "############ PASS G: bit-exactness against $REF ############"
    for cfg in rect_jacobi refined_yp82_rect_jacobi; do
        for side in ref new; do
            run="$RES/gate_${cfg}_${side}"
            [ -f "$run/run.log" ] && continue
            stage "$cfg" "$run" "$NSTEPS_GATE" || { echo "MISSING CONFIG $cfg"; continue; }
            echo "=== gate $cfg $side ($(date '+%F %T'))"
            if [ "$side" = ref ]; then launch "$REF" 4 0 "$run"; else launch "$EXE" 4 0 "$run"; fi
        done
        a=$(ls "$RES/gate_${cfg}_ref"/overhead_*.h5 2>/dev/null | tail -1)
        b=$(ls "$RES/gate_${cfg}_new"/overhead_*.h5 2>/dev/null | tail -1)
        if [ -n "$a" ] && [ -n "$b" ]; then
            echo "--- h5maxdiff $cfg"
            "$H5MAXDIFF" "$a" "$b" | tee "$RES/gate_${cfg}.txt"
        else
            echo "--- gate $cfg: NO SNAPSHOT PAIR (a='$a' b='$b')" | tee "$RES/gate_${cfg}.txt"
        fi
        rm -f "$RES/gate_${cfg}_ref"/overhead_*.h5 "$RES/gate_${cfg}_new"/overhead_*.h5
    done
fi

##########################  PASS 1 -- the op split  ##########################
has_pass 1 && echo "############ PASS 1: exchange volume by op ############"
# config:rank-list, ranks joined by "+". NOT by a comma: sbatch --export splits
# its value on commas, so a comma here silently truncates OP_SPECS to its first
# entry -- which is exactly what happened in job 5139461.
# The two twins the handout asks for (rect / refined_yp82) plus
# a VOLUME SWEEP at fixed rank count: base_jacobi (unblocked, zero local copy,
# 7-11 peers), nb16_jacobi (nb = 16, ~4x rect's copy load) and refined_big. Time
# per exchange CALL against points per call over that spread separates a fixed
# per-round cost from a per-point one -- which no single pair of configs can.
# refined_big holds ~39 GB of field state, so 1 rank does not fit on A100-40.
for entry in ${OP_SPECS:-base_jacobi:4+16 rect_jacobi:1+4+16 nb16_jacobi:4+16 \
                        refined_yp82_rect_jacobi:1+4+16 refined_big_rect_jacobi:4+16}; do
    has_pass 1 || continue
    cfg="${entry%%:*}"; ranks="${entry#*:}"
    for n in ${ranks//+/ }; do
        run="$RES/op_${cfg}_n${n}"
        [ -f "$run/run.log" ] && { echo "--- skip $cfg n=$n (done)"; continue; }
        stage "$cfg" "$run" "$NSTEPS_P1" || { echo "MISSING CONFIG $cfg"; continue; }
        echo "=== op $cfg n=$n ($(date '+%F %T'))"
        launch "$EXE" "$n" 0 "$run"
        grep -E "exchange sizes|exchange by op" "$run/run.log" | sed 's/^/    /'
        rm -f "$run"/overhead_*.h5
    done
done

###########################  PASS 2 -- the A0 probe  #########################
has_pass 2 && echo "############ PASS 2: A0 overlap probe ############"
for spec in ${A0_SPECS:-refined_yp82_rect_jacobi:8:2 refined_yp82_rect_jacobi:16:4 rect_jacobi:16:4}; do
    has_pass 2 || continue
    cfg="${spec%%:*}"; rest="${spec#*:}"; ranks="${rest%%:*}"; nodes="${rest##*:}"
    for side in base probe; do
        run="$RES/a0_${cfg}_r${ranks}N${nodes}_${side}"
        [ -f "$run/run.log" ] && { echo "--- skip $spec $side (done)"; continue; }
        stage "$cfg" "$run" "$NSTEPS_P2" || { echo "MISSING CONFIG $cfg"; continue; }
        echo "=== a0 $spec $side ($(date '+%F %T'))"
        if [ "$side" = probe ]; then
            export MOBY_A0_PROBE="$A0_ITER"
        else
            unset MOBY_A0_PROBE
        fi
        launch "$EXE" "$ranks" "$nodes" "$run"
        unset MOBY_A0_PROBE
        grep -E "exch_timing: (mpi_wait|a0_probe)" "$run/run.log" | sed 's/^/    /'
        rm -f "$run"/overhead_*.h5
    done
done

echo "=== exchange passes done $(date '+%F %T')"
