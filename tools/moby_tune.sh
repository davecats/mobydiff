#!/usr/bin/env bash
# moby_tune -- find this machine's best rank-to-GPU mapping by MEASURING it.
#
#     tools/moby_tune.sh <moby_solve> <case.ini> [outdir]
#     RANKS=16 NODES=4 tools/moby_tune.sh ./build_gpu/moby_solve case.ini
#
# WHY MEASURE INSTEAD OF MODEL. The two ends of a cross-node link must sit in the
# same GPU affinity class or the exchange collapses (8x on the wait, 20-27 % of
# the step -- overheadTest/results_horeka_2026-09-10.md). WHICH devices share a
# class is a property of the node's PCIe/NUMA layout, and deducing it needs a
# vendor topology API and a model of how the fabric reacts. Both were got wrong
# twice in the session that found the effect. Timing four short runs is not.
#
# WHAT IT SWEEPS, and why that is the whole space. comm.f90's select_target_device
# serves a node's FIRST local rank from MOBY_GPU_ORDER[0] and its LAST from
# MOBY_GPU_ORDER[1], identically on every node -- so every cross-node link joins
# those two devices. The mapping question is therefore exactly "which unordered
# pair of devices should lead the order", which is C(ndev,2) candidates: 6 on a
# 4-GPU node. Devices not named follow in index order and only affect intra-node
# pairs, which go over NVLink and do not care.
#
# WHAT IT DOES NOT TUNE. Not `nb`, not the rank count, not the decomposition.
# Those interact with the case rather than the machine, change memory footprint
# and leaf tables, and their effects sit in the few-percent band where a tuner
# would be fitting node-to-node variation (~6 % here) instead of an optimum.
#
# TUNE ONCE PER MACHINE. The answer is a topology property: it transfers across
# cases and runs. Put the printed line in your submit script -- do NOT re-tune
# inside a production run, or its performance stops being reproducible.
set -uo pipefail

EXE="${1:?usage: moby_tune.sh <moby_solve> <case.ini> [outdir]}"
CASE="${2:?usage: moby_tune.sh <moby_solve> <case.ini> [outdir]}"
OUT="${3:-moby_tune_out}"

RANKS="${RANKS:-8}"
NODES="${NODES:-2}"
NSTEPS="${NSTEPS:-40}"          # ratios within the sweep; not comparable to production
# Percent. Below this, a candidate is not better than the incumbent -- it is a
# different draw. Set from MEASURED node-to-node variation on the target machine
# (~6 % on HoreKa), not from within-allocation repeatability, which is much
# tighter and would let a tuner recommend a winner that does not survive the next
# allocation. Calibration evidence: at 16 ranks pair 2,3 led by 4.07 %, while at
# 8 ranks the SAME pair lost to 0,1 -- the ranking between two equally-matched
# classes is not a stable node property, so 3 % was too permissive.
NOISE="${NOISE:-5.0}"
MPIRUN_EXTRA="${MPIRUN_EXTRA:-}"

command -v mpirun >/dev/null || { echo "ERROR: mpirun not on PATH" >&2; exit 1; }
[ -x "$EXE" ]  || { echo "ERROR: not executable: $EXE" >&2; exit 1; }
[ -f "$CASE" ] || { echo "ERROR: no such case file: $CASE" >&2; exit 1; }
if [ $(( RANKS % NODES )) -ne 0 ]; then
    echo "ERROR: $RANKS ranks do not divide over $NODES nodes" >&2; exit 1
fi
PER_NODE=$(( RANKS / NODES ))
mkdir -p "$OUT"
CASE_ABS="$(cd "$(dirname "$CASE")" && pwd)/$(basename "$CASE")"
EXE_ABS="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"
OUT_ABS="$(cd "$OUT" && pwd)"

# One trial. Echoes "<s/step> <wait_us_per_round> <ndev>", or nothing on failure.
trial() {
    # Separate statements on purpose: `local a=$1 b=$a` evaluates b before a is
    # visible under `set -u` in this bash, which cost a whole calibration job.
    local tag="$1"
    local order="${2:-}"
    local d="$OUT_ABS/$tag"
    rm -rf "$d"; mkdir -p "$d"
    cp "$CASE_ABS" "$d/config.ini"
    if grep -qE '^[[:space:]]*profile[[:space:]]*=' "$d/config.ini"; then
        sed -i 's/^profile *=.*/profile = true/' "$d/config.ini"
    else
        printf '\n[output]\nprofile = true\n' >> "$d/config.ini"
    fi
    sed -i -e "s/^nsteps *=.*/nsteps = $NSTEPS/" \
           -e "s/^t_final *=.*/t_final = 0.0/" \
           -e "s/^field_interval *=.*/field_interval = 0/" \
           -e "s/^runtime_interval *=.*/runtime_interval = $NSTEPS/" "$d/config.ini"
    ( cd "$d" && env ${order:+MOBY_GPU_ORDER=$order} \
        mpirun -n "$RANKS" --map-by "ppr:${PER_NODE}:node" --bind-to core \
        -x MOBY_GPU_ORDER $MPIRUN_EXTRA "$EXE_ABS" config.ini > run.log 2>&1 < /dev/null )
    [ -f "$d/run.log" ] || return 1
    local sps wait ndev
    sps=$(grep -m1 '^timing: nsteps' "$d/run.log" | awk '{print $NF}')
    wait=$(grep -m1 '^exch_timing: mpi_wait' "$d/run.log" | awk '{printf "%.1f", $8/$4*1e6}')
    ndev=$(grep -m1 'devices available' "$d/run.log" | awk '{print $NF}')
    [ -n "$sps" ] || return 1
    echo "$sps ${wait:-0} ${ndev:-1}"
    rm -f "$d"/*.h5
}

echo "=== moby_tune: $RANKS ranks on $NODES node(s) = $PER_NODE/node, $NSTEPS steps"
echo "    case $CASE_ABS"
echo

# Probe with the built-in order, both to get a baseline and to learn ndev.
read -r base_sps base_wait NDEV <<< "$(trial probe "")"
if [ -z "${base_sps:-}" ]; then
    echo "ERROR: the probe run failed -- see $OUT_ABS/probe/run.log" >&2; exit 1
fi
printf "  %-14s %-18s s/step %s   wait/round %s us   (%s device(s) visible)\n" \
    "built-in" "(no MOBY_GPU_ORDER)" "$base_sps" "$base_wait" "$NDEV"

if [ "$NDEV" -lt 2 ]; then
    echo
    echo "Only $NDEV device visible per rank: the mapping is not this process's to"
    echo "choose (CUDA_VISIBLE_DEVICES is already pinned externally). Nothing to tune."
    exit 0
fi
if [ "$NODES" -lt 2 ]; then
    echo
    echo "Single node: there are no cross-node links to protect, and"
    echo "select_target_device keeps the identity mapping. Nothing to tune."
    exit 0
fi

names=("built-in"); orders=(""); spss=("$base_sps"); waits=("$base_wait")
for ((a = 0; a < NDEV; a++)); do
    for ((b = a + 1; b < NDEV; b++)); do
        read -r s w _ <<< "$(trial "pair_${a}_${b}" "$a,$b")"
        [ -n "${s:-}" ] || { echo "  pair $a,$b FAILED"; continue; }
        printf "  %-14s %-18s s/step %s   wait/round %s us\n" "pair $a,$b" "MOBY_GPU_ORDER=$a,$b" "$s" "$w"
        names+=("pair $a,$b"); orders+=("$a,$b"); spss+=("$s"); waits+=("$w")
    done
done

# Rank by s/step -- the quantity that matters -- and confirm the leader against
# the runner-up with a repeat, because a single trial cannot tell a real margin
# from node-to-node drift.
best=0; second=-1
for i in "${!spss[@]}"; do
    awk "BEGIN{exit !(${spss[$i]} < ${spss[$best]})}" && { second=$best; best=$i; continue; }
    if [ $second -lt 0 ] || awk "BEGIN{exit !(${spss[$i]} < ${spss[$second]})}"; then second=$i; fi
done

echo
echo "--- repeat of the incumbent and the leader, to separate a margin from drift ---"
# The decision is BEST vs the INCUMBENT (the built-in order), not best vs
# runner-up: several candidates tying with each other says nothing about whether
# any of them beats what the solver already does. An early version got this
# backwards and recommended keeping a built-in order that both leaders beat.
for i in 0 $best; do
    [ $i -eq 0 ] && [ $best -eq 0 ] && { [ $i -eq 0 ] || continue; }
    read -r s w _ <<< "$(trial "repeat_$i" "${orders[$i]}")"
    printf "  %-14s repeat  s/step %s   wait/round %s us\n" "${names[$i]}" "${s:-FAILED}" "${w:-}"
    [ -n "${s:-}" ] && spss[$i]=$(awk "BEGIN{printf \"%.8e\", (${spss[$i]}+$s)/2}")
    [ $best -eq 0 ] && break
done

gain=$(awk "BEGIN{printf \"%.2f\", 100*(${spss[0]}-${spss[$best]})/${spss[0]}}")
margin="n/a"
[ $second -ge 0 ] && margin=$(awk "BEGIN{printf \"%.2f\", 100*(${spss[$second]}-${spss[$best]})/${spss[$best]}}")

echo
echo "=== result"
printf "  incumbent   : %s  (s/step %s)\n" "${names[0]}" "${spss[0]}"
printf "  best        : %s  (s/step %s)\n" "${names[$best]}" "${spss[$best]}"
[ $second -ge 0 ] && printf "  runner-up   : %s  (s/step %s)\n" "${names[$second]}" "${spss[$second]}"
printf "  gain vs incumbent : %s %%   (tie threshold %s %%)\n" "$gain" "$NOISE"
printf "  lead over runner-up: %s %%\n" "$margin"
echo
if [ "$best" -eq 0 ]; then
    echo "  The built-in order is already the best measured. Set nothing."
elif awk "BEGIN{exit !($gain < $NOISE)}"; then
    echo "  Nothing to gain: the best candidate is within noise of what the solver"
    echo "  already does. Keep the built-in order -- a deterministic default beats"
    echo "  a coin flip."
    echo "      (set nothing)"
else
    echo "  Put this in your submit script, once, for this machine:"
    echo "      export MOBY_GPU_ORDER=${orders[$best]}"
    if awk "BEGIN{exit !($gain < 2*$NOISE)}"; then
        echo
        echo "  MARGINAL ($gain %, threshold $NOISE %). Confirm it in a SEPARATE"
        echo "  allocation before adopting it: this margin is the size of node-to-node"
        echo "  variation, and a gain that does not reproduce on other nodes is a draw."
    fi
    if [ "$margin" != "n/a" ] && awk "BEGIN{exit !($margin < $NOISE)}"; then
        echo
        echo "  (Several candidates are equivalent within noise -- ${names[$second]} would"
        echo "   do as well. They are presumably the same affinity class.)"
    fi
fi
