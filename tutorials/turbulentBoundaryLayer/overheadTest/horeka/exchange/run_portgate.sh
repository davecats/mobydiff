#!/bin/bash
# The jacobi-interface PORT gate (docs/next_session_port_finish.md).
#
#   run_portgate.sh <new_exe> <ref_exe> <results_dir>
#
# WHY THIS IS NOT run_mapgate.sh. The port's headline change is the S3 skew
# lockdown, which is a PHYSICS change: `main` defaulted to divergence-form
# momentum convection, HEAD hardwires skew. So the two binaries do NOT agree on
# a bare ini, and a plain A/B would just confirm that. What must be bit-exact is
# the lockdown against the OLD binary RUNNING `[flow] convection = skew` -- the
# configuration the lockdown makes unconditional. This driver therefore injects
# per-SIDE ini fragments, which run_mapgate.sh cannot do.
#
# GROUP A (no scalars): ref gets `[flow] convection = skew`, new gets nothing.
# GROUP B (scalars): the old key drove BOTH the momentum and the scalar
#   operator; the port split them, because they are different operators (the
#   scalar's was the ADVECTIVE form, f = 1, not the momentum kernel's f = 1/2).
#   So ref gets `[flow] convection = skew` and new gets `[scalar] convection =
#   advective`, which is the same pair of operators either side.
#
# DATASETS ARE NAMED for every scalar case. tools/h5maxdiff's default list is
# velocity, pressure and the RANS scalars ONLY -- a passive-scalar case gated
# without arguments compares four datasets and never looks at the scalar. That
# happened during this port and is why the argument is explicit here.
set -uo pipefail

NEW="${1:?usage: run_portgate.sh <new_exe> <ref_exe> <results_dir>}"
REF="${2:?usage: run_portgate.sh <new_exe> <ref_exe> <results_dir>}"
RES="${3:?usage: run_portgate.sh <new_exe> <ref_exe> <results_dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${MOBY_ROOT:-$(cd "$HERE/../../../../.." && pwd)}"
H5="${H5MAXDIFF:-$ROOT/tools/h5maxdiff}"
NSTEPS="${NSTEPS:-200}"

SKEW=$'\n[flow]\nconvection = skew\n'
ADV=$'\n[scalar]\nconvection = advective\n'

mkdir -p "$RES"
pass=0; fail=0

leg() {   # leg <tag> <ini> <ranks> <ref-extra> <new-extra> <datasets...>
    local tag="$1" ini="$2" ranks="$3" rx="$4" nx="$5"; shift 5
    local ds="$*" side bin extra d pfx a b out
    local src="$ROOT/$ini"
    [ -f "$src" ] || { echo "  MISSING CASE $ini"; fail=$((fail+1)); return; }
    for side in ref new; do
        bin="$REF"; extra="$rx"
        [ "$side" = new ] && { bin="$NEW"; extra="$nx"; }
        d="$RES/${tag}_${side}"
        [ -f "$d/run.log" ] && continue
        rm -rf "$d"; cp -a "$(dirname "$src")" "$d"
        sed -e "s/^nsteps *=.*/nsteps = $NSTEPS/" -e "s/^t_final *=.*/t_final = 0.0/" \
            "$src" > "$d/cfg.ini"
        printf '%s' "$extra" >> "$d/cfg.ini"
        ( cd "$d" && mpirun -n "$ranks" --bind-to core --map-by numa "$bin" cfg.ini \
              > run.log 2>&1 )
        [ $? -ne 0 ] && { echo "  $tag $side FAILED -- $d/run.log"; fail=$((fail+1)); return; }
    done
    pfx=$(awk -F= '/^field_prefix/ {gsub(/ /,"",$2); print $2}' "$RES/${tag}_new/cfg.ini" | tail -1)
    a=$(ls "$RES/${tag}_ref/${pfx}"_[0-9]*.h5 2>/dev/null | tail -1)
    b=$(ls "$RES/${tag}_new/${pfx}"_[0-9]*.h5 2>/dev/null | tail -1)
    if [ -z "$a" ] || [ -z "$b" ]; then
        echo "  $tag: NO SNAPSHOT PAIR"; fail=$((fail+1)); return
    fi
    out=$("$H5" "$a" "$b" $ds | tail -1)
    printf "  %-16s %s\n" "$tag" "$out" | tee -a "$RES/summary.txt"
    case "$out" in OK*) pass=$((pass+1));; *) fail=$((fail+1));; esac
}

echo "=== port gate start $(date '+%F %T')  nsteps=$NSTEPS"
echo "--- GROUP A: no scalars (ref runs convection = skew)"
leg min_channel    tutorials/min_channel/input_gpu.ini                      4 "$SKEW" ""
leg min_channel1   tutorials/min_channel/input_gpu.ini                      1 "$SKEW" ""
leg beltrami_slaby validation/beltrami/slab_y.ini                           1 "$SKEW" ""
leg les_ibm        validation/channel_interface/les_ibm/channel_ibm.ini     1 "$SKEW" ""
leg turb180        validation/rans_sst/turb180.ini                          1 "$SKEW" ""
leg wf180_y30      validation/rans_sst/wf180_y30.ini                        1 "$SKEW" ""
leg lam30t         validation/rans_sst/lam30t.ini                           1 "$SKEW" ""

echo "--- GROUP B: scalars (ref [flow] skew == new [scalar] advective)"
leg conduction     validation/scalar/conduction.ini  1 "$SKEW" "$ADV" un vn wn pn s1
leg prsweep        validation/scalar/prsweep.ini     1 "$SKEW" "$ADV" un vn wn pn s1
leg wave           validation/scalar/wave.ini        1 "$SKEW" "$ADV" un vn wn pn s1
leg ibmwavy        validation/scalar/ibmwavy.ini     1 "$SKEW" "$ADV" un vn wn pn theta phi
leg ibmwavyr       validation/scalar/ibmwavyr.ini    1 "$SKEW" "$ADV" un vn wn pn theta phi

echo "=== port gate done $(date '+%F %T'):  $pass OK, $fail NOT OK"
exit $(( fail > 0 ))
