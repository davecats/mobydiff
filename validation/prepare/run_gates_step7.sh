#!/usr/bin/env bash
# Step-7 gates (docs/next_session_prepare_everything.md, increments 7-1/7-2):
# a case PREPARED by moby_prepare and SOLVED from the case file must
# reproduce the reference binary's inline solve at max_abs 0, at production
# flags (the file path copies values, nothing is recomputed), for
#
#   body-free cases with an explicit [blocks] nb   (7-1: prepare accepts them;
#       prepare 1 == 4 ranks writes IDENTICAL files),
#   body-free cases with NO [blocks] nb            (the nb rule: prepare on P
#       ranks stores nb = grid/dims(P); DIFFERENT blocks tables, SAME fields),
#   an analytic-body case                          (7-2: the solver takes the
#       leaf table from the file and never rebuilds it).
#
#   REF=~/step7_ref_binaries/moby_solve_cpu NEW=../../build_cpu ./run_gates_step7.sh
#
# Environment: REF (reference SOLVER binary, required), NEW (build dir holding
# moby_solve + moby_prepare, default ../../build_cpu), MODE (label, cpu|gpu),
# PY (python with h5py), SOLVE_RANKS (default "1 4"). GPU: REF is the GPU
# reference solver, PREP the CPU prepare (canonical) and SOLVE_RANKS=1 --
# the GPU solve from the CPU-prepared file vs the GPU inline reference.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
REF=${REF:?set REF to the reference moby_solve}
NEW=${NEW:-$ROOT/build_cpu}
MODE=${MODE:-cpu}
PY=${PY:-python3}
PREP="${PREP:-$NEW/moby_prepare}"
SOLVE="$NEW/moby_solve"
CMP="$PY $ROOT/tools/compare_fields.py --tolerance 0"
tag="s7_${MODE}"
sel=${1:-all}
pass=0; fail=0

check() { local name="$1"; shift
    if "$@" >/dev/null 2>&1; then echo "PASS  $name"; pass=$((pass+1))
    else echo "FAIL  $name"; fail=$((fail+1)); fi; }
expect_fail() { local name="$1"; shift
    if "$@" >/dev/null 2>&1; then echo "FAIL  $name (ran, expected a stop)"; fail=$((fail+1))
    else echo "PASS  $name"; pass=$((pass+1)); fi; }

# short_ini <src> <dst> <steps> <prefix> [case.h5]: a fixed-step run with its
# own output prefix; with a case file, [ibm] coeff_file points at it (the
# section is appended when the ini has none).
short_ini() {
    sed -e "s/^nsteps.*/nsteps = $3/" -e "s/^t_final.*/t_final = 0.0/" \
        -e "s/^field_interval.*/field_interval = $3/" -e "s/^field_prefix.*/field_prefix = $4/" \
        -e "s/^stats_sample_interval.*/stats_sample_interval = 0/" \
        -e "s/^stats_write_interval.*/stats_write_interval = 0/" "$1" > "$2"
    if [ -n "${5:-}" ]; then
        if grep -q '^\[ibm\]' "$2"; then
            sed -i "s|^\[ibm\]|[ibm]\ncoeff_file = $5|" "$2"
        else
            printf '\n[ibm]\ncoeff_file = %s\n' "$5" >> "$2"
        fi
    fi
}

# gate <name> <dir> <ini> <steps> <datasets> <explicit_nb 0|1> [solve ranks]
# An ini that pins [mpi] dims = 1 1 1 (the P0 wavy cases) solves on 1 rank
# only; prepare still runs on 1 and 4 (the dims describe the SOLVE there).
gate() {
    local name="$1" dir="$2" ini="$3" steps="$4" ds="$5" explicit="$6" ranks="${7:-${SOLVE_RANKS:-1 4}}"
    [ "$sel" = all ] || [ "$sel" = "$name" ] || return
    echo "== $name ($MODE)"
    ( cd "$dir" || exit 1
      rm -f "${tag}_${name}_"*.h5 ".${tag}_${name}_"*.ini
      # reference: inline solve on each solve rank count
      for r in $ranks; do
          short_ini "$ini" ".${tag}_${name}_ref${r}.ini" "$steps" "${tag}_${name}_ref${r}"
          mpirun -n $r --oversubscribe "$REF" ".${tag}_${name}_ref${r}.ini" > ".${tag}_${name}_ref${r}.log" 2>&1 \
              || { echo "   REF run failed ($r ranks)"; exit 2; }
      done
      # prepare on 1 and 4 ranks (no coeff_file in the prepare input)
      for p in 1 4; do
          short_ini "$ini" ".${tag}_${name}_prep${p}.ini" "$steps" "${tag}_${name}_x"
          mpirun -n $p --oversubscribe "$PREP" ".${tag}_${name}_prep${p}.ini" "${tag}_${name}_case${p}.h5" \
              > ".${tag}_${name}_prep${p}.log" 2>&1 || { echo "   PREPARE failed ($p ranks)"; exit 2; }
      done ) || { fail=$((fail+1)); return; }
    local d="$dir"
    if [ "$explicit" = 1 ]; then
        check "$name: prepare 1 == 4 ranks (identical files)" $PY h5same.py \
            "$d/${tag}_${name}_case1.h5" "$d/${tag}_${name}_case4.h5"
    else
        expect_fail "$name: nb-less -> prepare 1 != 4 ranks (different tables)" $PY h5same.py \
            "$d/${tag}_${name}_case1.h5" "$d/${tag}_${name}_case4.h5"
    fi
    # solve from each prepared file on 1 and 4 ranks vs the matching reference
    for p in 1 4; do for r in $ranks; do
        local pfx="${tag}_${name}_p${p}r${r}"
        ( cd "$d" && short_ini "$ini" ".${pfx}.ini" "$steps" "$pfx" "${tag}_${name}_case${p}.h5" \
            && mpirun -n $r --oversubscribe "$SOLVE" ".${pfx}.ini" > ".${pfx}.log" 2>&1 )
        local ok=$?
        if [ "$explicit" = 0 ] && [ $p = 1 ] && [ $r = 4 ]; then
            # a 1-rank-derived nb is the whole grid as ONE block: 4 ranks cannot own it
            if [ $ok -ne 0 ] && grep -q "rank owns no blocks" "$d/.${pfx}.log"; then
                echo "PASS  $name: 1-rank auto file on 4 ranks stops (rank owns no blocks)"; pass=$((pass+1))
            else
                echo "FAIL  $name: 1-rank auto file on 4 ranks did not stop as documented"; fail=$((fail+1))
            fi
            continue
        fi
        if [ $ok -ne 0 ]; then echo "FAIL  $name: solve from case${p} on $r ranks (run failed)"; fail=$((fail+1)); continue; fi
        check "$name: prepared(${p}) + solved(${r}) == inline ref(${r})" $CMP \
            "$d/${tag}_${name}_ref${r}_${steps}.h5" "$d/${pfx}_${steps}.h5" $ds
    done; done
}

gate min_channel "$ROOT/tutorials/min_channel"      input.ini     20 "un vn wn pn" 1
gate beltrami    "$ROOT/validation/beltrami"        slab_y.ini    5  "un vn wn pn" 1
gate wf180_y30   "$ROOT/validation/rans_sst"        wf180_y30.ini 20 "un vn wn pn nut k omega" 0
gate conduction  "$ROOT/validation/scalar"          conduction.ini 50 "un vn wn pn s1" 0
gate wavy        "$ROOT/validation/prepare"         wavy.ini      1  "un vn wn pn" 1 "1"
# Step 7-4: the live analytic-IBM inis that used to run the solver's inline
# coefficient kernel now run as prepare+solve pairs (the solver prepares
# in-process; here the file is prepared explicitly and gated the same way).
# On the GPU compare these against the CPU reference (REF=moby_solve_cpu):
# the reference GPU binary computed the analytic coefficients on the device
# (libm ulps), the case file carries the host kernel's.
gate conj_wavy   "$ROOT/validation/conjugate"       wavy.ini      20 "un vn wn pn theta" 1
gate ibmwavy     "$ROOT/validation/scalar"          ibmwavy.ini   20 "un vn wn pn theta phi" 1
gate ibmwavyr    "$ROOT/validation/scalar"          ibmwavyr.ini  20 "un vn wn pn theta phi" 1
gate rg_wavy     "$ROOT/validation/rans_geometry"   wavy.ini      1  "un vn wn pn" 1 "1"
gate rg_wavyref  "$ROOT/validation/rans_geometry"   wavy_refine.ini 1 "un vn wn pn" 1 "1"

echo
echo "passed: $pass  failed: $fail"
[ "$fail" -eq 0 ]
