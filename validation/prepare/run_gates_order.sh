#!/usr/bin/env bash
# Gates for the leaf-key bit order (blocks.f90 leaf_key / min_surface_key_order,
# docs/next_session_after_step9.md item 2).
#
# A case file built by this tree records its block order (the block_key_order
# attribute, the minimum-surface rule); a file without the record -- every file
# written before 2026-10-01 -- is read with the legacy interleave. The order
# only decides WHO owns a block, so:
#
#   unit    keyorder_test: the rule's bit orders + the partition count
#   order   a LEGACY-order case file (prepared by OLD) and a re-prepared one
#           (NEW) solve to IDENTICAL fields with the new solver, on every rank
#           count, and identical to OLD's own solve; the new file carries the
#           record and its rows are in the order tools/partition_analysis.py
#           --order minsurface gives for the old file's leaves; a body case's
#           coefficient / wall-distance tiles are the same numbers in both
#   print   the solver's "partition: face cells shared across ranks" line
#           equals tools/partition_analysis.py --order file on its case file
#   restart a snapshot written under one order restarts under the other and
#           continues exactly like a restart under its own order
#   tool    tools/h5maxdiff matches rows on (origin, level) across the orders
#
#   OLD=~/step9b_ref_binaries/build_cpu NEW=../../build_cpu ./run_gates_order.sh [group]
#
# Environment: OLD (build dir of a pre-order tree: moby_prepare + moby_solve,
# required), NEW (build dir, default ../../build_cpu), PY (python with h5py),
# H5MAXDIFF (default ../../tools/h5maxdiff; the tool group is skipped when it
# is not built). Cases: the min_channel tutorial (xyz, two refine boxes), the
# analytic wavy wall under refine_body, and an nb-less box at 8 ranks (one
# block per rank, the layout of the 2026-09-30 node-boundary finding).
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
OLD=${OLD:?set OLD to a build dir of a tree that predates the order record}
OLD=$(cd "$OLD" && pwd)
NEW=${NEW:-$ROOT/build_cpu}
NEW=$(cd "$NEW" && pwd)
PY=${PY:-python3}
H5MAXDIFF=${H5MAXDIFF:-$ROOT/tools/h5maxdiff}
CMP="$PY $ROOT/tools/compare_fields.py --tolerance 0"
sel=${1:-all}
pass=0; fail=0
want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
check() { local name="$1"; shift
    if "$@" > .check.out 2>&1; then echo "PASS  $name"; pass=$((pass+1))
    else echo "FAIL  $name"; sed 's/^/        /' .check.out | tail -5; fail=$((fail+1)); fi; }

mkdir -p order && cd order
export OMPI_MCA_hwloc_base_binding_policy=none PRTE_MCA_hwloc_default_binding_policy=none

if want unit; then
    check "keyorder_test" mpirun -n 1 "$NEW/keyorder_test"
fi

cat > box.ini <<'INI'
[case]
name = generic

[grid]
nx = 128
ny = 8
nz = 16
lx = 16.0
ly = 1.0
lz = 2.0

[grid.x]
distribution = uniform
[grid.y]
distribution = uniform
[grid.z]
distribution = uniform

[flow]
re = 100.0
initial_u = 1.0
initial_v = 0.5
initial_w = 0.25

[time]
dt = 0.001
dtmax = 0.001
nsteps = 5
t_final = 0.0
cflmax = 0.8
pecletmax = 0.5

[pressure]
niter = 4

[boundary]
periodic_x = true
periodic_y = true
periodic_z = true

[output]
field_interval = 5
field_prefix = box
INI

# name | base ini | prepare ranks | solve ranks | steps | datasets | grid | nb | periodic
CASES=(
  "minch|$ROOT/tutorials/min_channel/input.ini|1|1 4|20|un vn wn pn|128 64 8|8 8 8|1 0 1"
  "wavyr|$ROOT/validation/prepare/wavy_refine.ini|1|1|1|un vn wn pn|64 64 16|8 8 8|1 0 1"
  "box|box.ini|8|8|5|un vn wn pn|128 8 16|64 4 8|1 1 1"
)

# derive <out.ini> from <base>: prefix, case file, steps, optional restart file
derive() {
    local base=$1 out=$2 pfx=$3 casef=$4 steps=$5 restart=${6:-}
    sed -e "s/^nsteps.*/nsteps = $steps/" -e "s/^t_final.*/t_final = 0.0/" \
        -e "s/^field_interval.*/field_interval = $steps/" \
        -e "s/^field_prefix.*/field_prefix = $pfx/" \
        -e "s/^stats_sample_interval.*/stats_sample_interval = 0/" \
        -e "s/^stats_write_interval.*/stats_write_interval = 0/" "$base" > "$out"
    "$ROOT/tools/set_case_file.sh" "$out" "$casef"
    [ -n "$restart" ] && printf '\n[restart]\nfile = %s\n' "$restart" >> "$out"
    return 0
}
last() { ls -t "$1"_[0-9]*.h5 2>/dev/null | head -1; }
run() {  # ranks binary ini log
    mpirun -n "$1" "$2" "$3" > "$4" 2>&1
}

for entry in "${CASES[@]}"; do
    IFS='|' read -r name base pranks sranks steps datasets grid nb per <<< "$entry"
    echo "== $name"
    rm -f ${name}_*.h5 ${name}_*.ini ${name}_*.log
    # the two case files: the same ini, the two builders
    for side in old new; do
        dir=$OLD; [ $side = new ] && dir=$NEW
        derive "$base" ${name}_prep_$side.ini ${name}_$side ${name}_$side.case.h5 "$steps"
        mpirun -n "$pranks" "$dir/moby_prepare" ${name}_prep_$side.ini > ${name}_prep_$side.log 2>&1 \
            || { echo "FAIL  $name prepare ($side) -- order/${name}_prep_$side.log"; fail=$((fail+1)); continue 2; }
    done

    if want order; then
        check "$name: the new case file records its order, the old one does not" $PY - "$name" <<'PYEOF'
import sys, h5py
n = sys.argv[1]
with h5py.File(f"{n}_new.case.h5") as a, h5py.File(f"{n}_old.case.h5") as b:
    assert "block_key_order" in a.attrs and "block_key_order" not in b.attrs
    ra, rb = a["blocks"][...], b["blocks"][...]
    assert sorted(map(tuple, ra)) == sorted(map(tuple, rb)), "the leaf SETS differ"
    print("order", "".join("xyz"[d-1] for d in a.attrs["block_key_order"]),
          " rows moved:", int((ra != rb).any(axis=1).sum()), "of", len(ra))
PYEOF
        sed 's/^/        /' .check.out | tail -1
        check "$name: the builder's row order is the mirror's (partition_analysis.py minsurface)" $PY - "$name" "$ROOT" "$grid" "$nb" "$per" <<'PYEOF'
import sys, h5py
sys.path.insert(0, sys.argv[2] + "/tools")
import partition_analysis as pa
n = sys.argv[1]
grid, nb, per = ([int(v) for v in sys.argv[i].split()] for i in (3, 4, 5))
with h5py.File(f"{n}_new.case.h5") as a, h5py.File(f"{n}_old.case.h5") as b:
    new = [tuple(int(v) for v in r) for r in a["blocks"][...]]
    old = [tuple(int(v) for v in r) for r in b["blocks"][...]]
    rec = [int(d) - 1 for d in a.attrs["block_key_order"]]
nlev = max(r[3] for r in old) + 1
seq = pa.min_surface_order(grid, nb, per, nlev)
assert seq == rec, f"recorded {rec}, mirror {seq}"
assert pa.sort_rows(old, nb, [1, 1, 1], nlev, seq) == new, "row order differs"
assert pa.sort_rows(new, nb, [1, 1, 1], nlev, pa.legacy_order(grid, nb, nlev)) == old
PYEOF
        # a body case: the per-leaf tiles of the two files are the same numbers
        if $PY -c "import h5py,sys; sys.exit(0 if 'coef_blocks' in h5py.File('${name}_new.case.h5') else 1)"; then
            check "$name: coefficient and wall-distance tiles identical once rows are matched" \
                $PY ../compare_case.py ${name}_new.case.h5 ${name}_old.case.h5 --coef-tol 0 --dwall-tol 0
        fi
    fi

    # solves: the new solver on both files at every rank count, the old solver on its own
    first=
    for side in old new; do
        for r in $sranks; do
            pfx=${name}_${side}_r$r
            derive "$base" $pfx.ini $pfx ${name}_$side.case.h5 "$steps"
            run "$r" "$NEW/moby_solve" $pfx.ini $pfx.log || { echo "FAIL  $name solve ($side, $r ranks)"; fail=$((fail+1)); continue; }
            [ -z "$first" ] && first=$(last $pfx)
            if want order; then
                check "$name: $side-order file, $r rank(s) == the first run (max_abs 0)" $CMP "$first" "$(last $pfx)" $datasets
            fi
            if want print; then
                check "$name: partition print == partition_analysis.py ($side-order file, $r ranks)" $PY - "$pfx.log" "$ROOT" "${name}_$side.case.h5" "$r" "$grid" "$nb" "$per" <<'PYEOF'
import re, subprocess, sys
log, root, casef, ranks = sys.argv[1:5]
grid, nb, per = sys.argv[5].split(), sys.argv[6].split(), sys.argv[7].split()
m = re.search(r"partition: face cells shared across ranks (\d+) \((\d+) ranks\), across nodes (\d+)", open(log).read())
assert m, "no partition line in " + log
out = subprocess.run([sys.executable, root + "/tools/partition_analysis.py", casef, "--grid", *grid,
                      "--nb", *nb, "--periodic", *per, "--ranks", ranks, "--per-node", ranks,
                      "--order", "file"], capture_output=True, text=True, check=True).stdout
t = re.search(r"across ranks\s+([\d,]+)\s+across nodes\s+([\d,]+)", out)
tool = int(t.group(1).replace(",", ""))
print("solver", m.group(1), "tool", tool)
assert int(m.group(1)) == tool and int(m.group(2)) == int(ranks) and int(m.group(3)) == 0
PYEOF
            fi
        done
    done
    if want order; then
        r=${sranks%% *}
        pfx=${name}_oldbin
        derive "$base" $pfx.ini $pfx ${name}_old.case.h5 "$steps"
        run "$r" "$OLD/moby_solve" $pfx.ini $pfx.log \
            && check "$name: == the old solver on its own (legacy) file" $CMP "$first" "$(last $pfx)" $datasets
    fi

    if want restart; then
        r=${sranks##* }
        more=$((2*steps))
        for from in old new; do
            snap=$(last ${name}_${from}_r$r)
            for under in old new; do
                pfx=${name}_rs_${from}_under_${under}
                derive "$base" $pfx.ini $pfx ${name}_$under.case.h5 "$more" "$snap"
                run "$r" "$NEW/moby_solve" $pfx.ini $pfx.log || { echo "FAIL  $name restart ($from under $under)"; fail=$((fail+1)); }
            done
            check "$name: $from-order snapshot restarts identically under both orders" \
                $CMP "$(last ${name}_rs_${from}_under_old)" "$(last ${name}_rs_${from}_under_new)" $datasets
        done
        check "$name: the cross-order restart says so" grep -q "different block order" ${name}_rs_old_under_new.log
        check "$name: an own-order restart does not" bash -c "! grep -q 'different block order' ${name}_rs_new_under_new.log"
    fi

    if want tool && [ -x "$H5MAXDIFF" ]; then
        r=${sranks%% *}
        check "$name: h5maxdiff matches the rows of the two orders" bash -c \
            "'$H5MAXDIFF' '$(last ${name}_old_r$r)' '$(last ${name}_new_r$r)' $datasets | grep -q '^OK:'"
    fi
done
rm -f .check.out

echo
echo "order gates: $pass passed, $fail failed"
[ $fail -eq 0 ]
