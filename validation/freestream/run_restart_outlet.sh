#!/usr/bin/env bash
# Restart gate for the outlet face (docs/next_session_outlet.md, O2): a run
# stopped and restarted at mid-length must EQUAL the continuous run.
#
# The normal velocity of an outlet face is state. On a low face it is index 1
# of un/vn/wn; on a high face it is the un_xmax / vn_ymax / wn_zmax dataset
# (io.f90 write_outlet_planes). Each case runs N1 steps (a), restarts from
# that snapshot for N2 more (b), and runs N1+N2 in one go (c); b and c are
# compared with tools/compare_fields.py at tolerance 0. A last leg strips the
# plane dataset from the snapshot and restarts: the solver must say so and
# fall back to the zero-gradient value (the difference is then reported, not
# gated).
#
#   SOLVER=/abs/build_cpu/moby_solve [RANKS=1] ./run_restart_outlet.sh [case]
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
SOLVER=${SOLVER:-$ROOT/build_cpu/moby_solve}
RANKS=${RANKS:-1}
sel=${1:-all}
status=0

# case | dir | ini | N1 | N2 | ranks (0 = $RANKS)
# blasius2d sets no [blocks] nb, so its block layout -- and with it the
# initial condition IC_blasius.h5 -- belongs to the rank count that minted
# it (run_blasius.sh: one rank): it is pinned.
CASES=(
  "outlet_high|$ROOT/validation/freestream|outlet_high.ini|100|100|0"
  "outlet_low|$ROOT/validation/freestream|outlet_low.ini|100|100|0"
  "pois_io|$ROOT/validation/freestream|pois_io.ini|50|50|0"
  "lamboseen|$ROOT/validation/freestream|lamboseen.ini|60|60|0"
  "blasius2d|$ROOT/validation/blasius|blasius2d.ini|30|30|1"
)

# variant of an ini: steps, prefix, and (optionally) the restart file
variant() {  # in out steps prefix [restart-file]
    awk -v rst="${5:-}" '
        /^[ \t]*\[/ {sec = $0; gsub(/[ \t]/, "", sec)}
        sec == "[case]" && /^[ \t]*file[ \t]*=/ {next}
        rst != "" && sec == "[restart]" && /^[ \t]*file[ \t]*=/ {print "file = " rst; done = 1; next}
        {print}
        END {if (rst != "" && !done) print "[restart]\nfile = " rst}' "$1" |
    sed -e "s/^nsteps.*/nsteps = $3/" -e "s/^t_final.*/t_final = 0.0/" \
        -e "s/^field_interval.*/field_interval = 0/" \
        -e "s/^field_prefix.*/field_prefix = $4/" > "$2"
}

run() {  # ini log   (on $nr ranks, set per case)
    "$ROOT/tools/prepare_if_missing.sh" "$nr" "$SOLVER" "$1" "$2.prep" || return 1
    mpirun -n "$nr" "$SOLVER" "$1" > "$2" 2>&1
}

for entry in "${CASES[@]}"; do
    IFS='|' read -r name dir ini n1 n2 nr <<< "$entry"
    [ "$sel" = all ] || [ "$sel" = "$name" ] || continue
    [ "$nr" = 0 ] && nr=$RANKS
    echo "== $name ($nr rank(s), $n1 + $n2 steps)"
    ( cd "$dir" || exit 1
      p="rst_${name}"
      rm -f "${p}"_[abcd]_*.h5 "${p}"_[abcd].case.h5 "${p}_stripped.h5"
      variant "$ini" ".${p}_a.ini" "$n1" "${p}_a"
      variant "$ini" ".${p}_c.ini" "$((n1 + n2))" "${p}_c"
      run ".${p}_a.ini" ".${p}_a.log" || { echo "   RUN FAILED (a)"; exit 2; }
      run ".${p}_c.ini" ".${p}_c.log" || { echo "   RUN FAILED (c)"; exit 2; }
      snap=$(ls -t "${p}_a_"*.h5 | head -1)
      variant "$ini" ".${p}_b.ini" "$n2" "${p}_b" "$snap"
      run ".${p}_b.ini" ".${p}_b.log" || { echo "   RUN FAILED (b)"; exit 2; }
      b=$(ls -t "${p}_b_"*.h5 | head -1); c=$(ls -t "${p}_c_"*.h5 | head -1)
      if python3 "$ROOT/tools/compare_fields.py" "$c" "$b" un vn wn pn --tolerance 0 | sed 's/^/   /'; then
          echo "   PASS (restart == continuous, max_abs 0)"
      else
          echo "   FAIL"; exit 3
      fi
      # the fallback: the same snapshot without its outlet planes
      python3 - "$snap" "${p}_stripped.h5" <<'PY'
import sys, shutil, h5py
shutil.copy(sys.argv[1], sys.argv[2])
with h5py.File(sys.argv[2], "a") as h:
    gone = [k for k in ("un_xmax", "vn_ymax", "wn_zmax") if k in h]
    for k in gone:
        del h[k]
print("   stripped:", " ".join(gone) if gone else "(no plane dataset: a low-side outlet)")
PY
      variant "$ini" ".${p}_d.ini" "$n2" "${p}_d" "${p}_stripped.h5"
      run ".${p}_d.ini" ".${p}_d.log" || { echo "   RUN FAILED (d, stripped)"; exit 2; }
      grep -h "note: restart file has no" ".${p}_d.log" | sed 's/^/  /'
      d=$(ls -t "${p}_d_"*.h5 | head -1)
      python3 "$ROOT/tools/compare_fields.py" "$c" "$d" un vn wn pn --tolerance 1e300 \
          | sed 's/^/   fallback: /'
    ) || status=1
done

echo
[ $status -eq 0 ] && echo "outlet restart gate: ALL PASS" || echo "outlet restart gate: FAILURES"
exit $status
