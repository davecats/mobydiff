#!/usr/bin/env bash
# A3 INCREMENT 0 gates — multi-level (3-level) refine_body prerequisite
# (docs/next_session_airfoil.md, A3): the leaf table, exchange entries and
# dwall tiles at levels > 1 before any airfoil physics.
#
#   ./run_gates.sh              # all gates, sequentially (one job at a time)
#   ./run_gates.sh uniform      # one group: uniform | dwall
#
# Environment:  BIN=<solver>   (default ../../build_cpu/main)
#               BIN_GPU=<gpu>  (default ../../build_gpu/main)
set -o pipefail
cd "$(dirname "$0")"
set -u

BIN=${BIN:-../../build_cpu/main}
BIN_GPU=${BIN_GPU:-../../build_gpu/main}
CMP="python3 ../../tools/compare_fields.py"
U0=0.9396926207859084
V0=0.3420201433256687
W0=0.2
sel=${1:-all}
status=0

run() { # bin ini ranks log
    local bin=$1 ini=$2 ranks=$3 log=$4
    echo "== $log (ranks $ranks) =="
    if ! mpirun -n "$ranks" "$bin" "$ini" > "$log.log" 2>&1; then
        echo "   FAILED — see $log.log"; status=1; return 1
    fi
    tail -n 2 "$log.log" | sed 's/^/   /'
}

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }

# --- gate (a): uniform oblique flow across every interface level, EXACT ---
if want uniform; then
    for r in 1 4; do
        sed -e "s/^field_prefix.*/field_prefix = uni_r${r}/" uniform.ini > .r.ini
        rm -f uni_r${r}_*.h5
        run "$BIN" .r.ini $r uni_r$r && \
        python3 check_uniform.py "$(ls -t uni_r${r}_*.h5 | head -1)" \
            --u0 $U0 --v0 $V0 --w0 $W0 --levels 3 || status=1
        rm -f .r.ini
    done
    if [ -n "$(ls -t uni_r1_*.h5 2>/dev/null)" ] && [ -n "$(ls -t uni_r4_*.h5 2>/dev/null)" ]; then
        $CMP "$(ls -t uni_r1_*.h5 | head -1)" "$(ls -t uni_r4_*.h5 | head -1)" --tolerance 0 \
            && echo "uniform 1==4 ranks: EXACT" || { echo "uniform ranks MISMATCH"; status=1; }
    fi
    if [ -x "$BIN_GPU" ]; then
        sed -e 's/^field_prefix.*/field_prefix = uni_gpu/' uniform.ini > .g.ini
        rm -f uni_gpu_*.h5
        run "$BIN_GPU" .g.ini 1 uni_gpu && \
        python3 check_uniform.py "$(ls -t uni_gpu_*.h5 | head -1)" \
            --u0 $U0 --v0 $V0 --w0 $W0 --levels 3 || status=1
        rm -f .g.ini
    fi
fi

# --- gate (b): per-level dwall vs the exact STL-polygon distance ---
if want dwall; then
    rm -f mlb_dwall_*.h5 mlb_dwall_ransgeom.h5
    run "$BIN" dwall.ini 1 dwall && \
    python3 check_dwall_cylinder.py mlb_dwall_ransgeom.h5 \
        ../cylinder/cylinder.stl || status=1
    # rank independence of the dump (T1 gate pattern)
    mv mlb_dwall_ransgeom.h5 .dwall_r1.h5 2>/dev/null
    rm -f mlb_dwall_*.h5
    run "$BIN" dwall.ini 4 dwall_r4 && {
        python3 - <<'EOF' && echo "ransgeom 1==4 ranks: identical" \
            || { echo "ransgeom rank MISMATCH"; status=1; }
import sys
import h5py
import numpy as np
with h5py.File(".dwall_r1.h5") as a, h5py.File("mlb_dwall_ransgeom.h5") as b:
    ok = all(np.array_equal(a[n][...], b[n][...]) for n in a)
sys.exit(0 if ok else 1)
EOF
    }
fi

echo
[ $status -eq 0 ] && echo "ALL SELECTED GATES PASS" || echo "GATE FAILURES (see above)"
exit $status
