#!/usr/bin/env bash
# RANS T2 gate runner (docs/next_session_iddes.md, phase T2) — the LONG
# physics runs, meant for a big machine:
#
#   rsync the repo over, build (see below), then:
#       ./run_gates.sh            # all gates, sequentially
#       ./run_gates.sh turb180    # one gate
#   rsync the produced *_<step>.h5 + *.log back and run ./check_gates.sh
#   (needs python3 + h5py/numpy only).
#
# Build on the target machine (NVHPC or any MPI Fortran stack):
#       ./compile.sh cpu
# Environment:
#       RANKS=<n>   MPI ranks (default 4; the small channels have 12-32
#                   blocks, so more than ~8 ranks does not help them —
#                   ibm180 is the big one and takes any count)
#       BIN=<path>  solver binary (default ../../build_cpu/main)
#
# Prerequisite files (already in the tree; regenerate only if missing):
#   ../rans_geometry/ibm_coeff_blocks_l1.h5   (../rans_geometry/setup.sh)
#   ../channel_interface/les_ibm/IC.h5        (committed)
set -uo pipefail
cd "$(dirname "$0")"

BIN=${BIN:-../../build_cpu/main}
RANKS=${RANKS:-4}
sel=${1:-all}
status=0

run_case() {
    local ini=$1 ranks=$2
    local name=${ini%.ini}
    echo "== $name (ranks $ranks) =="
    if ! mpirun -n "$ranks" "$BIN" "$ini" > "$name.log" 2>&1; then
        echo "   FAILED — see $name.log"
        status=1
        return
    fi
    tail -n 2 "$name.log" | sed 's/^/   /'
}

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }

# Gate (a): laminar channel, SST on — parabola + k decay.        (~minutes)
want laminar  && run_case laminar.ini  "$RANKS"

# Gate (b): developed turbulent channel vs the log law.          (~hours CPU)
want turb180  && run_case turb180.ini  "$RANKS"
want turb395  && run_case turb395.ini  "$RANKS"

# Gate (c): the les_ibm off-grid IBM channel through the IBM wall
# treatment — THE key IBM gate.                                  (heaviest)
if want ibm180; then
    if [ ! -f ../rans_geometry/ibm_coeff_blocks_l1.h5 ]; then
        echo "== ibm180 SKIPPED: ../rans_geometry/ibm_coeff_blocks_l1.h5 missing"
        echo "   (generate with ../rans_geometry/setup.sh — needs python+trimesh)"
        status=1
    else
        run_case ibm180.ini "$RANKS"
    fi
fi

# Gate (d): 2:1 band-refined RANS channel + its single-level twin.
want base180u  && run_case base180u.ini  "$RANKS"
want refine180 && run_case refine180.ini "$RANKS"

echo
echo "runs done (status $status); rsync *_*.h5 and *.log back, then ./check_gates.sh"
exit $status
