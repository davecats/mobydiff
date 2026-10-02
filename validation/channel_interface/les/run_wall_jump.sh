#!/usr/bin/env bash
# The SGS model at a level jump that touches a wall (wall_jump.ini,
# check_wall_jump.py): one step, CPU.
#
#     [SOLVER=../../../build_cpu/moby_solve] ./run_wall_jump.sh
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../../.. && pwd)
SOLVER=${SOLVER:-$ROOT/build_cpu/moby_solve}
rm -f wall_jump.case.h5 wall_jump_1.h5
"$ROOT/tools/prepare_if_missing.sh" 1 "$SOLVER" wall_jump.ini wall_jump.prep.log || { echo "PREPARE FAILED"; exit 1; }
mpirun -n 1 "$SOLVER" wall_jump.ini > wall_jump.log 2>&1 || { echo "RUN FAILED: wall_jump.log"; exit 1; }
python3 check_wall_jump.py wall_jump_1.h5
