#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: NX=... NY=... NZ=... DIMS='dx dy dz' $0 output.ini" >&2
    exit 1
fi

OUT="$1"
mkdir -p "$(dirname "$OUT")"

NX="${NX:-200}"
NY="${NY:-$NX}"
NZ="${NZ:-$NX}"
LX="${LX:-1.0}"
LY="${LY:-1.0}"
LZ="${LZ:-1.0}"
DIMS="${DIMS:-0 0 0}"
RE="${RE:-100.0}"
DT="${DT:-1.0e-4}"
NSTEPS="${NSTEPS:-20}"
CFLMAX="${CFLMAX:-0.0}"
DTMAX="${DTMAX:-1.0e-3}"
NITER="${NITER:-3}"
SOR="${SOR:-1.5}"
FORCING_X="${FORCING_X:-1.0}"
FORCING_Y="${FORCING_Y:-0.0}"
FORCING_Z="${FORCING_Z:-0.0}"
IBM_ENABLED="${IBM_ENABLED:-true}"
FIELD_INTERVAL="${FIELD_INTERVAL:-0}"
FIELD_PREFIX="${FIELD_PREFIX:-field}"
GRID_X_DISTRIBUTION="${GRID_X_DISTRIBUTION:-uniform}"
GRID_Y_DISTRIBUTION="${GRID_Y_DISTRIBUTION:-uniform}"
GRID_Z_DISTRIBUTION="${GRID_Z_DISTRIBUTION:-uniform}"
GRID_X_STRETCH="${GRID_X_STRETCH:-0.0}"
GRID_Y_STRETCH="${GRID_Y_STRETCH:-0.0}"
GRID_Z_STRETCH="${GRID_Z_STRETCH:-0.0}"
PERIODIC_X="${PERIODIC_X:-true}"
PERIODIC_Y="${PERIODIC_Y:-false}"
PERIODIC_Z="${PERIODIC_Z:-true}"

cat > "$OUT" <<EOF
[grid]
nx = $NX
ny = $NY
nz = $NZ
lx = $LX
ly = $LY
lz = $LZ

[grid.x]
distribution = $GRID_X_DISTRIBUTION
stretch = $GRID_X_STRETCH

[grid.y]
distribution = $GRID_Y_DISTRIBUTION
stretch = $GRID_Y_STRETCH

[grid.z]
distribution = $GRID_Z_DISTRIBUTION
stretch = $GRID_Z_STRETCH

[mpi]
dims = $DIMS

[flow]
re = $RE
forcing_x = $FORCING_X
forcing_y = $FORCING_Y
forcing_z = $FORCING_Z

[time]
dt = $DT
nsteps = $NSTEPS
cflmax = $CFLMAX
dtmax = $DTMAX

[pressure]
niter = $NITER
sor = $SOR

[ibm]
enabled = $IBM_ENABLED

[boundary]
periodic_x = $PERIODIC_X
periodic_y = $PERIODIC_Y
periodic_z = $PERIODIC_Z

y_min_p_type = neumann
y_max_p_type = neumann

[output]
field_interval = $FIELD_INTERVAL
field_prefix = $FIELD_PREFIX
EOF
