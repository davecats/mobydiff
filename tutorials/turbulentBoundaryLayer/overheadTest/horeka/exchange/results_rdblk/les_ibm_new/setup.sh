#!/usr/bin/env bash
# Regenerate ALL prerequisite data files for the LES<->IBM coupling validation.
# Needed only if you don't have the committed files / want to rebuild from scratch.
# Requires: a built mobygrid (build_cpu/mobygrid) + the geometry venv with
# trimesh + libigl + h5py (here: /home/davide/ibmc/bin/python). On another host,
# point PY/MOBYGRID/MPIRUN at the local equivalents, or just rsync this whole
# directory (the .h5 files travel with it) and skip setup entirely.
#
# Usage:  ./setup.sh
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-/home/davide/ibmc/bin/python}"
ROOT=../../..
MOBYGRID="${MOBYGRID:-$ROOT/build_cpu/mobygrid}"
MAIN="${MAIN:-$ROOT/build_gpu/main}"
MPIRUN="${MPIRUN:-/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/comm_libs/13.1/hpcx/hpcx-2.25.1/ompi/bin/mpirun}"

echo "== 1. grid.h5 (mobygrid, restart-stripped copy so it doesn't follow [restart])"
grep -v '^file = IC.h5' channel_ibm.ini | sed '/^\[restart\]/d' > _grid_only.ini
$MPIRUN -n 1 "$MOBYGRID" _grid_only.ini grid.h5

echo "== 2. wall STLs (two solid slabs, inside=solid convention)"
$PY make_walls_stl.py

echo "== 3. single-level coefficient file (ibm_coeff.h5)"
$PY "$ROOT/tools/mobygeom.py" stl-ibm-coeff \
    --geometry wall_lo.stl wall_hi.stl --grid-file grid.h5 --re 180 \
    --output ibm_coeff.h5 --no-tiled-output

echo "== 4. block-table coefficient file for refine_body (ibm_coeff_blocks.h5)"
$PY "$ROOT/tools/mobygeom.py" block-table \
    --geometry wall_lo.stl wall_hi.stl --grid-file grid.h5 --re 180 \
    --block-nb 8 --levels 2 --output ibm_coeff_blocks.h5

echo "== 5. cold-start a 1-step run to mint a field with correct attrs (cs_1.h5)"
sed -e '/^\[restart\]/,$d' \
    -e 's/^large_disturbance_amplitude = .*/large_disturbance_amplitude = 1.0e-2/' \
    -e 's/^small_noise_amplitude = .*/small_noise_amplitude = 1.0e-3/' \
    -e 's/^field_interval = .*/field_interval = 1/' \
    -e 's/^field_prefix = .*/field_prefix = cs/' \
    -e 's/^t_final = .*/t_final = 1.0e-4/' -e 's/^nsteps = .*/nsteps = 1/' \
    channel_ibm.ini > _coldstart.ini
$MPIRUN -n 1 "$MAIN" _coldstart.ini

echo "== 6. ICs: KMM180 developed field mapped into the fluid gap"
$PY make_ibm_ic.py                                           # IC.h5 (single level)
$PY make_ibm_ic.py --leaves ibm_coeff_blocks.h5 --out IC_refine.h5   # IC_refine.h5

rm -f _grid_only.ini _coldstart.ini cs_1.h5
echo "== done. Prerequisites ready. Now: python3 run_ibm_les.py --mpirun \"$MPIRUN\""
