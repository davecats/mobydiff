#!/usr/bin/env bash
# Regenerate the prerequisite data files for the LES<->IBM coupling validation.
# Needed only if you don't have the committed files / want to rebuild from scratch.
# Requires: the CPU build (moby_prepare), a solver build for the cold start, and
# the geometry venv with trimesh + h5py (here: /home/davide/ibmc/bin/python) for
# the STL and IC generators. On another host, point PY/PREP/SOLVE/MPIRUN at the
# local equivalents, or just rsync this whole directory (the .h5 files travel
# with it) and skip setup entirely.
#
# UPDATED 2026-09-26. The original built grid.h5 with `build_cpu/mobygrid` and
# the coefficient files with `tools/mobygeom.py`; mobygrid was DELETED in the
# prepare/solve split P3 and mobygeom's geometry subcommands retired, so this
# script could not run. The block-table file for case c is now written by
# moby_prepare (validation/prepare flat_refine gate: blocks + masks identical to
# the mobygeom file). NOT regenerated any more, because they are committed:
#   grid.h5      -- the old mobygrid grid; still read as --grid-file by the
#                   mobygeom reference in ../../rans_geometry/setup.sh
#   ibm_coeff_case.h5 -- the retired mobygeom's single-level legacy-format
#                   coefficients converted once, number for number, into the
#                   case-file layout (`moby_prepare --convert-legacy` at
#                   numerics-review step 7-3, gated max_abs 0 against a solve
#                   from the legacy file, which left the tree on 2026-09-30).
#                   The standard suite's les_ibm case and measure_nut.py read
#                   it. A moby_prepare file of channel_ibm.ini with stl_file
#                   is its recomputed equivalent, ~1e-10 apart.
#
# Usage:  ./setup.sh
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-/home/davide/ibmc/bin/python}"
ROOT=../../..
PREP="${PREP:-$ROOT/build_cpu/moby_prepare}"
SOLVE="${SOLVE:-${MAIN:-$ROOT/build_gpu/moby_solve}}"
MPIRUN="${MPIRUN:-/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/comm_libs/13.1/hpcx/hpcx-2.25.1/ompi/bin/mpirun}"

echo "== 1. wall STLs (two solid slabs, inside=solid convention)"
$PY make_walls_stl.py

echo "== 2. block-table coefficient file for refine_body (ibm_coeff_blocks.h5, moby_prepare)"
# The case ini with the coefficient file swapped for the STLs and [restart]
# stripped (the IC does not exist yet and must not drive the grid).
sed -e '/^\[restart\]/,$d' \
    -e 's|^\[ibm\]|[ibm]\nstl_file = wall_lo.stl\nstl_file = wall_hi.stl|' \
    channel_ibm_refine.ini > _prep.ini
$MPIRUN -n 1 "$PREP" _prep.ini ibm_coeff_blocks.h5

echo "== 3. cold-start a 1-step run to mint a field with correct attrs (cs_1.h5)"
sed -e '/^\[restart\]/,$d' \
    -e 's/^large_disturbance_amplitude = .*/large_disturbance_amplitude = 1.0e-2/' \
    -e 's/^small_noise_amplitude = .*/small_noise_amplitude = 1.0e-3/' \
    -e 's/^field_interval = .*/field_interval = 1/' \
    -e 's/^field_prefix = .*/field_prefix = cs/' \
    -e 's/^t_final = .*/t_final = 1.0e-4/' -e 's/^nsteps = .*/nsteps = 1/' \
    channel_ibm.ini > _coldstart.ini
$MPIRUN -n 1 "$SOLVE" _coldstart.ini

echo "== 4. ICs: KMM180 developed field mapped into the fluid gap"
$PY make_ibm_ic.py                                           # IC.h5 (single level)
$PY make_ibm_ic.py --leaves ibm_coeff_blocks.h5 --out IC_refine.h5   # IC_refine.h5

rm -f _prep.ini _coldstart.ini cs_1.h5
echo "== done. Prerequisites ready. Now: python3 run_ibm_les.py --mpirun \"$MPIRUN\""
