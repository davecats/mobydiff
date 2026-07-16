#!/usr/bin/env python3
"""Shift-invariance gate for the STL minimum-image logic: the sphere
translated by exactly -lx/2 (an EXACT float32 translation of the same
mesh) STRADDLES the x-periodic boundary, so its case file must be the
x-rolled copy of the centred sphere's -- masks and the leaf-table SET
exactly, matching leaves' coef/dwall tiles to float64 round-off (the
absolute coordinates entering the arithmetic differ).

usage: shift_check.py centred.h5 shifted.h5 <shift_cells_level0>
"""
import sys

import h5py
import numpy as np


def main():
    cen_path, shf_path, shift0 = sys.argv[1], sys.argv[2], int(sys.argv[3])
    ok = True
    with h5py.File(cen_path, "r") as cen, h5py.File(shf_path, "r") as shf:
        nb = int(cen.attrs["block_nb"])
        gs = np.array([cen.attrs["nx"], cen.attrs["ny"], cen.attrs["nz"]], dtype=int)

        # masks: roll level-l rasters by the shift expressed in level-l blocks.
        level = 0
        while f"block_touch_l{level}" in cen:
            lat = (gs // nb) * 2**level
            sblk = (shift0 * 2**level) // nb
            for name in (f"block_touch_l{level}", f"block_buried_l{level}"):
                a = cen[name][...].reshape(lat[2], lat[1], lat[0])
                b = shf[name][...].reshape(lat[2], lat[1], lat[0])
                if not np.array_equal(np.roll(a, -sblk, axis=2), b):
                    print(f"{name}: shifted raster differs")
                    ok = False
                else:
                    print(f"{name}: shift-identical")
            level += 1

        # leaf tables: identical as SETS after shifting origins (Morton
        # order differs, so match leaves through an origin lookup).
        cb, sb = cen["blocks"][...], shf["blocks"][...]
        if cb.shape != sb.shape:
            print(f"leaf counts differ: {cb.shape} vs {sb.shape}")
            sys.exit(1)
        index = {}
        for row in range(sb.shape[0]):
            index[tuple(sb[row])] = row
        cc, cd = cen["coef_blocks"], shf["coef_blocks"]
        ca = cen["dwall_blocks"][...] if "dwall_blocks" in cen else None
        da = shf["dwall_blocks"][...] if "dwall_blocks" in shf else None
        max_coef = 0.0
        max_dwall = 0.0
        unmatched = 0
        for row in range(cb.shape[0]):
            ox, oy, oz, lev = cb[row]
            nxl = gs[0] * 2**lev
            key = ((ox + shift0 * 2**lev) % nxl, oy, oz, lev)
            srow = index.get(key)
            if srow is None:
                unmatched += 1
                continue
            max_coef = max(max_coef, float(np.max(np.abs(cc[row] - cd[srow]))))
            if ca is not None and da is not None:
                max_dwall = max(max_dwall, float(np.max(np.abs(ca[row] - da[srow]))))
        print(f"leaves: {cb.shape[0]}, unmatched after shift: {unmatched}")
        print(f"coef tiles: max |dc| = {max_coef:.3e}   dwall tiles: max |dd| = {max_dwall:.3e}")
        if unmatched > 0 or max_coef > 1.0e-9 or max_dwall > 1.0e-11:
            ok = False
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
