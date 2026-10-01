#!/usr/bin/env python3
"""Compare a moby_prepare case file against a mobygeom block-table
reference: blocks / per-level masks / block_active must be IDENTICAL
(integer classification -- parity and winding number must agree on clean
geometry); coef_blocks and dwall_blocks match to a tolerance (bisection
crossings and indicator-polished distances vs mobygeom's exact ray/igl
queries).

The blocks table must hold the same LEAVES; its row order may differ (a
prepared xyz file carries the minimum-surface key-bit order, a mobygeom
reference the legacy interleave), in which case the per-leaf tiles are
matched on (origin, level) before they are compared (h5rows.py).

Reference masks may be WINDOWED (mask_win_lo_l{l}/mask_win_dims_l{l},
level-l block coords, x-fastest within the window); prepared files carry
full rasters. Both are expanded to the full lattice before comparing.

usage: compare_case.py candidate.h5 reference.h5 --coef-tol T --dwall-tol T
exit 0 when everything passes.
"""
import argparse
import sys

import h5py
import numpy as np

from h5rows import match_rows


def full_raster(f, name, level, lattice_dims):
    """Return dataset `name` expanded onto the full level lattice
    (z-slowest/x-fastest raster, stored as a numpy array indexed
    [z, y, x])."""
    raw = f[name][...]
    nx, ny, nz = lattice_dims
    lo = f.attrs.get(f"mask_win_lo_l{level}")
    dims = f.attrs.get(f"mask_win_dims_l{level}")
    if lo is None:
        return raw.reshape(nz, ny, nx)
    wx, wy, wz = int(dims[0]), int(dims[1]), int(dims[2])
    out = np.zeros((nz, ny, nx), dtype=raw.dtype)
    out[lo[2]:lo[2]+wz, lo[1]:lo[1]+wy, lo[0]:lo[0]+wx] = raw.reshape(wz, wy, wx)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidate")
    ap.add_argument("reference")
    ap.add_argument("--coef-tol", type=float, default=1.0e-6,
                    help="max |dc|/(|c|+1) on coef_blocks")
    ap.add_argument("--dwall-tol", type=float, default=2.0e-9,
                    help="max |dd| on dwall_blocks")
    args = ap.parse_args()

    ok = True
    with h5py.File(args.candidate, "r") as cand, h5py.File(args.reference, "r") as ref:
        # Per-direction block size; the scalar block_nb is what cubic and
        # mobygeom-written files carry.
        def block_nb(f):
            if "block_nb_xyz" in f.attrs:
                return np.array(f.attrs["block_nb_xyz"], dtype=int)
            return np.full(3, int(f.attrs["block_nb"]), dtype=int)
        nb = block_nb(ref)
        assert np.array_equal(nb, block_nb(cand))
        gs = np.array([ref.attrs["nx"], ref.attrs["ny"], ref.attrs["nz"]], dtype=int)
        mask = np.ones(3, dtype=int)
        if "refine_dims" in ref.attrs:
            mask = np.array(ref.attrs["refine_dims"], dtype=int)

        # blocks leaf table: identical.
        cb, rb = cand["blocks"][...], ref["blocks"][...]
        perm = match_rows(cb, rb)
        if perm is None:
            print(f"blocks tables differ: {cb.shape} vs {rb.shape}")
            sys.exit(1)
        if np.array_equal(cb, rb):
            print(f"blocks: {cb.shape[0]} leaves identical")
        else:
            print(f"blocks: {cb.shape[0]} leaves identical as a set; row orders"
                  " differ, tiles matched on (origin, level)")

        # per-level masks: identical after window expansion. mobygeom
        # always writes them; prepare only for refine_body cases -- a
        # candidate without any masks skips the check.
        level = 0
        while f"block_touch_l{level}" in ref:
            if "block_touch_l0" not in cand:
                print("masks: not in candidate (non-refine_body prepare), skipped")
                break
            lattice = (gs // nb) * mask * 2**level + (gs // nb) * (1 - mask)
            for name in (f"block_touch_l{level}", f"block_buried_l{level}"):
                if name not in cand:
                    print(f"{name}: missing in candidate")
                    ok = False
                    continue
                a = full_raster(cand, name, level, lattice)
                b = full_raster(ref, name, level, lattice)
                if not np.array_equal(a, b):
                    print(f"{name}: {int(np.sum(a != b))} cells differ")
                    ok = False
                else:
                    print(f"{name}: identical ({int(b.sum())} set)")
            level += 1

        if "block_active" in ref:
            if not np.array_equal(cand["block_active"][...], ref["block_active"][...]):
                print("block_active differs")
                ok = False
            else:
                print("block_active: identical")

        # coefficient tiles. The solid/fluid CLASSIFICATION must agree
        # exactly (parity vs winding number on a clean mesh); the graded
        # values agree to the bisection tolerance except at near-grazing
        # crossings, where ((d0-d)/d)/d0^2 has exploding relative
        # sensitivity while the absolute value stays negligible -- their
        # count is reported.
        c, r = cand["coef_blocks"][...][perm], ref["coef_blocks"][...]
        solid = 0.5e30 / float(ref.attrs["re"])
        flips = int(np.sum((c >= solid) != (r >= solid)))
        rel_all = np.abs(c - r) / (np.abs(r) + 1.0)
        graded = ~((c >= solid) | (r >= solid))
        rel = float(np.max(np.where(graded, rel_all, 0.0)))
        outliers = int(np.sum(np.where(graded, rel_all, 0.0) > 1.0e-6))
        print(f"coef_blocks: classification flips = {flips}; graded max"
              f" |dc|/(|c|+1) = {rel:.3e} (tol {args.coef_tol:.1e},"
              f" {outliers} cells > 1e-6 of {c.size})")
        if flips > 0 or rel > args.coef_tol:
            ok = False

        if "dwall_blocks" in ref and "dwall_blocks" in cand:
            # Interior cells only: at ghost coordinates beyond a periodic
            # boundary prepare stores the PERIODIC (minimum-image)
            # distance -- the solver's own analytic-walldist convention --
            # while mobygeom stores the base-mesh distance. The ghost gap
            # is reported informationally.
            d, s = cand["dwall_blocks"][...][perm], ref["dwall_blocks"][...]
            diff = np.abs(d - s)
            dd = np.max(diff[:, 1:-1, 1:-1, 1:-1])
            ghost = diff.copy()
            ghost[:, 1:-1, 1:-1, 1:-1] = 0.0
            print(f"dwall_blocks: interior max |dd| = {dd:.3e} (tol {args.dwall_tol:.1e});"
                  f" ghost convention gap {np.max(ghost):.3e}")
            if dd > args.dwall_tol:
                ok = False

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
