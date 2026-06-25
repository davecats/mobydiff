#!/usr/bin/env python3
"""Coarse-side interface gate for a mobydiff Beltrami slab field (block-table).

The block-averaged roughness in interface_diagnostics.py dilutes a single-cell-row
artifact over the whole block, so it cannot see the spurious band that sits on the
ONE coarse interior cell row touching a 2:1 interface. This tool reassembles the
LEVEL-0 (coarse) blocks onto the coarse global grid and measures, PER coarse y-row,
the tangential (x-z) discrete-Laplacian roughness of the ERROR vs the exact Beltrami
solution (a high-k / checkerboard proxy). A clean interface keeps the interface-
adjacent coarse rows ~ the coarse interior rows; the band shows up as a localized
spike at exactly the rows bordering the refined band.

Assumes a y-normal band (the channel-relevant geometry). Data axes are (k,j,i)=(z,y,x).

Usage: python3 tools/interface_coarse_gate.py FIELD.h5 [--var u|v|w|p]
"""
from __future__ import annotations
import argparse, h5py, numpy as np


def exact_coarse(xc, yc, zc, F):
    """Exact Beltrami at coarse cell centres (broadcast to (z,y,x))."""
    sx, cx = np.sin(xc), np.cos(xc)
    sy, cy = np.sin(yc), np.cos(yc)
    sz, cz = np.sin(zc), np.cos(zc)
    nz, ny, nx = zc.size, yc.size, xc.size
    one = np.ones((nz, ny, nx))
    ex = {}
    ex["u"] = F * (sz[:, None, None] + cy[None, :, None]) * one
    ex["v"] = F * (cz[:, None, None] + sx[None, None, :]) * one
    ex["w"] = F * (sy[None, :, None] + cx[None, None, :]) * one
    ex["p"] = np.zeros((nz, ny, nx))  # not checked
    return ex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("field")
    ap.add_argument("--vars", default="uvw")
    a = ap.parse_args()
    h5 = h5py.File(a.field, "r")
    lx = float(h5.attrs["lx"]); re = float(h5.attrs["re"]); t = float(h5.attrs["t_current"])
    nx = int(h5.attrs["nx"]); ny = int(h5.attrs["ny"]); nz = int(h5.attrs["nz"])
    nb = int(h5.attrs["block_nb_x"])
    nu = 1.0 / re; k0 = 2.0 * np.pi / lx
    F = np.exp(-nu * k0 * k0 * t)
    blocks = h5["blocks"][...]
    data = {v: h5[{"u": "un", "v": "vn", "w": "wn", "p": "pn"}[v]][...] for v in a.vars}
    h5.close()

    h = lx / nx  # coarse cell size
    # Dense coarse global grid (z,y,x); NaN where no level-0 block (refined band).
    glob = {v: np.full((nz, ny, nx), np.nan) for v in a.vars}
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        if lev != 0:
            continue
        for v in a.vars:
            glob[v][oz:oz+nb, oy:oy+nb, ox:ox+nb] = data[v][bid]

    # Coarse cell-centre coords.
    xc = (np.arange(nx) + 0.5) * h * k0
    yc = (np.arange(ny) + 0.5) * h * k0
    zc = (np.arange(nz) + 0.5) * h * k0
    ex = exact_coarse(xc, yc, zc, F)

    # Which coarse y-rows are present (not all-NaN)?
    present = np.array([not np.all(np.isnan(glob[a.vars[0]][:, j, :])) for j in range(ny)])

    print(f"{a.field}")
    print(f"  Re={re:.0f} nx={nx} t={t:.4f} F={F:.4f}")
    print(f"  coarse y-rows present: {np.where(present)[0].tolist()}")
    print(f"  {'j':>3} {'y':>6}  " + "  ".join(f"{v}:roughRMS" for v in a.vars))
    for j in range(ny):
        if not present[j]:
            continue
        cells = []
        for v in a.vars:
            e = glob[v][:, j, :] - ex[v][:, j, :]   # (z,x) slice of the error
            # tangential discrete Laplacian (periodic in x and z)
            lap = (np.roll(e, 1, 0) + np.roll(e, -1, 0)
                   + np.roll(e, 1, 1) + np.roll(e, -1, 1) - 4 * e)
            cells.append(np.sqrt(np.nanmean(lap ** 2)))
        ystar = (j + 0.5) * h
        flags = "  <-- interface" if (j + 1 < ny and not present[j+1]) or (j - 1 >= 0 and not present[j-1]) else ""
        print(f"  {j:3d} {ystar:6.3f}  " + "  ".join(f"{c:.3e}" for c in cells) + flags)


if __name__ == "__main__":
    main()
