#!/usr/bin/env python3
"""Independent level-0 touch classification of the egg-carton case file.

    check_eggclass.py <eggclass.case.h5> <default-wall twin.case.h5>

A level-0 block is TOUCHED when its one-cell-dilated box holds cell centres
on both sides of the surface y = offset + amp sin(kx x + px) sin(kz z + pz)
(periodic in x and z, the wall at y = 0). The solver's block_touch_l0 must
equal this raster exactly, and must DIFFER from the default-wall twin's --
the twin is what the pre-fix moby_prepare produced for ANY analytic wall.
Parameters are read from the input echo in the case file (case_inputs).
"""
import re
import sys

import h5py
import numpy as np


def echo_values(echo, key):
    """`key = v1 v2 ...` of the case-file input echo, as floats."""
    m = re.search(rf"^\s*{key}\s*=\s*(.+)$", echo, re.M)
    return [float(v) for v in m.group(1).split()]


def main():
    egg, twin = sys.argv[1], sys.argv[2]
    with h5py.File(egg, "r") as f:
        touch = f["block_touch_l0"][...].reshape(-1)
        echo = f.attrs["case_inputs"]
        echo = echo.decode() if isinstance(echo, bytes) else str(echo)
        nx, ny, nz = (int(f.attrs[k]) for k in ("nx", "ny", "nz"))
        lx, ly, lz = (float(f.attrs[k]) for k in ("lx", "ly", "lz"))
        nb = [int(v) for v in np.atleast_1d(f.attrs["block_nb"])]
        nb = nb*3 if len(nb) == 1 else nb
    with h5py.File(twin, "r") as f:
        touch_twin = f["block_touch_l0"][...].reshape(-1)
    amp = echo_values(echo, "wall_amp")[0]; off = echo_values(echo, "wall_offset")[0]
    nwx, nwz = echo_values(echo, "wall_n_wave"); px, pz = echo_values(echo, "wall_phase")

    xc = (np.arange(nx) + 0.5)*lx/nx; yc = (np.arange(ny) + 0.5)*ly/ny; zc = (np.arange(nz) + 0.5)*lz/nz
    X, Z = np.meshgrid(xc, zc, indexing="ij")
    h = off + amp*np.sin(2*np.pi*nwx*X/lx + px)*np.sin(2*np.pi*nwz*Z/lz + pz)
    solid = yc[None, :, None] < h[:, None, :]
    bx, by, bz = nx//nb[0], ny//nb[1], nz//nb[2]
    mine = np.zeros((bx, by, bz), bool)
    for bi in range(bx):
        for bj in range(by):
            for bk in range(bz):
                ii = np.arange(bi*nb[0] - 1, (bi + 1)*nb[0] + 1) % nx
                jj = np.arange(max(bj*nb[1] - 1, 0), min((bj + 1)*nb[1] + 1, ny))
                kk = np.arange(bk*nb[2] - 1, (bk + 1)*nb[2] + 1) % nz
                s = solid[np.ix_(ii, jj, kk)]
                mine[bi, bj, bk] = s.any() and (~s).any()
    mine = mine.transpose(2, 1, 0).reshape(-1)          # x-fastest raster
    same = np.array_equal(mine.astype(touch.dtype), touch)
    differs = not np.array_equal(touch, touch_twin)
    print(f"touched level-0 blocks: solver {int(touch.sum())}, independent {int(mine.sum())}, "
          f"default-wall twin {int(touch_twin.sum())}")
    print(f"solver == independent: {same};  solver != twin: {differs}")
    sys.exit(0 if (same and differs) else 1)


if __name__ == "__main__":
    main()
