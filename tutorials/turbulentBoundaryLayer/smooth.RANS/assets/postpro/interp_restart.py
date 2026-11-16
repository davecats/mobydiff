#!/usr/bin/env python3
"""Interpolate a converged smooth.RANS snapshot onto another grid of the same
domain, as the initial condition of a finer (or coarser) run:

    python3 interp_restart.py SOURCE.h5 TARGET.case.h5 OUT.h5

TARGET.case.h5 is the moby_prepare case file of the target run: its node lines
and block table define the output; the run attributes (time, step, boundary
conditions, run control) are the SOURCE snapshot's, with the grid and block
attributes replaced by the target's and dt cut to DT_START so the time-step
controller starts below the finer grid's diffusive limit. Every
field is z-averaged (the RANS is spanwise uniform), interpolated bilinearly
in (x, y) at its own staggered location (un on the low x faces, vn on the low
y faces, the rest at cell centres; nearest value beyond the source points)
and broadcast over z. The outlet-plane state (un_xmax, vn_ymax) is DELETED --
an IC must write the plane or drop it (CLAUDE.md landmines), and the solver
rebuilds it from the field. The projection removes the O(h^2) divergence of
the interpolated velocity in the first steps; nut is recomputed from k and
omega at the first assembly.
"""
import sys

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator

DT_START = 1.0e-3
FIELDS = ("un", "vn", "wn", "pn", "k", "omega", "nut", "gamma", "rethetat")


def read2d(f, name):
    """z-averaged (ny, nx) array of a single-level block-layout dataset."""
    blocks = f["blocks"][...]
    nb = [int(f.attrs[f"block_nb_{d}"]) for d in "xyz"]
    ny, nx = int(f.attrs["ny"]), int(f.attrs["nx"])
    out = np.empty((ny, nx))
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        assert lev == 0, "single-level layouts only"
        out[oy:oy + nb[1], ox:ox + nb[0]] = f[name][bid].mean(axis=0)
    return out


def coords(x, y, name):
    """(y, x) coordinates of a dataset's points: low faces for un / vn."""
    xc, yc = 0.5 * (x[1:] + x[:-1]), 0.5 * (y[1:] + y[:-1])
    return (yc if name != "vn" else y[:-1]), (xc if name != "un" else x[:-1])


def main():
    src, case, out = sys.argv[1:4]
    with h5py.File(src, "r") as fs, h5py.File(case, "r") as fc, h5py.File(out, "w") as fo:
        for key, val in fs.attrs.items():
            fo.attrs[key] = val
        blocks = fc["blocks"][...]
        nb = [int(v) for v in fc.attrs["block_nb_xyz"]]
        for d, n in zip("xyz", nb):
            fo.attrs[f"block_nb_{d}"] = np.int32(n)
        for key in ("nx", "ny", "nz", "lx", "ly", "lz", "re", "periodic",
                    "grid_distribution", "grid_stretch", "grid_natural_dyw_plus"):
            fo.attrs[key] = fc.attrs[key]
        fo.attrs["n_blocks"] = np.int32(len(blocks))
        fo.attrs["dt"] = DT_START
        fo.create_dataset("blocks", data=blocks)
        xt, yt, zt = (fc[f"{d}_nodes"][...] for d in "xyz")
        for d, line in zip("xyz", (xt, yt, zt)):
            fo.create_dataset(d, data=line)
        nz = int(fc.attrs["nz"])
        xs_n, ys_n = fs["x"][...], fs["y"][...]
        for name in FIELDS:
            if name not in fs:
                continue
            ys, xs = coords(xs_n, ys_n, name)
            interp = RegularGridInterpolator((ys, xs), read2d(fs, name),
                                             bounds_error=False, fill_value=None)
            yq, xq = coords(xt, yt, name)
            yq = np.clip(yq, ys[0], ys[-1])          # nearest beyond the source points
            xq = np.clip(xq, xs[0], xs[-1])
            Y, X = np.meshgrid(yq, xq, indexing="ij")
            full = interp(np.stack([Y, X], axis=-1))   # (ny, nx)
            data = np.empty((len(blocks), nb[2], nb[1], nb[0]))
            for bid, (ox, oy, oz, lev) in enumerate(blocks):
                data[bid] = full[None, oy:oy + nb[1], ox:ox + nb[0]]
            fo.create_dataset(name, data=data)
        assert nz % nb[2] == 0
        print(f"wrote {out}: {', '.join(n for n in FIELDS if n in fo)} from {src} "
              f"(step {int(fs.attrs['step'])}, t {float(fs.attrs['t_current']):.1f}) "
              f"onto {int(fc.attrs['nx'])} x {int(fc.attrs['ny'])} x {nz}; outlet planes dropped")


if __name__ == "__main__":
    main()
