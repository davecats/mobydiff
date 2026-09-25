#!/usr/bin/env python3
"""Interpolate a turbulentBoundaryLayer restart field onto a finer wall-normal
grid (same nx, nz; more ny), for the Dy+_max=4 case.

  make_finewall_restart.py <source.h5> <new_y_nodes.npy> <output.h5> <dyw_plus>

The BL field is single-block: datasets (1, nz, ny, nx) for un/vn/wn/pn plus the
x/y/z node lines. On the staggered grid u/w/p live at y cell CENTRES and v at
the lower y FACES (nodes[:-1]) -- so u/w/p interpolate center->center and v
face->face (the convention in tools/.../interpolate_channel_restart_dyw.py).
The new y node line comes from the solver's own grid build for the target ny
(a cold start), so the interpolated field lands exactly on the grid the restart
will rebuild from config. Everything else (blocks table, x/z, attrs) is copied;
only ny/block_nb_y/y change.
"""
import sys

import h5py
import numpy as np


def centers(nodes):
    return 0.5 * (nodes[:-1] + nodes[1:])


def interp_y(data, old_y, new_y):
    # data is (nz, ny, nx); interpolate every z-x line along y (linear, clamped).
    moved = np.moveaxis(data, 1, 0)               # (ny, nz, nx)
    flat = moved.reshape(moved.shape[0], -1)
    hi = np.clip(np.searchsorted(old_y, new_y, side="left"), 1, len(old_y) - 1)
    lo = hi - 1
    w = (new_y - old_y[lo]) / (old_y[hi] - old_y[lo])
    out = (1.0 - w[:, None]) * flat[lo, :] + w[:, None] * flat[hi, :]
    return np.moveaxis(out.reshape((len(new_y),) + moved.shape[1:]), 0, 1)


def main():
    source, ynpy, output, dyw = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
    new_y_nodes = np.load(ynpy).astype(np.float64)
    ny_new = len(new_y_nodes) - 1

    with h5py.File(source, "r") as src, h5py.File(output, "w") as dst:
        old_y_nodes = np.array(src["y"], dtype=np.float64)
        oc, nc = centers(old_y_nodes), centers(new_y_nodes)
        ovf, nvf = old_y_nodes[:-1], new_y_nodes[:-1]     # interior v faces

        for k, v in src.attrs.items():
            dst.attrs[k] = v
        dst.attrs["ny"] = np.int32(ny_new)
        dst.attrs["block_nb_y"] = np.int32(ny_new)
        # CRITICAL: the wall-normal grid is rebuilt from the stored grid params
        # on restart, so grid_natural_dyw_plus[1] must match the new y-line's
        # wall clustering -- otherwise the solver builds the OLD grid under the
        # interpolated field (mismatch -> blow-up, and the Peclet dt of the old
        # spacing). Copied blindly it would keep the source's value.
        gd = np.array(src.attrs["grid_natural_dyw_plus"], dtype=np.float64)
        gd[1] = dyw
        dst.attrs["grid_natural_dyw_plus"] = gd

        dst.create_dataset("blocks", data=np.array(src["blocks"]))
        dst.create_dataset("x", data=np.array(src["x"], dtype=np.float64))
        dst.create_dataset("y", data=new_y_nodes)
        dst.create_dataset("z", data=np.array(src["z"], dtype=np.float64))

        for name, yold, ynew in [("un", oc, nc), ("wn", oc, nc), ("pn", oc, nc),
                                 ("vn", ovf, nvf)]:
            d = np.array(src[name], dtype=np.float64)       # (1, nz, ny, nx)
            out = interp_y(d[0], yold, ynew)[None, ...]     # (1, nz, ny_new, nx)
            dst.create_dataset(name, data=out)
            assert dst[name].shape == (1, d.shape[1], ny_new, d.shape[3]), dst[name].shape

    print(f"wrote {output}: ny {len(old_y_nodes)-1} -> {ny_new}")
    print(f"  old dy_wall={old_y_nodes[1]-old_y_nodes[0]:.6e}  new dy_wall={new_y_nodes[1]-new_y_nodes[0]:.6e}")


if __name__ == "__main__":
    main()
