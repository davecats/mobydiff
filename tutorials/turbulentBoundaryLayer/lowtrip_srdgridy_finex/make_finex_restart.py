#!/usr/bin/env python3
"""Interpolate a boundaryLayer restart field onto a finer STREAMWISE grid (more
nx; same ny, nz), for the fine-x case.

  make_finex_restart.py <source.h5> <new_x_nodes.npy> <output.h5>

Single-block field (1, nz, ny, nx). On the staggered grid the streamwise
velocity u lives at x FACES and v/w/p at x CENTRES (mirror of the wall-normal
case, where v was the face-staggered one). So u interpolates face->face
(nodes[:-1]) and v/w/p centre->centre, along the last (x) axis. The x grid is
GEOMETRIC, rebuilt on restart from nx + grid_stretch[0] alone (no dyw-style
attr), so only nx / block_nb_x / the x dataset change -- no attr to resync.
"""
import sys
import h5py
import numpy as np


def centers(n):
    return 0.5 * (n[:-1] + n[1:])


def interp_x(data, ox, nx_):
    # data (nz, ny, nx) -> interpolate along the last axis onto nx_
    hi = np.clip(np.searchsorted(ox, nx_, side="left"), 1, len(ox) - 1)
    lo = hi - 1
    w = (nx_ - ox[lo]) / (ox[hi] - ox[lo])
    return (1.0 - w) * data[..., lo] + w * data[..., hi]


def main():
    source, xnpy, output = sys.argv[1], sys.argv[2], sys.argv[3]
    nxn = np.load(xnpy).astype(np.float64)
    nx_new = len(nxn) - 1

    with h5py.File(source, "r") as src, h5py.File(output, "w") as dst:
        oxn = np.array(src["x"], dtype=np.float64)
        oc, nc = centers(oxn), centers(nxn)
        oxf, nxf = oxn[:-1], nxn[:-1]              # interior x faces (u DOFs)

        for k, v in src.attrs.items():
            dst.attrs[k] = v
        dst.attrs["nx"] = np.int32(nx_new)
        dst.attrs["block_nb_x"] = np.int32(nx_new)

        dst.create_dataset("blocks", data=np.array(src["blocks"]))
        dst.create_dataset("x", data=nxn)
        dst.create_dataset("y", data=np.array(src["y"], dtype=np.float64))
        dst.create_dataset("z", data=np.array(src["z"], dtype=np.float64))

        for name, ox, nxx in [("un", oxf, nxf), ("vn", oc, nc),
                              ("wn", oc, nc), ("pn", oc, nc)]:
            d = np.array(src[name], dtype=np.float64)      # (1, nz, ny, nx)
            out = interp_x(d[0], ox, nxx)[None, ...]       # (1, nz, ny, nx_new)
            dst.create_dataset(name, data=out)
            assert dst[name].shape == (1, d.shape[1], d.shape[2], nx_new), dst[name].shape

    print(f"wrote {output}: nx {len(oxn)-1} -> {nx_new}")
    print(f"  old dx0={oxn[1]-oxn[0]:.6e}  new dx0={nxn[1]-nxn[0]:.6e}")


if __name__ == "__main__":
    main()
