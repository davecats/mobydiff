#!/usr/bin/env python3
"""Error of a mobydiff field against the exact 3D Beltrami / ABC flow.

The Beltrami flow (k = 2*pi/Lx, nu = 1/Re) is a fully 3D exact incompressible-NS
solution on a 2*pi-periodic CUBE; curl(u)=k u, so it decays self-similarly:

    u = sin(k z) + cos(k y)        v = sin(k x) + cos(k z)
    w = sin(k y) + cos(k x)        u(t) = u(0) F(t),  F(t) = exp(-nu k^2 t)
    p = -1/2 |u|^2

Each component depends on two coordinates and is independent of its OWN
staggered direction, so unlike the 2D Taylor-Green vortex every direction is
exercised (in particular z and w). Reads a solver field (uniform or
block-refined), evaluates the exact u, v, w at each cell's staggered coordinate
and the field's own time, and reports volume-weighted L2 / Linf velocity errors.
A uniform-resolution sweep should give L2 ~ h^2.

Usage:
  python3 tools/check_beltrami.py FIELD.h5
"""

from __future__ import annotations

import argparse

import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("field")
    args = parser.parse_args()

    h5 = h5py.File(args.field, "r")
    lx = float(h5.attrs["lx"]); re = float(h5.attrs["re"]); t = float(h5.attrs["t_current"])
    nx = int(h5.attrs["nx"]); ny = int(h5.attrs["ny"]); nz = int(h5.attrs["nz"])
    k0 = 2.0*np.pi/lx
    nu = 1.0/re
    F = np.exp(-nu*k0*k0*t)
    blocks = h5["blocks"][...]
    nbx = int(h5.attrs["block_nb_x"]); nby = int(h5.attrs["block_nb_y"]); nbz = int(h5.attrs["block_nb_z"])
    # Per-direction refinement mask (absent = xyz octree): an unrefined
    # direction keeps the level-0 spacing at every level (refine_dims = xz).
    mask = np.asarray(h5.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
    U = h5["un"]; V = h5["vn"]; W = h5["wn"]

    s2 = {"u": 0.0, "v": 0.0, "w": 0.0}; vol = 0.0
    linf = {"u": 0.0, "v": 0.0, "w": 0.0}
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        hx = lx/(nx*2**(lev*mask[0]))
        hy = lx/(ny*2**(lev*mask[1]))
        hz = lx/(nz*2**(lev*mask[2]))
        gi = ox + np.arange(nbx); gj = oy + np.arange(nby); gk = oz + np.arange(nbz)
        xc = (gi + 0.5)*hx; yc = (gj + 0.5)*hy; zc = (gk + 0.5)*hz   # cell centres
        # Each component is independent of its own staggered direction, so only
        # the cell-centred coordinates of the two relevant directions enter.
        sx, cx = np.sin(k0*xc), np.cos(k0*xc)   # (i,)
        sy, cy = np.sin(k0*yc), np.cos(k0*yc)   # (j,)
        sz, cz = np.sin(k0*zc), np.cos(k0*zc)   # (k,)
        exu = F*(sz[:, None, None] + cy[None, :, None])   # u = sin(kz)+cos(ky), broadcast over i
        exv = F*(cz[:, None, None] + sx[None, None, :])   # v = sin(kx)+cos(kz), broadcast over j
        exw = F*(sy[None, :, None] + cx[None, None, :])   # w = sin(ky)+cos(kx), broadcast over k
        eu = U[bid] - exu; ev = V[bid] - exv; ew = W[bid] - exw
        w = hx*hy*hz
        s2["u"] += np.sum(eu**2)*w; s2["v"] += np.sum(ev**2)*w; s2["w"] += np.sum(ew**2)*w
        vol += eu.size*w
        linf["u"] = max(linf["u"], np.abs(eu).max())
        linf["v"] = max(linf["v"], np.abs(ev).max())
        linf["w"] = max(linf["w"], np.abs(ew).max())
    h5.close()

    l2 = {c: np.sqrt(s2[c]/vol) for c in s2}
    l2vel = np.sqrt(sum(s2.values())/vol)
    print(f"{args.field}")
    print(f"  t={t:.4f}  F=exp(-nu k^2 t)={F:.6f}  Re={re:.0f}  nx={nx}")
    print(f"  L2:  u={l2['u']:.3e}  v={l2['v']:.3e}  w={l2['w']:.3e}  vel={l2vel:.3e}")
    print(f"  Linf u={linf['u']:.3e}  v={linf['v']:.3e}  w={linf['w']:.3e}")
    return l2vel


if __name__ == "__main__":
    main()
