#!/usr/bin/env python3
"""Error of a mobydiff field against the exact 2D Taylor-Green vortex.

The decaying TGV (k = 2*pi/Lx, nu = 1/Re) is an exact incompressible-NS
solution on a 2*pi-periodic square, z-independent with w = 0:

    u = -cos(k x) sin(k y) F(t),   v = sin(k x) cos(k y) F(t),
    p = -1/4 (cos 2kx + cos 2ky) F(t)^2,   F(t) = exp(-2 nu k^2 t).

This reads a solver field (single-block/uniform or block-refined; both stored
as per-block datasets), evaluates the exact u, v at each cell's staggered
coordinate and the field's own time, and reports volume-weighted L2 and Linf
velocity errors. With several uniform resolutions the L2 error should fall as
h^2; with a refinement patch, --error-map localises any interface artifact.

Usage:
  python3 tools/check_tgv.py FIELD.h5 [--error-map out.png]
"""

from __future__ import annotations

import argparse

import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("field")
    parser.add_argument("--error-map", default=None,
                        help="write an x-y |velocity error| map (PNG) to localise artifacts")
    args = parser.parse_args()

    h5 = h5py.File(args.field, "r")
    lx = float(h5.attrs["lx"]); re = float(h5.attrs["re"]); t = float(h5.attrs["t_current"])
    nx = int(h5.attrs["nx"])
    k0 = 2.0*np.pi/lx
    nu = 1.0/re
    F = np.exp(-2.0*nu*k0*k0*t)
    blocks = h5["blocks"][...]
    nbx = int(h5.attrs["block_nb_x"]); nby = int(h5.attrs["block_nb_y"]); nbz = int(h5.attrs["block_nb_z"])
    U = h5["un"]; V = h5["vn"]

    s2 = {"u": 0.0, "v": 0.0}; vol = 0.0; linf = {"u": 0.0, "v": 0.0}
    map_pts = []   # (x, y, |err|) at one z-plane, for the error map
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        h = lx/(nx*2**lev)                      # isotropic cell size at this level
        gi = ox + np.arange(nbx); gj = oy + np.arange(nby)
        xc = (gi + 0.5)*h; yc = (gj + 0.5)*h    # cell centres
        xf = gi*h; yf = gj*h                    # west / south faces (staggered locations)
        # exact, [j, i]; z-independent so broadcast over k
        exu = -np.outer(np.sin(k0*yc), np.cos(k0*xf))*F     # u at (x-face, y-centre)
        exv = np.outer(np.cos(k0*yf), np.sin(k0*xc))*F      # v at (x-centre, y-face)
        eu = U[bid] - exu[None, :, :]           # (k, j, i)
        ev = V[bid] - exv[None, :, :]
        w = h**3
        s2["u"] += np.sum(eu**2)*w; s2["v"] += np.sum(ev**2)*w; vol += eu.size*w
        linf["u"] = max(linf["u"], np.abs(eu).max()); linf["v"] = max(linf["v"], np.abs(ev).max())
        if args.error_map is not None:
            err = np.sqrt(eu[0]**2 + ev[0]**2)  # z = first plane
            XX, YY = np.meshgrid(xc, yc)
            map_pts.append((XX.ravel(), YY.ravel(), err.ravel(), h))
    h5.close()

    l2u = np.sqrt(s2["u"]/vol); l2v = np.sqrt(s2["v"]/vol)
    l2 = np.sqrt((s2["u"] + s2["v"])/vol)
    print(f"{args.field}")
    print(f"  t={t:.4f}  F=exp(-2 nu k^2 t)={F:.6f}  Re={re:.0f}  nx={nx}")
    print(f"  L2:  u={l2u:.3e}  v={l2v:.3e}  vel={l2:.3e}")
    print(f"  Linf u={linf['u']:.3e}  v={linf['v']:.3e}")

    if args.error_map is not None:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 6))
        vmax = max(p[2].max() for p in map_pts) or 1.0
        for XX, YY, err, h in map_pts:
            sc = ax.scatter(XX, YY, c=err, s=(6*h/(lx/nx))**2, cmap="inferno", vmin=0, vmax=vmax)
        fig.colorbar(sc, ax=ax, label="|velocity error|")
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect("equal")
        ax.set_title(f"TGV |velocity error|  L2={l2:.2e}  (marker size ~ cell)")
        fig.tight_layout(); fig.savefig(args.error_map, dpi=150)
        print(f"  wrote {args.error_map}")
    return l2


if __name__ == "__main__":
    main()
