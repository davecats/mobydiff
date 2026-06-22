#!/usr/bin/env python3
"""Error of a mobydiff field against the exact Poiseuille + blowing/suction flow.

A channel (walls in one direction) with a uniform wall-normal crossflow V (blown
in at the min wall, sucked out at the max wall, so global mass is conserved) and
a body force G in a tangential direction is an exact STEADY incompressible-NS
solution:

    nu u'' - V u' = -G,   u(y) = (G/V)[ y - L (e^{V y/nu}-1)/(e^{V L/nu}-1) ],
    v = V (wall-normal),  third component = 0,   y = wall-normal coordinate.

It exercises momentum convection (V du/dy) and diffusion (nu u'') AND the
pressure projection (sustaining the crossflow divergence-free). The wall
direction, flow direction, V and G are all read from the field metadata, so the
check works for walls in x, y or z. A wall-normal resolution sweep should give
L2 ~ h^2.

Usage:
  python3 tools/check_poiseuille.py FIELD.h5
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
    re = float(h5.attrs["re"]); nu = 1.0/re
    leng = np.array([h5.attrs["lx"], h5.attrs["ly"], h5.attrs["lz"]], float)
    nx = int(h5.attrs["nx"]); t = float(h5.attrs["t_current"])
    periodic = np.array(h5.attrs["periodic"]).astype(int)
    forcing = np.array([h5.attrs["forcing_x"], h5.attrs["forcing_y"], h5.attrs["forcing_z"]], float)
    bcv = np.array(h5.attrs["bc_value"], float)

    dwall = int(np.argwhere(periodic == 0)[0][0])     # 0-based wall direction
    dflow = int(np.argmax(np.abs(forcing)))           # 0-based flow direction
    G = forcing[dflow]
    V = bcv[np.argmax(np.abs(bcv))]                    # the (only) non-zero wall velocity
    L = leng[dwall]
    if V == 0.0:
        raise SystemExit("no crossflow (V=0): not a blowing/suction case")
    eVL = np.exp(V*L/nu)

    def profile(yc):
        return (G/V)*(yc - L*(np.exp(V*yc/nu) - 1.0)/(eVL - 1.0))

    blocks = h5["blocks"][...]
    nb = [int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]), int(h5.attrs["block_nb_z"])]
    comps = [("un", h5["un"]), ("vn", h5["vn"]), ("wn", h5["wn"])]
    s2 = [0.0, 0.0, 0.0]; vol = 0.0; linf = [0.0, 0.0, 0.0]
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        h = leng[0]/(nx*2**lev)   # isotropic cubic cell at this level (lx/nx == ly/ny == lz/nz)
        org = [ox, oy, oz]
        # cell-centred wall-normal coordinate, broadcast over the field axes (k,j,i)
        wc = (org[dwall] + np.arange(nb[dwall]) + 0.5)*h
        prof = profile(wc)
        # axis of the wall-normal direction in a (k,j,i)=(z,y,x) array
        kji_axis = {0: 2, 1: 1, 2: 0}[dwall]
        shp = [1, 1, 1]; shp[kji_axis] = nb[dwall]
        prof_b = prof.reshape(shp)
        w3 = h**3
        for c, (name, dset) in enumerate(comps):
            if c == dflow:
                ex = np.broadcast_to(prof_b, dset[bid].shape)
            elif c == dwall:
                ex = V
            else:
                ex = 0.0
            e = dset[bid] - ex
            s2[c] += np.sum(e**2)*w3
            linf[c] = max(linf[c], np.abs(e).max())
        vol += comps[0][1][bid].size*w3
    h5.close()

    l2 = [np.sqrt(x/vol) for x in s2]
    l2vel = np.sqrt(sum(s2)/vol)
    wall = "xyz"[dwall]; flow = "xyz"[dflow]
    print(f"{args.field}")
    print(f"  walls in {wall}, flow in {flow}, V={V:.4g}, G={G:.4g}, Re={re:.0f}, "
          f"Re_cross=VL/nu={V*L/nu:.3g}, t={t:.3f}")
    print(f"  L2:  u={l2[0]:.3e}  v={l2[1]:.3e}  w={l2[2]:.3e}  vel={l2vel:.3e}")
    print(f"  Linf u={linf[0]:.3e}  v={linf[1]:.3e}  w={linf[2]:.3e}")
    return l2vel


if __name__ == "__main__":
    main()
