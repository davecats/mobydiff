#!/usr/bin/env python3
"""Map a developed field onto the FLAGEUL-MATCHED grid, and seed the scalars.

    ./make_flageul_ic.py kstat_80000.h5 fmint_1.h5 IC_flageul.h5

The source run and the target differ in EVERY direction -- 160x384x224 uniform
on ly = 3.2 with the fluid at [0.6, 2.6], against 224x224x224 with a stretched
y line on ly = 4.0 with the fluid at [1.0, 3.0] -- so the velocity is
interpolated rather than copied. The fluid GAP is the same size (2.0), so the
map is a shift, y_new = y_old + 0.4, plus trilinear interpolation; x and z are
periodic and span the same lengths.

THE SCALARS ARE SEEDED, NOT INTERPOLATED, and the level is chosen freely:
with an imposed FLUX at both outer faces the problem is pure Neumann, so the
absolute temperature is undetermined. Putting theta = 0 AT THE INTERFACE keeps
the fluid O(15) and the solid O(-100) instead of the O(600) means the Dirichlet
version carried for kappa_s = 0.1 -- the same physics, far better conditioned
for a variance computed as <s^2> - <s>^2.

Solid: constant flux q = 1 gives the exact linear profile
    theta(y) = -(d - |y - y_interface|) / (kappa_s D_f)   (0 at the interface)
Fluid: the Reynolds analogy on the target run's own mean velocity.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "scalar"))
from scalar_tools import BlockGeometry                        # noqa: E402

RE, PR = 149.0, 0.71
SRC = dict(ylo=0.6, yhi=2.6)
DST = dict(ylo=1.0, yhi=3.0)
SWEEP = [("k1", 1.0), ("k2", 1.0), ("k3", 1.0),
         ("a1", 0.1), ("a2", 10.0), ("a3", 100.0)]


def assemble(f, name):
    """Blocks -> one global array (nx, ny, nz).

    The block datasets are stored (nbz, nby, nbx) -- scalar_tools.mesh says so
    -- so each tile is transposed on the way in. Without that, x and z are
    swapped WITHIN each block while the blocks themselves land correctly,
    which yields a field that looks plausible but is tiled.
    """
    geo = BlockGeometry(f)
    nb = geo.nb
    bl = f["blocks"][...]
    nx = int(bl[:, 0].max()) + nb[0]
    ny = int(bl[:, 1].max()) + nb[1]
    nz = int(bl[:, 2].max()) + nb[2]
    g = np.zeros((nx, ny, nz))
    d = f[name]
    for b, (ox, oy, oz, _) in enumerate(bl):
        g[ox:ox + nb[0], oy:oy + nb[1], oz:oz + nb[2]] = d[b].transpose(2, 1, 0)
    return g


def scatter(f, name, g):
    geo = BlockGeometry(f)
    nb = geo.nb
    bl = f["blocks"][...]
    out = np.zeros(f[name].shape)
    for b, (ox, oy, oz, _) in enumerate(bl):
        out[b] = g[ox:ox + nb[0], oy:oy + nb[1], oz:oz + nb[2]].transpose(2, 1, 0)
    f[name][...] = out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("template", help="a 1-step snapshot on the TARGET grid")
    ap.add_argument("out")
    ap.add_argument("--source", type=float, default=0.064631)
    a = ap.parse_args()

    shutil.copyfile(a.template, a.out)
    fs = h5py.File(a.src, "r")
    xs, ys, zs = fs["x"][...], fs["y"][...], fs["z"][...]
    xsc, ysc, zsc = (0.5 * (v[:-1] + v[1:]) for v in (xs, ys, zs))

    with h5py.File(a.out, "r+") as fd:
        xd, yd, zd = fd["x"][...], fd["y"][...], fd["z"][...]
        xdc, ydc, zdc = (0.5 * (v[:-1] + v[1:]) for v in (xd, yd, zd))
        shift = DST["ylo"] - SRC["ylo"]

        # ---- velocity + pressure: trilinear, y shifted into the new gap ----
        for name in ("un", "vn", "wn", "pn"):
            g = assemble(fs, name)
            # separable 1-D interpolation, x and z periodic
            def interp_axis(arr, axis, src, dst, period=None):
                if period is not None:
                    src = np.concatenate([src - period, src, src + period])
                    arr = np.concatenate([arr] * 3, axis=axis)
                return np.apply_along_axis(
                    lambda col: np.interp(dst, src, col), axis, arr)
            g = interp_axis(g, 0, xsc, xdc, period=xs[-1])
            g = interp_axis(g, 2, zsc, zdc, period=zs[-1])
            g = interp_axis(g, 1, ysc + shift, ydc)
            scatter(fd, name, g)
            print(f"   {name}: interpolated  ({g.shape})")

        # ---- the mean velocity of the TARGET field, for the analogy --------
        u = assemble(fd, "un")
        umean = u.mean(axis=(0, 2))

        df = 1.0 / (RE * PR)
        d = DST["ylo"]
        fl = (ydc > DST["ylo"]) & (ydc < DST["yhi"])
        j = np.clip(np.minimum(ydc - DST["ylo"], DST["yhi"] - ydc), None, 0.0)
        depth = np.where(ydc < DST["ylo"], DST["ylo"] - ydc,
                         np.where(ydc > DST["yhi"], ydc - DST["yhi"], 0.0))
        for nm, ks in SWEEP:
            if nm not in fd:
                continue
            th = np.empty(ydc.size)
            th[~fl] = -depth[~fl] / (ks * df)      # exact steady solid profile
            th[fl] = umean[fl] * 1.0               # analogy, theta_tau = 1
            scatter(fd, nm, np.broadcast_to(th[None, :, None],
                                            (xdc.size, ydc.size, zdc.size)))
            print(f"   {nm}: kappa_s {ks:<6g} outer face {th[0]:10.3f}  "
                  f"interface 0.0  centre {th[fl].max():8.3f}")
        fd.attrs["step"] = np.int32(0)
        fd.attrs["t_current"] = 0.0
        fd.attrs["re"] = float(RE)
    print(f"{a.out}: written")


if __name__ == "__main__":
    main()
