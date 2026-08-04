#!/usr/bin/env python3
"""Seed a developed channel snapshot with the Reynolds-analogy theta profile.

The S2 turbulent-channel gates restart from the ../channel_interface/les
campaign's DEVELOPED velocity fields, so only the scalar has to spin up.
From the linear (conduction) profile that takes ~3 relaxation times of the
resolved turbulent mixing -- pure wall-clock spent on a transient nobody
measures.

This writes instead the profile the REYNOLDS ANALOGY implies from the
snapshot's OWN mean velocity: with antisymmetric isothermal walls and no
source, theta and the streamwise velocity obey the same constant-flux
transport equation up to the molecular Prandtl number, so

    theta(y) = theta_w0 - (theta_w0 - theta_w1)/2 * <u>(y)/<u>(centre)

on the lower half, antisymmetric on the upper.

This IS NOT circular with the gate: the profile comes from the run's own
velocity field, not from Kader's correlation (the external reference the gate
measures theta+ against), and it is deliberately WRONG near the wall by
exactly the Pr-dependent offset (beta(0.71) = 3.83 vs the velocity's B = 5.5)
that the gate's sublayer and log-intercept checks resolve. The near-wall
correction relaxes on the fast diffusive time y^2/D; what the IC buys is the
SLOW core adjustment.

  make_theta_ic.py in.h5 out.h5 [--scalar theta] [--walls 1 -1]
                                [--zero-pressure]

Multi-level (2:1-refined) files are handled: every leaf row belongs to exactly
one level, so the mean-velocity profile is the union of the per-level rows
sorted by y, and each leaf is filled by interpolating it at its own cell
centres.
"""

from __future__ import annotations

import argparse
import shutil
import sys

import h5py
import numpy as np

from scalar_tools import BlockGeometry


def row_profile(f, geo, name):
    """(y, <q>) of every stored cell row, x/z-averaged, sorted by y."""
    dset = f[name]
    acc = {}
    for bid in range(geo.n_blocks):
        (_, _), (yc, _), (_, _) = geo.block_axes(bid)
        block = dset[bid]
        n = block.shape[0] * block.shape[2]
        for jj, yv in enumerate(yc):
            key = round(float(yv), 12)
            tot, cnt = acc.get(key, (0.0, 0))
            acc[key] = (tot + float(np.mean(block[:, jj, :])) * n, cnt + n)
    ys = np.array(sorted(acc))
    return ys, np.array([acc[y][0] / acc[y][1] for y in ys])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--scalar", default="theta")
    ap.add_argument("--walls", type=float, nargs=2, default=[1.0, -1.0])
    ap.add_argument("--zero-pressure", action="store_true")
    a = ap.parse_args()

    if a.infile != a.outfile:
        shutil.copyfile(a.infile, a.outfile)

    with h5py.File(a.outfile, "r+") as f:
        geo = BlockGeometry(f)
        ly = float(f.attrs["ly"])
        y, u = row_profile(f, geo, "un")
        # no-slip wall end points
        yy = np.concatenate([[0.0], y, [ly]])
        uu = np.concatenate([[0.0], u, [0.0]])
        uc = uu.max()
        half = 0.5 * (a.walls[0] - a.walls[1])
        theta = np.where(yy <= 0.5 * ly,
                         a.walls[0] - half * uu / uc,
                         a.walls[1] + half * uu / uc)

        if a.scalar not in f:
            src = f["pn"]
            f.create_dataset(a.scalar, shape=src.shape, dtype=src.dtype)
        dset = f[a.scalar]
        for bid in range(geo.n_blocks):
            _, ymesh, _, _ = geo.mesh(bid)
            dset[bid] = np.broadcast_to(np.interp(ymesh, yy, theta),
                                        dset.shape[1:]).astype(dset.dtype)

        if a.zero_pressure:
            f["pn"][...] = 0.0

    print(f"wrote {a.outfile}: {a.scalar}(y) = Reynolds analogy of the snapshot's "
          f"own <u>(y) (U_c = {uc:.4f}), walls {a.walls[0]:+g} / {a.walls[1]:+g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
