#!/usr/bin/env python3
"""Reduce Neuhauser's NekRS pipe DNS to the radial profiles this comparison
needs, and cache them as `neuhauser_profiles.npz`.

WHY THE CACHE EXISTS. The published dataset (DOI 10.35097/26za20q32xsz43yk)
is an 11.8 GB bag holding a 6.8 GB zarr store of (time, r, phi) moments for
every case, Prandtl number and boundary condition. Everything the comparison
uses is a handful of 226-point RADIAL profiles -- 1 MB once reduced -- so the
tutorial ships the reduction, not the archive, and stays runnable without the
download. Re-running this script against the archive REPRODUCES the cache;
that is the point of keeping it.

    # one-off, only if you want to rebuild the cache from the archive
    tar xf '10.35097-26za20q32xsz43yk(1).tar'
    tar xf .../data/dataset/cht_short/joined_datasets.interp.zarr.tar -C /somewhere
    NEUHAUSER_ZARR=/somewhere/joined_datasets.interp.zarr ./extract_neuhauser.py

The reduction is exactly what `compare_neuhauser.polar_moments` + `batch` do
(the cartesian-to-polar rotation, second moments formed BEFORE the azimuthal
average, then the batch mean and standard error over the `time` batch axis),
plus the <theta'^2> budget terms `variance_budget.py` needs for case c0.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_neuhauser import CASE_MAP, ZARR, batch, open_case, polar_moments

PR = 0.71
NU = 1.0 / 5300.0
ALPHA = NU / PR
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "neuhauser_profiles.npz")


def budget(d):
    """The <theta'^2> budget terms, azimuthally and batch averaged.

    Kept separate from `polar_moments` because it needs the TRIPLE products
    and the gradient-squared moments, which only the reference carries and
    only case c0 is compared on.
    """
    r = d["r"].values
    c, s = np.cos(d["phi"].values), np.sin(d["phi"].values)
    A = lambda k: np.asarray(d[k])                       # (time, r, phi)
    t, tt = A("t"), A("t*t")
    ur = A("u") * c + A("v") * s
    urt = A("t*u") * c + A("t*v") * s
    urtt = A("t*t*u") * c + A("t*t*v") * s
    var = np.mean(tt - t ** 2, axis=(0, 2))
    flux = np.mean(urt - ur * t, axis=(0, 2))
    trip = np.mean(urtt - 2 * t * urt - ur * tt + 2 * ur * t ** 2, axis=(0, 2))
    grad = np.mean(A("d(t)dx*d(t)dx") + A("d(t)dy*d(t)dy") + A("d(t)dz*d(t)dz")
                   - A("d(t)dx") ** 2 - A("d(t)dy") ** 2 - A("d(t)dz") ** 2,
                   axis=(0, 2))
    tmean = np.mean(t, axis=(0, 2))
    return dict(
        r=r,
        P=-2 * flux * np.gradient(tmean, r),
        eps=2 * ALPHA * grad,
        Tt=-np.gradient(r * trip, r) / np.maximum(r, 1e-12),
        Tn=ALPHA * np.gradient(r * np.gradient(var, r), r) / np.maximum(r, 1e-12),
        var=var, flux=flux, trip=trip,
        # the Kasagi source is proportional to u_z, so it FLUCTUATES and feeds
        # the variance directly: 2 S <u_z' t'> with S = 1/u_b.
        tuz=np.mean(A("t*w") - t * A("w"), axis=(0, 2)),
        uz=np.mean(A("w"), axis=(0, 2)),
    )


def main():
    print(f"reading {ZARR}")
    out = {}
    for nm in CASE_MAP:
        d = open_case(nm, PR)
        p = batch(polar_moments(d, thermal=True))
        for k, v in p.items():
            out[f"{nm}/{k}"] = np.asarray(v, dtype=np.float64)
        print(f"  {nm:5s} {len(p['r'])} radii")
        if nm == "c0":
            for k, v in budget(d).items():
                out[f"budget/{k}"] = np.asarray(v, dtype=np.float64)
            print("        + the <theta'^2> budget terms")
    np.savez_compressed(OUT, **out)
    print(f"wrote {OUT}  ({os.path.getsize(OUT) / 1e6:.2f} MB,"
          f" {len(out)} arrays)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
