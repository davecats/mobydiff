#!/usr/bin/env python3
"""Re-seed the MEAN scalar profile of a developed run from its own transport.

The thermal field of the S2 turbulent-channel gates relaxes on the slowest
mode of the effective wall-normal diffusivity, tau = L^2/(pi^2 D_eff) ~ 6
t.u. at Re_tau 180; reaching a stationary mean from any IC therefore costs
~3 tau of wall clock in which nothing is measured.

This shortcuts it WITHOUT touching the physics: from the time-averaged
profiles of an already-running case it measures the TOTAL wall-normal
diffusivity the run itself realises,

    D_tot(y) = J(y) / (-d<theta>/dy),
    J(y) = <v'theta'> - (D + nu_t/Pr_t) d<theta>/dy     (all measured)

solves the constant-flux profile that diffusivity implies,

    d(theta_target)/dy = -J_new/D_tot(y),   J_new fixed by the wall values,

and writes a restart in which every cell's FLUCTUATION is kept exactly and
only the (x,z)-mean is replaced:

    theta_new(x,y,z) = theta(x,y,z) - <theta>(y) + theta_target(y).

The developed fluctuating field survives, so only the residual mean
adjustment remains -- and whether the result is genuinely stationary is then
a MEASUREMENT (check_scalar_turb.py's constant-J check plus a first-half /
second-half comparison of the averaging window), not an assumption.

  retarget_theta.py field.h5 out.h5 --stats 'case_*.h5' [--skip N]
                    [--scalar theta] [--pr] [--prt] [--re] [--walls 1 -1]
"""

from __future__ import annotations

import argparse
import shutil
import sys

import h5py
import numpy as np

from scalar_tools import BlockGeometry
from check_scalar_turb import sorted_snapshots, snapshot_stats, level_bands


def merged_face_profile(paths, levels, scalar, pr, prt, re, walls):
    """Time-averaged (y_face, J, d<theta>/dy) merged over levels."""
    face = {}
    for lev, band in levels:
        acc = None
        for p in paths:
            s = snapshot_stats(p, lev, scalar, pr, prt, re, walls, band)
            acc = s if acc is None else {k: (acc[k] + s[k] if k in ("flux", "theta")
                                             else acc[k]) for k in acc}
        n = len(paths)
        thm = acc["theta"] / n
        flux = acc["flux"] / n
        yc, yf = acc["y"], acc["yface"]
        for j in range(1, yc.size):
            if not np.isfinite(flux[j]):
                continue
            grad = (thm[j] - thm[j - 1]) / (yc[j] - yc[j - 1])
            face[round(float(yf[j]), 12)] = (flux[j], grad)
    ys = np.array(sorted(face))
    return ys, np.array([face[y][0] for y in ys]), np.array([face[y][1] for y in ys])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--stats", nargs="+", required=True)
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--scalar", default="theta")
    ap.add_argument("--pr", type=float, default=0.71)
    ap.add_argument("--prt", type=float, default=0.85)
    ap.add_argument("--re", type=float, default=180.0)
    ap.add_argument("--walls", type=float, nargs=2, default=[1.0, -1.0])
    a = ap.parse_args()

    paths = sorted_snapshots(a.stats, a.skip)
    with h5py.File(paths[0], "r") as f:
        tab = f["blocks"][...]
        nb = (int(f.attrs["block_nb_x"]), int(f.attrs["block_nb_y"]),
              int(f.attrs["block_nb_z"]))
        levels = [(lev, b) for lev in sorted(set(int(v) for v in tab[:, 3]))
                  for b in range(len(level_bands(tab, nb, lev)))]
        ly = float(f.attrs["ly"])

    yf, flux, grad = merged_face_profile(paths, levels, a.scalar, a.pr, a.prt,
                                         a.re, a.walls)
    # D_tot from the measured flux and gradient; the walls are pure molecular.
    dmol = 1.0 / (a.re * a.pr)
    dtot = np.where(np.abs(grad) > 1e-12, flux / np.maximum(-grad, 1e-30), dmol)
    dtot = np.maximum(dtot, dmol)
    yy = np.concatenate([[0.0], yf, [ly]])
    dd = np.concatenate([[dmol], dtot, [dmol]])
    g = 1.0 / dd
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(yy))])
    jnew = (a.walls[0] - a.walls[1]) / integ[-1]
    target = a.walls[0] - jnew * integ

    if a.infile != a.outfile:
        shutil.copyfile(a.infile, a.outfile)
    with h5py.File(a.outfile, "r+") as f:
        geo = BlockGeometry(f)
        dset = f[a.scalar]
        data = dset[...]
        # per-row (x,z) mean of the snapshot itself
        rows = {}
        for bid in range(geo.n_blocks):
            (_, _), (yc, _), (_, _) = geo.block_axes(bid)
            for jj, yv in enumerate(yc):
                key = round(float(yv), 12)
                s, n = rows.get(key, (0.0, 0))
                plane = data[bid][:, jj, :]
                rows[key] = (s + float(plane.sum()), n + plane.size)
        mean = {k: v[0] / v[1] for k, v in rows.items()}
        for bid in range(geo.n_blocks):
            (_, _), (yc, _), (_, _) = geo.block_axes(bid)
            for jj, yv in enumerate(yc):
                key = round(float(yv), 12)
                data[bid][:, jj, :] += np.interp(yv, yy, target) - mean[key]
        dset[...] = data

    print(f"wrote {a.outfile}: {a.scalar} mean re-targeted from {len(paths)} "
          f"snapshots (t = {min(snapshot_time(p) for p in paths):.2f} .. "
          f"{max(snapshot_time(p) for p in paths):.2f})")
    print(f"  measured D_tot: {dtot.min():.4e} .. {dtot.max():.4e} "
          f"(molecular {dmol:.4e}), implied J = {jnew:.6f} "
          f"-> theta+ at the centre = {1.0/jnew:.2f}")
    return 0


def snapshot_time(p):
    with h5py.File(p, "r") as f:
        return float(f.attrs["t_current"])


if __name__ == "__main__":
    sys.exit(main())
