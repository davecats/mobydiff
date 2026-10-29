#!/usr/bin/env python3
"""Upper-surface maps of the wing from a strip_loads pieces file: Cp and the
chordwise skin-friction sign (the separation line) in plan view.
Usage: surface_plot.py PIECES.npz OUT.png [--title T] [--smooth SIGMA]

--smooth SIGMA (chords; default 0 = off) replaces each piece's Cp and Cf,x by
the area-weighted Gaussian mean over the upper-surface pieces within 3 SIGMA
(in the plan-view distance). The pieces are STL triangles (~0.003 c) sampled
from the c/512 field (h = 0.002 c): the raw Cf,x carries piece-to-piece noise
at the cell scale, so SIGMA of a few cells (0.005-0.01) shows the mean
separation line without moving it."""
import sys, argparse, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument("pieces"); ap.add_argument("out")
ap.add_argument("--title"); ap.add_argument("--smooth", type=float, default=0.0)
a = ap.parse_args()
d = np.load(a.pieces); cen, nrm, cp, cf, tan, area = d["cen"], d["nrm"], d["cp"], d["cf_fit"], d["tangent"], d["area"]
up = (nrm[:, 2] > 0.2) & (cen[:, 1] > 0.3) & (cen[:, 0] < 8.6) & np.isfinite(cp)
x, y, w = cen[up, 0], cen[up, 1], area[up]
cp = cp[up]
cfx = cf[up] * np.sign(tan[up, 0])          # chordwise friction: positive downstream

def smooth(f, sigma):
    """Area-weighted Gaussian mean of f over the pieces within 3 sigma (plan
    view). Pieces without a value (a failed friction fit is NaN) get weight 0
    and stay NaN themselves, so the raw plot's gaps neither spread nor fill."""
    from scipy.spatial import cKDTree
    tree = cKDTree(np.column_stack([x, y]))
    ok = np.isfinite(f); wf = np.where(ok, w, 0.0); f0 = np.where(ok, f, 0.0)
    out = np.full_like(f, np.nan)
    for i, nb in enumerate(tree.query_ball_point(np.column_stack([x, y]), 3 * sigma)):
        if not ok[i]: continue
        nb = np.asarray(nb)
        k = wf[nb] * np.exp(-((x[nb] - x[i])**2 + (y[nb] - y[i])**2) / (2 * sigma**2))
        out[i] = np.sum(k * f0[nb]) / np.sum(k)
    return out

if a.smooth > 0:
    cp, cfx = smooth(cp, a.smooth), smooth(cfx, a.smooth)
tag = f" (Gaussian smoothing sigma = {a.smooth:g} c)" if a.smooth > 0 else ""
fig, ax = plt.subplots(2, 1, figsize=(12, 11), sharex=True)
s0 = ax[0].scatter(x, y, c=cp, s=0.6, cmap="RdBu_r", vmin=-1.5, vmax=1.0, rasterized=True)
ax[0].set_title((a.title or a.pieces) + tag); ax[0].set_ylabel("y / c"); ax[0].set_aspect("equal")
fig.colorbar(s0, ax=ax[0], label="Cp (upper surface)", shrink=0.8)
s1 = ax[1].scatter(x, y, c=np.clip(cfx, -0.004, 0.008), s=0.6, cmap="RdYlBu_r", vmin=-0.004, vmax=0.008, rasterized=True)
ax[1].set_ylabel("y / c"); ax[1].set_xlabel("x / c"); ax[1].set_aspect("equal")
fig.colorbar(s1, ax=ax[1], label="Cf,x (upper surface; blue = reversed flow)", shrink=0.8)
fig.tight_layout(); fig.savefig(a.out, dpi=130); print("wrote", a.out)
