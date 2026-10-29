#!/usr/bin/env python3
"""Upper-surface maps of the wing from a strip_loads pieces file: Cp and the
chordwise skin-friction sign (the separation line) in plan view.
Usage: surface_plot.py PIECES.npz OUT.png [--title T]"""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = np.load(sys.argv[1]); cen, nrm, cp, cf, tan = d["cen"], d["nrm"], d["cp"], d["cf_fit"], d["tangent"]
up = (nrm[:, 2] > 0.2) & (cen[:, 1] > 0.3) & (cen[:, 0] < 8.6) & np.isfinite(cp)
x, y = cen[up, 0], cen[up, 1]
cfx = cf[up] * np.sign(tan[up, 0])          # chordwise friction: positive downstream
fig, ax = plt.subplots(2, 1, figsize=(12, 11), sharex=True)
s0 = ax[0].scatter(x, y, c=cp[up], s=0.6, cmap="RdBu_r", vmin=-1.5, vmax=1.0, rasterized=True)
ax[0].set_title(sys.argv[3] if len(sys.argv) > 3 else sys.argv[1]); ax[0].set_ylabel("y / c"); ax[0].set_aspect("equal")
fig.colorbar(s0, ax=ax[0], label="Cp (upper surface)", shrink=0.8)
s1 = ax[1].scatter(x, y, c=np.clip(cfx, -0.004, 0.008), s=0.6, cmap="RdYlBu_r", vmin=-0.004, vmax=0.008, rasterized=True)
ax[1].set_ylabel("y / c"); ax[1].set_xlabel("x / c"); ax[1].set_aspect("equal")
fig.colorbar(s1, ax=ax[1], label="Cf,x (upper surface; blue = reversed flow)", shrink=0.8)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=130); print("wrote", sys.argv[2])
