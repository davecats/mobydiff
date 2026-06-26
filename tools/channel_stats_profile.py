#!/usr/bin/env python3
"""Single-snapshot channel statistics: mean velocity profile and rms of the
velocity fluctuations, from ONE block-table field (step), for one or more runs.

For each wall-normal cell-row (per refinement level, since the 2:1 wall bands are
finer) it averages over the homogeneous x,z directions across all blocks sharing
that row:  U(y) = <u>_{x,z},  u'_rms(y) = sqrt(<(u-<u>)^2>_{x,z})  (and v,w).
The field is a single instantaneous snapshot, so these are x,z-ensemble means at
that one time -- not time-averaged statistics.

Usage:
  python3 tools/channel_stats_profile.py OUT.png FIELD1.h5[:LABEL] [FIELD2.h5[:LABEL] ...]
"""
from __future__ import annotations
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def level_line(base, lev):
    line = base.copy()
    for _ in range(lev):
        mid = 0.5 * (line[:-1] + line[1:])
        new = np.empty(2 * len(line) - 1)
        new[0::2] = line
        new[1::2] = mid
        line = new
    return line


def profiles(path):
    """Return dict y -> (yc, mean[u,v,w], rms[u,v,w]) sorted by yc.
    h5py block axes are [block, k(z), j(y), i(x)]."""
    with h5py.File(path, "r") as f:
        nb = int(f.attrs["block_nb_x"])
        yb = f["y"][...]
        blocks = f["blocks"][...]
        D = {v: f[{"u": "un", "v": "vn", "w": "wn"}[v]][...] for v in "uvw"}
    rows = {}
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        for jj in range(nb):
            key = (int(lev), int(oy) + jj)
            # plane over (z,x) at this y-row jj (axis1 = y)
            rows.setdefault(key, []).append({v: D[v][bid, :, jj, :].ravel() for v in "uvw"})
    lines = {L: level_line(yb, L) for L in set(k[0] for k in rows)}
    recs = []
    for (lev, gj), planes in rows.items():
        yc = 0.5 * (lines[lev][gj] + lines[lev][gj + 1])
        cells = {v: np.concatenate([pl[v] for pl in planes]) for v in "uvw"}
        mean = {v: cells[v].mean() for v in "uvw"}
        rms = {v: np.sqrt(np.mean((cells[v] - mean[v]) ** 2)) for v in "uvw"}
        recs.append((yc, mean, rms))
    recs.sort(key=lambda r: r[0])
    return recs


def main():
    out = sys.argv[1]
    runs = []
    for sp in sys.argv[2:]:
        path, label = (sp.rsplit(":", 1) if ":" in sp[2:] else (sp, sp))
        runs.append((path, label))

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2))
    colors = ["tab:blue", "tab:red", "tab:green"]
    titles = ["mean U(y)", "u'_rms(y)", "v'_rms(y)", "w'_rms(y)"]
    for ri, (path, label) in enumerate(runs):
        recs = profiles(path)
        y = np.array([r[0] for r in recs])
        U = np.array([r[1]["u"] for r in recs])
        urms = np.array([r[2]["u"] for r in recs])
        vrms = np.array([r[2]["v"] for r in recs])
        wrms = np.array([r[2]["w"] for r in recs])
        c = colors[ri % len(colors)]
        for ax, data in zip(axes, [U, urms, vrms, wrms]):
            ax.plot(data, y, "-o", ms=2.5, lw=1.2, color=c, label=label)
        # print the near-wall interface rows for a numeric check
        print(f"\n== {label} ==  (y, U, u'rms, v'rms, w'rms)")
        for r in recs:
            if 0.55 < r[0] < 0.75 or 1.25 < r[0] < 1.45:  # around the interfaces
                print(f"  y={r[0]:.4f}  U={r[1]['u']:8.4f}  "
                      f"u'={r[2]['u']:.4f} v'={r[2]['v']:.4f} w'={r[2]['w']:.4f}")
    for ax, t in zip(axes, titles):
        ax.set_title(t)
        ax.set_ylabel("y")
        ax.grid(alpha=0.3)
        ax.axhline(0.643, color="grey", ls="--", lw=0.6)
        ax.axhline(1.357, color="grey", ls="--", lw=0.6)
    axes[0].legend(fontsize=9)
    fig.suptitle("Channel x,z-averaged profiles at step 250 (single snapshot; dashed = 2:1 interface)")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=120)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
