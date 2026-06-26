#!/usr/bin/env python3
"""Plot TIME-AVERAGED channel statistics from a developed-flow run's stats files
(channel_stats.h5 [+ _l1, _l2 ...] per refinement level), one or more runs.

For each run it reads every level file, keeps the rows that actually carried
cells (count > 0), combines them across levels and sorts by y -- giving the full
wall-normal profile at its native (fine near the walls, coarse in the core)
resolution. From profile[:, s] = time+space mean of stat s it forms:
  mean U(y);  u'/v'/w' rms = sqrt(<qq> - <q>^2);  Reynolds shear -<u'v'>.
u_tau = sqrt(forcing_x * h), h = half the wall-to-wall extent (from coord), so
the mean is shown as U+ vs y+ (semilog) with the law-of-the-wall references.

Usage: python3 tools/plot_channel_stats.py OUT.png STATS1.h5[:LABEL] [STATS2.h5[:LABEL] ...]
  STATS*.h5 is the LEVEL-0 stats file (channel_stats.h5); _l1 etc. are found
  automatically next to it.
"""
from __future__ import annotations
import os
import sys
import glob
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# profile column indices (Fortran STAT_* minus 1)
U, V, W, UU, VV, WW, UV = 0, 1, 2, 3, 4, 5, 6


def level_files(path0):
    """channel_stats.h5 -> [channel_stats.h5, channel_stats_l1.h5, ...]."""
    root, ext = os.path.splitext(path0)
    files = [path0]
    for lf in sorted(glob.glob(f"{root}_l*{ext}")):
        files.append(lf)
    return files


def read_run(path0):
    ys, prof, re, fx = [], [], None, None
    for lf in level_files(path0):
        with h5py.File(lf, "r") as f:
            coord = f["coord"][...]
            P = f["profile"][...]
            count = f["count"][...]
            re = float(f.attrs["re"]); fx = float(f.attrs["forcing_x"])
        m = count > 0
        ys.append(coord[m]); prof.append(P[m])
    y = np.concatenate(ys); P = np.concatenate(prof, axis=0)
    order = np.argsort(y)
    y, P = y[order], P[order]
    h = 0.5 * (y.min() + y.max())          # half wall-to-wall extent
    utau = np.sqrt(fx * h)
    out = dict(
        y=y, h=h, re=re, utau=utau,
        Umean=P[:, U],
        urms=np.sqrt(np.clip(P[:, UU] - P[:, U] ** 2, 0, None)),
        vrms=np.sqrt(np.clip(P[:, VV] - P[:, V] ** 2, 0, None)),
        wrms=np.sqrt(np.clip(P[:, WW] - P[:, W] ** 2, 0, None)),
        uv=-(P[:, UV] - P[:, U] * P[:, V]),
    )
    return out


def main():
    out = sys.argv[1]
    runs = [(sp.rsplit(":", 1) if ":" in sp[2:] else (sp, sp)) for sp in sys.argv[2:]]
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    colors = ["tab:blue", "tab:red", "tab:green"]
    yint = 0.643
    for ri, (path, label) in enumerate(runs):
        d = read_run(path)
        c = colors[ri % len(colors)]
        y, h, Re, ut = d["y"], d["h"], d["re"], d["utau"]
        wall = np.minimum(y, 2 * h - y)        # distance to nearest wall
        yp = wall * Re * ut
        # mean U+ vs y+ (fold both walls)
        o = np.argsort(yp)
        ax[0, 0].semilogx(yp[o], d["Umean"][o] / ut, "-", color=c, lw=1.3, label=label)
        # rms and Reynolds stress vs y (full span)
        ax[0, 1].plot(d["urms"], y, "-", color=c, lw=1.3, label=label)
        ax[0, 2].plot(d["vrms"], y, "-", color=c, lw=1.3, label=label)
        ax[1, 0].plot(d["wrms"], y, "-", color=c, lw=1.3, label=label)
        ax[1, 1].plot(d["uv"], y, "-", color=c, lw=1.3, label=label)
        ax[1, 2].plot(d["urms"] / ut, yp, "-", color=c, lw=1.3, label=label)  # u'+ vs y+
        print(f"{label}: Re_tau={Re*ut:.1f}, u_tau={ut:.4f}, h={h:.4f}, "
              f"U+_core={d['Umean'].max()/ut:.2f}")
    yp_ref = np.logspace(0, np.log10(180), 50)
    ax[0, 0].semilogx(yp_ref, yp_ref, "k:", lw=1)
    ax[0, 0].semilogx(yp_ref, np.log(yp_ref) / 0.41 + 5.2, "k--", lw=1)
    ax[0, 0].set(xlabel="y+", ylabel="U+", title="mean U+ vs y+ (log law)", xlim=(0.8, 200), ylim=(0, 20))
    for a, t in [(ax[0, 1], "u'_rms(y)"), (ax[0, 2], "v'_rms(y)"), (ax[1, 0], "w'_rms(y)"),
                 (ax[1, 1], "-<u'v'>(y)")]:
        a.set(title=t, ylabel="y"); a.axhline(yint, color="grey", ls="--", lw=0.6)
        a.axhline(2 - yint, color="grey", ls="--", lw=0.6); a.grid(alpha=0.3)
    ax[1, 2].set(title="u'+ vs y+", xlabel="u'+", ylabel="y+"); ax[1, 2].grid(alpha=0.3)
    for a in ax.ravel():
        a.legend(fontsize=8)
    fig.suptitle("Developed channel, time-averaged statistics (dashed = 2:1 interface)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out, dpi=115)
    print("wrote", out)


if __name__ == "__main__":
    main()
