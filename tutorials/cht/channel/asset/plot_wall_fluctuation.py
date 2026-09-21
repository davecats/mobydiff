#!/usr/bin/env python3
"""Instantaneous temperature fluctuation at the conjugate interface.

    ./plot_wall_fluctuation.py [--snap Cstat_715000.h5] [--out wall_fluctuation.png]

A plan (x-z) view of theta' in the first fluid cell, for three walls that share
ONE velocity field and differ only in effusivity K = sqrt(kappa_s C_s). That is
the whole point of the sweep made visible: the same turbulence writes the same
thermal footprint on every wall, and the wall's own effusivity decides how much
of it survives.

TWO ROWS, because one scale cannot show both things:
  top     a SHARED colour scale -- the honest amplitude comparison. The
          near-isoflux wall (K = 0.1) fluctuates freely; the near-isothermal
          one (K = 100) is almost blank, which IS the result.
  bottom  each panel normalised by its OWN rms -- the structure comparison.
          The streaks do not change shape as K rises; only their amplitude
          collapses. A single scale would hide that, and per-panel scales
          alone would hide the collapse.

Colour: a DIVERGING map, two hues with a neutral midpoint, because theta' is
signed about zero. Never a rainbow, and never a hue at the midpoint.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "scalar"))
from scalar_tools import BlockGeometry                        # noqa: E402

RE = 149.0
YLO, YHI = 1.0, 3.0
LX, LZ = 4.0 * np.pi, 2.0 * np.pi
SHOW = [("a1", 0.1, "near-isoflux"), ("k1", 1.0, "matched"),
        ("a3", 100.0, "near-isothermal")]

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8b8a85"
SURFACE = "#fcfcfb"
# validated categorical slots 1 (blue) and 2 (orange) as the two poles, with a
# neutral midpoint -- the diverging rule from references/palette.md
DIVERGING = LinearSegmentedColormap.from_list(
    "cool_warm", ["#123f70", "#2a78d6", "#a8c8ee", "#efeeea",
                  "#f4b49a", "#eb6834", "#8c3413"])


def plane(f, name, jt, bl, nb, nx, nz):
    """An x-z plane assembled from the blocks, as (nx, nz).

    THE BLOCK DATASETS ARE STORED (nbz, nby, nbx) -- see
    scalar_tools.BlockGeometry.mesh, which says so. Slicing the middle (y)
    index therefore yields a (z, x) tile that must be TRANSPOSED before it is
    placed. Getting this wrong does not fail loudly: it produces a plausible
    field with a checkerboard of block-sized tiles, because each tile is
    internally transposed while the tiles themselves sit in the right places.
    The check that catches it is continuity -- the mean |d/dx| across a block
    boundary should match the mean inside a block, and it was 5-8x larger.
    """
    p = np.zeros((nx, nz))
    for b, (ox, oy, oz, _) in enumerate(bl):
        ox, oy, oz = int(ox), int(oy), int(oz)
        if oy <= jt < oy + nb[1]:
            p[ox:ox + nb[0], oz:oz + nb[2]] = f[name][b][:, jt - oy, :].T
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    snaps = sorted([q for q in glob.glob(os.path.join(here, "Cstat_*.h5"))
                    if re.fullmatch(r"Cstat_\d+\.h5", os.path.basename(q))],
                   key=lambda q: int(re.findall(r"\d+", os.path.basename(q))[0]))
    ap.add_argument("--snap", default=snaps[-1] if snaps else None)
    ap.add_argument("--stats", default=os.path.join(here, "flageul_clean.h5"))
    ap.add_argument("--out", default=os.path.join(here, "wall_fluctuation.png"))
    a = ap.parse_args()

    with h5py.File(a.stats, "r") as g:
        P, yy = g["profile"][...], g["coord"][...]
    m = (yy > YLO) & (yy < YHI)
    i0 = int(np.argmax(m))
    i1 = int(len(m) - 1 - np.argmax(m[::-1]))
    ttau = 0.5 * (abs(P[i0, 4]) + abs(P[i1, 6]))

    with h5py.File(a.snap, "r") as f:
        y = f["y"][...]
        yc = 0.5 * (y[:-1] + y[1:])
        bl = f["blocks"][...]
        nb = BlockGeometry(f).nb
        jt = int(np.argmax(yc > YLO))
        nx = int(bl[:, 0].max()) + nb[0]
        nz = int(bl[:, 2].max()) + nb[2]
        t = float(f.attrs["t_current"])
        fields = [( (plane(f, nm, jt, bl, nb, nx, nz)) , nm, K, lab)
                  for nm, K, lab in SHOW]
    yp_plane = (yc[jt] - YLO) * RE

    flucts = []
    for p, nm, K, lab in fields:
        flucts.append(((p - p.mean()) / ttau, nm, K, lab))
    shared = max(abs(np.percentile(fl, 0.5)) for fl, _, _, _ in flucts)
    shared = max(abs(np.percentile(flucts[0][0], [0.2, 99.8])).max(), 1e-9)

    xw = np.linspace(0, LX * RE, nx)
    zw = np.linspace(0, LZ * RE, nz)
    fig, axg = plt.subplots(2, 3, figsize=(15.2, 6.6), facecolor=SURFACE)
    ims = {}
    for row in range(2):
        for col, (fl, nm, K, lab) in enumerate(flucts):
            A = axg[row, col]
            if row == 0:
                dat, lim = fl, shared
                ttl = f"$K={K:g}$  ({lab})"
            else:
                dat, lim = fl / fl.std(), 3.0
            im = A.pcolormesh(xw, zw, dat.T, cmap=DIVERGING, vmin=-lim, vmax=lim,
                              shading="auto", rasterized=True)
            A.set_aspect("equal")
            A.set_facecolor(SURFACE)
            A.tick_params(colors=INK2, labelsize=8)
            for sp in A.spines.values():
                sp.set_color(MUTED)
            if row == 0:
                A.set_title(ttl, color=INK, fontsize=11, pad=6)
                A.annotate(f"rms $= {fl.std():.3f}\\,\\theta_\\tau$",
                           xy=(0.975, 0.055), xycoords="axes fraction",
                           ha="right", va="bottom", color=INK, fontsize=9,
                           fontweight="bold",
                           bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE,
                                     ec="none", alpha=0.85))
                A.tick_params(labelbottom=False)
            else:
                A.set_xlabel("$x^+$", color=INK2, fontsize=9.5)
            if col == 0:
                A.set_ylabel("$z^+$", color=INK2, fontsize=9.5)
            else:
                A.tick_params(labelleft=False)
            if col == 2:
                ims[row] = im

    axg[0, 0].annotate("SHARED scale — the amplitude collapses", xy=(0.015, 0.94),
                       xycoords="axes fraction", ha="left", va="top",
                       color=INK, fontsize=9.5, fontweight="bold",
                       bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE,
                                 ec=MUTED, lw=0.6, alpha=0.9))
    axg[1, 0].annotate("each by its OWN rms — the structure does not",
                       xy=(0.015, 0.94), xycoords="axes fraction", ha="left",
                       va="top", color=INK, fontsize=9.5, fontweight="bold",
                       bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE,
                                 ec=MUTED, lw=0.6, alpha=0.9))

    fig.suptitle("Instantaneous temperature fluctuation at the conjugate "
                 "interface — three walls, ONE velocity field",
                 color=INK, fontsize=12.5, fontweight="bold", x=0.008,
                 ha="left", y=0.985)
    fig.text(0.008, 0.012,
             f"Plan view of $\\theta'$ in the first fluid cell ($y^+={yp_plane:.2f}$) "
             f"at $t={t:.0f}$, $Re_\\tau=149$, $Pr=0.71$. The three walls differ ONLY "
             f"in effusivity $K=\\sqrt{{\\kappa_s C_s}}$: the same turbulence writes the "
             f"same footprint on each,\nand the wall decides how much survives — rms "
             f"{flucts[0][0].std():.2f}, {flucts[1][0].std():.2f}, "
             f"{flucts[2][0].std():.3f}$\\,\\theta_\\tau$, a factor "
             f"{flucts[0][0].std()/flucts[2][0].std():.0f} across the sweep.",
             color=MUTED, fontsize=8.5, ha="left", va="bottom")
    fig.subplots_adjust(left=0.045, right=0.905, top=0.885, bottom=0.135,
                        hspace=0.06, wspace=0.05)
    for row, lab in ((0, r"$\theta'/\theta_\tau$"),
                     (1, r"$\theta'\,/\,\mathrm{rms}(\theta')$")):
        pos = axg[row, 2].get_position()
        cax = fig.add_axes([0.918, pos.y0, 0.011, pos.height])
        cb = fig.colorbar(ims[row], cax=cax)
        cb.ax.tick_params(colors=INK2, labelsize=8)
        cb.outline.set_edgecolor(MUTED)
        cb.set_label(lab, color=INK2, fontsize=9.5)
    fig.savefig(a.out, dpi=150, facecolor=SURFACE)
    print(f"{a.out} written  (snapshot {os.path.basename(a.snap)}, t={t:.1f}, "
          f"plane y+ {yp_plane:.2f})")
    for fl, nm, K, lab in flucts:
        print(f"   {nm}  K={K:<6g} rms {fl.std():7.4f}  min {fl.min():+7.3f}  "
              f"max {fl.max():+7.3f}")


if __name__ == "__main__":
    main()
