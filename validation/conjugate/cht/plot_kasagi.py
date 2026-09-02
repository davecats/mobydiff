#!/usr/bin/env python3
"""The like-for-like conjugate channel against Flageul et al. (2015).

    ./plot_kasagi.py [--stats kasagi_stats.h5] [--out kasagi_vs_flageul.png]

Two panels, and the first is the PRECONDITION for the second:

  (a) THE FLUX PROFILE. A variance comparison only means something if both
      runs solve the same thermal problem, and the flux profile is what that
      means in practice. Ours is driven by Kasagi's beta*u_x source, so the
      mean flux must follow the CUMULATIVE FLOW RATE -- and the panel checks
      exactly that, against A*int(<u> dy) built from the run's OWN measured
      mean velocity. The two lie on top of each other by construction, which
      is the licence to call panel (b) one-to-one.
  (b) THE VARIANCE, on a linear axis so a 5 % agreement is actually visible
      (the log-axis version, plus the kappa_s sweep and the deviation, is in
      cht_validation.png).

The earlier constant-flux and bulk-heating campaigns BRACKETED the reference
and are what motivated this run; they are documented in README.md and are
deliberately not drawn here.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "scalar"))

RE, PR = 149.0, 0.71                 # the like-for-like run
YLO, YHI = 0.6, 2.6
NS = 7
S, SS, US, CLO, JLO, CHI, JHI = range(NS)
SOURCE = 0.064134                    # [scalar.N] source of cht149_kasagi.ini

HERE = os.path.dirname(os.path.abspath(__file__))
REF = {k: os.path.join(HERE, f"flageul_fig5a_{k}.dat")
       for k in ("conjug", "isoq", "isot")}

C_US, C_REF = "#2a78d6", "#eb6834"   # validated categorical slots 1 and 2
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8b8a85"
SURFACE, GRID = "#fcfcfb", "#e4e3df"


def load(path, isc=0):
    with h5py.File(path, "r") as f:
        P, y = f["profile"][...], f["coord"][...]
    m = (y > YLO) & (y < YHI)
    i0 = int(np.argmax(m))
    i1 = int(len(m) - 1 - np.argmax(m[::-1]))
    mean = P[:, NS * isc + S]
    var = np.maximum(P[:, NS * isc + SS] - mean ** 2, 0.0)
    # |.| per face: the two wall fluxes point opposite ways when the interior
    # is heated, where abs(j_lo + j_hi) would read zero.
    tt = 0.5 * (abs(P[i0, NS * isc + JLO]) + abs(P[i1, NS * isc + JHI]))
    yp = np.minimum(y - YLO, YHI - y) * RE
    low = m & (y < 0.5 * (YLO + YHI))
    o = np.argsort(yp[low])
    return dict(yp=yp[low][o], var=(var / tt ** 2)[low][o],
                flux=(P[:, NS * isc + JLO] / tt)[low][o], ttau=tt)


def kasagi_shape(pattern="kstat_*0000.h5"):
    """A*int(<u> dy) from the run's OWN snapshots -- what the source must give."""
    files = sorted(glob.glob(os.path.join(HERE, pattern)))
    if not files:
        return None, None
    from scalar_tools import BlockGeometry
    acc, ny, ynode = None, None, None
    for p in files:
        with h5py.File(p, "r") as f:
            ynode = f["y"][...]
            yc = 0.5 * (ynode[:-1] + ynode[1:])
            bl = f["blocks"][...]
            nb = BlockGeometry(f).nb[1]
            ny = yc.size
            u = np.zeros(ny)
            c = np.zeros(ny)
            for b, (ox, oy, oz, _) in enumerate(bl):
                u[oy:oy + nb] += f["un"][b].mean(axis=(0, 2))
                c[oy:oy + nb] += 1
            u /= c
        acc = u if acc is None else acc + u
    u = acc / len(files)
    yc = 0.5 * (ynode[:-1] + ynode[1:])
    low = (yc > YLO) & (yc < 0.5 * (YLO + YHI))
    dy = np.diff(ynode)[low]
    cum = np.cumsum((u[low] * dy)[::-1])[::-1]      # int_y^centre <u> dy'
    return (yc[low] - YLO) * RE, SOURCE * cum


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats", default=os.path.join(HERE, "kasagi_stats.h5"))
    ap.add_argument("--out", default=os.path.join(HERE, "kasagi_vs_flageul.png"))
    a = ap.parse_args()

    d = load(a.stats)
    ref = np.loadtxt(REF["conjug"])
    isoq, isot = np.loadtxt(REF["isoq"]), np.loadtxt(REF["isot"])

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.6), facecolor=SURFACE)
    for x in ax:
        x.set_facecolor(SURFACE)
        for s in ("top", "right"):
            x.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            x.spines[s].set_color(GRID)
        x.tick_params(colors=INK2, labelsize=9)
        x.grid(True, color=GRID, lw=0.6, zorder=0)
        x.set_axisbelow(True)
        x.set_xlabel(r"$y^+$", color=INK2, fontsize=10)

    # ---- (a) the flux profile: the precondition ----------------------------
    ky, kf = kasagi_shape()
    ax[0].semilogx(d["yp"], np.abs(d["flux"]), "-", color=C_US, lw=2.4, zorder=4)
    if ky is not None:
        ax[0].semilogx(ky, kf, ls=(0, (6, 2.5)), color=C_REF, lw=2.2, zorder=5)
        dev = np.abs(np.interp(ky, d["yp"], np.abs(d["flux"])) - kf).max()
        ax[0].annotate("Kasagi's law: $A\\int_y^h\\langle u\\rangle\\,dy$,\n"
                       "from this run's OWN mean velocity\n"
                       f"max departure {dev:.3f}",
                       xy=(0.03, 0.06), xycoords="axes fraction", ha="left",
                       va="bottom", color=C_REF, fontsize=8.5, fontweight="bold")
    ax[0].annotate("measured flux", xy=(0.60, 0.62), xycoords="axes fraction",
                   ha="left", color=C_US, fontsize=9.5, fontweight="bold")
    ax[0].set_xlim(0.45, 160)
    ax[0].set_ylim(-0.03, 1.12)
    ax[0].set_ylabel(r"$|J(y)|\,/\,J_{\rm wall}$", color=INK2, fontsize=10)
    ax[0].set_title("(a)  the flux profile IS theirs — the precondition",
                    color=INK, fontsize=11, loc="left", fontweight="bold")

    # ---- (b) the variance, linear axis so 5 % is visible -------------------
    gg = np.logspace(np.log10(max(isoq[0, 0], isot[0, 0])),
                     np.log10(min(isoq[-1, 0], isot[-1, 0])), 400)
    qv, tv = (np.interp(gg, isoq[:, 0], isoq[:, 1]),
              np.interp(gg, isot[:, 0], isot[:, 1]))
    ok = tv >= 0.3           # below this their isoT symbols are axis-clipped
    ax[1].fill_between(gg[ok], tv[ok], qv[ok], color=C_REF, alpha=0.10, lw=0, zorder=1)
    ax[1].semilogx(gg, qv, ls=":", color=C_REF, lw=1.3, zorder=3)
    ax[1].semilogx(gg[ok], tv[ok], ls=":", color=C_REF, lw=1.3, zorder=3)
    ax[1].semilogx(ref[:, 0], ref[:, 1], ls=(0, (6, 2.5)), color=C_REF, lw=2.4, zorder=5)
    ax[1].semilogx(d["yp"], d["var"], "-", color=C_US, lw=2.4, zorder=6)
    b = (d["yp"] > 5) & (d["yp"] < 40)
    lo, hi = max(ref[0, 0], d["yp"][0]), min(ref[-1, 0], d["yp"][-1])
    g = np.logspace(np.log10(lo), np.log10(hi), 300)
    rms = np.sqrt(np.mean((np.interp(g, d["yp"], d["var"])
                           - np.interp(g, ref[:, 0], ref[:, 1])) ** 2))
    ax[1].annotate(f"this work, $K=1$\npeak {d['var'][b].max():.2f}",
                   xy=(0.04, 0.985), xycoords="axes fraction", ha="left",
                   va="top", color=C_US, fontsize=9.5, fontweight="bold")
    ax[1].annotate(f"Flageul 2015, digitised\nconjugate, peak {ref[:,1].max():.2f}\n"
                   f"shaded: isoT-isoQ brackets",
                   xy=(0.04, 0.80), xycoords="axes fraction", ha="left",
                   va="top", color=C_REF, fontsize=9, fontweight="bold")
    ax[1].annotate(f"full-profile rms {rms:.2f} = {100*rms/ref[:,1].max():.1f} % of their peak",
                   xy=(0.5, 0.03), xycoords="axes fraction", ha="center",
                   va="bottom", color=INK2, fontsize=9, fontweight="bold")
    ax[1].set_xlim(0.45, 160)
    ax[1].set_ylim(0, 7.0)
    ax[1].set_ylabel(r"$\langle\theta'^2\rangle/\theta_\tau^2$", color=INK2, fontsize=10)
    ax[1].set_title("(b)  temperature variance, $\\kappa_s=\\alpha_s=1$",
                    color=INK, fontsize=11, loc="left", fontweight="bold")

    fig.suptitle("Conjugate heat transfer at an immersed interface, "
                 "$Re_\\tau=149$, $Pr=0.71$ — one-to-one against "
                 "Flageul et al. (2015)",
                 color=INK, fontsize=12, fontweight="bold", x=0.012, ha="left",
                 y=0.985)
    fig.text(0.012, 0.006,
             "Same thermal problem (Kasagi's beta*u_x source) and same Reynolds number as the "
             "reference, so this is a one-to-one comparison.\nReference DIGITISED from their "
             "fig. 5 (they tabulate nothing), cross-validated against their second panel to 0.015.",
             color=MUTED, fontsize=8.5, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.085, 1, 0.945))
    fig.savefig(a.out, dpi=170, facecolor=SURFACE)
    print(f"{a.out} written   (rms {rms:.3f}, peak {d['var'][b].max():.3f} "
          f"vs {ref[:,1].max():.3f})")


if __name__ == "__main__":
    main()
