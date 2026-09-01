#!/usr/bin/env python3
"""The two thermal problems against Flageul et al. 2015.

    ./plot_bulk.py [--out bulk_vs_constflux.png]

WHAT THIS FIGURE IS FOR. The conjugate scheme has been run at Re_tau 180 under
two different thermal problems on the SAME flow, the same grid and the same six
conjugate walls:

  CONSTANT FLUX   antisymmetric Dirichlet at the outer solid faces, no source.
                  The wall-normal flux is then the same at every height, so
                  thermal production never switches off and the variance rises
                  all the way to the centreline.
  BULK HEATING    a uniform volumetric source in the fluid, both outer faces at
                  zero. Now dJ/dy = S, the flux falls to zero at the centre and
                  the variance peaks near the wall -- the STRUCTURE of
                  Flageul's case, whose source is Kasagi's f_T u_x.

Panel (b) is the reason the two differ and is drawn from the runs' own data: it
puts our two flux profiles beside the Kasagi shape, computed here by
integrating THIS run's measured mean velocity. Flageul's flux lies between
ours, and so does their variance -- which is the honest way to report the
comparison.
"""

from __future__ import annotations

import argparse
import os

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RE, PR = 180.0, 0.71
YLO, YHI = 0.6, 2.6
NS = 7
S, SS, US, CLO, JLO, CHI, JHI = range(NS)

# Flageul et al. 2015 fig. 5, conjugate case G = G_2 = 1 (= our k1), read off
# the plot, so +-0.1 rather than digits.
FL_PEAK, FL_PEAK_YP, FL_WALL, FL_ERR = 6.3, 18.0, 1.1, 0.1

# validated categorical slots 1/3/2 (see references/palette.md); all three
# series are direct-labelled, which is also what the aqua's contrast WARN asks
# for.
C_CONST, C_BULK, C_REF = "#2a78d6", "#1baf7a", "#eb6834"
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
    # |.| per face: the two wall fluxes point opposite ways under bulk heating
    tt = 0.5 * (abs(P[i0, NS * isc + JLO]) + abs(P[i1, NS * isc + JHI]))
    yp = np.minimum(y - YLO, YHI - y) * RE
    lower = m & (y < 0.5 * (YLO + YHI))
    o = np.argsort(yp[lower])
    return dict(yp=yp[lower][o], var=(var / tt ** 2)[lower][o],
                flux=(P[:, NS * isc + JLO] / tt)[lower][o], ttau=tt,
                wall=0.5 * (var[i0] + var[i1]) / tt ** 2)


def kasagi_flux(vel_path):
    """The Kasagi flux shape from THIS run's mean velocity: with a source
    proportional to u, the flux at y is what has been generated below it."""
    with h5py.File(vel_path, "r") as f:
        P, y = f["profile"][...], f["coord"][...]
    m = (y > YLO) & (y < YHI)
    lower = m & (y < 0.5 * (YLO + YHI))
    u = np.maximum(P[lower, 0], 0.0)
    yy = y[lower]
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (u[1:] + u[:-1]) * np.diff(yy))])
    return (yy - YLO) * RE, 1.0 - cum / cum[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--const", default="cht_stats.h5")
    ap.add_argument("--bulk", default="bulk_stats.h5")
    ap.add_argument("--vel", default="cht_vel_stats.h5")
    ap.add_argument("--out", default="bulk_vs_constflux.png")
    a = ap.parse_args()

    cf, bk = load(a.const), load(a.bulk)

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.5), facecolor=SURFACE)
    for x in ax:
        x.set_facecolor(SURFACE)
        for s in ("top", "right"):
            x.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            x.spines[s].set_color(GRID)
        x.tick_params(colors=INK2, labelsize=9)
        x.grid(True, color=GRID, lw=0.6, zorder=0)
        x.set_axisbelow(True)
        x.set_xscale("log")
        x.set_xlim(0.7, 180)
        x.set_xlabel(r"$y^+$", color=INK2, fontsize=10)

    # ---- (a) the temperature variance, THEIR conjugate case = our k1 -------
    ax[0].plot(cf["yp"], cf["var"], color=C_CONST, lw=2.0, zorder=3)
    ax[0].plot(bk["yp"], bk["var"], color=C_BULK, lw=2.0, zorder=3)
    ax[0].errorbar([FL_PEAK_YP], [FL_PEAK], yerr=[FL_ERR], color=C_REF, marker="o",
                   ms=8, mew=0, lw=0, elinewidth=1.6, capsize=3, zorder=5)
    ax[0].errorbar([0.8], [FL_WALL], yerr=[FL_ERR], color=C_REF, marker="o",
                   ms=8, mew=0, lw=0, elinewidth=1.6, capsize=3, zorder=5)
    ax[0].annotate("constant flux\npeak 7.07", (60, 7.6), color=C_CONST,
                   fontsize=9, ha="center", fontweight="bold")
    ax[0].annotate("bulk heating\npeak 4.90", (60, 3.0), color=C_BULK,
                   fontsize=9, ha="center", fontweight="bold")
    ax[0].annotate("Flageul 2015\n(conjugate, $K=1$)", (18, 6.3),
                   xytext=(6.0, 8.6), color=C_REF, fontsize=9, ha="center",
                   fontweight="bold",
                   arrowprops=dict(arrowstyle="-", color=C_REF, lw=1.2))
    ax[0].set_ylim(0, 9.6)
    ax[0].set_ylabel(r"$\langle\theta'^2\rangle/\theta_\tau^2$", color=INK2, fontsize=10)
    ax[0].set_title("(a)  temperature variance, $\\kappa_s=\\alpha_s=1$",
                    color=INK, fontsize=11, loc="left", fontweight="bold")

    # ---- (b) WHY: the flux profile each problem imposes --------------------
    ax[1].plot(cf["yp"], np.abs(cf["flux"]), color=C_CONST, lw=2.0, zorder=3)
    ax[1].plot(bk["yp"], np.abs(bk["flux"]), color=C_BULK, lw=2.0, zorder=3)
    if os.path.exists(a.vel):
        kyp, kfl = kasagi_flux(a.vel)
        ax[1].plot(kyp, kfl, color=C_REF, lw=2.0, ls=(0, (5, 2)), zorder=4)
        ax[1].annotate("Kasagi source $\\propto u(y)$\n(Flageul's problem)",
                       (30, 0.86), xytext=(3.0, 0.42), color=C_REF, fontsize=9,
                       ha="center", fontweight="bold",
                       arrowprops=dict(arrowstyle="-", color=C_REF, lw=1.2))
    ax[1].annotate("constant flux", (14, 1.02), color=C_CONST, fontsize=9,
                   ha="center", fontweight="bold")
    ax[1].annotate("bulk heating\n$dJ/dy=S$", (110, 0.30), color=C_BULK,
                   fontsize=9, ha="center", fontweight="bold")
    ax[1].set_ylim(-0.03, 1.13)
    ax[1].set_ylabel(r"$|J(y)|\,/\,J_{\rm wall}$", color=INK2, fontsize=10)
    ax[1].set_title("(b)  the flux profile each problem imposes",
                    color=INK, fontsize=11, loc="left", fontweight="bold")

    fig.suptitle("Conjugate heat transfer at an immersed interface, "
                 "$Re_\\tau=180$, $Pr=0.71$ — the thermal problem sets the "
                 "variance level", color=INK, fontsize=12, fontweight="bold",
                 x=0.012, ha="left", y=0.985)
    fig.text(0.012, 0.005,
             "Both curves are the same solver, flow, grid and conjugate wall; "
             "only the thermal problem differs. Flageul's flux profile lies "
             "between ours, and so does their variance.",
             color=MUTED, fontsize=8.5, ha="left")
    fig.tight_layout(rect=(0, 0.035, 1, 0.945))
    fig.savefig(a.out, dpi=170, facecolor=SURFACE)
    print(f"{a.out} written")
    for nm, d in (("constant flux", cf), ("bulk heating ", bk)):
        # the NEAR-WALL peak and the max over the fluid are the same number
        # only when the flux vanishes at the centreline. Print both, always.
        b = (d["yp"] > 5.0) & (d["yp"] < 40.0)
        print(f"   {nm}: wall {d['wall']:.3f}   near-wall peak "
              f"{d['var'][b].max():.3f} at y+ {d['yp'][b][int(np.argmax(d['var'][b]))]:.1f}"
              f"   max over the fluid {d['var'].max():.3f}")


if __name__ == "__main__":
    main()
