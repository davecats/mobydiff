#!/usr/bin/env python3
"""The conjugate-channel validation figure: mean profiles and variances
against Flageul et al. (2015).

    ./plot_cht.py [--out cht_validation.png]

THREE PANELS, one question each:

  (a) theta+(y+)  -- the MEAN, against Kader's correlation. Kader is a
      published analytic curve, so it is plotted as a curve; it describes the
      constant-flux wall layer, so it is drawn solid where it applies and
      dashed beyond.
  (b) <theta'^2>/theta_tau^2 (y+) -- the VARIANCE, for the K sweep, with
      Flageul's values as discrete markers.
  (c) the interface coupling: the wall/near-wall-peak ratio against the
      wall-normal resolution, which is the quantity that converges to them.

HONESTY ABOUT THE REFERENCE. Flageul et al. publish no table; their figure 5
is a plot. The four values used here -- isoQ wall 4.2, conjugate wall 1.1,
isoT wall 0, peak 6.3 at y+ ~ 18 -- were read off that figure and are good to
about +-0.1, NOT to the digits printed. They are therefore drawn as MARKERS
WITH AN UNCERTAINTY BAR and never as a curve: drawing a reference line through
four eyeballed points would dress up a reading as data.

AND THE COMPARISON TRAP THIS FIGURE MAKES VISIBLE. Their thermal problem is
Kasagi's wall-flux formulation, whose flux falls linearly to zero at the
centreline; ours holds the flux CONSTANT across the channel. So their variance
decays toward the centre and ours does not -- panel (b) shows both, and the
near-wall peak (shaded band) is the only place the two are comparable. The
first pass at this comparison took max() over the whole channel and so
compared our centreline against their near-wall peak.
"""

from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RE, PR = 180.0, 0.71
YLO, YHI = 0.6, 2.6
NS = 7
S, SS, US, CLO, JLO, CHI, JHI = range(NS)

# the sweep, in the ini's order
SWEEP = [("k1", 1.0, 1.0), ("k2", 1.0, 100.0), ("k3", 1.0, 1.0e4),
         ("a1", 0.1, 0.1), ("a2", 10.0, 10.0), ("a3", 100.0, 100.0)]

# Flageul et al. 2015, figure 5, READ FROM THE PLOT (+-0.1, not digits)
FL = dict(wall_isoQ=4.2, wall_conjug=1.1, wall_isoT=0.0, peak=6.3, peak_yp=18.0,
          err=0.1)

# validated ordinal ramp (one hue, light->dark with K): see references/palette.md
RAMP = ["#86b6ef", "#3987e5", "#256abf", "#104281"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8b8a85"
SURFACE = "#fcfcfb"
ACCENT = "#eb6834"          # the reference, a different job -> a different hue
GRID = "#e4e3df"


def kader(yp, pr=PR):
    yp = np.asarray(yp, dtype=float)
    b = (3.85 * pr ** (1.0 / 3.0) - 1.3) ** 2 + 2.12 * np.log(pr)
    g = 0.01 * (pr * yp) ** 4 / (1.0 + 5.0 * pr ** 3 * yp)
    return pr * yp * np.exp(-g) + (2.12 * np.log(1.0 + yp) + b) * np.exp(-1.0 / np.maximum(g, 1e-300))


def load(path):
    with h5py.File(path, "r") as f:
        P, y = f["profile"][...], f["coord"][...]
    m = (y > YLO) & (y < YHI)
    i0 = int(np.argmax(m))
    i1 = int(len(m) - 1 - np.argmax(m[::-1]))
    yp = np.minimum(y - YLO, YHI - y) * RE          # distance to the NEARER wall
    out = []
    for i in range(P.shape[1] // NS):
        mean = P[:, NS * i + S]
        var = np.maximum(P[:, NS * i + SS] - mean ** 2, 0.0)
        tt = 0.5 * abs(P[i0, NS * i + JLO] + P[i1, NS * i + JHI])
        band = m & (yp > 5.0) & (yp < 40.0)
        out.append(dict(y=y, m=m, yp=yp, mean=mean, var=var / tt ** 2, ttau=tt,
                        i0=i0, i1=i1,
                        wall=0.5 * (var[i0] + var[i1]) / tt ** 2,
                        peak=float((var / tt ** 2)[band].max())))
    return out


def half(d, arr):
    """The lower half of the fluid, sorted by y+ — the two halves are
    statistically identical, so plotting both would only double the ink."""
    sel = d["m"] & (d["y"] < 0.5 * (YLO + YHI))
    o = np.argsort(d["yp"][sel])
    return d["yp"][sel][o], arr[sel][o]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats", default="cht_stats.h5")
    ap.add_argument("--yconv", nargs="*", default=["yc1p5_stats.h5", "yc1p0_stats.h5", "yc0p75_stats.h5"])
    ap.add_argument("--dyp", nargs="*", type=float, default=[1.5, 1.0, 0.75])
    ap.add_argument("--out", default="cht_validation.png")
    a = ap.parse_args()

    sweep = load(a.stats)
    # the four K values with alpha_s = 1 (a1, k1, a2, a3): one clean decade
    # ladder, and every one of them converged (see the README's reseed note)
    show = [("a1", 3, 0.1), ("k1", 0, 1.0), ("a2", 4, 10.0), ("a3", 5, 100.0)]

    fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.3), facecolor=SURFACE)
    for x in ax:
        x.set_facecolor(SURFACE)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            x.spines[sp].set_color(MUTED)
        x.tick_params(colors=INK2, labelsize=9, length=3)
        x.grid(True, color=GRID, lw=0.7, zorder=0)
        x.set_axisbelow(True)

    # ---- (a) the mean, against Kader ---------------------------------------
    A = ax[0]
    d = sweep[0]
    yp, mean = half(d, d["mean"])
    thp = (d["mean"][d["i0"]] - mean) / d["ttau"]
    kd = kader(yp) - kader(np.array([yp[0]]))[0]
    A.semilogx(yp, kd, "-", color=ACCENT, lw=2.0, zorder=3)
    A.semilogx(yp[yp > 40], kd[yp > 40], "--", color=ACCENT, lw=2.0, zorder=3)
    A.semilogx(yp, thp, "-", color=RAMP[1], lw=2.0, zorder=4)
    A.annotate("this work, K = 1", xy=(yp[-1], thp[-1]), xytext=(-4, 6),
               textcoords="offset points", ha="right", color=RAMP[2], fontsize=9.5)
    A.annotate("Kader (1981)", xy=(30, kader(np.array([30.0]))[0] - kader(np.array([yp[0]]))[0]),
               xytext=(4, -14), textcoords="offset points", color=ACCENT, fontsize=9.5)
    A.annotate("dashed: Kader beyond the\nconstant-flux layer it describes",
               xy=(0.97, 0.05), xycoords="axes fraction", color=MUTED,
               fontsize=8, ha="right", va="bottom")
    A.set_xlabel("$y^+$", color=INK2, fontsize=10)
    A.set_ylabel(r"$\theta^+$", color=INK2, fontsize=10)
    A.set_title("(a)  mean temperature", color=INK, fontsize=11, loc="left", pad=8)

    # ---- (b) the variance, with Flageul's read-off values -------------------
    B = ax[1]
    B.axvspan(5, 40, color=GRID, alpha=0.55, zorder=0)
    B.annotate("near-wall peak:\nthe only comparable region", xy=(0.5, 0.985),
               xycoords="axes fraction", color=MUTED, fontsize=8,
               ha="center", va="top")
    for (nm, idx, K), c in zip(show, RAMP):
        d = sweep[idx]
        yp, v = half(d, d["var"])
        B.loglog(yp, v, "-", color=c, lw=2.0, zorder=4)
        B.annotate(f"K = {K:g}", xy=(yp[0], v[0]), xytext=(7, -1),
                   textcoords="offset points", ha="left", va="center",
                   color=c, fontsize=9.5, zorder=7)
    # the reference: markers with an uncertainty bar, never a curve
    B.errorbar([0.52], [FL["wall_isoQ"]], yerr=FL["err"], fmt="s", ms=7,
               color=ACCENT, mfc=SURFACE, mew=1.8, capsize=3, zorder=6)
    B.errorbar([0.52], [FL["wall_conjug"]], yerr=FL["err"], fmt="o", ms=7,
               color=ACCENT, mfc=ACCENT, mew=1.8, capsize=3, zorder=6)
    B.errorbar([FL["peak_yp"]], [FL["peak"]], yerr=FL["err"], fmt="D", ms=7,
               color=ACCENT, mfc=SURFACE, mew=1.8, capsize=3, zorder=6)
    B.annotate("Flageul isoQ wall", xy=(0.52, FL["wall_isoQ"]), xytext=(11, 9),
               textcoords="offset points", ha="left", color=ACCENT, fontsize=8.5)
    B.annotate("Flageul conjugate wall", xy=(0.52, FL["wall_conjug"]),
               xytext=(11, -12), textcoords="offset points", ha="left",
               color=ACCENT, fontsize=8.5)
    B.annotate("Flageul peak", xy=(FL["peak_yp"], FL["peak"]), xytext=(6, -20),
               textcoords="offset points", ha="left", color=ACCENT, fontsize=8.5)
    B.annotate("their isoT wall is 0\n(off a log axis)", xy=(0.97, 0.04),
               xycoords="axes fraction", ha="right", va="bottom",
               color=ACCENT, fontsize=8)
    B.set_xlabel("$y^+$", color=INK2, fontsize=10)
    B.set_ylabel(r"$\langle\theta'^2\rangle/\theta_\tau^2$", color=INK2, fontsize=10)
    B.set_title("(b)  temperature variance", color=INK, fontsize=11, loc="left", pad=8)
    B.set_ylim(0.03, 16)
    B.set_xlim(0.40, 300)

    # ---- (c) the interface coupling vs wall-normal resolution --------------
    C = ax[2]
    dyp, ratio = [], []
    for path, d_ in zip(a.yconv, a.dyp):
        if not os.path.exists(path):
            continue
        s = load(path)[0]
        dyp.append(d_)
        ratio.append(s["wall"] / s["peak"])
    if len(dyp):
        o = np.argsort(dyp)[::-1]
        dyp = np.array(dyp)[o]; ratio = np.array(ratio)[o]
        C.plot(dyp, ratio, "-o", color=RAMP[2], lw=2.0, ms=8, zorder=4)
        for x_, y_ in zip(dyp, ratio):
            C.annotate(f"{y_:.4f}", xy=(x_, y_), xytext=(0, 9), textcoords="offset points",
                       ha="center", color=RAMP[3], fontsize=9)
    ref = FL["wall_conjug"] / FL["peak"]
    C.axhline(ref, color=ACCENT, lw=2.0, zorder=3)
    C.axhspan(ref * (1 - 0.09), ref * (1 + 0.09), color=ACCENT, alpha=0.12, zorder=1)
    C.annotate(f"Flageul  {ref:.4f}", xy=(0.03, ref), xycoords=("axes fraction", "data"),
               xytext=(0, 6), textcoords="offset points", ha="left",
               color=ACCENT, fontsize=9.5)
    C.annotate("shaded: their $\\pm0.1$ figure read-off", xy=(0.97, 0.04),
               xycoords="axes fraction", ha="right", va="bottom",
               color=ACCENT, fontsize=8)
    C.set_xticks(list(dyp))
    C.set_xticklabels([f"{x:g}" for x in dyp])
    C.set_xlim(1.62, 0.63)
    C.set_xlabel(r"$\Delta y^+$ at the wall   (refining $\rightarrow$)", color=INK2, fontsize=10)
    C.set_ylabel(r"$\langle\theta'^2\rangle_{wall}\,/\,\langle\theta'^2\rangle_{peak}$",
                 color=INK2, fontsize=10)
    C.set_title("(c)  interface coupling, K = 1", color=INK, fontsize=11, loc="left", pad=8)
    C.set_ylim(0.15, 0.205)

    fig.suptitle("Conjugate heat transfer in a turbulent channel, $Re_\\tau=180$, $Pr=0.71$ "
                 "— against Flageul et al. (2015), $Re_\\tau=149$",
                 color=INK, fontsize=12, x=0.008, ha="left", y=0.985)
    fig.text(0.008, 0.015,
             "Curves: this work (6 conjugate scalars on one velocity field). Markers: Flageul et al. fig. 5, "
             "read off the plot to $\\pm0.1$ — not tabulated.\n"
             "Their flux falls linearly to zero at the centreline (Kasagi wall-flux problem); ours is constant "
             "across the channel, so only the shaded near-wall band is comparable.",
             color=MUTED, fontsize=8, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.075, 1, 0.945))
    fig.savefig(a.out, dpi=170, facecolor=SURFACE)
    print(f"{a.out} written")

    # the table view the contrast WARN obliges
    print("\n  K      <t'2>_wall   <t'2>_nwpeak   wall/peak     (Flageul: 1.10 / 6.3 / 0.1746)")
    for nm, idx, K in show:
        d = sweep[idx]
        print(f"  {K:6g}   {d['wall']:10.3f}   {d['peak']:12.3f}   {d['wall']/d['peak']:9.4f}")
    if len(dyp):
        print("\n  dy+ at the wall:  " + "  ".join(f"{x:g}" for x in dyp))
        print("  wall/peak      :  " + "  ".join(f"{r:.4f}" for r in ratio))
    return 0


if __name__ == "__main__":
    sys.exit(main())
