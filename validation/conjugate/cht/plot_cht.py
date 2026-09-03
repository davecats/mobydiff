#!/usr/bin/env python3
"""The conjugate-channel validation figure: OUR LIKE-FOR-LIKE RUN vs Flageul.

    ./plot_cht.py [--stats flageul_stats.h5] [--re 149] [--interfaces 1.0 3.0]

This shows ONE case of ours -- the FLAGEUL-MATCHED run: Kasagi's beta*u_x
source, Re_tau = 149, their grid resolution (dy+ 0.49 -> 4.8, d+ = 149) and
their outer-wall boundary condition (imposed heat flux). The constant-flux and
bulk-heating campaigns BRACKETED the reference and the Kasagi/Re_tau-149 run
closed most of the gap; all three are documented in README.md and deliberately
not drawn here, so that every curve belongs to a single one-to-one comparison.

The interfaces and Re_tau are OPTIONS, not constants: the earlier campaigns
sit at 0.6/2.6 and Re_tau 180, and a wrong pair silently shifts every profile
(they set both the fluid mask and the y+ origin).

FOUR PANELS, each a comparison:

  (a) theta+(y+) against KADER's correlation -- a published analytic curve, so
      it is drawn as a curve, solid over the wall layer it describes and dashed
      beyond.
  (b) the VARIANCE for the kappa_s sweep, against their DIGITISED conjugate
      curve and the band spanned by their two IDEAL brackets. The headline.
  (c) the same comparison as a DEVIATION, which is the only way to read a 5 %
      agreement off a log plot.
  (d) the interface response against the effusivity K, with their conjugate
      value and their isoQ asymptote marked -- the sweep is what shows the
      agreement is not a K = 1 coincidence.

HONESTY ABOUT THE REFERENCE. Flageul et al. tabulate nothing; every reference
number here is DIGITISED from their figure 5 by digitize_flageul.py and
cross-validated against the same series in their second panel (0.015 for the
conjugate line, 0.2-0.4 for the two symbol series). isoT is unusable below
<T'^2> ~ 0.3, where their symbols are clipped by their own axis -- so that band
edge is drawn only where it is data, and the exact analytic limit (isoT -> 0 at
the wall) is used instead of the digitised floor.
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

# Defaults for the Flageul-MATCHED case (cht149_flageul.ini). The earlier
# campaigns sat at Re_tau 180 with the interfaces at 0.6/2.6 and a uniform y
# line; --re and --interfaces cover them.
RE, PR = 149.0, 0.71
YLO, YHI = 1.0, 3.0
NS = 7
S, SS, US, CLO, JLO, CHI, JHI = range(NS)

# the sweep, in the ini's order
SWEEP = [("k1", 1.0, 1.0), ("k2", 1.0, 100.0), ("k3", 1.0, 1.0e4),
         ("a1", 0.1, 0.1), ("a2", 10.0, 10.0), ("a3", 100.0, 100.0)]

# All three reference series are DIGITISED (digitize_flageul.py): the conjugate
# case from
# its solid line, and the two IDEAL BRACKETS from their symbol series -- isoQ
# (green x, the K -> 0 limit) and isoT (blue +, the K -> infinity limit). The
# brackets are good to ~0.2-0.4 in <T'^2> against ~0.015 for the line, because
# a symbol centroid is a coarser estimator; both numbers are in the .dat
# headers, measured by digitising each series off BOTH panels of their fig. 5.
HERE = os.path.dirname(os.path.abspath(__file__))
REF_DAT = os.path.join(HERE, "flageul_fig5a_conjug.dat")
REF_ISOQ = os.path.join(HERE, "flageul_fig5a_isoq.dat")
REF_ISOT = os.path.join(HERE, "flageul_fig5a_isot.dat")

# validated ordinal ramp (one hue, light->dark with K): see references/palette.md
RAMP = ["#86b6ef", "#3987e5", "#256abf", "#104281"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8b8a85"
SURFACE = "#fcfcfb"
ACCENT = "#eb6834"          # the reference, a different job -> a different hue
C_US   = "#2a78d6"          # this work, where one curve is drawn
OTHER  = "#1baf7a"          # the OTHER thermal problem: a third job, third hue
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
        # |.| of each face, not of the sum: identical where one flux crosses
        # the whole channel, but in the bulk-heating problem the two walls
        # are fed from the interior and their fluxes point OPPOSITE ways.
        tt = 0.5 * (abs(P[i0, NS * i + JLO]) + abs(P[i1, NS * i + JHI]))
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
    ap.add_argument("--stats", default="flageul_stats.h5",
                    help="the FLAGEUL-MATCHED run: Kasagi source, Re_tau = 149, "
                         "their grid resolution and their outer-wall BC")
    ap.add_argument("--re", type=float, default=RE)
    ap.add_argument("--interfaces", type=float, nargs=2, default=(YLO, YHI),
                    help="the two grid-aligned fluid/solid interfaces. They set "
                         "which rows are fluid AND the y+ origin, so a wrong "
                         "pair silently shifts every profile.")
    ap.add_argument("--out", default="cht_validation.png")
    a = ap.parse_args()
    globals()["RE"] = a.re
    globals()["YLO"], globals()["YHI"] = a.interfaces

    sweep = load(a.stats)
    ref = np.loadtxt(REF_DAT)
    isoq, isot = np.loadtxt(REF_ISOQ), np.loadtxt(REF_ISOT)
    # the four kappa_s = alpha_s cases: one clean decade ladder in K
    show = [("a1", 3, 0.1), ("k1", 0, 1.0), ("a2", 4, 10.0), ("a3", 5, 100.0)]

    fig, axg = plt.subplots(2, 2, figsize=(11.6, 8.0), facecolor=SURFACE)
    ax = axg.ravel()
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
    thp = np.abs(mean - d["mean"][d["i0"]]) / d["ttau"]
    kd = kader(yp) - kader(np.array([yp[0]]))[0]
    A.semilogx(yp, kd, "-", color=ACCENT, lw=2.0, zorder=3)
    A.semilogx(yp[yp > 40], kd[yp > 40], "--", color=ACCENT, lw=2.0, zorder=3)
    A.semilogx(yp, thp, "-", color=C_US, lw=2.2, zorder=4)
    A.annotate("this work, $K=1$", xy=(yp[-1], thp[-1]), xytext=(-4, 7),
               textcoords="offset points", ha="right", color=C_US, fontsize=9.5)
    A.annotate("Kader (1981)", xy=(30, kader(np.array([30.0]))[0] - kader(np.array([yp[0]]))[0]),
               xytext=(4, -15), textcoords="offset points", color=ACCENT, fontsize=9.5)
    dev = np.abs(thp - kd).max() / kd.max()
    A.annotate(f"Kader dashed beyond the wall layer it describes.\n"
               f"max departure {100*dev:.1f} % of the centreline value.",
               xy=(0.98, 0.03), xycoords="axes fraction", color=MUTED,
               fontsize=8, ha="right", va="bottom")
    A.set_xlabel("$y^+$", color=INK2, fontsize=10)
    A.set_ylabel(r"$\theta^+$", color=INK2, fontsize=10)
    A.set_title("(a)  mean temperature", color=INK, fontsize=11, loc="left", pad=8)

    # ---- (b) THE HEADLINE: the variance sweep vs their curve + brackets -----
    B = ax[1]
    gg = np.logspace(np.log10(max(isoq[0, 0], isot[0, 0])),
                     np.log10(min(isoq[-1, 0], isot[-1, 0])), 400)
    qv = np.interp(gg, isoq[:, 0], isoq[:, 1])
    tv = np.interp(gg, isot[:, 0], isot[:, 1])
    ok = tv >= 0.3          # below this their isoT symbols are axis-clipped
    B.fill_between(gg[ok], tv[ok], qv[ok], color=ACCENT, alpha=0.10, lw=0, zorder=1)
    B.loglog(gg, qv, ls=":", color=ACCENT, lw=1.3, zorder=5)
    B.loglog(gg[ok], tv[ok], ls=":", color=ACCENT, lw=1.3, zorder=5)
    B.loglog(ref[:, 0], ref[:, 1], ls=(0, (6, 2.5)), color=ACCENT, lw=2.4, zorder=6)
    for (nm, idx, K), c in zip(show, RAMP):
        dd = sweep[idx]
        ypc, v = half(dd, dd["var"])
        B.loglog(ypc, v, "-", color=c, lw=2.0, zorder=4)
        B.annotate(f"K = {K:g}", xy=(ypc[0], v[0]), xytext=(8, 0),
                   textcoords="offset points", ha="left", va="center",
                   color=c, fontsize=9.5, zorder=7,
                   bbox=dict(boxstyle="round,pad=0.12", fc=SURFACE, ec="none",
                             alpha=0.85))
    B.annotate("Flageul 2015, digitised:\ndashed = their CONJUGATE case\n"
               "dotted band = isoT to isoQ,\ntheir two IDEAL brackets",
               xy=(0.03, 0.97), xycoords="axes fraction", ha="left",
               va="top", color=ACCENT, fontsize=8.5, zorder=8)
    # K = 0.1 leaves the band in the CORE. That is documented physics, not a
    # discrepancy to hide: a near-insulating wall gives the bulk-heated core no
    # sink for large-scale thermal structure (README, campaign 3). Its WALL
    # value is the one that compares, and it matches their isoQ to 3 % -- see
    # panel (d).
    B.annotate("$K=0.1$ leaves the band in the CORE: a\n"
               "near-insulating wall gives the core no sink\n"
               "for large-scale structure. Its WALL value\n"
               "matches their isoQ to 3 % (panel d).",
               xy=(0.97, 0.04), xycoords="axes fraction", ha="right",
               va="bottom", color=MUTED, fontsize=7.8, zorder=8)
    B.set_xlabel("$y^+$", color=INK2, fontsize=10)
    B.set_ylabel(r"$\langle\theta'^2\rangle/\theta_\tau^2$", color=INK2, fontsize=10)
    B.set_title("(b)  temperature variance", color=INK, fontsize=11, loc="left", pad=8)
    B.set_ylim(0.015, 26)
    B.set_xlim(0.40, 300)

    # ---- (c) the same comparison as a DEVIATION ----------------------------
    # A 5 % agreement is invisible on a log plot; this is where it is read.
    C = ax[2]
    d1 = sweep[0]
    ypc, v = half(d1, d1["var"])
    lo, hi = max(ref[0, 0], ypc[0]), min(ref[-1, 0], ypc[-1])
    g = np.logspace(np.log10(lo), np.log10(hi), 300)
    ours, theirs = np.interp(g, ypc, v), np.interp(g, ref[:, 0], ref[:, 1])
    rms = np.sqrt(np.mean((ours - theirs)**2))
    band = 0.05*ref[:, 1].max()
    C.axhline(0.0, color=ACCENT, lw=2.0, zorder=3)
    C.axhspan(-band, band, color=ACCENT, alpha=0.10, lw=0, zorder=1)
    C.semilogx(g, ours - theirs, "-", color=C_US, lw=2.2, zorder=4)
    C.annotate(f"full-profile rms {rms:.2f}\n= {100*rms/ref[:,1].max():.1f} % of their peak",
               xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
               color=C_US, fontsize=9, fontweight="bold")
    C.annotate("shaded: $\\pm5$ % of their peak", xy=(0.03, 0.05),
               xycoords="axes fraction", ha="left", va="bottom",
               color=ACCENT, fontsize=8)
    C.annotate("the residual is a slightly FLATTER profile:\n"
               "high at the wall, low at the peak", xy=(0.97, 0.05),
               xycoords="axes fraction", ha="right", va="bottom",
               color=MUTED, fontsize=8)
    C.set_xlabel("$y^+$", color=INK2, fontsize=10)
    C.set_ylabel(r"$\langle\theta'^2\rangle - {\rm Flageul}$", color=INK2, fontsize=10)
    C.set_title("(c)  deviation from Flageul, $K=1$", color=INK, fontsize=11,
                loc="left", pad=8)
    C.set_ylim(-1.05, 1.05)

    # ---- (d) the interface response against the effusivity K ---------------
    D = ax[3]
    yq = float(np.minimum(sweep[0]["y"] - YLO, YHI - sweep[0]["y"])[sweep[0]["i0"]] * RE)
    # The LADDER is the four alpha_s = 1 cases; k2/k3 share K = 10/100 with
    # a2/a3 at a different alpha_s, so they are the EFFUSIVITY-COLLAPSE
    # partners and are drawn as separate markers, not joined into the line --
    # a connected line through a same-K pair would read as a vertical jump.
    lad = [(np.sqrt(k*c), sweep[i]["wall"]) for nm, i, K in show
           for _, k, c in [SWEEP[i]]]
    lad.sort()
    D.loglog([x for x, _ in lad], [y for _, y in lad], "-o", color=RAMP[2],
             lw=2.0, ms=8, zorder=4)
    for nm, i in (("k2", 1), ("k3", 2)):
        _, k, c = SWEEP[i]
        D.loglog([np.sqrt(k*c)], [sweep[i]["wall"]], "o", color=RAMP[2],
                 mfc=SURFACE, mew=2.0, ms=8, zorder=5)
    D.annotate("hollow: same $K$, different $\\alpha_s$\n"
               "(the effusivity-collapse spread)", xy=(0.97, 0.05),
               xycoords="axes fraction", ha="right", va="bottom",
               color=MUTED, fontsize=8)
    D.axhline(np.interp(yq, isoq[:, 0], isoq[:, 1]), color=ACCENT, lw=2.0,
              ls=":", zorder=3)
    D.annotate("their isoQ ($K\\rightarrow0$ limit)",
               xy=(0.03, np.interp(yq, isoq[:, 0], isoq[:, 1])),
               xycoords=("axes fraction", "data"), xytext=(0, 6),
               textcoords="offset points", color=ACCENT, fontsize=8.5)
    fc = float(np.interp(yq, ref[:, 0], ref[:, 1]))
    D.errorbar([1.0], [fc], yerr=[0.05*fc], fmt="D", ms=9, color=ACCENT,
               mfc=SURFACE, mew=2.0, capsize=3, zorder=6)
    D.annotate(f"their CONJUGATE case\n$K=1$:  {fc:.2f}", xy=(1.0, fc),
               xytext=(12, -4), textcoords="offset points", ha="left",
               color=ACCENT, fontsize=8.5)
    ours1 = float(sweep[0]["wall"])
    D.annotate(f"ours {ours1:.2f}  ({100*(ours1/fc-1):+.0f} %)", xy=(1.0, ours1),
               xytext=(-10, 10), textcoords="offset points", ha="right",
               color=RAMP[3], fontsize=8.5, fontweight="bold")
    D.set_xlabel(r"effusivity  $K=\sqrt{\kappa_s C_s}$", color=INK2, fontsize=10)
    D.set_ylabel(r"$\langle\theta'^2\rangle$ at $y^+=%.2f$" % yq, color=INK2, fontsize=10)
    D.set_title("(d)  interface response vs effusivity", color=INK, fontsize=11,
                loc="left", pad=8)

    fig.suptitle("Conjugate heat transfer at an immersed interface, $Re_\\tau=149$, "
                 "$Pr=0.71$ - one-to-one against Flageul et al. (2015)",
                 color=INK, fontsize=12, x=0.008, ha="left", y=0.985)
    fig.text(0.008, 0.010,
             "Ours: Kasagi's beta*u_x source at Re_tau = 149 - the SAME thermal problem and the SAME Reynolds number as the reference, so "
             "every curve here is a\none-to-one comparison. Reference DIGITISED from their fig. 5 (they tabulate nothing), cross-validated against "
             "their second panel to 0.015.",
             color=MUTED, fontsize=8, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.045, 1, 0.955))
    fig.savefig(a.out, dpi=170, facecolor=SURFACE)
    print(f"{a.out} written")
    print(f"\n  K       <t'2>_wall (y+ {yq:.2f})   near-wall peak    vs Flageul peak {ref[:,1].max():.3f}")
    for nm, idx, K in show:
        s = sweep[idx]
        print(f"  {K:6g}   {s['wall']:16.3f} {s['peak']:16.3f}"
              f"      ({100*(s['peak']/ref[:,1].max()-1):+.1f} %)")


if __name__ == "__main__":
    sys.exit(main())
