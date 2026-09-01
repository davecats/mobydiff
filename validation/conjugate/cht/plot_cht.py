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
Kasagi's wall-flux formulation, whose flux falls to zero at the centreline. Two
of ours are drawn: the CONSTANT-FLUX campaign, whose flux is the same at every
height so its variance never stops rising, and the BULK-HEATING one (dJ/dy = S),
whose flux vanishes at the centre as theirs does and which therefore peaks near
the wall and decays. Flageul lies between them -- see the caption. The first
pass at this comparison took max() over the whole channel and so compared our
centreline against their near-wall peak; the shaded band is where the
constant-flux curves may be read.
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
# The CONJUGATE case is no longer read off by eye: flageul_fig5a_conjug.dat is
# the digitised black solid curve (digitize_flageul.py, validated against the
# same curve in their panel 5b to 0.02). Its peak is 6.208 at y+ 17.6 -- and
# its value at OUR first cell, y+ = 0.75, is 1.270, NOT the 1.1 that reading
# the curve where it meets the axis suggested: their first point is at
# y+ = 0.49 and the curve is still climbing there.
# isoQ/isoT remain eyeballed brackets, so they stay markers with a bar.
# All three are digitised now (digitize_flageul.py): the conjugate case from
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
    ap.add_argument("--stats", default="cht_stats.h5")
    ap.add_argument("--bulk", default="bulk_stats.h5",
                    help="the BULK-HEATING campaign (dJ/dy = S). Drawn beside "
                         "the constant-flux one wherever the two differ, "
                         "because Flageul's flux profile lies between them.")
    ap.add_argument("--yconv", nargs="*", default=["yc1p5_stats.h5", "yc1p0_stats.h5", "yc0p75_stats.h5"])
    ap.add_argument("--dyp", nargs="*", type=float, default=[1.5, 1.0, 0.75])
    ap.add_argument("--out", default="cht_validation.png")
    a = ap.parse_args()

    sweep = load(a.stats)
    bulk = load(a.bulk) if os.path.exists(a.bulk) else None
    # the four K values with alpha_s = 1 (a1, k1, a2, a3): one clean decade
    # ladder, and every one of them converged (see the README's reseed note)
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
    thp = (d["mean"][d["i0"]] - mean) / d["ttau"]
    kd = kader(yp) - kader(np.array([yp[0]]))[0]
    A.semilogx(yp, kd, "-", color=ACCENT, lw=2.0, zorder=3)
    A.semilogx(yp[yp > 40], kd[yp > 40], "--", color=ACCENT, lw=2.0, zorder=3)
    A.semilogx(yp, thp, "-", color=RAMP[1], lw=2.0, zorder=4)
    A.annotate("constant flux, K = 1", xy=(yp[-1], thp[-1]), xytext=(-4, 6),
               textcoords="offset points", ha="right", color=RAMP[2], fontsize=9.5)
    if bulk is not None:
        db = bulk[0]
        ypb, mb = half(db, db["mean"])
        # the bulk-heated fluid is HOTTER than its wall, the opposite
        # sense; abs() measures both along the direction heat flows and
        # is exact here because each profile is monotone from its wall.
        thb = np.abs(mb - db["mean"][db["i0"]]) / db["ttau"]
        A.semilogx(ypb, thb, "-", color=OTHER, lw=2.0, zorder=4)
        A.annotate("bulk heating, K = 1", xy=(ypb[-1], thb[-1]), xytext=(-4, -16),
                   textcoords="offset points", ha="right", color=OTHER, fontsize=9.5)
    A.annotate("Kader (1981)", xy=(30, kader(np.array([30.0]))[0] - kader(np.array([yp[0]]))[0]),
               xytext=(4, -14), textcoords="offset points", color=ACCENT, fontsize=9.5)
    # measured, not asserted: max|theta+ - Kader| over the whole profile.
    # Kader is FITTED to channel/pipe data, whose flux decays in the outer
    # region -- so it matches the bulk-heated channel, not the literally
    # constant-flux one, which is the reverse of what its derivation suggests.
    A.annotate("Kader dashed beyond the layer it describes.\n"
               "Bulk heating follows it to within 4.3 % of the\n"
               "centreline value; constant flux departs by 34 % — (b) is why.",
               xy=(0.98, 0.02), xycoords="axes fraction", color=MUTED,
               fontsize=8, ha="right", va="bottom")
    A.set_xlabel("$y^+$", color=INK2, fontsize=10)
    A.set_ylabel(r"$\theta^+$", color=INK2, fontsize=10)
    A.set_title("(a)  mean temperature", color=INK, fontsize=11, loc="left", pad=8)

    # ---- (b) the LOCAL SLOPE: is the profile logarithmic? -------------------
    # It is not, and this panel is here because the eye cannot tell on a
    # semilog plot. In this thermal problem the TOTAL flux is constant by
    # construction, so dtheta+/dy+ = Pr/(1 + Pr D_t/nu): it equals Pr in the
    # conduction sublayer, follows 1/(Pr_t kappa y+) only through the log
    # layer, and then TURNS BACK UP toward the centreline because the eddy
    # diffusivity falls off there while the flux may not. A logarithmic
    # profile would keep decreasing as 1/y+.
    Sx = ax[1]
    d = sweep[0]
    yp, mean = half(d, d["mean"])
    thp_ = (d["mean"][d["i0"]] - mean) / d["ttau"]
    slope = np.gradient(thp_, yp)
    ypc, conv = half(d, np.zeros_like(d["mean"]))     # placeholder for shape
    Sx.loglog(yp, slope, "-", color=RAMP[1], lw=2.0, zorder=5)
    Sx.axhline(PR, color=MUTED, lw=1.4, ls=":", zorder=3)
    Sx.loglog(yp[yp > 8], 1.0 / (0.85 * 0.41 * yp[yp > 8]), "--", color=ACCENT,
              lw=1.8, zorder=4)
    jmin = int(np.argmin(slope[yp > 8])) + int((yp <= 8).sum())
    Sx.plot([yp[jmin]], [slope[jmin]], "o", ms=8, color=RAMP[3], mfc=SURFACE,
            mew=2.0, zorder=6)
    Sx.annotate(f"minimum at $y^+$ = {yp[jmin]:.0f},\nthen it turns back UP",
                xy=(yp[jmin], slope[jmin]), xytext=(-6, -34),
                textcoords="offset points", ha="right", color=RAMP[3], fontsize=8.5)
    Sx.annotate("$Pr$ = 0.71 (conduction)", xy=(0.55, PR), xytext=(0, 6),
                textcoords="offset points", color=MUTED, fontsize=8.5)
    Sx.annotate("$1/(Pr_t\\kappa y^+)$  (logarithmic)", xy=(22, 1/(0.85*0.41*22)),
                xytext=(-6, -20), textcoords="offset points", ha="right",
                color=ACCENT, fontsize=8.5)
    Sx.annotate("constant flux", xy=(2.0, slope[int(np.argmin(np.abs(yp-2.0)))]),
                xytext=(6, -14), textcoords="offset points", color=RAMP[2], fontsize=9)
    if bulk is not None:
        db = bulk[0]
        ypb, mb = half(db, db["mean"])
        sb = np.abs(np.gradient(np.abs(mb - db["mean"][db["i0"]]) / db["ttau"], ypb))
        Sx.loglog(ypb, np.maximum(sb, 2e-3), "-", color=OTHER, lw=2.0, zorder=5)
        Sx.set_ylim(1.5e-3, 1.4)
        Sx.annotate("bulk heating: keeps FALLING,\n"
                    "because $J\\rightarrow0$ at the centre",
                    xy=(ypb[int(0.82*len(ypb))],
                        max(sb[int(0.82*len(ypb))], 2e-3)), xytext=(-14, -30),
                    textcoords="offset points", ha="right", color=OTHER,
                    fontsize=8.5)
    Sx.set_xlabel("$y^+$", color=INK2, fontsize=10)
    Sx.set_ylabel(r"$d\theta^+/dy^+$", color=INK2, fontsize=10)
    Sx.set_title("(b)  the profile is NOT logarithmic outside the log layer",
                 color=INK, fontsize=11, loc="left", pad=8)
    Sx.set_xlim(0.5, 260)

    # ---- (c) the variance, with Flageul's read-off values -------------------
    B = ax[2]
    # THE BRACKETS, as a band: their four boundary conditions span isoT (ideal
    # Dirichlet, the wall kills the fluctuation) to isoQ (ideal Neumann, the
    # wall follows it freely). Any conjugate wall must lie inside, and the band
    # is the useful statement -- the two edges separately are not.
    if os.path.exists(REF_ISOQ) and os.path.exists(REF_ISOT):
        q, t = np.loadtxt(REF_ISOQ), np.loadtxt(REF_ISOT)
        gg = np.logspace(np.log10(max(q[0, 0], t[0, 0])),
                         np.log10(min(q[-1, 0], t[-1, 0])), 400)
        qv, tv = np.interp(gg, q[:, 0], q[:, 1]), np.interp(gg, t[:, 0], t[:, 1])
        # Draw the isoT edge only where it is DATA. Below <T\'^2> ~ 0.3 their
        # '+' symbols are clipped by their own axis, so the digitised values
        # there are a floor, not a measurement (see digitize_flageul.py).
        ok = tv >= 0.3
        B.fill_between(gg[ok], tv[ok], qv[ok], color=ACCENT, alpha=0.10, lw=0, zorder=1)
        B.loglog(gg, qv, ls=":", color=ACCENT, lw=1.3, zorder=5)
        B.loglog(gg[ok], tv[ok], ls=":", color=ACCENT, lw=1.3, zorder=5)
        B.annotate("isoQ", xy=(gg[int(0.62*gg.size)], qv[int(0.62*gg.size)]),
                   xytext=(0, 7), textcoords="offset points", ha="center",
                   color=ACCENT, fontsize=8.5, zorder=8)
        B.annotate("isoT", xy=(gg[int(0.30*gg.size)], tv[int(0.30*gg.size)]),
                   xytext=(2, -12), textcoords="offset points", ha="center",
                   color=ACCENT, fontsize=8.5, zorder=8)
    for (nm, idx, K), c in zip(show, RAMP):
        d = sweep[idx]
        yp, v = half(d, d["var"])
        B.loglog(yp, v, "-", color=c, lw=2.0, zorder=4)
        B.annotate(f"K = {K:g}", xy=(yp[0], v[0]), xytext=(7, -1),
                   textcoords="offset points", ha="left", va="center",
                   color=c, fontsize=9.5, zorder=7)
    # The SAME wall (K = 1) under the other thermal problem. Flageul's peak
    # falls between the two curves, which is the whole comparison: their flux
    # profile lies between a constant one and our linear one.
    if bulk is not None:
        db = bulk[0]
        ypb, vb = half(db, db["var"])
        B.loglog(ypb, vb, "-", color=OTHER, lw=2.0, zorder=5)
        B.annotate("bulk heating, K = 1\n(peaks near the wall and\ndecays,"
                   " as theirs does)",
                   xy=(150, vb[int(np.argmin(np.abs(ypb - 150)))]),
                   xytext=(-4, -46), textcoords="offset points", ha="right",
                   color=OTHER, fontsize=8.5, zorder=7)
    # THEIR CONJUGATE CASE, as the digitised curve rather than three points
    if os.path.exists(REF_DAT):
        r = np.loadtxt(REF_DAT)
        B.loglog(r[:, 0], r[:, 1], ls=(0, (6, 2.5)), color=ACCENT, lw=2.2, zorder=6)
        B.annotate("Flageul 2015, digitised:\n"
                   "dashed = their CONJUGATE case\n"
                   "dotted band = their IDEAL brackets",
                   xy=(0.03, 0.97), xycoords="axes fraction", ha="left",
                   va="top", color=ACCENT, fontsize=8.5, zorder=8)
    B.set_xlabel("$y^+$", color=INK2, fontsize=10)
    B.set_ylabel(r"$\langle\theta'^2\rangle/\theta_\tau^2$", color=INK2, fontsize=10)
    B.set_title("(c)  temperature variance", color=INK, fontsize=11, loc="left", pad=8)
    B.set_ylim(0.03, 16)
    B.set_xlim(0.40, 300)

    # ---- (c) the interface coupling vs wall-normal resolution --------------
    # THE RATIO MUST BE READ AT A MATCHED HEIGHT. Taking each study's own
    # FIRST CELL compares different y+ -- ours moves with the resolution
    # (0.75 / 0.5 / 0.375) while theirs sits at 0.49 -- so the ratio appeared
    # to "converge" onto the reference simply because the sampling point slid
    # down the rising near-wall curve. Both resolve y+ = 0.75, so read there.
    C = ax[3]
    YQ = 0.75
    dyp, ratio, raw = [], [], []
    for path, d_ in zip(a.yconv, a.dyp):
        if not os.path.exists(path):
            continue
        s = load(path)[0]
        ypq, vq = half(s, s["var"])
        dyp.append(d_)
        ratio.append(float(np.interp(YQ, ypq, vq)) / s["peak"])
        raw.append(s["wall"] / s["peak"])
    if len(dyp):
        o = np.argsort(dyp)[::-1]
        dyp = np.array(dyp)[o]; ratio = np.array(ratio)[o]; raw = np.array(raw)[o]
        C.plot(dyp, ratio, "-o", color=RAMP[2], lw=2.0, ms=8, zorder=4)
        C.plot(dyp, raw, "--o", color=MUTED, lw=1.4, ms=5, zorder=3)
        for x_, y_ in zip(dyp, ratio):
            C.annotate(f"{y_:.4f}", xy=(x_, y_), xytext=(0, 9), textcoords="offset points",
                       ha="center", color=RAMP[3], fontsize=9)
        C.annotate("grey: each run's own FIRST CELL — not a fixed\n"
                   "height (0.75/0.5/0.375), so not comparable",
                   xy=(0.5, 0.055), xycoords="axes fraction",
                   ha="center", va="bottom", color=MUTED, fontsize=8)
    r = np.loadtxt(REF_DAT)
    ref = float(np.interp(YQ, r[:, 0], r[:, 1])) / r[:, 1].max()
    C.axhline(ref, color=ACCENT, lw=2.0, zorder=3)
    C.axhspan(ref * (1 - 0.03), ref * (1 + 0.03), color=ACCENT, alpha=0.12, zorder=1)
    C.annotate(f"Flageul  {ref:.4f}  (digitised)", xy=(0.03, ref), xycoords=("axes fraction", "data"),
               xytext=(0, 6), textcoords="offset points", ha="left",
               color=ACCENT, fontsize=9.5)
    C.annotate("shaded: $\\pm3$ %", xy=(0.985, 0.90),
               xycoords="axes fraction", ha="right", va="top",
               color=ACCENT, fontsize=8)
    C.set_xticks(list(dyp))
    C.set_xticklabels([f"{x:g}" for x in dyp])
    C.set_xlim(1.62, 0.63)
    C.set_xlabel(r"$\Delta y^+$ at the wall   (refining $\rightarrow$)", color=INK2, fontsize=10)
    C.set_ylabel(r"$\langle\theta'^2\rangle_{wall}\,/\,\langle\theta'^2\rangle_{peak}$",
                 color=INK2, fontsize=10)
    C.set_title("(d)  interface coupling at a MATCHED $y^+=0.75$ (constant flux)",
                color=INK, fontsize=11, loc="left", pad=8)
    C.set_ylim(0.168, 0.215)

    fig.suptitle("Conjugate heat transfer in a turbulent channel, $Re_\\tau=180$, $Pr=0.71$ "
                 "— against Flageul et al. (2015), $Re_\\tau=149$",
                 color=INK, fontsize=12, x=0.008, ha="left", y=0.985)
    fig.text(0.008, 0.010,
             "Curves: this work. Dashed orange: Flageul et al.'s conjugate case DIGITISED from their fig. 5 (their paper\n"
             "tabulates nothing); cross-validated against the same curve in their panel 5b to 0.02. Blue ramp = the K sweep\n"
             "under CONSTANT FLUX; green = the same K = 1 wall under BULK HEATING. Their Kasagi flux profile lies BETWEEN\n"
             "the two, and so does their variance: near-wall peak 7.07 / 6.21 (Flageul) / 4.89. The reference is BRACKETED.",
             color=MUTED, fontsize=8, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.085, 1, 0.955))
    fig.savefig(a.out, dpi=170, facecolor=SURFACE)
    print(f"{a.out} written")

    # the table view the contrast WARN obliges
    rr = np.loadtxt(REF_DAT)
    print(f"\n  K      <t'2>_wall   <t'2>_nwpeak   wall/peak     "
          f"(Flageul digitised: peak {rr[:,1].max():.3f} at y+ "
          f"{rr[int(np.argmax(rr[:,1])),0]:.1f}; at y+0.75 {np.interp(0.75,rr[:,0],rr[:,1]):.3f})")
    for nm, idx, K in show:
        d = sweep[idx]
        print(f"  {K:6g}   {d['wall']:10.3f}   {d['peak']:12.3f}   {d['wall']/d['peak']:9.4f}")
    if len(dyp):
        print("\n  dy+ at the wall:  " + "  ".join(f"{x:g}" for x in dyp))
        print("  wall/peak      :  " + "  ".join(f"{r:.4f}" for r in ratio))
    return 0


if __name__ == "__main__":
    sys.exit(main())
