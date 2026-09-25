#!/usr/bin/env python3
"""Publication figures for the conjugate pipe vs Neuhauser's NekRS DNS.

Every panel is built from the accumulators pipe_stats.py writes and from the
reference zarr through compare_neuhauser.py, so the figures and the tables in
pipe_report.md come from the same numbers.

CONVENTIONS, identical on both sides of every comparison:
  * wall units from THEIR definitions (u_tau from d<w>/dr at r -> 0.5,
    y+ = (0.5-r)/delta_v, delta_v = nu/u_tau);
  * our theta_tau = q_w/u_tau with q_w from the INTERFACE-HEAT BALANCE
    |Q|/(2 pi R L), never from the near-wall slope -- the slope is differenced
    across the cut cells and reads 6.6 % low on the coarse grid;
  * the mean temperature is referenced to the interface value on both sides,
    because their `t` is wall-referenced and ours is absolute.

    ./plot_pipe.py            # writes fig1..fig6 png + pipe_report tables
"""

from __future__ import annotations

import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))      # pipe_stats.py, beside the inis
from pipe_stats import profiles, azimuthal_mean
from compare_neuhauser import ref_case

NU = 1.0 / 5300.0
U_TAU = 0.068379
DELTA_V = 2.759318e-3
R, LZ = 0.5, 12.5
PR = 0.71
NAMES = ["c0", "c1", "c2", "c3", "c4", "mbc", "isof"]
LABEL = {"c0": r"CHT $K{=}1,\ \lambda_{sf}{=}1$", "c1": r"CHT $\lambda_{sf}{=}2$",
         "c2": r"CHT $\lambda_{sf}{=}0.5$", "c3": r"CHT $K{=}4$",
         "c4": r"CHT $K{=}0.25$", "mbc": "MBC", "isof": "IF"}
COLOR = {"c0": "#1f77b4", "c1": "#2ca02c", "c2": "#17becf", "c3": "#d62728",
         "c4": "#9467bd", "mbc": "#7f7f7f", "isof": "#ff7f0e"}
# interface heat of each run -> q_w = |Q|/(2 pi R L)
QW = {"coarse": 9.474 / (2 * np.pi * R * LZ), "fine": 9.649 / (2 * np.pi * R * LZ),
      "z-fine": 0.24659, "production": 0.24854}
# The campaign ran four grids; the tutorial keeps the PRODUCTION one, whose
# accumulators are shipped in this directory. The other three survive as the
# resolution study of Fig. 6 -- measured error levels (hard-coded there, they
# are numbers not fields) and the wall-residual maps in grid_residual.npz.
GRIDS = {"coarse": dict(dxp=2.94, ls="-", mk="o"),
         "fine": dict(dxp=1.84, ls="--", mk="s"),
         "z-fine": dict(dxp=2.94, ls="-.", mk="^"),
         "production": dict(stats="pipe_prod_statsD.npz",
                            snaps="pipe_prod_snaps.npz",
                            dxp=1.84, ls="-", mk="D")}

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5,
    "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "lines.linewidth": 1.3, "figure.dpi": 150, "savefig.dpi": 200,
    "savefig.bbox": "tight", "axes.axisbelow": True,
})


# The profile figures show the finest grid alone, which is how the comparison
# is quoted: the coarser grids are a resolution study (Fig. 6), not a
# validation statement. The flag is kept as a constant rather than deleted
# because every panel below still reads it to choose its labels and limits.
FINAL_ONLY = True


def asset(name):
    """A file shipped beside this script, whatever the working directory."""
    return os.path.join(HERE, name)


def load(tag):
    g = dict(GRIDS[tag])
    g["P"] = profiles(dict(np.load(asset(g["stats"]), allow_pickle=True)))
    try:
        g["S"] = profiles(dict(np.load(asset(g["snaps"]), allow_pickle=True)))
    except FileNotFoundError:
        g["S"] = None
    g["qw"] = QW[tag]
    g["ttau"] = g["qw"] / U_TAU
    return g


def wall_value(r, t, interface=0.5):
    m = (r < interface) & np.isfinite(t)
    a, b = np.polyfit(r[m][-6:], t[m][-6:], 1)
    return a * interface + b


def fluid(P):
    return (P["count"][:, 0] > 0) & (P["r"] < 0.5)


def solid(P, depth=0.1):
    return (P["count"][:, 1] > 0) & (P["r"] > 0.5) & (P["r"] < 0.5 + 0.95 * depth)


def yplus(r):
    return (0.5 - r) / DELTA_V


def main():
    keys = ["production"] if FINAL_ONLY else list(GRIDS)
    G = {k: load(k) for k in keys}
    REF = {nm: ref_case(nm, PR) for nm in NAMES}
    r_ref = REF["c0"]["r"]
    yp_ref = yplus(r_ref)
    ttau_ref = 1.00330 / U_TAU
    tw_ref = {nm: wall_value(r_ref, REF[nm]["t"]) for nm in NAMES}
    mf = r_ref < 0.5

    # ---------------- Fig 1: hydrodynamics -----------------------------
    fig, ax = plt.subplots(1, 4, figsize=(13.5, 3.1))
    for tag, g in G.items():
        S = g["S"]
        if S is None:
            continue
        f = fluid(S)
        yp = yplus(S["r"][f])
        lab = (r"present, $\Delta^+$ 1.84 / $\Delta z^+$ 5.1" if FINAL_ONLY
               else rf"present, $\Delta^+={g['dxp']}$")
        ax[0].semilogx(yp, S["uz"][f, 0] / U_TAU, g["ls"], color="C0" if tag == "coarse" else "C3", label=lab)
        for k, (key, c) in enumerate((("uz_rms", "C0"), ("ur_rms", "C2"), ("up_rms", "C1"))):
            ax[1].semilogx(yp, S[key][f, 0] / U_TAU, g["ls"], color=c,
                           label=None if tag == "fine" else ["$u_z'$", "$u_r'$", r"$u_\theta'$"][k])
        ax[2].plot(S["r"][f] / R, S["urz_cov"][f, 0] / U_TAU**2, g["ls"], color="C0" if tag == "coarse" else "C3")
        tau = (-NU * np.gradient(S["uz"][f, 0], S["r"][f]) + S["urz_cov"][f, 0]) / U_TAU**2
        ax[2].plot(S["r"][f] / R, tau, g["ls"], color="C2", lw=1.0,
                   label="total stress" if tag == "coarse" else None)
        if tag == "coarse":
            ax[2].plot([], [], "-", color="C0", label=r"$\langle u_r'u_z'\rangle$")
    ax[0].semilogx(yp_ref[mf], REF["c0"]["uz"][mf] / U_TAU, "k-", lw=1.0, label="Neuhauser (NekRS)")
    yv = np.logspace(-0.3, 0.8, 20)
    ax[0].semilogx(yv, yv, ":", color="0.4", lw=0.9)
    yl = np.logspace(1.1, 2.3, 20)
    ax[0].semilogx(yl, np.log(yl) / 0.41 + 5.2, ":", color="0.4", lw=0.9)
    ax[0].text(2.0, 5.5, r"$U^+=y^+$", fontsize=7, color="0.35")
    ax[0].text(30, 11.5, r"$\frac{1}{\kappa}\ln y^+ + B$", fontsize=7, color="0.35")
    for key, c in (("uz_rms", "C0"), ("ur_rms", "C2"), ("ut_rms", "C1")):
        ax[1].semilogx(yp_ref[mf], REF["c0"][key][mf] / U_TAU, "k-", lw=0.9)
    ax[2].plot(r_ref[mf] / R, REF["c0"]["urz"][mf] / U_TAU**2, "k-", lw=0.9)
    ax[2].plot([0, 1], [0, 1], ":", color="0.4", lw=0.9)
    ax[2].text(0.45, 0.75, r"$u_\tau^2 r/R$ (exact)", fontsize=7, color="0.35", rotation=32)
    for tag, g in G.items():
        S = g["S"]
        if S is None:
            continue
        f = fluid(S)
        q = np.interp(S["r"][f], r_ref, REF["c0"]["uz"])
        ax[3].semilogx(yplus(S["r"][f]), 100 * (S["uz"][f, 0] - q) / q, g["ls"],
                       color="C0" if tag == "coarse" else "C3", label=rf"$\langle u_z\rangle$, $\Delta^+={g['dxp']}$")
        q = np.interp(S["r"][f], r_ref, REF["c0"]["uz_rms"])
        ax[3].semilogx(yplus(S["r"][f]), 100 * (S["uz_rms"][f, 0] - q) / q, g["ls"],
                       color="C2", label=rf"$u_z'$, $\Delta^+={g['dxp']}$")
    ax[3].axhline(0, color="k", lw=0.8)
    ax[3].set_ylim(-25, 15)
    for a, xl, yl_ in zip(ax, [r"$y^+$", r"$y^+$", "$r/R$", r"$y^+$"],
                          [r"$\langle u_z\rangle^+$", r"$u_i'/u_\tau$",
                           r"$\langle u_r'u_z'\rangle/u_\tau^2$", "deviation [%]"]):
        a.set_xlabel(xl); a.set_ylabel(yl_)
    ax[0].set_xlim(0.5, 200); ax[1].set_xlim(0.5, 200); ax[3].set_xlim(0.5, 200)
    ax[0].set_title("(a) mean velocity"); ax[1].set_title("(b) turbulence intensities")
    ax[2].set_title("(c) Reynolds and total stress"); ax[3].set_title("(d) deviation from the DNS")
    ax[0].legend(loc="upper left"); ax[1].legend(loc="upper left")
    ax[2].legend(loc="upper left"); ax[3].legend(loc="lower left")
    fig.tight_layout(); fig.savefig(asset("figures/fig1_hydrodynamics.png")); plt.close(fig)
    print("   fig1_hydrodynamics.png")

    # ---------------- Fig 2/3: mean temperature and fluctuation --------
    for fignum, (key, rkey, ylab, title) in enumerate(
            (("_mean", "t", r"$\langle\theta\rangle^+$", "mean temperature"),
             ("_rms", "t_rms", r"$\theta'_{rms}/\theta_\tau$", "temperature fluctuation")), start=2):
        fig, ax = plt.subplots(1, 1 + len(keys),
                               figsize=(4.4 * (1 + len(keys)), 3.4))
        for nm in NAMES:
            q = REF[nm][rkey][mf] / ttau_ref
            if key == "_mean":
                q = (REF[nm]["t"][mf] - tw_ref[nm]) / ttau_ref
            ax[0].semilogx(yp_ref[mf], q, "-", color=COLOR[nm], lw=0.9, label=LABEL[nm])
            g0 = keys[-1]
            P = G[g0]["P"]; f = fluid(P)
            v = P[f"{nm}{key}"][f, 0] / G[g0]["ttau"]
            if key == "_mean":
                v = (P[f"{nm}{key}"][f, 0] - wall_value(P["r"][f], P[f"{nm}{key}"][f, 0])) / G[g0]["ttau"]
            ax[0].semilogx(yplus(P["r"][f]), v, "o", color=COLOR[nm], ms=2.0, mfc="none", mew=0.6)
        ax[0].set_xlim(0.5, 200); ax[0].set_xlabel(r"$y^+$"); ax[0].set_ylabel(ylab)
        ax[0].set_title(f"(a) {title}: lines DNS, symbols present"
                        + (r" ($\Delta^+$ 1.84, $\Delta z^+$ 5.1)" if FINAL_ONLY
                           else r" ($\Delta^+=2.94$)"))
        ax[0].legend(ncol=2, loc="upper left" if key == "_mean" else "lower right",
                     framealpha=0.92, fontsize=7)
        for j, (tag, g) in enumerate((k, GRIDS[k]) for k in keys):
            P = G[tag]["P"]; f = fluid(P)
            for nm in NAMES:
                v = P[f"{nm}{key}"][f, 0] / G[tag]["ttau"]
                qq = np.interp(P["r"][f], r_ref, REF[nm][rkey]) / ttau_ref
                if key == "_mean":
                    v = (P[f"{nm}{key}"][f, 0] - wall_value(P["r"][f], P[f"{nm}{key}"][f, 0])) / G[tag]["ttau"]
                    qq = np.interp(P["r"][f], r_ref, REF[nm]["t"] - tw_ref[nm]) / ttau_ref
                ax[1 + j].semilogx(yplus(P["r"][f]), 100 * (v - qq) / qq, "-",
                                   color=COLOR[nm], lw=1.0, label=LABEL[nm])
            ax[1 + j].axhline(0, color="k", lw=0.8)
            ax[1 + j].set_xlim(0.5, 200)
            ax[1 + j].set_ylim(*((-8, 8) if key == "_mean" else (-12, 22))
                               if FINAL_ONLY else
                               (-30, 25 if key == "_mean" else 60))
            ax[1 + j].set_xlabel(r"$y^+$"); ax[1 + j].set_ylabel("deviation [%]")
            ax[1 + j].set_title(f"({'bcde'[j]}) $\\Delta^+={g['dxp']}$, "
                                f"$\\Delta z^+$="
                                f"{dict(coarse=14.2, fine=10.1)
                                   .get(tag, 5.1)}")
        ax[1].legend(ncol=2, fontsize=6.5, loc="lower left")
        fig.tight_layout()
        fig.savefig(asset(f"figures/fig{fignum}_"
                          f"{'mean' if key == '_mean' else 'fluctuation'}.png"))
        plt.close(fig)
        print(f"   fig{fignum}_{'mean' if key=='_mean' else 'fluctuation'}.png")

    # ---------------- Fig 4: the conjugate signature -------------------
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.2))
    # NORMALISED AT (r-R)/d = 0.1, NOT AT THE INTERFACE. The interface bin is
    # a cut cell whose penalization blend under-represents theta', so using it
    # as the denominator inflates the whole profile by 5-9 % -- on every grid,
    # which is what gave that away. One cell deeper the decay is reproduced to
    # 1-2 %.
    DEPTH0 = 0.1
    for nm in ["c0", "c1", "c2", "c3", "c4", "isof"]:
        rs = REF[nm]["r"]; ms = (rs > 0.5) & (rs < 0.6)
        iface = np.interp(0.5 + DEPTH0 * 0.1, rs, REF[nm]["t_rms"])
        ax[0].plot((rs[ms] - 0.5) / 0.1, REF[nm]["t_rms"][ms] / iface, "-",
                   color=COLOR[nm], lw=0.9, label=LABEL[nm])
        for tg, mk, ms_ in ((keys[-1], "o", 2.4),) if FINAL_ONLY else \
                (("coarse", "o", 2.2), ("fine", "s", 2.0)):
            P = G[tg]["P"]; sm = solid(P)
            xs = (P["r"][sm] - 0.5) / 0.1
            i0 = np.argmin(np.abs(xs - DEPTH0))
            ax[0].plot(xs, P[f"{nm}_rms"][sm, 1] / P[f"{nm}_rms"][sm, 1][i0],
                       mk, color=COLOR[nm], ms=ms_, mfc="none", mew=0.6)
    ax[0].set_xlabel(r"$(r-R)/d$")
    ax[0].set_ylabel(r"$\theta'_{rms}\,/\,\theta'_{rms}(0.1d)$")
    ax[0].set_title("(a) through the solid shell"
                    + ("" if FINAL_ONLY else r" (o: $\Delta^+$ 2.94, s: 1.84)"))
    ax[0].legend(fontsize=6.8)
    ax[0].set_ylim(0.25, 1.15); ax[0].set_xlim(0.05, 1.0)

    Ks = {"mbc": 0, "c4": 1, "c0": 2, "c3": 3, "isof": 4}   # categorical
    # Evaluated just OUTSIDE the cut cells: the outermost fluid bin holds a
    # penalization blend over a cell that straddles the wall, which is not
    # like-for-like against a body-fitted DNS (and is where the profiles are
    # steepest). At y+ ~ 7 every case agrees to 0.2-1 %.
    P = G[keys[-1]]["P"]; f = fluid(P)
    fl = np.flatnonzero(f)[np.argmin(np.abs(P["r"][f] - 0.48))]; rr = P["r"][fl]
    Pf = G[keys[-1]]["P"]; ff = fluid(Pf)
    ifl = np.argmin(np.abs(Pf["r"][ff] - rr))
    pres = [("prod", P, fl, "o", r"present, $\Delta^+$ 1.84 / $\Delta z^+$ 5.1")] \
        if FINAL_ONLY else \
        [("coarse", P, fl, "o", r"present $\Delta^+=2.94$"),
         ("fine", Pf, np.flatnonzero(ff)[ifl], "s", r"present $\Delta^+=1.84$")]
    for tag, P_, idx, mk, lab in pres:
        b = P_[f"c0_rms"][idx, 0]
        ax[1].plot([Ks[n] for n in Ks], [P_[f"{n}_rms"][idx, 0] / b for n in Ks],
                   mk, ms=4, mfc="none", label=lab)
    b = np.interp(rr, r_ref, REF["c0"]["t_rms"])
    ax[1].plot([Ks[n] for n in Ks],
               [np.interp(rr, r_ref, REF[n]["t_rms"]) / b for n in Ks], "k^-", ms=4, lw=0.9,
               label="Neuhauser")
    ax[1].set_xticks(list(Ks.values()))
    ax[1].set_xticklabels(["MBC\n$K\\to0$", "0.25", "1", "4", "IF\n$K\\to\\infty$"])
    ax[1].set_xlabel("effusivity ratio $K$")
    ax[1].set_ylabel(r"$\theta'_{rms}(R)/\theta'_{rms,c0}(R)$")
    ax[1].set_title(rf"(b) the effusivity bracket, $y^+={(0.5-rr)/DELTA_V:.1f}$")
    ax[1].legend()

    lam = {"c2": 0.5, "c0": 1.0, "c1": 2.0}
    for tag, P_, idx, mk, lab in pres:
        b = P_["c0_rms"][idx, 0]
        ax[2].plot(list(lam.values()), [P_[f"{n}_rms"][idx, 0] / b for n in lam],
                   mk, ms=4, mfc="none", label=lab)
    b = np.interp(rr, r_ref, REF["c0"]["t_rms"])
    ax[2].plot(list(lam.values()),
               [np.interp(rr, r_ref, REF[n]["t_rms"]) / b for n in lam], "k^-", ms=4, lw=0.9,
               label="Neuhauser")
    ax[2].set_xscale("log")
    ax[2].set_xticks([0.5, 1, 2]); ax[2].set_xticklabels(["0.5", "1", "2"])
    ax[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[2].set_xlabel(r"conductivity ratio $\lambda_{sf}=\kappa_s/\kappa_f$")
    ax[2].set_ylabel(r"$\theta'_{rms}(R)/\theta'_{rms,c0}(R)$")
    ax[2].set_title("(c) the conductivity sweep"); ax[2].legend()
    fig.tight_layout(); fig.savefig(asset("figures/fig4_conjugate_signature.png")); plt.close(fig)
    print("   fig4_conjugate_signature.png")

    # ---------------- Fig 5: turbulent heat fluxes ---------------------
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.2))
    for tag, g in G.items():
        S = g["S"]
        if S is None:
            continue
        f = fluid(S)
        col = "C0" if tag == "coarse" else "C3"
        ax[0].semilogx(yplus(S["r"][f]), S["c0_flux_ur"][f, 0] / (g["ttau"] * U_TAU),
                       g["ls"], color=col, label=rf"present $\Delta^+={g['dxp']}$")
        ax[1].semilogx(yplus(S["r"][f]), S["c0_flux_uz"][f, 0] / (g["ttau"] * U_TAU),
                       g["ls"], color=col)
    ax[0].semilogx(yp_ref[mf], REF["c0"]["tur"][mf] / (ttau_ref * U_TAU), "k-", lw=0.9,
                   label="Neuhauser")
    ax[1].semilogx(yp_ref[mf], REF["c0"]["tuz"][mf] / (ttau_ref * U_TAU), "k-", lw=0.9)
    ax[0].set_xlabel(r"$y^+$"); ax[0].set_ylabel(r"$\langle u_r'\theta'\rangle^+$")
    ax[1].set_xlabel(r"$y^+$"); ax[1].set_ylabel(r"$\langle u_z'\theta'\rangle^+$")
    ax[0].set_title("(a) radial turbulent heat flux"); ax[1].set_title("(b) axial turbulent heat flux")
    ax[0].set_xlim(0.5, 200); ax[1].set_xlim(0.5, 200); ax[0].legend()
    fig.tight_layout(); fig.savefig(asset("figures/fig5_heat_flux.png")); plt.close(fig)
    print("   fig5_heat_flux.png")

    # ---------------- Fig 6: grid convergence and the artifact ---------
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.3))
    # AGAINST THE AXIAL SPACING, which is what these quantities actually
    # converge in: the z-fine grid keeps the COARSE radial spacing and still
    # collapses onto the trend, while the "fine" grid -- refined radially but
    # barely axially -- sits with the coarse one.
    dxp = np.array([14.16, 10.11, 5.06])
    series = {r"bulk velocity deficit": np.array([2.42, 0.69, 0.18]),
              r"$\langle\theta\rangle^+$ deviation": np.array([6.0, 2.4, 0.3]),
              r"$\theta'^+$ deviation (axis)": np.array([16.6, 14.5, 6.0])}
    # the production grid sits on top of the z-fine point in Dz+, with half
    # the radial spacing -- the fluid-side errors barely move, which is the
    # cleanest statement that Dz+ is the controlling direction
    ax[0].loglog([5.06], [6.4], "D", ms=6, mfc="none", color="0.25",
                 label=r"production ($\Delta^+$ 1.84)")
    for (lab, v), mk in zip(series.items(), "os^d"):
        ax[0].loglog(dxp, v, mk + "-", ms=5, mfc="none", label=lab)
    for pwr, st in ((1, ":"), (2, "--")):
        ax[0].loglog(dxp, 6.0 * (dxp / dxp[0]) ** pwr, st, color="0.5", lw=0.9)
    ax[0].text(7.0, 3.4, r"$h^1$", color="0.4", fontsize=8)
    ax[0].text(7.0, 1.5, r"$h^2$", color="0.4", fontsize=8)
    ax[0].set_xticks(dxp); ax[0].set_xticklabels(["14.2", "10.1", "5.1"])
    ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[0].set_xlim(4.4, 16)
    ax[0].set_xlabel(r"$\Delta z^+$   (radial spacing: 2.94, 1.84, 2.94)")
    ax[0].set_ylabel("error measure [%]")
    ax[0].set_title("(a) convergence in the AXIAL spacing"); ax[0].legend(fontsize=6.8)

    NSEC = 64
    # The z-averaged cross-sections, reduced from the 52 MB solver statistics
    # files by make_caches.py -- the spectrum needs the 2D field, nothing else.
    RES = np.load(asset("grid_residual.npz"))
    # COARSE AND FINE ONLY, deliberately. grid_residual.npz also carries the
    # production grid, and it is NOT drawn here: the measurement band
    # 0.49 <= r < 0.4999 holds a different set of cell layers on each grid
    # (about 1.2 layers at Delta = 1.3/160, about 1.9 at 1.3/256), so the
    # amplitudes are not a convergence sequence and must not be read as one.
    # What the panel shows is that the Cartesian signature EXISTS and sits at
    # ~1 % of Delta T -- see the README for why its rate was withdrawn.
    for tag, col, off in (("coarse", "C0", -0.2), ("fine", "C3", 0.2)):
        x, y = RES[f"{tag}/x"], RES[f"{tag}/y"]
        X, Y = np.meshgrid(x, y, indexing="ij")
        rr_ = np.hypot(X - 0.65, Y - 0.65).ravel()
        ph = np.arctan2(Y - 0.65, X - 0.65).ravel()
        # A ring that straddles the shell/jacket boundary mixes shell cells
        # with INERT jacket cells held at a constant value, which makes its
        # "azimuthal mean" meaningless and its residual enormous -- 667 units
        # against 4-6 inside the isoflux shell. Excluded everywhere it is
        # measured or drawn.
        dxg = 1.3 / int(RES[f"{tag}/nx"])
        sv = RES[f"{tag}/sv"]
        res = (sv - azimuthal_mean(sv, np.hypot(X - 0.65, Y - 0.65), dxg)).ravel()
        dT = 17.2 * (QW[tag] / U_TAU)
        m = (rr_ >= 0.49) & (rr_ < 0.4999)
        sec = np.minimum(((ph[m] + np.pi) / (2 * np.pi) * NSEC).astype(int), NSEC - 1)
        v = np.array([np.nanmean(res[m][sec == i]) for i in range(NSEC)])
        A = np.abs(np.fft.rfft(v - np.nanmean(v))) / NSEC * 2 / dT * 100
        ax[1].bar(np.arange(1, 17) + off, A[1:17], width=0.4, color=col,
                  label=rf"$\Delta^+={GRIDS[tag]['dxp']}$")
        floor = np.mean([A[k] for k in range(1, 25) if k % 4])
        ax[1].axhline(floor, color=col, ls=":", lw=1.0)
        if tag == "coarse":
            cmap_data = (res, rr_, X, Y, dT)
    ax[1].set_xlabel("azimuthal mode $m$"); ax[1].set_ylabel(r"amplitude [% of $\Delta T$]")
    ax[1].set_title("(b) azimuthal spectrum at the wall\n"
                    "$m=4,8,12$: the Cartesian grid; dotted: sampling floor",
                    fontsize=8.5)
    ax[1].set_xticks([4, 8, 12, 16]); ax[1].legend()

    res, rr_, X, Y, dT = cmap_data
    n = int(np.sqrt(res.size))
    dxc = 1.3 / 160
    fld = np.where(rr_ < 0.6 - 1.5 * dxc, res, np.nan).reshape(n, n) / dT * 100
    lim = np.nanpercentile(np.abs(fld), 99)
    pc = ax[2].pcolormesh(X, Y, fld, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="nearest")
    th = np.linspace(0, 2 * np.pi, 400)
    for rad, st in ((0.5, "-"), (0.6, "--")):
        ax[2].plot(0.65 + rad * np.cos(th), 0.65 + rad * np.sin(th), st, color="k", lw=0.8)
    ax[2].set_aspect("equal"); ax[2].set_xlim(0.05, 1.25); ax[2].set_ylim(0.05, 1.25)
    ax[2].set_xlabel("$x$"); ax[2].set_ylabel("$y$")
    ax[2].set_title("(c) residual about the axisymmetric mean\n"
                    r"$\Delta^+=2.94$", fontsize=8.5)
    fig.colorbar(pc, ax=ax[2], label=r"% of $\Delta T$")
    fig.tight_layout(); fig.savefig(asset("figures/fig6_convergence.png")); plt.close(fig)
    print("   fig6_convergence.png")

    # ---------------- Fig 7: instantaneous cross-section ---------------
    # One z-plane of one production snapshot, cached by make_caches.py: a
    # snapshot is 5.6 GB and this figure draws 2 MB of it.
    PL = np.load(asset("pipe_prod_plane.npz"))
    fields = {n: PL[n].astype(np.float64) for n in
              ("un", "vn", "wn", "pn") + tuple(NAMES)}
    tcur = float(PL["t_current"])
    zpl = float(PL["z_plane"])
    cx = 0.5 * (PL["x"][:-1] + PL["x"][1:])
    cy = 0.5 * (PL["y"][:-1] + PL["y"][1:])
    X, Y = np.meshgrid(cx, cy)
    rr_ = np.hypot(X - 0.65, Y - 0.65)
    cs, sn = (X - 0.65) / rr_, (Y - 0.65) / rr_
    uc = 0.5 * (fields["un"] + np.roll(fields["un"], -1, axis=1))
    vc = 0.5 * (fields["vn"] + np.roll(fields["vn"], -1, axis=0))
    ur, ut = uc * cs + vc * sn, -uc * sn + vc * cs
    fluid_m = rr_ < 0.5

    def azimuthal_fluct(f):
        """The field minus its azimuthal mean at the same radius. This is what
        makes the conjugate signature visible: how far the fluctuation reaches
        into the shell IS the measurement."""
        # The ring width ADAPTS (pipe_stats.azimuthal_mean): a fixed width
        # narrower than a cell samples an aliased subset of azimuths, which
        # empties the bins at the axis -- the speckle this used to draw there
        # -- and lands on m = 4, 8, 12 where a Cartesian artifact would.
        return f - azimuthal_mean(f, rr_, float(cx[1] - cx[0]))

    fig, axs = plt.subplots(3, 4, figsize=(12.4, 9.6))
    th = np.linspace(0, 2 * np.pi, 400)

    def panel(a, fld, title, cmap, mask_solid=True, sym=True, pct=99.0):
        f = np.where(fluid_m, fld, np.nan) if mask_solid else fld
        if sym:
            lim = np.nanpercentile(np.abs(f), pct)
            kw = dict(vmin=-lim, vmax=lim)
        else:
            kw = dict(vmin=np.nanpercentile(f, 1), vmax=np.nanpercentile(f, pct))
        pc = a.pcolormesh(X, Y, f, cmap=cmap, shading="nearest", **kw)
        for rad, st, lw in ((0.5, "-", 0.9), (0.6, "--", 0.7)):
            a.plot(0.65 + rad * np.cos(th), 0.65 + rad * np.sin(th), st,
                   color="k", lw=lw)
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        a.set_xlim(0.03, 1.27); a.set_ylim(0.03, 1.27)
        a.set_title(title, fontsize=9)
        cb = fig.colorbar(pc, ax=a, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=6.5)

    panel(axs[0, 0], fields["wn"] / U_TAU, r"$u_z/u_\tau$ (axial)", "viridis",
          sym=False)
    panel(axs[0, 1], ur / U_TAU, r"$u_r/u_\tau$ (radial)", "RdBu_r")
    panel(axs[0, 2], ut / U_TAU, r"$u_\theta/u_\tau$ (azimuthal)", "RdBu_r")
    panel(axs[0, 3], fields["pn"] / U_TAU**2, r"$p/u_\tau^2$", "PuOr_r")
    # Each case is normalised by ITS OWN interface fluctuation, so the panels
    # are directly comparable in what they are meant to show -- how far the
    # fluctuation reaches into the shell. (A common scale would be dominated
    # by the isoflux case, whose near-insulating shell carries a temperature
    # drop ~100x the fluid's.) The inert jacket is masked: it holds a constant
    # value, so a fluctuation there is meaningless.
    # Normalised by each case's fluctuation at MID-RADIUS, where all seven
    # collapse onto one another (Fig. 3a): that puts the fluid structure on a
    # common scale and leaves the near-wall and shell differences -- which is
    # what the panels exist to show. Normalising by the INTERFACE value
    # instead would saturate the core of MBC and K = 0.25, whose interface
    # fluctuation is small by construction.
    Pf_ = G[keys[-1]]["P"]
    fsel = (Pf_["count"][:, 0] > 0) & (Pf_["r"] < 0.5)
    ffl = np.flatnonzero(fsel)[np.argmin(np.abs(Pf_["r"][fsel] - 0.25))]
    # ...and the shell/jacket boundary ring is excluded for the same reason
    # it is excluded from the artifact measurement.
    shown = rr_ < 0.6 - 1.5 * float(cx[1] - cx[0])
    for k, nm in enumerate(["c0", "c1", "c2", "c3", "c4", "mbc", "isof"]):
        scale_ = Pf_[f"{nm}_rms"][ffl, 0]
        fld = np.where(shown, azimuthal_fluct(fields[nm]) / scale_, np.nan)
        a = axs[1 + k // 4, k % 4]
        pc = a.pcolormesh(X, Y, fld, cmap="RdBu_r", shading="nearest",
                          vmin=-3, vmax=3)
        for rad, st, lw in ((0.5, "-", 0.9), (0.6, "--", 0.7)):
            a.plot(0.65 + rad * np.cos(th), 0.65 + rad * np.sin(th), st,
                   color="k", lw=lw)
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        a.set_xlim(0.03, 1.27); a.set_ylim(0.03, 1.27)
        a.set_title(rf"$\theta^\prime/\theta^\prime_{{rms}}(r{{=}}R/2)$   {LABEL[nm]}",
                    fontsize=9)
        cb = fig.colorbar(pc, ax=a, fraction=0.046, pad=0.02, extend="both")
        cb.ax.tick_params(labelsize=6.5)
    axs[2, 3].axis("off")
    axs[2, 3].text(0.0, 0.66,
                   "Instantaneous cross-section\n"
                   r"$\Delta^+ = 1.84$, " + f"z/D = {zpl:.2f}, "
                   + f"t = {tcur:.0f}" + r" $D/u_b$"
                   "\n\nsolid line: pipe wall, r = R\n"
                   "dashed: shell outer surface\n\n"
                   "velocity panels are masked in the\n"
                   "solid; scalar panels show fluid and\n"
                   "shell, each normalised by its own\n"
                   "mid-radius fluctuation, where all\n"
                   "seven collapse -- so what differs is\n"
                   "the wall and the shell: penetration\n"
                   "deep for K = 4, shallow for K = 0.25,\n"
                   "and exactly zero for MBC.\n\n"
                   "The outermost shell ring is cut:\n"
                   "it mixes shell with inert jacket.",
                   fontsize=8, va="top", color="0.35",
                   transform=axs[2, 3].transAxes)
    fig.tight_layout(); fig.savefig(asset("figures/fig7_instantaneous.png")); plt.close(fig)
    print("   fig7_instantaneous.png")
    return 0




def figure_blocks():
    """The BLOCK DECOMPOSITION of the production grid.

    It is SINGLE LEVEL -- no 2:1 refinement anywhere (1792 leaves, all level
    0) -- so these are the leaves as they stand, which are also the objects
    refinement would subdivide. Why refinement was not used: the conjugate
    cut-face coefficient is a SAME-LEVEL arm, so both the pipe wall and the
    band boundary 0.1 deeper would have to lie inside the finest level, which
    is nearly the whole annulus. See the README.
    """
    from matplotlib.collections import LineCollection
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    BK = np.load(asset("pipe_prod_blocks.npz"))
    blocks = BK["blocks"]
    nb = int(BK["block_nb"])
    node = {d: BK[d] for d in "xyz"}
    cx, cy = 0.65, 0.65

    def extent(b):
        ox, oy, oz, _l = b
        return ((node["x"][ox], node["x"][ox + nb]),
                (node["y"][oy], node["y"][oy + nb]),
                (node["z"][oz], node["z"][oz + nb]))

    def classify(xr, yr):
        """What the block holds, from the radii of its corners."""
        corners = [(x - cx, y - cy) for x in xr for y in yr]
        rs = [np.hypot(a, b) for a, b in corners]
        # nearest point of the rectangle to the axis
        dx0 = max(xr[0] - cx, 0, cx - xr[1]); dy0 = max(yr[0] - cy, 0, cy - yr[1])
        rmin, rmax = np.hypot(dx0, dy0), max(rs)
        if rmax <= 0.5:
            return "fluid"
        if rmin >= 0.6:
            return "jacket"
        if rmin < 0.5 < rmax:
            return "cut"
        return "shell"

    COL = {"fluid": "#0e6f7c", "cut": "#b06a12", "shell": "#d9a441",
           "jacket": "#9aa6b2"}
    LAB = {"fluid": "entirely fluid", "cut": "cut by the pipe wall",
           "shell": "shell and jacket", "jacket": "entirely inert jacket"}

    # two block layers deep, which is what makes the cubes read as cubes
    ZTOP = 2.0 * float(node["z"][nb] - node["z"][0])

    fig = plt.figure(figsize=(11.5, 4.6))
    ax0 = fig.add_subplot(1, 2, 1)
    seen = set()
    counts = {k: 0 for k in COL}
    for b in blocks:
        xr, yr, zr = extent(b)
        k = classify(xr, yr)
        if b[2] == 0:                      # one z layer for the cross-section
            ax0.add_patch(plt.Rectangle((xr[0], yr[0]), xr[1] - xr[0], yr[1] - yr[0],
                                        facecolor=COL[k], edgecolor="k", lw=0.8,
                                        alpha=0.45,
                                        label=LAB[k] if k not in seen else None))
            seen.add(k)
        counts[k] += 1
    th = np.linspace(0, 2 * np.pi, 400)
    for rad, st in ((0.5, "-"), (0.6, "--")):
        ax0.plot(cx + rad * np.cos(th), cy + rad * np.sin(th), st, color="k", lw=1.4)
    # the cell grid inside one block, to show what a block contains
    xr, yr, _ = extent(blocks[0])
    for t in np.linspace(xr[0], xr[1], nb + 1):
        ax0.plot([t, t], [yr[0], yr[1]], color="k", lw=0.15, alpha=0.5)
        ax0.plot([xr[0], xr[1]], [t, t], color="k", lw=0.15, alpha=0.5)
    ax0.set_aspect("equal"); ax0.set_xlim(0, 1.3); ax0.set_ylim(0, 1.3)
    ax0.set_xlabel("$x$"); ax0.set_ylabel("$y$")
    nbx = (len(node["x"]) - 1) // nb
    nbz = (len(node["z"]) - 1) // nb
    ax0.set_title(f"(a) block lattice, $\\Delta^+=1.84$: "
                  f"{nbx}×{nbx}×{nbz} blocks of ${nb}^3$\n"
                  "(one block drawn with its cells)")
    ax0.legend(loc="upper right", fontsize=7, framealpha=0.95)

    ax1 = fig.add_subplot(1, 2, 2, projection="3d")
    segs, cols = [], []
    for b in blocks:
        xr, yr, zr = extent(b)
        k = classify(xr, yr)
        if k == "jacket" or zr[0] > 1.5 * ZTOP / 2:  # two block layers, so cubes read as cubes
            continue
        c = [(xr[i], yr[j], zr[m]) for i in (0, 1) for j in (0, 1) for m in (0, 1)]
        for a, bb in ((0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
                      (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)):
            segs.append((c[a], c[bb])); cols.append(COL[k])
    ax1.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=0.45,
                                          alpha=0.8))
    zc = np.linspace(0, ZTOP, 2)
    TH, ZC = np.meshgrid(np.linspace(0, 2 * np.pi, 120), zc)
    ax1.plot_surface(cx + 0.5 * np.cos(TH), cy + 0.5 * np.sin(TH), ZC,
                     color="#0e6f7c", alpha=0.22, linewidth=0, shade=False)
    ax1.set_xlim(0, 1.3); ax1.set_ylim(0, 1.3); ax1.set_zlim(0, ZTOP)
    ax1.set_box_aspect((1, 1, max(0.8, 1.9 * ZTOP / 2.5)))
    ax1.set_xlabel("$x$", labelpad=-6); ax1.set_ylabel("$y$", labelpad=-6)
    ax1.set_zlabel("$z$", labelpad=-2)
    ax1.set_xticks([0, 0.65, 1.3]); ax1.set_yticks([0, 0.65, 1.3])
    ax1.set_zticks([0, ZTOP / 2, ZTOP])
    ax1.set_zticklabels([f"{v:.2f}" for v in (0, ZTOP / 2, ZTOP)])
    ax1.tick_params(labelsize=6.5, pad=-1)
    ax1.view_init(elev=20, azim=-62)
    ax1.set_title("(b) two block layers in 3D (jacket-only blocks removed):\n"
                  r"how blocks of $32^3$ cells tile a curved wall", fontsize=9)
    fig.tight_layout()
    fig.savefig(asset("figures/fig8_blocks.png"))
    plt.close(fig)
    print("   fig8_blocks.png  " + ", ".join(f"{LAB[k]}: {v}" for k, v in counts.items()))


if __name__ == "__main__":
    if "--blocks" in sys.argv:
        sys.exit(figure_blocks())
    sys.exit(main())
