#!/usr/bin/env python3
"""LES<->IBM coupling validation -- developed-statistics analysis (gates 3 & 4).

Reads the stats-leg channel_stats + nut snapshots (run_ibm_les.py) and produces:

  ibm_les_profiles.png:
    * Gate 3 -- mean U+(y+) law of the wall (a_wale vs b_none control), with the
      grid-aligned LES channel (../les/, no IBM) and the linear/log lines. The IBM
      wall sits at y=0.259375 (bottom) / 2.259375 (top); y+ is measured from it.
    * resolved stresses u'/v'/w'/-<u'v'>(y) for a_wale vs b_none.
    * nut(y) for a_wale and c_refine -- Gate 4: across the 2:1 interface at the wall
      (case c) nut must STEP by ~the filter-width ratio (delta^2, ~4x for a 2:1 face)
      with NO spurious band (a smooth step, not an overshoot ridge).

channel_stats here is per-level (x,z homogeneous), reassembled to the finest y
lattice; the solid rows (velocity penalized to ~0) are masked out by count/coef.

Usage:  python3 ibm_les_stats.py
"""
from __future__ import annotations
import glob
import os

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
LESDIR = os.path.join(HERE, "..", "les", "runs")     # grid-aligned LES reference
RETAU = 180.0
Y_LO, Y_HI = 0.259375, 2.259375                      # IBM wall locations
U, V, W, UU, VV, WW, UV = range(7)
COL = {"a_wale": "tab:red", "b_none": "tab:gray", "c_refine": "tab:blue",
       "uniform": "tab:green"}


def read_stats(run_dir, prefix):
    """channel_stats profile (level 0 + _l* concatenated), sorted by y."""
    p = os.path.join(run_dir, "stats", f"{prefix}.h5")
    if not os.path.exists(p):
        return None
    ys, pr = [], []
    for lf in [p] + sorted(glob.glob(p.replace(".h5", "") + "_l*.h5")):
        with h5py.File(lf) as f:
            c = f["coord"][...]; P = f["profile"][...]; n = f["count"][...]
        m = n > 0
        ys.append(c[m]); pr.append(P[m])
    if not ys:
        return None
    y = np.concatenate(ys); P = np.concatenate(pr, 0); o = np.argsort(y)
    return y[o], P[o]


def comps(P):
    urms = np.sqrt(np.clip(P[:, UU] - P[:, U] ** 2, 0, None))
    vrms = np.sqrt(np.clip(P[:, VV] - P[:, V] ** 2, 0, None))
    wrms = np.sqrt(np.clip(P[:, WW] - P[:, W] ** 2, 0, None))
    uv = -(P[:, UV] - P[:, U] * P[:, V])
    return P[:, U], urms, vrms, wrms, uv


def nut_profile(run_dir, prefix):
    """time+plane-mean nut(y) reassembled to the finest y-lattice."""
    snaps = [s for s in sorted(glob.glob(os.path.join(run_dir, "stats", f"{prefix}_*.h5")))
             if "stats" not in os.path.basename(s)]
    if not snaps:
        return None
    with h5py.File(snaps[0]) as f:
        if "nut" not in f:
            return None
        nb = f["nut"].shape[1]; bl = f["blocks"][...]; y = f["y"][...]
    lmax = int(bl[:, 3].max()); nyf = (len(y) - 1) * (2 ** lmax)
    acc = np.zeros(nyf); cnt = np.zeros(nyf)
    for s in snaps:
        with h5py.File(s) as f:
            n = f["nut"][...]
        for bid, (ox, oy, oz, lev) in enumerate(bl):
            ff = 2 ** (lmax - lev)
            for jj in range(nb):
                g0 = (oy + jj) * ff
                acc[g0:g0 + ff] += n[bid, :, jj, :].sum(); cnt[g0:g0 + ff] += nb * nb
    prof = acc / np.maximum(cnt, 1)
    yf = np.interp(np.linspace(0, len(y) - 1, nyf + 1), np.arange(len(y)), y)
    return 0.5 * (yf[:-1] + yf[1:]), prof


def yplus_lower(y):
    """y+ from the nearer IBM wall (lower half measured from Y_LO)."""
    return (y - Y_LO) * RETAU


def interface_ys(coef_blocks=os.path.join(HERE, "ibm_coeff_blocks.h5")):
    """y of the two 2:1 interfaces (coarse-core / fine-band boundaries), read
    from the refine_body leaf table: the coarse (level-0) leaves span a central
    y band; the interfaces are at its lower and upper physical edges."""
    if not os.path.exists(coef_blocks):
        return (0.75, 1.75)
    with h5py.File(coef_blocks) as f:
        bl = f["blocks"][...]; ly = float(f.attrs["ly"]); ny = int(f.attrs["ny"])
    dy = ly / ny
    coarse = bl[bl[:, 3] == 0]
    nb = 8
    ylo = coarse[:, 1].min() * dy            # origin (base cells) -> physical y
    yhi = (coarse[:, 1].max() + nb) * dy
    return (ylo, yhi)


def main():
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    IF_LO, IF_HI = interface_ys()        # 2:1 interface y (coarse-core edges)

    # --- Gate 3: mean U+(y+) law of the wall ---
    a = ax[0, 0]
    yp = np.logspace(-0.5, 2.4, 50)
    a.plot(yp, yp, "k:", lw=0.8, label="U+=y+")
    a.plot(yp[yp > 8], 2.44 * np.log(yp[yp > 8]) + 5.0, "k--", lw=0.8,
           label="2.44 ln y+ +5")
    for case, prefix in (("a_wale", "channel_ibm_stats"),
                         ("b_none", "channel_ibm_stats")):
        st = read_stats(os.path.join(RUNS, case), prefix)
        if st is None:
            continue
        y, P = st
        lo = y < 0.5 * (Y_LO + Y_HI)
        a.semilogx(yplus_lower(y[lo]), P[lo, U], color=COL[case], lw=1.5, label=case)
    # grid-aligned LES channel (no IBM) reference if present
    ref = read_stats(os.path.join(LESDIR, "uniform"), "channel_stats")
    if ref is not None:
        y, P = ref
        lo = y < 1.0
        a.semilogx(y[lo] * RETAU, P[lo, U], color=COL["uniform"], lw=1.2, ls="-.",
                   label="grid-aligned LES (../les)")
    a.set_xlabel("y+ (from IBM wall)"); a.set_ylabel("U+")
    a.set_title("Gate 3: mean velocity / law of the wall"); a.legend(fontsize=8)
    a.set_xlim(0.5, 250)

    # --- resolved stresses ---
    a = ax[0, 1]
    for case, prefix in (("a_wale", "channel_ibm_stats"), ("b_none", "channel_ibm_stats")):
        st = read_stats(os.path.join(RUNS, case), prefix)
        if st is None:
            continue
        y, P = st
        u, ur, vr, wr, uv = comps(P)
        a.plot(y, ur, color=COL[case], lw=1.4, label=f"{case} u'")
        a.plot(y, vr, color=COL[case], lw=1.0, ls="--")
        a.plot(y, wr, color=COL[case], lw=1.0, ls=":")
        a.plot(y, uv, color=COL[case], lw=1.0, ls="-.")
    for x in (Y_LO, Y_HI):
        a.axvline(x, color="k", lw=0.5, alpha=0.3)
    a.set_xlabel("y"); a.set_ylabel("rms (solid u', dash v', dot w', dashdot -<u'v'>)")
    a.set_title("resolved stresses (a_wale vs b_none)"); a.legend(fontsize=8)

    # --- Gate 4: nut(y) with interface step ---
    a = ax[1, 0]
    for case in ("a_wale", "c_refine"):
        npf = nut_profile(os.path.join(RUNS, case), "channel_ibm")
        if npf is None:
            continue
        yc, nu = npf
        a.plot(nu / (1.0 / RETAU), yc, color=COL[case], lw=1.4, label=case)
    for x in (Y_LO, Y_HI):
        a.axhline(x, color="k", lw=0.5, alpha=0.3)
    # mark the 2:1 interfaces of the refine_body case (block-row 2 and 7 boundaries)
    a.axhline(IF_LO, color="tab:blue", lw=0.6, ls=":", alpha=0.6)
    a.axhline(IF_HI, color="tab:blue", lw=0.6, ls=":", alpha=0.6)
    a.set_xlabel("nut / nu_mol"); a.set_ylabel("y")
    a.set_title("Gate 4: nut(y) -- step (not band) across the 2:1 interface")
    a.legend(fontsize=8); a.set_ylim(0, 2.5)

    # --- nut zoom near the bottom wall+interface ---
    a = ax[1, 1]
    for case in ("a_wale", "c_refine"):
        npf = nut_profile(os.path.join(RUNS, case), "channel_ibm")
        if npf is None:
            continue
        yc, nu = npf
        m = yc < 0.8
        a.plot(nu[m] / (1.0 / RETAU), yc[m], color=COL[case], lw=1.4, marker=".",
               ms=3, label=case)
    a.axhline(Y_LO, color="k", lw=0.6, alpha=0.4, label="IBM wall")
    a.axhline(IF_LO, color="tab:blue", lw=0.6, ls=":", alpha=0.6, label="2:1 interface")
    a.set_xlabel("nut / nu_mol"); a.set_ylabel("y")
    a.set_title("nut zoom: wall + interface (no band)"); a.legend(fontsize=8)
    a.set_ylim(0.2, IF_LO + 0.25)

    out = os.path.join(HERE, "ibm_les_profiles.png")
    fig.tight_layout(); fig.savefig(out, dpi=120)
    print(f"wrote {out}")

    # printed gate-4 metric: nut step ratio across the interface for c_refine
    npf = nut_profile(os.path.join(RUNS, "c_refine"), "channel_ibm")
    if npf is not None:
        yc, nu = npf
        fine = (yc > IF_LO - 0.2) & (yc < IF_LO - 0.01)
        coarse = (yc > IF_LO + 0.01) & (yc < IF_LO + 0.2)
        if fine.any() and coarse.any():
            rf, rc = nu[fine].max(), nu[coarse].max()
            print(f"Gate 4: nut fine-side max {rf/(1/RETAU):.3f}, coarse-side max "
                  f"{rc/(1/RETAU):.3f}, step ratio coarse/fine = {rc/max(rf,1e-30):.2f} "
                  f"(physical delta^2 ~ 4x; a smooth step, NOT a band)")


if __name__ == "__main__":
    main()
