#!/usr/bin/env python3
"""LES validation analysis (run after run_les.py).

Homogeneous cases (reference / uniform / slab): read the solver's time-averaged
per-y channel_stats (mean U + resolved stresses, per level, reassembled) and plot
the LES cases against the filtered-DNS reference (128^3 no-LES):
  les_profiles.png : U+ vs y+ (log law), u'/v'/w' rms, -<u'v'> vs y
  + a printed CORE ratio table LES/reference.
nut(y): time-averaged per-level from the stats-leg field snapshots; overlaid with
the 2:1 interface marked -- gate: nut steps by ~the physical filter-width ratio
(delta^2 ~ 4x across a 2:1 face), NO spike/band beyond it.
Patch (edges+corners): delegates to tools/patch_interface_stats.py (nut-aware),
patch run vs the uniform LES run as the matched base control.

Reference is UNFILTERED DNS, so LES resolved stresses sit slightly below the
reference peak by the SGS contribution (expected); the mean U and -<u'v'> (total
momentum balance) should match closely.
"""
import os, glob, subprocess, sys
import numpy as np
import h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RUNS = os.path.join(HERE, "runs")
RETAU = 180.0
U, V, W, UU, VV, WW, UV = range(7)
COL = {"reference": "k", "uniform": "tab:green", "slab": "tab:blue", "patch": "tab:red"}


def read_stats(case):
    """channel_stats per-y profile (level 0 + _l* concatenated), sorted by y."""
    p = os.path.join(RUNS, case, "stats", "channel_stats.h5")
    if not os.path.exists(p):
        return None
    ys, pr = [], []
    for lf in [p] + sorted(glob.glob(p.replace(".h5", "") + "_l*.h5")):
        with h5py.File(lf) as f:
            c = f["coord"][...]; P = f["profile"][...]; n = f["count"][...]
        m = n > 0
        ys.append(c[m]); pr.append(P[m])
    y = np.concatenate(ys); P = np.concatenate(pr, 0); o = np.argsort(y)
    return y[o], P[o]


def comps(P):
    urms = np.sqrt(np.clip(P[:, UU] - P[:, U] ** 2, 0, None))
    vrms = np.sqrt(np.clip(P[:, VV] - P[:, V] ** 2, 0, None))
    wrms = np.sqrt(np.clip(P[:, WW] - P[:, W] ** 2, 0, None))
    uv = -(P[:, UV] - P[:, U] * P[:, V])
    return P[:, U], urms, vrms, wrms, uv


def nut_profile(case):
    """time-averaged nut(y) over the stats-leg snapshots, reassembled to the
    finest y-lattice (coarse leaves repeated). Returns (yc_finest, nut_y) or None."""
    snaps = sorted(glob.glob(os.path.join(RUNS, case, "stats", f"{case}_*.h5")))
    if not snaps:
        return None
    with h5py.File(snaps[0]) as f:
        if "nut" not in f:
            return None
        nb = int(f.attrs["block_nb_x"]); bl = f["blocks"][...]; y = f["y"][...]
    lmax = int(bl[:, 3].max()); nyf = (len(y) - 1) * (2 ** lmax)
    acc = np.zeros(nyf); cnt = np.zeros(nyf)
    for s in snaps:
        with h5py.File(s) as f:
            n = f["nut"][...]
        for bid, (ox, oy, oz, lev) in enumerate(bl):
            ff = 2 ** (lmax - lev)
            for jj in range(nb):
                g0 = (oy + jj) * ff
                val = n[bid, :, jj, :].sum()
                acc[g0:g0 + ff] += val; cnt[g0:g0 + ff] += nb * nb
    prof = acc / np.maximum(cnt, 1)
    yf = np.interp(np.linspace(0, len(y) - 1, nyf + 1), np.arange(len(y)), y)
    return 0.5 * (yf[:-1] + yf[1:]), prof


def homogeneous_figure():
    data = {c: read_stats(c) for c in ("reference", "uniform", "slab")}
    data = {k: v for k, v in data.items() if v is not None}
    if not data:
        print("no channel_stats found -- run run_les.py first"); return
    C = {k: comps(data[k][1]) for k in data}
    Y = {k: data[k][0] for k in data}

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    # U+ vs y+ (log law), lower half only
    a = ax[0, 0]
    yp = np.logspace(0, np.log10(RETAU), 100)
    a.semilogx(yp, yp, "k:", lw=0.8, label="U+=y+")
    a.semilogx(yp[yp > 11], np.log(yp[yp > 11]) / 0.41 + 5.2, "k--", lw=0.8, label="log law")
    for k in data:
        y = Y[k]; lo = y <= 1.0
        a.semilogx(y[lo] * RETAU, C[k][0][lo], color=COL[k], lw=1.4, label=k)
    a.set_xlabel("y+"); a.set_ylabel("U+"); a.set_title("mean velocity (log law)"); a.legend(fontsize=8)

    for axi, (ci, t) in zip([ax[0, 1], ax[0, 2], ax[1, 0], ax[1, 1]],
                            [(1, "u' rms"), (2, "v' rms"), (3, "w' rms"), (4, "-<u'v'>")]):
        for k in data:
            axi.plot(C[k][ci], Y[k], color=COL[k], lw=1.4, label=k)
        axi.set_xlabel(t); axi.set_ylabel("y"); axi.set_ylim(0, 2); axi.set_title(t)
    ax[0, 1].legend(fontsize=8)

    # nut(y) overlay
    a = ax[1, 2]
    for k in ("uniform", "slab"):
        np_ = nut_profile(k)
        if np_ is None:
            continue
        yc, nu = np_
        a.plot(nu / (1 / RETAU), yc, color=COL[k], lw=1.4, label=f"{k} nut/nu_mol")
    a.set_xlabel("nut / nu_mol"); a.set_ylabel("y"); a.set_ylim(0, 2)
    a.set_title("eddy viscosity (interface step ~ delta^2)"); a.legend(fontsize=8)

    fig.suptitle("LES vs filtered-DNS reference (Re_tau 180, coarse 64^3 + WALE)")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(HERE, "les_profiles.png"); fig.savefig(out, dpi=120)
    print("wrote", out)

    # CORE ratio table LES/reference (interp ref to LES y), 0.3<y<1.7
    if "reference" in data:
        yr, Pr = data["reference"]; Cr = comps(Pr)
        print("\nCORE (0.3<y<1.7) mean ratio LES/reference:    U      u'     v'     w'    -<u'v'>")
        for k in ("uniform", "slab"):
            if k not in data:
                continue
            y = Y[k]; m = (y > 0.3) & (y < 1.7)
            r = []
            for ci in (0, 1, 2, 3, 4):
                num = C[k][ci][m]
                den = np.interp(y[m], yr, Cr[ci])
                r.append(np.mean(num / np.where(den == 0, np.nan, den)))
            print(f"  {k:9s}  {r[0]:.3f}  {r[1]:.3f}  {r[2]:.3f}  {r[3]:.3f}  {r[4]:.3f}")


def patch_band():
    pglob = os.path.join(RUNS, "patch", "stats", "patch_*.h5")
    bglob = os.path.join(RUNS, "uniform", "stats", "uniform_*.h5")
    if not glob.glob(pglob) or not glob.glob(bglob):
        print("\npatch band metric: need patch + uniform snapshots -- skipped")
        return
    print("\n== patch edge/corner band metric (patch vs uniform LES control) ==")
    subprocess.run([sys.executable, os.path.join(ROOT, "tools/patch_interface_stats.py"),
                    "--patch", pglob, "--base", bglob], check=False)


if __name__ == "__main__":
    homogeneous_figure()
    patch_band()
