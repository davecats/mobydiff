#!/usr/bin/env python3
"""2:1-interface turbulent-validation figures, from the developed/runs stats +
final fields. Produces, into this directory:
  * uniform128_4way.png      -- u'/v'/w'/-<u'v'> vs y for reflux_on, reflux_off,
                                uniform-256 (reference), uniform-128 (coarse control)
  * uniform128_isolation.png -- core-zoom u'/v'/w' for reference/uniform128/reflux_off
                                (does reflux_off track uniform128 -> resolution, or
                                drop below -> interface loss?), + a printed ratio table
  * xsection_xy.png / xsection_zy.png -- fluctuation cross-sections, uniform vs
                                reflux_off, all variables, interfaces marked

Run names are the runs/<name>/stats subdirs; missing runs are skipped. Reusable
for the edge/corner + LES validation (point at the new runs)."""
import os, glob
import numpy as np
import h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "runs")
IFACE = [0.6205, 1.3795]
U, V, W, UU, VV, WW, UV = range(7)
VARS = ["un", "vn", "wn", "pn"]; LAB = ["u", "v", "w", "p"]


def read_stats(name):
    p = os.path.join(BASE, name, "stats", "channel_stats.h5")
    if not os.path.exists(p):
        return None
    ys, pr = [], []
    for lf in [p] + sorted(glob.glob(p.replace(".h5", "") + "_l*.h5")):
        with h5py.File(lf) as f:
            c = f["coord"][...]; P = f["profile"][...]; n = f["count"][...]
        m = n > 0; ys.append(c[m]); pr.append(P[m])
    y = np.concatenate(ys); P = np.concatenate(pr, 0); o = np.argsort(y)
    return y[o], P[o]


def comp(P):
    return (np.sqrt(np.clip(P[:, UU] - P[:, U] ** 2, 0, None)),
            np.sqrt(np.clip(P[:, VV] - P[:, V] ** 2, 0, None)),
            np.sqrt(np.clip(P[:, WW] - P[:, W] ** 2, 0, None)),
            -(P[:, UV] - P[:, U] * P[:, V]))


def stats_figures():
    runs = {n: read_stats(n) for n in ["reflux_on", "reflux_off", "reference", "uniform128"]}
    runs = {k: v for k, v in runs.items() if v is not None}
    C = {k: comp(runs[k][1]) for k in runs}
    Y = {k: runs[k][0] for k in runs}

    if {"reflux_off", "reference", "uniform128"} <= set(runs):
        yctrl = Y["uniform128"]; core = (yctrl > 0.7) & (yctrl < 1.3)
        def ratio(num, den, ci):
            a = np.interp(yctrl, Y[num], C[num][ci]); b = np.interp(yctrl, Y[den], C[den][ci])
            return np.mean(a[core] / b[core])
        print("CORE (0.7<y<1.3) mean fluctuation ratios:        u'     v'     w'")
        for num, den in [("reflux_off", "reference"), ("uniform128", "reference"),
                         ("reflux_off", "uniform128")]:
            r = [ratio(num, den, ci) for ci in range(3)]
            tag = "  <-- ISOLATION (==1 => interface clean)" if den == "uniform128" else ""
            print(f"  {num+' / '+den:24s}  {r[0]:.3f}  {r[1]:.3f}  {r[2]:.3f}{tag}")

    col = {"reflux_on": "tab:red", "reflux_off": "tab:blue", "reference": "k", "uniform128": "tab:green"}
    order = [k for k in ["reference", "uniform128", "reflux_off", "reflux_on"] if k in runs]
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    for a, (t, idx) in zip(ax.flat, [("u' rms", 0), ("v' rms", 1), ("w' rms", 2), ("-<u'v'>", 3)]):
        for k in order:
            a.plot(C[k][idx], Y[k], color=col[k], lw=1.3, label=k, alpha=0.85)
        for yi in IFACE: a.axhline(yi, color="gray", ls="--", lw=0.6, alpha=0.6)
        a.set_title(t); a.set_ylabel("y"); a.set_ylim(0, 2)
    ax[0, 0].legend(fontsize=8)
    fig.suptitle("4-way: reflux on/off vs uniform-256 (ref) vs uniform-128 (coarse-core control)")
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(os.path.join(HERE, "uniform128_4way.png"), dpi=120)
    print("wrote uniform128_4way.png")

    if {"reflux_off", "reference", "uniform128"} <= set(runs):
        fig, ax = plt.subplots(1, 3, figsize=(13, 5), sharey=True)
        for a, (t, idx) in zip(ax, [("u' rms", 0), ("v' rms", 1), ("w' rms", 2)]):
            for k in ["reference", "uniform128", "reflux_off"]:
                a.plot(C[k][idx], Y[k], color=col[k], lw=1.6, marker='.', ms=3, label=k, alpha=0.85)
            for yi in IFACE: a.axhline(yi, color="gray", ls="--", lw=0.6)
            a.set_title(t); a.set_ylim(0.55, 1.45)
        ax[0].set_ylabel("y (core zoom)"); ax[0].legend(fontsize=9)
        fig.suptitle("Isolation: reflux_off CORE tracks uniform-128 (resolution) or drops below (interface loss)?")
        fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(os.path.join(HERE, "uniform128_isolation.png"), dpi=120)
        print("wrote uniform128_isolation.png")


def reassemble(fn):
    """block-table field -> global [z,y,x]. origins are LEVEL-l cell units;
    fine start = origin * 2^(maxlev-l). Returns (dict var->[z,y,x], y-nodes)."""
    with h5py.File(fn) as f:
        nb = int(f.attrs["block_nb_x"]); blocks = f["blocks"][...]; yb = f["y"][...]
        D = {v: f[v][...] for v in VARS}
    maxlev = int(blocks[:, 3].max()); nx, ny, nz = 256, 128, 256
    G = {v: np.full((nz, ny, nx), np.nan) for v in VARS}
    yl = yb.copy()
    for _ in range(maxlev):
        mid = 0.5 * (yl[:-1] + yl[1:]); new = np.empty(2 * len(yl) - 1)
        new[0::2] = yl; new[1::2] = mid; yl = new
    for bid, (ox, oy, oz, l) in enumerate(blocks):
        s = 2 ** (maxlev - int(l)); x0, y0, z0 = int(ox) * s, int(oy) * s, int(oz) * s
        for v in VARS:
            ex = np.repeat(np.repeat(np.repeat(D[v][bid], s, 0), s, 1), s, 2)
            G[v][z0:z0 + nb * s, y0:y0 + nb * s, x0:x0 + nb * s] = ex
    return G, yl


def xsection_figures():
    uni = os.path.join(BASE, "reference/stats/channel_field_80001.h5")
    ref = os.path.join(BASE, "reflux_off/stats/channel_field_80001.h5")
    if not (os.path.exists(uni) and os.path.exists(ref)):
        print("xsection: need reference + reflux_off final fields -- skipped"); return
    with h5py.File(uni) as f:
        UU_ = {v: f[v][0] for v in VARS}; yun = f["y"][...]
    GG, yrn = reassemble(ref)
    for cut, (name, title, xn) in enumerate([("xy", "x-y (mid z)", "x"), ("zy", "z-y (mid x)", "z")]):
        fig, ax = plt.subplots(4, 2, figsize=(11, 13), sharex=True, sharey=True)
        for r, (v, lab) in enumerate(zip(VARS, LAB)):
            if cut == 0: pu = UU_[v][128]; pr = GG[v][128]
            else:        pu = UU_[v][:, :, 128].T; pr = GG[v][:, :, 128].T
            fu = pu - pu.mean(1, keepdims=True); fr = pr - pr.mean(1, keepdims=True)
            lim = np.percentile(np.abs(np.concatenate([fu, fr])), 99)
            for c, (dat, yn, tag) in enumerate([(fu, yun, "uniform"), (fr, yrn, "reflux-off")]):
                im = ax[r, c].pcolormesh(np.arange(dat.shape[1] + 1), yn, dat,
                                         cmap="RdBu_r", vmin=-lim, vmax=lim, shading="flat")
                for yi in IFACE: ax[r, c].axhline(yi, color="k", lw=0.7, ls="--", alpha=0.7)
                if r == 0: ax[r, c].set_title(tag)
                if c == 0: ax[r, c].set_ylabel(f"{lab}'   y")
                if r == 3: ax[r, c].set_xlabel(f"{xn} index")
                plt.colorbar(im, ax=ax[r, c], fraction=0.046, pad=0.02)
        fig.suptitle(f"{title} cross-section: fluctuations, uniform vs reflux-off (dashed = 2:1 interface)")
        fig.tight_layout(rect=[0, 0, 1, 0.985]); fig.savefig(os.path.join(HERE, f"xsection_{name}.png"), dpi=110)
        print(f"wrote xsection_{name}.png")


if __name__ == "__main__":
    stats_figures()
    xsection_figures()
