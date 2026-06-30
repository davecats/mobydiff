#!/usr/bin/env python3
"""Time-averaged fluctuation-rms + mean-nut cross-sections (z-y, x-normal),
reassembled onto the finest lattice, for a refined case vs the uniform LES
control. A spurious interface band shows as a localised ridge in u'/v'/w' rms or
a nut spike AT the interface beyond the physical resolution/filter-width step.
Usage: python3 fig_interface_rms.py slab   (or patch)"""
import os, sys, glob
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
from compare_fields import block_geometry           # noqa: E402
VARS = ["un", "vn", "wn", "nut"]


def fine_y(yb, lmax):
    yl = yb.copy()
    for _ in range(lmax):
        mid = 0.5 * (yl[:-1] + yl[1:]); new = np.empty(2 * len(yl) - 1)
        new[0::2] = yl; new[1::2] = mid; yl = new
    return yl


def accumulate(case):
    snaps = sorted(glob.glob(os.path.join(HERE, "runs", case, "stats", f"{case}_*.h5")))
    if not snaps:
        return None
    with h5py.File(snaps[0]) as f:
        bl, nb, lmax, shape = block_geometry(f); yb = f["y"][...]
        has = {v: (v in f) for v in VARS}
    nz, ny, nx = shape
    s1 = {v: np.zeros(shape) for v in VARS if has[v]}
    s2 = {v: np.zeros(shape) for v in VARS if has[v]}
    n = 0
    for s in snaps:
        with h5py.File(s) as f:
            for v in s1:
                D = f[v][...]
                G = np.full(shape, np.nan)
                for bid, (ox, oy, oz, l) in enumerate(bl):
                    sc = 2 ** (lmax - int(l)); x0, y0, z0 = ox * sc, oy * sc, oz * sc
                    ex = np.repeat(np.repeat(np.repeat(D[bid], sc, 0), sc, 1), sc, 2)
                    G[z0:z0 + nb[2] * sc, y0:y0 + nb[1] * sc, x0:x0 + nb[0] * sc] = ex
                s1[v] += G; s2[v] += G * G
        n += 1
    mean = {v: s1[v] / n for v in s1}
    rms = {v: np.sqrt(np.maximum(s2[v] / n - mean[v] ** 2, 0.0)) for v in s1}
    return mean, rms, fine_y(yb, lmax), bl, nb, lmax, n


def interfaces(case, bl, nb, lmax, yfn):
    """y-coordinates of level changes (refined-case interface planes)."""
    if case != "slab":
        return []
    lev = bl[:, 3]
    fine_rows = sorted(set((bl[lev > 0][:, 1] // nb[1]).tolist()))  # fine base rows
    # interface = top of a fine row that borders a coarse row (and mirror)
    ys = []
    gnb = 48  # base ny
    for r in range(6):
        if r in fine_rows and (r + 1) not in fine_rows:
            ys.append(yfn[(r + 1) * nb[1] * (2 ** lmax) // (2 ** lmax)] if False else None)
    # simpler: interfaces at base nodes bounding the fine band
    base_y = yfn[::2 ** lmax]
    out = []
    for r in range(1, 6):
        if (r in fine_rows) != ((r - 1) in fine_rows):
            out.append(base_y[r * nb[1]])
    return out


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else "slab"
    R = accumulate(case); U = accumulate("uniform")
    if R is None or U is None:
        sys.exit("need both <case> and uniform snapshots")
    rm, rr, yfn, bl, nb, lmax, nR = R
    um, ur, yfu, blu, nbu, lmu, nU = U
    iy = interfaces(case, bl, nb, lmax, yfn)
    # each column on its own grid (uniform 64^3, refined finest 128x96x128)
    zfnR = np.arange(rr["un"].shape[0] + 1); xR = rr["un"].shape[2] // 2
    zfnU = np.arange(ur["un"].shape[0] + 1); xU = ur["un"].shape[2] // 2

    panels = [("un", rr, ur, "u' rms"), ("vn", rr, ur, "v' rms"),
              ("wn", rr, ur, "w' rms"), ("nut", rm, um, "mean nut")]
    fig, ax = plt.subplots(len(panels), 2, figsize=(10, 4 * len(panels)),
                           sharex=True)
    for r, (v, Rd, Ud, t) in enumerate(panels):
        if v not in Rd:
            continue
        pr = Rd[v][:, :, xR]; pu = Ud[v][:, :, xU]
        vmax = np.nanpercentile(np.concatenate([pr.ravel(), pu.ravel()]), 99)
        vmin = 0.0
        for c, (dat, yn, zn, tag) in enumerate([(pu, yfu, zfnU, "uniform"),
                                                (pr, yfn, zfnR, case)]):
            m = ax[r, c].pcolormesh(yn, zn, dat, cmap="viridis", vmin=vmin, vmax=vmax,
                                    shading="flat")
            for yi in iy:
                ax[r, c].axvline(yi, color="w", ls="--", lw=0.8)
            if r == 0:
                ax[r, c].set_title(tag)
            if c == 0:
                ax[r, c].set_ylabel(f"{t}\n z-index")
            if r == len(panels) - 1:
                ax[r, c].set_xlabel("y")
            plt.colorbar(m, ax=ax[r, c], fraction=0.025, pad=0.02)
    fig.suptitle(f"{case} vs uniform: time-avg fluctuation rms + mean nut "
                 f"(x-normal; dashed = 2:1 interface; {nR} snaps)")
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    out = os.path.join(HERE, f"fig_{case}_rms_slice.png")
    fig.savefig(out, dpi=120); print("wrote", out)


if __name__ == "__main__":
    main()
