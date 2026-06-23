#!/usr/bin/env python3
"""Slice visualisation of a mobydiff Beltrami field vs the exact solution.

For a block-table field dump, draws a constant-z slice through the (refined)
domain: each block is rendered at its NATIVE resolution (pcolormesh per block),
so the 2:1 coarse/fine interface and the patch edges/corners are visible. For
every variable (u, v, w, p) three panels are shown: solver, exact, error
(solver-exact). The exact 3D Beltrami / ABC solution (see check_beltrami.py) is
  u = sin(kz)+cos(ky)  v = sin(kx)+cos(kz)  w = sin(ky)+cos(kx)  p = -1/2|u|^2
decaying by F = exp(-nu k^2 t); pressure is compared after removing the mean
(defined up to a constant). The error panels share a symmetric colour scale so
variants are comparable side by side.

Usage:
  python3 tools/plot_beltrami_fields.py FIELD.h5 [--z 3.1416] [--out fig.png]
"""
from __future__ import annotations

import argparse
import os

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def block_exact(ox, oy, oz, lev, lx, nx, nb, F):
    """Exact (u,v,w,p) at the cell centres of one block, shape (k,j,i)=(z,y,x)."""
    h = lx / (nx * 2 ** lev)
    k0 = 2.0 * np.pi / lx
    xc = (ox + np.arange(nb) + 0.5) * h
    yc = (oy + np.arange(nb) + 0.5) * h
    zc = (oz + np.arange(nb) + 0.5) * h
    sx, cx = np.sin(k0 * xc), np.cos(k0 * xc)
    sy, cy = np.sin(k0 * yc), np.cos(k0 * yc)
    sz, cz = np.sin(k0 * zc), np.cos(k0 * zc)
    # Each component is independent of its OWN direction, so the raw sums are
    # size-1 in that axis; broadcast to the full (k,j,i) block so a slice is 2D.
    full = (nb, nb, nb)
    exu = np.broadcast_to(F * (sz[:, None, None] + cy[None, :, None]), full).copy()  # sin(kz)+cos(ky)
    exv = np.broadcast_to(F * (cz[:, None, None] + sx[None, None, :]), full).copy()  # sin(kx)+cos(kz)
    exw = np.broadcast_to(F * (sy[None, :, None] + cx[None, None, :]), full).copy()  # sin(ky)+cos(kx)
    exp = -0.5 * (exu ** 2 + exv ** 2 + exw ** 2)       # p=-1/2|u|^2
    return {"u": exu, "v": exv, "w": exw, "p": exp}, h


def read_divslice(path, z0, lx, nx, nb):
    """Read the 'pn' slot (= divergence) of a companion div dump and return the
    per-block z-slice quads (xe, ye, C). Returns [] if the file is absent."""
    import os
    if not os.path.exists(path):
        return []
    with h5py.File(path, "r") as h:
        blocks = h["blocks"][...]
        P = h["pn"][...]
    quads = []
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        h = lx / (nx * 2 ** lev)
        kc = int(np.floor(z0 / h)) - oz
        if kc < 0 or kc >= nb:
            continue
        xe = (ox + np.arange(nb + 1)) * h
        ye = (oy + np.arange(nb + 1)) * h
        quads.append((xe, ye, P[bid][kc].astype(float)))
    return quads


def main():
    import os
    ap = argparse.ArgumentParser()
    ap.add_argument("field")
    ap.add_argument("--z", type=float, default=None, help="slice z (default lz/2)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    h5 = h5py.File(a.field, "r")
    lx = float(h5.attrs["lx"]); re = float(h5.attrs["re"]); t = float(h5.attrs["t_current"])
    nx = int(h5.attrs["nx"]); nb = int(h5.attrs["block_nb_x"])
    nu = 1.0 / re; k0 = 2.0 * np.pi / lx
    F = np.exp(-nu * k0 * k0 * t)
    z0 = a.z if a.z is not None else lx / 2.0
    blocks = h5["blocks"][...]
    data = {v: h5[{"u": "un", "v": "vn", "w": "wn", "p": "pn"}[v]][...] for v in "uvwp"}
    h5.close()

    # Mean offset for pressure (defined up to a constant): remove the volume mean
    # of both solver and exact so the error panel shows structure, not the offset.
    pvol = pwt = epvol = 0.0
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        ex, h = block_exact(ox, oy, oz, lev, lx, nx, nb, F)
        wgt = h ** 3
        pvol += np.sum(data["p"][bid]) * wgt; pwt += data["p"][bid].size * wgt
        epvol += np.sum(ex["p"]) * wgt
    p_off = pvol / pwt; ep_off = epvol / pwt

    # Collect per-block slice quads: (Xedges, Yedges, solver, exact, err) per var.
    quads = {v: [] for v in "uvwp"}
    rng = {v: {"f": 0.0, "e": 0.0} for v in "uvwp"}
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        h = lx / (nx * 2 ** lev)
        kc = int(np.floor(z0 / h)) - oz
        if kc < 0 or kc >= nb:
            continue
        ex, _ = block_exact(ox, oy, oz, lev, lx, nx, nb, F)
        xe = (ox + np.arange(nb + 1)) * h
        ye = (oy + np.arange(nb + 1)) * h
        for v in "uvwp":
            sol = data[v][bid][kc].astype(float)        # (j,i)=(y,x)
            exa = ex[v][kc].astype(float)
            if v == "p":
                sol = sol - p_off; exa = exa - ep_off
            err = sol - exa
            quads[v].append((xe, ye, sol, exa, err))
            rng[v]["f"] = max(rng[v]["f"], np.abs(sol).max(), np.abs(exa).max())
            rng[v]["e"] = max(rng[v]["e"], np.abs(err).max())

    # Companion divergence dumps (MOBY_DIVDUMP): <prefix>_divpre_<step>.h5 etc.
    base = a.field
    dirn, fn = os.path.split(base)
    stem, _, rest = fn.rpartition("_")          # "tst", "_", "1.h5"
    pre_path = os.path.join(dirn, f"{stem}_divpre_{rest}")
    post_path = os.path.join(dirn, f"{stem}_divpost_{rest}")
    divq = {"pre": read_divslice(pre_path, z0, lx, nx, nb),
            "post": read_divslice(post_path, z0, lx, nx, nb)}
    have_div = bool(divq["pre"]) or bool(divq["post"])

    nrows = 5 if have_div else 4
    fig, axs = plt.subplots(nrows, 3, figsize=(13, 4 * nrows), constrained_layout=True)
    cols = ["solver", "exact", "error (solver-exact)"]
    for r, v in enumerate("uvwp"):
        fmax = rng[v]["f"] or 1.0
        emax = rng[v]["e"] or 1.0
        for c in range(3):
            ax = axs[r, c]
            cmap = "RdBu_r"
            vmin, vmax = (-emax, emax) if c == 2 else (-fmax, fmax)
            for (xe, ye, sol, exa, err) in quads[v]:
                C = [sol, exa, err][c]
                pcm = ax.pcolormesh(xe, ye, C, cmap=cmap, vmin=vmin, vmax=vmax,
                                    shading="flat")
            ax.set_aspect("equal"); ax.set_xlim(0, lx); ax.set_ylim(0, lx)
            if r == 0:
                ax.set_title(cols[c])
            if c == 0:
                ax.set_ylabel(f"{v}    y", fontsize=12)
            fig.colorbar(pcm, ax=ax, shrink=0.8)

    if have_div:
        # bottom row: divergence before / after correction (exact div = 0), shared
        # symmetric scale; col 2 unused.
        dmax = max([np.abs(C).max() for q in divq.values() for (_, _, C) in q] or [1.0])
        labels = {0: "div BEFORE correction (predictor)", 1: "div AFTER correction"}
        for c, key in ((0, "pre"), (1, "post")):
            ax = axs[4, c]
            for (xe, ye, C) in divq[key]:
                pcm = ax.pcolormesh(xe, ye, C, cmap="RdBu_r", vmin=-dmax, vmax=dmax, shading="flat")
            ax.set_aspect("equal"); ax.set_xlim(0, lx); ax.set_ylim(0, lx)
            ax.set_title(labels[c])
            if c == 0:
                ax.set_ylabel("divU   y", fontsize=12)
            fig.colorbar(pcm, ax=ax, shrink=0.8)
        axs[4, 2].axis("off")
        mpre = max([np.abs(C).max() for (_, _, C) in divq["pre"]] or [0.0])
        mpost = max([np.abs(C).max() for (_, _, C) in divq["post"]] or [0.0])
        axs[4, 2].text(0.05, 0.5, f"max|divU|\n pre  = {mpre:.3e}\n post = {mpost:.3e}",
                       fontsize=12, va="center", family="monospace")

    fig.suptitle(f"{os.path.basename(a.field)}   Re={re:.0f}  t={t:.2f}  "
                 f"z-slice={z0:.3f}  (block-native resolution)", fontsize=13)
    out = a.out or (os.path.splitext(a.field)[0] + "_fields.png")
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")
    # also print the per-variable Linf error on this slice for a quick read
    for v in "uvwp":
        print(f"  {v}: slice max|err| = {rng[v]['e']:.4e}")


if __name__ == "__main__":
    main()
