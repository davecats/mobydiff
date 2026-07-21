#!/usr/bin/env python3
"""Six-panel span-averaged field visualisation of a C10 snapshot:
u (chordwise), w (lift), p, and the SST variables k, omega, nut.

  plot_c10_turb_fields.py <field.h5> [--window x0 x1 z0 z1] [--out f.png]

Fields are painted on the finest lattice (coarse cells repeated, the
cv_forces convention); solid cells (penalization-damped |u| ~ 1e-26)
and removed buried blocks show as blanks. k, omega and nut/nu are shown
in log10.
"""
import argparse

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NAMES = ["un", "wn", "pn", "k", "omega", "nut"]


def paint(path, x0, x1, z0, z1):
    """Span-averaged fields on the finest lattice inside the window."""
    with h5py.File(path, "r") as f:
        blocks = f["blocks"][...]
        nb = int(f.attrs["block_nb_x"])
        nx = int(f.attrs["nx"]); nz = int(f.attrs["nz"])
        lx = float(f.attrs["lx"]); lz = float(f.attrs["lz"])
        mask = np.asarray(f.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        lmax = int(blocks[:, 3].max())
        hx = lx/(nx*2**(lmax*int(mask[0])))
        hz = lz/(nz*2**(lmax*int(mask[2])))
        i0, i1 = int(x0/hx), int(np.ceil(x1/hx))
        j0, j1 = int(z0/hz), int(np.ceil(z1/hz))
        out = {n: np.full((j1-j0, i1-i0), np.nan) for n in NAMES}
        data = {n: f[n] for n in NAMES if n in f}
        for bid, (ox, oy, oz, lev) in enumerate(blocks):
            fx = 2**((lmax-int(lev))*int(mask[0]))
            fz = 2**((lmax-int(lev))*int(mask[2]))
            bx0, bz0 = ox*fx, oz*fz
            if bx0 >= i1 or bx0 + nb*fx <= i0 or bz0 >= j1 or bz0 + nb*fz <= j0:
                continue
            si0 = max(i0, bx0); si1 = min(i1, bx0 + nb*fx)
            sj0 = max(j0, bz0); sj1 = min(j1, bz0 + nb*fz)
            for n in data:
                m2 = data[n][bid][...].mean(axis=1)          # (z, x)
                rep = m2.repeat(fz, axis=0).repeat(fx, axis=1)
                out[n][sj0-j0:sj1-j0, si0-i0:si1-i0] = \
                    rep[sj0-bz0:sj1-bz0, si0-bx0:si1-bx0]
        xc = (np.arange(i0, i1) + 0.5)*hx
        zc = (np.arange(j0, j1) + 0.5)*hz
    return xc, zc, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--window", type=float, nargs=4,
                    default=[48.75, 53.25, 46.4, 49.6],
                    help="x0 x1 z0 z1 (chord-lift plane)")
    ap.add_argument("--re", type=float, default=4.0e5)
    ap.add_argument("--bl", action="store_true",
                    help="boundary-layer colour ranges (tight around the "
                         "ambient levels) for near-wall zooms")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    x0, x1, z0, z1 = a.window

    xc, zc, F = paint(a.h5, x0, x1, z0, z1)
    solid = (np.abs(F["un"]) + np.abs(F["wn"])) < 1e-20
    for n in NAMES:
        F[n] = np.where(solid, np.nan, F[n])
    nu = 1.0/a.re

    if a.bl:
        # boundary-layer contrast: ambient k = 3.75e-3 (log -2.43),
        # ambient omega = 150 (log 2.18), ambient nut/nu = 10 (log 1)
        panels = [
            ("un",    "u (chordwise)",           dict(cmap="RdBu_r", vmin=-0.3, vmax=1.7)),
            ("wn",    "w (lift)",                dict(cmap="RdBu_r", vmin=-0.5, vmax=0.5)),
            ("pn",    "p",                       dict(cmap="RdBu_r", vmin=-0.5, vmax=0.5)),
            ("k",     r"log$_{10}$ k",           dict(cmap="magma", vmin=-3.0, vmax=-1.0)),
            ("omega", r"log$_{10}$ $\omega$",    dict(cmap="magma", vmin=1.5, vmax=4.5)),
            ("nut",   r"log$_{10}$ $\nu_t/\nu$", dict(cmap="magma", vmin=0.5, vmax=2.5)),
        ]
    else:
        panels = [
            ("un",    "u (chordwise)",           dict(cmap="RdBu_r", vmin=-0.3, vmax=1.7)),
            ("wn",    "w (lift)",                dict(cmap="RdBu_r", vmin=-0.5, vmax=0.5)),
            ("pn",    "p",                       dict(cmap="RdBu_r", vmin=-0.5, vmax=0.5)),
            ("k",     r"log$_{10}$ k",           dict(cmap="magma", vmin=-6.0, vmax=-1.0)),
            ("omega", r"log$_{10}$ $\omega$",    dict(cmap="magma", vmin=0.0, vmax=6.0)),
            ("nut",   r"log$_{10}$ $\nu_t/\nu$", dict(cmap="magma", vmin=-1.0, vmax=3.0)),
        ]
    fig, axes = plt.subplots(3, 2, figsize=(13, 12.5))
    for ax, (name, title, kw) in zip(axes.ravel(), panels):
        v = F[name]
        if name == "k":
            v = np.log10(np.maximum(v, 1e-30))
        elif name == "omega":
            v = np.log10(np.maximum(v, 1e-30))
        elif name == "nut":
            v = np.log10(np.maximum(v/nu, 1e-30))
        im = ax.pcolormesh(xc, zc, v, shading="nearest", **kw)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x/c"); ax.set_ylabel("z/c")
        fig.colorbar(im, ax=ax, shrink=0.85)
    tag = a.h5.replace(".h5", "")
    fig.suptitle(f"{tag} (span-averaged; body/buried blank)", y=0.995)
    out = a.out or f"fields_{tag}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
