#!/usr/bin/env python3
"""Interface-quality diagnostics for a mobydiff Beltrami field (block-table).

Measures the physical properties a good 2:1 interface treatment must have,
INSTEAD of the (chaotic, uninformative) long-time Beltrami stability:

  (i)  spurious high-k near interfaces -> error and a roughness proxy (RMS of the
       per-cell discrete Laplacian of the ERROR) split interface-band vs interior.
  (iii) momentum conservation -> total volume-weighted Sum u,v,w (exactly 0 for
       Beltrami); the drift is the spurious net momentum.
  (iv) accuracy / order -> volume-weighted L2 / Linf error vs the exact solution,
       split band vs interior; pass two files (coarse, fine) for the order.

  (ii) (mass / correct divU at the interface) is measured separately with the
       solver's own dump_divergence (MOBY_RKDIV), not here -- python lacks the
       ghost layers and ifGrad to reproduce the interface divergence operator.

The exact 3D Beltrami / ABC solution (k=2pi/Lx, F=exp(-nu k^2 t)):
  u=sin(kz)+cos(ky)  v=sin(kx)+cos(kz)  w=sin(ky)+cos(kx)  p=-1/2|u|^2.

The interface BAND is the shell of cells within `--band` (default 2 coarse cells)
of the level-1 patch surface (derived from the level-1 blocks), i.e. inside the
box grown by band but outside the box shrunk by band -- the cells the 2:1
transfer actually touches (faces, edges, corners).

Usage:
  python3 tools/interface_diagnostics.py FIELD.h5 [FIELD_FINE.h5] [--band 0.4]
"""
from __future__ import annotations

import argparse
import h5py
import numpy as np


def exact_block(ox, oy, oz, lev, lx, nx, nb, F):
    h = lx / (nx * 2 ** lev)
    k0 = 2.0 * np.pi / lx
    xc = (ox + np.arange(nb) + 0.5) * h
    yc = (oy + np.arange(nb) + 0.5) * h
    zc = (oz + np.arange(nb) + 0.5) * h
    sx, cx = np.sin(k0 * xc), np.cos(k0 * xc)
    sy, cy = np.sin(k0 * yc), np.cos(k0 * yc)
    sz, cz = np.sin(k0 * zc), np.cos(k0 * zc)
    full = (nb, nb, nb)
    exu = np.broadcast_to(F * (sz[:, None, None] + cy[None, :, None]), full)
    exv = np.broadcast_to(F * (cz[:, None, None] + sx[None, None, :]), full)
    exw = np.broadcast_to(F * (sy[None, :, None] + cx[None, None, :]), full)
    cen = (xc, yc, zc)
    return {"u": exu, "v": exv, "w": exw}, h, cen


def lap_error(e):
    """RMS of the 3D discrete Laplacian of the error over interior cells of a
    block (high-k / checkerboard proxy). Returns (rms, count) or (nan, 0) if the
    block is too small (nb<3 leaves no interior cell)."""
    if min(e.shape) < 3:
        return np.nan, 0
    d2 = (e[2:, 1:-1, 1:-1] - 2 * e[1:-1, 1:-1, 1:-1] + e[:-2, 1:-1, 1:-1]
          + e[1:-1, 2:, 1:-1] - 2 * e[1:-1, 1:-1, 1:-1] + e[1:-1, :-2, 1:-1]
          + e[1:-1, 1:-1, 2:] - 2 * e[1:-1, 1:-1, 1:-1] + e[1:-1, 1:-1, :-2])
    return float(np.sum(d2 ** 2)), d2.size


def analyze(field, band):
    h5 = h5py.File(field, "r")
    lx = float(h5.attrs["lx"]); re = float(h5.attrs["re"]); t = float(h5.attrs["t_current"])
    nx = int(h5.attrs["nx"]); nb = int(h5.attrs["block_nb_x"])
    nu = 1.0 / re; k0 = 2.0 * np.pi / lx
    F = np.exp(-nu * k0 * k0 * t)
    blocks = h5["blocks"][...]
    D = {v: h5[{"u": "un", "v": "vn", "w": "wn"}[v]][...] for v in "uvw"}
    h5.close()

    # Level-1 patch bounding box (physical) -> the interface surface.
    fine = blocks[blocks[:, 3] == 1]
    if len(fine):
        hf = lx / (nx * 2)
        lo = np.array([fine[:, d].min() * hf for d in range(3)])
        hi = np.array([(fine[:, d].max() + nb) * hf for d in range(3)])
    else:
        lo = hi = None

    # Accumulators per region (band / interior=not-band), per variable.
    s2 = {r: {v: 0.0 for v in "uvw"} for r in ("band", "int")}
    vol = {r: 0.0 for r in ("band", "int")}
    linf = {r: {v: 0.0 for v in "uvw"} for r in ("band", "int")}
    rough = {r: {v: [0.0, 0] for v in "uvw"} for r in ("band", "int")}
    mom = {v: 0.0 for v in "uvw"}

    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        ex, h, (xc, yc, zc) = exact_block(ox, oy, oz, lev, lx, nx, nb, F)
        w = h ** 3
        # per-cell band membership: within `band` of the patch surface
        if lo is None:
            inband = np.zeros((nb, nb, nb), bool)
        else:
            X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")  # but data is (k,j,i)=(z,y,x)
            # build a (k,j,i) mask: coords broadcast as z(k),y(j),x(i)
            inO = ((zc[:, None, None] >= lo[2] - band) & (zc[:, None, None] <= hi[2] + band)
                   & (yc[None, :, None] >= lo[1] - band) & (yc[None, :, None] <= hi[1] + band)
                   & (xc[None, None, :] >= lo[0] - band) & (xc[None, None, :] <= hi[0] + band))
            inI = ((zc[:, None, None] >= lo[2] + band) & (zc[:, None, None] <= hi[2] - band)
                   & (yc[None, :, None] >= lo[1] + band) & (yc[None, :, None] <= hi[1] - band)
                   & (xc[None, None, :] >= lo[0] + band) & (xc[None, None, :] <= hi[0] - band))
            inband = np.broadcast_to(inO & ~inI, (nb, nb, nb))
        for v in "uvw":
            e = D[v][bid] - ex[v]
            mom[v] += float(np.sum(D[v][bid])) * w
            for r, mask in (("band", inband), ("int", ~inband)):
                if mask.any():
                    s2[r][v] += float(np.sum((e[mask]) ** 2)) * w
                    linf[r][v] = max(linf[r][v], float(np.abs(e[mask]).max()))
            vol_b = float(inband.sum()) * w; vol["band"] += vol_b
            vol["int"] += (e.size - inband.sum()) * w
            # roughness on whole-block interior, attributed by block band-fraction
            ss, n = lap_error(e)
            if n:
                frac = inband.mean()
                rough["band"][v][0] += ss * frac; rough["band"][v][1] += n * frac
                rough["int"][v][0] += ss * (1 - frac); rough["int"][v][1] += n * (1 - frac)

    # vol double-counted over 3 vars -> divide
    for r in vol:
        vol[r] /= 3.0
    res = {"t": t, "re": re, "nx": nx, "F": F, "mom": mom}
    for r in ("band", "int"):
        l2v = np.sqrt(sum(s2[r].values()) / vol[r]) if vol[r] else np.nan
        lf = max(linf[r].values())
        rg = {v: (np.sqrt(rough[r][v][0] / rough[r][v][1]) if rough[r][v][1] else np.nan)
              for v in "uvw"}
        res[r] = {"L2": l2v, "Linf": lf, "rough": max(rg.values())}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fields", nargs="+")
    ap.add_argument("--band", type=float, default=0.4, help="band half-width (physical)")
    a = ap.parse_args()
    out = []
    for f in a.fields:
        r = analyze(f, a.band)
        out.append(r)
        print(f"\n{f}")
        print(f"  Re={r['re']:.0f} nx={r['nx']} t={r['t']:.3f}")
        print(f"  momentum drift |Pu,Pv,Pw| = "
              f"{abs(r['mom']['u']):.2e} {abs(r['mom']['v']):.2e} {abs(r['mom']['w']):.2e}  (exact 0)")
        for reg in ("band", "int"):
            d = r[reg]
            tag = "interface-band" if reg == "band" else "interior      "
            print(f"  {tag}:  L2={d['L2']:.4e}  Linf={d['Linf']:.4e}  roughRMS(lap err)={d['rough']:.4e}")
    if len(out) == 2:
        c, fdat = out
        if c["nx"] < fdat["nx"]:
            pass
        else:
            c, fdat = fdat, c
        print("\norder (log2 coarse/fine):")
        for reg in ("band", "int"):
            o = np.log2(c[reg]["L2"] / fdat[reg]["L2"])
            print(f"  {reg:14s} L2 order = {o:.2f}")


if __name__ == "__main__":
    main()
