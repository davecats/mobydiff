#!/usr/bin/env python3
"""Window slice of a multi-level block-table field at z = mid, painted onto
the finest lattice (coarse cells replicated) — full reassembly at this
grid would be ~69 GB/field, a window around the airfoil is MBs.

  slice_field.py <field.h5> [--window x0 x1 y0 y1] [--out slice.npz]

Prints boundary-layer / separation / turbulence diagnostics for the
NACA 0012 runs (LE at (4.5, 6.0), chord 1) and stores the window arrays.
"""
import argparse
import sys

import h5py
import numpy as np


def load_window(path, names, x0, x1, y0, y1, span="z"):
    """Window slice at mid-SPAN. span = 'z' (default: window in x-y) or
    'y' (the refine_dims = xz orientation: window in x-z, 'yc' returns
    the LIFT z coordinate). Per-direction replication factors honour the
    refine_dims attribute (absent = xyz octree)."""
    with h5py.File(path, "r") as f:
        blocks = f["blocks"][...]
        nb = int(f.attrs["block_nb_x"])
        lmax = int(blocks[:, 3].max())
        nx = int(f.attrs["nx"]); ny = int(f.attrs["ny"]); nz = int(f.attrs["nz"])
        lx = float(f.attrs.get("lx", 12.0))
        lift_l = float(f.attrs.get("ly", 12.0)) if span == "z" else float(f.attrs.get("lz", 12.0))
        lift_n = ny if span == "z" else nz
        span_n = nz if span == "z" else ny
        mask = np.asarray(f.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        m_lift = int(mask[1] if span == "z" else mask[2])
        m_span = int(mask[2] if span == "z" else mask[1])
        hx = lx / (nx * 2**(lmax*int(mask[0])))
        hy = lift_l / (lift_n * 2**(lmax*m_lift))
        i0, i1 = int(x0 / hx), int(np.ceil(x1 / hx))
        j0, j1 = int(y0 / hy), int(np.ceil(y1 / hy))
        out = {n: np.full((j1 - j0, i1 - i0), np.nan) for n in names}
        # per-pixel refinement level: coarse cells are painted as replicated
        # fine pixels, which makes smooth coarse gradients look blocky --
        # read images with this map in hand (LE "checkerboard" analysis).
        out["level"] = np.full((j1 - j0, i1 - i0), -1.0)
        ks_f = (span_n * 2**(lmax*m_span)) // 2   # mid-span, finest lattice
        data = {n: f[n] for n in names}
        for bid, (ox, oy, oz, lev) in enumerate(blocks):
            fx = 2 ** ((lmax - int(lev))*int(mask[0]))
            fy = 2 ** ((lmax - int(lev))*int(mask[1]))
            fz = 2 ** ((lmax - int(lev))*int(mask[2]))
            fl = fy if span == "z" else fz
            fs = fz if span == "z" else fy
            bx0 = ox * fx
            bl0 = (oy * fy) if span == "z" else (oz * fz)
            bs0 = (oz * fz) if span == "z" else (oy * fy)
            if bx0 >= i1 or bx0 + nb * fx <= i0:
                continue
            if bl0 >= j1 or bl0 + nb * fl <= j0:
                continue
            if not (bs0 <= ks_f < bs0 + nb * fs):
                continue
            kk = (ks_f - bs0) // fs
            si0 = max(i0, bx0); si1 = min(i1, bx0 + nb * fx)
            sj0 = max(j0, bl0); sj1 = min(j1, bl0 + nb * fl)
            out["level"][sj0 - j0:sj1 - j0, si0 - i0:si1 - i0] = float(lev)
            for n in names:
                r = data[n][bid]                     # (nbz, nby, nbx)
                row = r[kk] if span == "z" else r[:, kk, :]  # (lift, chord)
                rep = row.repeat(fl, axis=0).repeat(fx, axis=1)
                out[n][sj0 - j0:sj1 - j0, si0 - i0:si1 - i0] = \
                    rep[sj0 - bl0:sj1 - bl0, si0 - bx0:si1 - bx0]
        xc = (np.arange(i0, i1) + 0.5) * hx
        yc = (np.arange(j0, j1) + 0.5) * hy
    return out, xc, yc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--window", type=float, nargs=4, default=[4.0, 6.5, 5.4, 6.6])
    ap.add_argument("--out", default="slice.npz")
    ap.add_argument("--re", type=float, default=1.0e5)
    ap.add_argument("--span", choices=("z", "y"), default="z",
                    help="span axis (y = the refine_dims xz orientation; the "
                         "window's second pair and 'yc' are then the LIFT z)")
    a = ap.parse_args()

    names = ["un", "vn", "pn", "k", "omega", "nut"]
    if a.span == "y":
        # the lift-normal fluctuating component is w in this orientation;
        # keep the npz key 'vn' meaning "lift-direction velocity"
        names = ["un", "wn", "pn", "k", "omega", "nut"]
    with h5py.File(a.h5) as f:
        names = [n for n in names if n in f]
    x0, x1, y0, y1 = a.window
    out, xc, yc = load_window(a.h5, names, x0, x1, y0, y1, span=a.span)
    if a.span == "y" and "wn" in out:
        out["vn"] = out.pop("wn")
    np.savez_compressed(a.out, xc=xc, yc=yc, **out)
    print(f"{a.out}: window [{x0},{x1}]x[{y0},{y1}], {out[names[0]].shape}, fields {names}")

    u = out["un"]
    solid = np.isnan(u) | (np.abs(u) < 1e-30)
    # diagnostics on the suction side (above the chord line y = 6.0)
    if "k" in out:
        kf = out["k"]
        upstream = (xc > 4.0) & (xc < 4.4)
        band = (yc > 5.9) & (yc < 6.1)
        print(f"ambient k just upstream of LE: {np.nanmean(kf[np.ix_(band, upstream)]):.3e} "
              f"(inlet k_inf would be 1.5(tu/100)^2)")
        print(f"max k in window: {np.nanmax(kf):.3e}")
    if "nut" in out:
        print(f"max nut/nu in window: {np.nanmax(out['nut'])*a.re:.1f}")
    # reversed-flow (u < -0.02) fraction along the suction surface
    rev = (u < -0.02)
    chord = (xc > 4.5) & (xc < 5.5)
    above = (yc > 6.0) & (yc < 6.15)
    sub = rev[np.ix_(above, chord)]
    print(f"reversed-flow cells above suction side: {np.nansum(sub)} "
          f"({100.0*np.nansum(sub)/sub.size:.1f}% of the band)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
