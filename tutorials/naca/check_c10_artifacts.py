#!/usr/bin/env python3
"""C10 aoa0 artifact battery (docs/next_session_naca_re4e5.md checklist):

  1. ambient tu arriving at the body (k on level 0/1 upstream of the nose)
  2. LE fan-strip metric: rms of the non-smooth part of the lift velocity
     w in surface-parallel strips at the R1 distances
     (0.006/0.012/0.023/0.047/0.094c), x/c in [0, 0.5]
  3. interface striping at the DISTANT interfaces: w along a vertical
     line above mid-chord; second-difference rms in a +-0.15c band
     around each level interface vs the quiet bands between them

  check_c10_artifacts.py <field.h5>
"""
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

import surface_cp_cf as sc

NOSE = (50.0, 48.0)


def probe_k_upstream(path):
    """Mean k in 1x1c boxes at 5c and 2c upstream of the nose."""
    out = {}
    with h5py.File(path, "r") as f:
        blocks = f["blocks"][...]
        nb = int(f.attrs["block_nb_x"])
        nx = int(f.attrs["nx"]); nz = int(f.attrs["nz"])
        lx = float(f.attrs["lx"]); lz = float(f.attrs["lz"])
        K = f["k"]
        for tag, dx in (("5c", -5.0), ("2c", -2.0)):
            x0, x1 = NOSE[0] + dx - 0.5, NOSE[0] + dx + 0.5
            z0, z1 = NOSE[1] - 0.5, NOSE[1] + 0.5
            vals = []
            for bid, (ox, oy, oz, lev) in enumerate(blocks):
                sx = 2**int(lev); sz = 2**int(lev)
                hx = lx/(nx*sx); hz = lz/(nz*sz)
                xc = (ox + 0.5 + np.arange(nb))*hx
                zc = (oz + 0.5 + np.arange(nb))*hz
                mx = (xc > x0) & (xc < x1)
                mz = (zc > z0) & (zc < z1)
                if mx.any() and mz.any():
                    km = K[bid][...].mean(axis=1)   # (z, x)
                    vals.append(km[np.ix_(mz, mx)].ravel())
            out[tag] = float(np.concatenate(vals).mean())
    return out


def fan_strips(path):
    """R1 fan metric: per-strip rms of w minus its 5-station running
    mean, over x/c in [0, 0.5] on the upper side."""
    x, z, u, w, p, p_inf, h = sc.load_fine_plane(path)
    fluid = (np.abs(u) + np.abs(w)) > 1e-20
    pts, nrm, tng, xoc, upper = sc.naca0012_polyline()
    tree = cKDTree(pts)
    d_all, idx = tree.query(np.column_stack([x, z]),
                            distance_upper_bound=0.15)
    near = np.isfinite(d_all) & (idx < pts.shape[0]) & fluid
    ci = idx[near]
    dx = np.column_stack([x[near], z[near]]) - pts[ci]
    dn = np.einsum("ij,ij->i", dx, nrm[ci])
    wv = w[near]
    sel_side = upper[ci] & (xoc[ci] > 0.0) & (xoc[ci] < 0.5) & (dn > 0)
    strips = [0.006, 0.012, 0.023, 0.047, 0.094]
    print("\nfan strips (upper side, x/c in [0, 0.5]):")
    print("  d/c      rms(w_hf)   [R1a unfiltered hit ~0.05-0.1 here]")
    for d in strips:
        m = sel_side & (dn > d/1.4) & (dn < d*1.4)
        if m.sum() < 50:
            print(f"  {d:.3f}    (only {m.sum()} cells)")
            continue
        s = xoc[ci[m]]
        o = np.argsort(s)
        vv = wv[m][o]
        # 100 arc bins; high-pass = value minus bin mean, plus
        # bin-to-bin roughness of the means themselves
        bins = np.linspace(0.0, 0.5, 101)
        bi = np.clip(np.digitize(s[o], bins) - 1, 0, 99)
        mean_b = np.zeros(100); cnt = np.zeros(100)
        np.add.at(mean_b, bi, vv); np.add.at(cnt, bi, 1)
        mean_b = np.where(cnt > 0, mean_b/np.maximum(cnt, 1), np.nan)
        hf = vv - mean_b[bi]
        # smooth trend of the bin means (running 5) -> stripe part
        mb = mean_b.copy()
        ok = np.isfinite(mb)
        sm = np.convolve(np.where(ok, mb, 0.0), np.ones(5)/5, "same") / \
            np.maximum(np.convolve(ok.astype(float), np.ones(5)/5, "same"), 1e-12)
        stripe = mb - sm
        r = np.sqrt(np.nanmean(hf**2) + np.nanmean(stripe[ok]**2))
        print(f"  {d:.3f}    {r:.4e}")


def interface_striping(path):
    """w along the vertical line x = 50.75 (mid-chord), z above the
    airfoil; localized second-difference rms around each interface."""
    from cv_forces import load_plane
    xq = 50.75
    xc, zc, u2, w2, p2, nut2, hx, hz = load_plane(
        path, xq - 0.1, xq + 0.1, 48.3, 56.5)
    icol = int(np.argmin(np.abs(xc - xq)))
    line = w2[:, icol]
    good = np.isfinite(line)
    zc, line = zc[good], line[good]
    d2 = np.abs(np.diff(line, 2))
    zmid = zc[1:-1]
    # C10 upper interface z locations (from the refine boxes)
    ifz = [48.30, 48.55, 49.05, 50.05, 51.05, 53.05, 55.80]
    print("\ninterface striping (|d2 w| rms, vertical line x = 50.75):")
    quiet = np.ones(zmid.size, dtype=bool)
    for zf in ifz:
        band = np.abs(zmid - zf) < 0.15
        quiet &= ~band
        if band.sum() > 3:
            print(f"  z = {zf:6.2f}: {np.sqrt(np.mean(d2[band]**2)):.3e}")
    print(f"  quiet bands: {np.sqrt(np.mean(d2[quiet]**2)):.3e}"
          f"   (ratio > ~3 at an interface = striping)")


if __name__ == "__main__":
    path = sys.argv[1]
    k = probe_k_upstream(path)
    kinf = 3.75e-3
    print(f"ambient k upstream: 5c = {k['5c']:.3e}, 2c = {k['2c']:.3e} "
          f"(k_inf = {kinf:.3e}; tu arriving = "
          f"{100*np.sqrt(k['2c']/kinf*0.0025):.2f}% of U, target 5%)")
    fan_strips(path)
    interface_striping(path)
