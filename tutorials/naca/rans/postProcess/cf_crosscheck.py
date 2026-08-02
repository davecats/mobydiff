#!/usr/bin/env python3
"""Two INDEPENDENT Cf estimators for the immersed NACA 0012, as a check on
the penalization-band estimator that surface_cp_cf.py uses by default.

  cf_crosscheck.py <field.h5> --coef <case.h5> [--npz cpcf_*.npz]
                   [--L 5.0] [--dmax-fit 4.0] [--stations 300]
                   [--plot cf_crosscheck.png]

The default estimator (surface_cp_cf.py) reads the wall gradient off the
PENALIZATION BAND, using the fact that the IBM coefficient encodes the
true sub-cell wall distance. The two here are built so that they fail
differently -- one never touches the band, the other never differences
near the wall at all:

  FIT — the band cells are DISCARDED and the wall gradient is fitted from
  the clean fluid outside it: u_t = g d + c d^2 through the origin
  (exact no-slip on the analytic surface), weighted 1/d, over d <=
  dmax_fit cells. The quadratic term is the near-wall expansion
  mu u_t'' = dp/ds; around the nose dp/ds is large and a straight line
  reads the wall slope low. Its weakness is the fit depth: g moves 3-20 %
  over 3..6 h, since 4 h is already y+ 6-8 (leaving the viscous sublayer
  at turbulent stations) while the laminar BL at the nose is only 1-3
  cells.

  VISCOUS INTEGRAL — the effective viscous term the momentum equation
  sees in the fluid is V = nu lap_h(u) - coef u (the raw Laplacian, which
  reads the penalized ~0 in the solid, plus the coefficient correction
  that moves the wall to its true sub-cell position). It approximates
  nu grad^2 u, so integrating along the wall normal from the wall to a
  height L gives
      int_0^L V dn = nu du_t/dn|_L - tau_w
  and hence tau_w = nu du_t/dn|_L - int_0^L V dn: the wall gradient from
  a gradient taken at a COMFORTABLE height L in clean fluid plus a volume
  integral, with no near-wall differencing. tau_w must not depend on L,
  and does not (median ratio 0.899-0.903 over L = 3..8 h), which is the
  check on the implementation. It is the noisiest of the three -- a
  difference of two larger numbers -- and overshoots by ~30 % around
  x/c 0.02-0.03 where the surface is steepest against the grid and the
  station-lateral viscous flux it neglects is largest.

Rejected alternative: the wall shear as the force-density integral
int coef*u_t dn, with no wall location at all. It does NOT work on this
case -- the ini runs without keep_buried, so the buried blocks carrying
much of that force are removed and the per-station integral scatters
over 0.04-4.6x OpenFOAM. It would need a keep-buried case file.
"""
import argparse
import os

import numpy as np
from scipy.spatial import cKDTree

import surface_cp_cf as S


def load_openfoam(of_dir):
    """OpenFOAM wall-shear sampling -> Cf(x/c) per side (compare_openfoam.py
    convention: rho = 1, U = 1, attached TE-ward flow has tau_x < 0)."""
    samp = os.path.join(of_dir, "postProcessing", "surface_sampling")
    t = sorted(os.listdir(samp), key=float)[-1]
    of = {}
    for side in ("suction", "pressure"):
        tau = np.loadtxt(os.path.join(samp, t, f"wallShearStress_{side}_side.raw"),
                         comments="#")
        o = np.argsort(tau[:, 0]); tau = tau[o]
        mag = np.linalg.norm(tau[:, 3:6], axis=1)
        of[side] = (tau[:, 0], 2.0*mag*np.sign(-tau[:, 3]))
    return of


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--coef", required=True)
    ap.add_argument("--npz", default=None, help="surface_cp_cf.py output to compare against")
    ap.add_argument("--L", type=float, default=5.0, help="integration height in cells")
    ap.add_argument("--dmax-fit", type=float, default=4.0, help="fit depth in cells")
    ap.add_argument("--stations", type=int, default=300)
    ap.add_argument("--re", type=float, default=4.0e5)
    ap.add_argument("--nose", type=float, nargs=2, default=[50.0, 48.0])
    ap.add_argument("--of", default="assets/openfoam",
                    help="OpenFOAM case dir to overlay in --plot")
    ap.add_argument("--plot", default=None)
    a = ap.parse_args()
    S.XLE, S.ZLE = a.nose
    nu = 1.0/a.re

    # Single-level surfaces only: the viscous integral differences the field
    # on ONE lattice, so a surface spanning two levels (refine_body_box)
    # would silently be evaluated on the finest patch alone.
    levs = S.surface_levels(a.h5)
    if len(levs) > 1:
        raise SystemExit(
            f"{a.h5}: the surface spans levels {levs}. cf_crosscheck.py handles "
            "single-level surfaces only — run it on the level-11 baseline "
            "(surface_cp_cf.py --levels handles the multi-level case).")

    D = S.build_lattice(a.h5, a.coef)
    h = D["hx"]; U, W, CU, CW = D["U"], D["W"], D["CU"], D["CW"]
    clean = S.clean_mask(D)

    def lap(F):
        L = np.full_like(F, np.nan)
        L[1:-1, 1:-1] = ((F[1:-1, 2:] - 2*F[1:-1, 1:-1] + F[1:-1, :-2])/D["hx"]**2
                         + (F[2:, 1:-1] - 2*F[1:-1, 1:-1] + F[:-2, 1:-1])/D["hz"]**2)
        return L

    effU = nu*lap(U) - CU*U     # coef already carries 1/Re
    effW = nu*lap(W) - CW*W
    # cell-centred velocities for the FIT (staggered faces averaged in)
    uc = np.full_like(U, np.nan); wc = np.full_like(W, np.nan)
    uc[:, :-1] = 0.5*(U[:, :-1] + U[:, 1:])
    wc[:-1, :] = 0.5*(W[:-1, :] + W[1:, :])

    kk, ii = np.mgrid[0:U.shape[0], 0:U.shape[1]]
    xc = (D["ox0"] + ii + 0.5)*D["hx"]; zc = (D["oz0"] + kk + 0.5)*D["hz"]

    # same polyline resolution as the default estimator, so the ratios below
    # compare estimators and not two discretizations of the surface
    pts, nrm, tng, xoc, upper = S.naca0012_polyline()
    arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    tree = cKDTree(pts)
    reach = 30*h
    dnU, ciU = S.surface_geometry(np.stack([xc - 0.5*D["hx"], zc], -1),
                                  np.isfinite(U), tree, pts, nrm, reach)
    dnW, ciW = S.surface_geometry(np.stack([xc, zc - 0.5*D["hz"]], -1),
                                  np.isfinite(W), tree, pts, nrm, reach)
    dnC, ciC = S.surface_geometry(np.stack([xc, zc], -1),
                                  clean & np.isfinite(uc) & np.isfinite(wc),
                                  tree, pts, nrm, reach)
    uu, ww = U.ravel(), W.ravel()
    eu, ew = effU.ravel(), effW.ravel()
    cu, cw = CU.ravel(), CW.ravel()
    ucr, wcr = uc.ravel(), wc.ravel()
    tU = tng[ciU]*np.where(tng[ciU, 0:1] >= 0, 1.0, -1.0)
    tW = tng[ciW]*np.where(tng[ciW, 0:1] >= 0, 1.0, -1.0)
    tC = tng[ciC]*np.where(tng[ciC, 0:1] >= 0, 1.0, -1.0)
    utC = ucr*tC[:, 0] + wcr*tC[:, 1]

    nst = a.stations
    edges = np.linspace(0.0, arc[-1], nst + 1)
    ds = arc[-1]/nst

    def station(ci):
        s = np.full(ci.shape, -1)
        g = ci >= 0
        s[g] = np.clip(np.searchsorted(edges, arc[ci[g]], side="right") - 1, 0, nst - 1)
        return s

    sU, sW, sC = station(ciU), station(ciW), station(ciC)
    mids = np.clip(np.searchsorted(arc, 0.5*(edges[:-1] + edges[1:])), 0, pts.shape[0] - 1)
    L = a.L*h

    fit = np.full(nst, np.nan); visc = np.full(nst, np.nan)
    for s in range(nst):
        # ---- FIT: anchored quadratic over the cells OUTSIDE the band ----
        m = (sC == s) & (dnC > 0) & (dnC <= a.dmax_fit*h)
        if m.sum() >= 4:
            d, v = dnC[m], utC[m]
            wq = 1.0/d
            S2 = np.sum(wq*d**2); S3 = np.sum(wq*d**3); S4 = np.sum(wq*d**4)
            b1 = np.sum(wq*d*v); b2 = np.sum(wq*d**2*v)
            det = S2*S4 - S3*S3
            if det != 0.0:
                fit[s] = 2.0*nu*(b1*S4 - b2*S3)/det

        # ---- VISCOUS INTEGRAL ----
        iu = (sU == s) & (dnU > 0) & (dnU < L) & np.isfinite(eu)
        iw = (sW == s) & (dnW > 0) & (dnW < L) & np.isfinite(ew)
        I = (np.sum(eu[iu]*tU[iu, 0]) + np.sum(ew[iw]*tW[iw, 1]))*h*h/ds
        # local slope at height L, from clean points only
        ou = (sU == s) & (cu == 0) & (np.abs(dnU - L) < 1.5*h)
        ow = (sW == s) & (cw == 0) & (np.abs(dnW - L) < 1.5*h)
        tt = np.r_[tU[ou, 0], tW[ow, 1]]
        aa = np.r_[dnU[ou], dnW[ow]]
        bb = np.r_[uu[ou], ww[ow]]
        if tt.size >= 6:
            M = np.column_stack([tt, tt*(aa - L)])
            try:
                c, *_ = np.linalg.lstsq(M, bb, rcond=None)
                visc[s] = 2.0*(nu*c[1] - I)
            except np.linalg.LinAlgError:
                pass

    out = dict(xoc=xoc[mids], up=upper[mids].astype(float), cf_fit=fit, cf_visc=visc)
    ref = np.load(a.npz) if a.npz else None
    if ref is not None:
        band = np.array([ref["cf"][np.argmin(np.abs(ref["xoc"] - x) + 1e3*(ref["up"] != u))]
                         for x, u in zip(out["xoc"], out["up"])])
        out["cf_band"] = band
        for name in ("cf_fit", "cf_visc"):
            r = out[name]/band
            m = np.isfinite(r) & (np.abs(band) > 1e-4)
            print(f"{name} / band (default) : median {np.median(r[m]):.3f}, "
                  f"5-95% [{np.percentile(r[m], 5):.3f}, {np.percentile(r[m], 95):.3f}] "
                  f"over {m.sum()} stations")
    tag = os.path.basename(a.h5).replace(".h5", "")
    np.savez(f"cf_crosscheck_{tag}.npz", **out)
    print(f"wrote cf_crosscheck_{tag}.npz")

    if a.plot:
        plot(out, a)


def plot(out, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    of = load_openfoam(a.of) if a.of and os.path.isdir(a.of) else None
    # Identity is carried by colour AND dash pattern, so the four curves stay
    # separable in greyscale and for colour-vision deficiency. The three
    # estimator hues are Okabe-Ito steps, validated (all-pairs CVD dE 11.0
    # deutan / 8.6 tritan, normal-vision 18.7, contrast >= 3:1 on white).
    series = [("cf_band", "penalization band (default)", "#0072B2", "-", 2.0),
              ("cf_fit", "wall-gradient fit", "#D55E00", "-.", 1.7),
              ("cf_visc", "viscous integral", "#009E73", ":", 1.9)]
    INK, MUTED = "#1a1a1a", "#5c5c5c"
    fig, axs = plt.subplots(2, 2, figsize=(12.0, 7.0))
    panels = ((0, 1.0, "suction side", (0.0, 1.0), (-0.004, 0.016)),
              (1, 1.0, "suction side — leading edge", (0.0, 0.12), (0.0, 0.045)),
              (2, 0.0, "pressure side", (0.0, 1.0), (-0.004, 0.016)),
              (3, 0.0, "pressure side — leading edge", (0.0, 0.12), (-0.030, 0.015)))
    for n, up, title, xl, yl in panels:
        ax = axs.ravel()[n]
        side = "suction" if up == 1.0 else "pressure"
        if of is not None:
            ax.plot(of[side][0], of[side][1], "--", lw=2.0, color=INK,
                    label="OpenFOAM (body-fitted)", zorder=5)
        m = out["up"] == up
        o = np.argsort(out["xoc"][m])
        for key, lab, col, ls, lw in series:
            if key not in out:
                continue
            ax.plot(out["xoc"][m][o], out[key][m][o], ls, lw=lw, color=col,
                    label=lab, zorder=4)
        ax.set_title(title, fontsize=10, color=INK)
        ax.set_xlim(*xl); ax.set_ylim(*yl)
        ax.axhline(0.0, lw=0.8, color=MUTED, alpha=0.5, zorder=1)
        ax.grid(alpha=0.25, lw=0.6)
        ax.tick_params(colors=MUTED, labelsize=9)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(MUTED)
        if n >= 2:
            ax.set_xlabel("x/c", color=INK)
        if n % 2 == 0:
            ax.set_ylabel(r"$C_f$", color=INK)
    axs[0, 0].legend(fontsize=8.5, framealpha=0.9, labelcolor=INK, loc="upper right")
    # full-chord panels clip the LE peak; park the note clear of both the
    # spike (x/c < 0.1) and the legend (x/c > 0.35)
    for ax in (axs[0, 0], axs[1, 0]):
        ax.annotate("LE off-scale —\nsee panel at right", xycoords="axes fraction",
                    xy=(0.145, 0.86), fontsize=7.5, color=MUTED, va="top")
    fig.suptitle(r"NACA 0012, Re = 4e5, $\alpha$ = 5$^\circ$: three independent "
                 "mobydiff $C_f$ estimators vs OpenFOAM", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(a.plot, dpi=150)
    print(f"wrote {a.plot}")


if __name__ == "__main__":
    main()
