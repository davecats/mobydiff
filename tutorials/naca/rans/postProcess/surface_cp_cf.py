#!/usr/bin/env python3
"""Surface Cp and Cf on the immersed NACA 0012 from a block-table field
snapshot (span-y / refine_dims = xz orientation).

  surface_cp_cf.py <field.h5> --coef <case.h5> [--out cpcf_aoaX.npz]
                   [--stations 300] [--dmax-cp 12.0] [--plot cpcf.png]

Method (per surface station on the ANALYTIC section polyline):

  Cf — PENALIZATION-BAND estimator. It reads the wall gradient that the
  solver itself imposes, so it needs no fit window and no fit depth.

  The IBM coefficient is assembled (ibm.f90 add_neighbor_coeff) as

      coef*Re = sum_neighbours ((d0-d)/d)/d0^2 = sum [1/(d d0) - 1/d0^2]

  over the neighbours B (spacing d0) of a fluid point A that lie inside
  the body, d being the distance from A to the SURFACE along AB. That is
  exactly the term converting "u = 0 at the solid neighbour, distance
  d0" into "u = 0 at the true surface, distance d": the scheme holds
  no-slip by an implicit linear extrapolation through the wall at its
  true sub-cell position. So the band cells -- the fluid points with
  coef > 0, reaching ~1 cell beyond the surface -- carry the wall
  gradient directly. Near a no-slip wall u ~ G dn t, so a through-origin
  least squares of the STAGGERED (un, wn) against (dn t_x, dn t_z) over
  a station's band points gives G = du_t/dn, and Cf = 2 nu G / U_inf^2,
  signed along the TE-ward tangent (attached flow -> positive).

  Each component is used at its own staggered location with its own
  normal distance, so nothing is interpolated and the half-cell offsets
  (un at x-h/2, wn at z-h/2) never enter as a wall-normal error. The
  identity above was verified against the geometry: predicted vs stored
  coef, median ratio 1.0003 over the band; and the scheme's own wall
  distance d = d0/(1 + coef Re d0^2) matches the geometric one.

  Two independent estimators -- a wall-gradient fit over the cells
  OUTSIDE the band, and a viscous-integral estimator that differences
  nowhere near the wall -- agree with this one to ~10 % over 300
  stations. See postProcess/cf_crosscheck.py.

  Cp — least-squares WALL EXTRAPOLATION: p of the near-wall CLEAN cells
  (those outside the penalization band entirely) is fitted linearly in d
  (p = p_w + b d, d <= dmax_cp cells) and read at d = 0. This keeps the
  thin-BL dp/dn ~ 0 behaviour where it holds and still captures the
  finite normal gradient near the curved LE. p_inf is the mean level-0
  pressure in a far-upstream box; Cp = (p_w - p_inf)/(0.5 U^2).

Geometry: the analytic closed-TE NACA 0012 section (the make_airfoil_stl
formula) as a dense polyline -- it matches the case's own STL to 0.013 h
in the solver's BVH wall distance. Stations are equal ARC-LENGTH bins,
not equal polyline-index bins: the polyline is cosine clustered, so
index bins are sub-cell narrow at the LE and span dozens of cells at
mid-chord.

Writes an npz (x/c, Cp, Cf, per side + metadata) and XFOIL-comparable
plain tables cp_<tag>.dat / cf_<tag>.dat (x/c, value; upper then lower
block, '#' headers).
"""
import argparse
import os
import re

import h5py
import numpy as np
from scipy.spatial import cKDTree

XLE, ZLE = 50.0, 48.0        # LE position in the chord-lift plane (--nose)
CHORD = 1.0
# far-upstream p_inf box relative to the nose (must lie on level 0)
PBOX_DX = (-12.0, -10.0)
PBOX_DZ = (-0.5, 0.5)
SOLID_MIN = 1.0e10           # fully solid points carry SOLID/Re, not a graded value


def naca0012_polyline(n=4096):
    """Closed-TE NACA 0012 section as (points, outward normals, tangents,
    x/c, upper-side flag), ordered TE -> upper -> LE -> lower -> TE (CCW
    in the x-z plane)."""
    tt = 0.12
    beta = np.linspace(0.0, np.pi, n//2)
    x = 0.5*(1.0 - np.cos(beta))
    yt = 5.0*tt*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2
                 + 0.2843*x**3 - 0.1036*x**4)
    xs = np.concatenate([x[::-1], x[1:-1]])
    zs = np.concatenate([yt[::-1], -yt[1:-1]])
    upper = np.zeros(xs.size, dtype=bool)
    upper[:n//2] = True
    pts = np.column_stack([XLE + CHORD*xs, ZLE + CHORD*zs])
    # outward normal of the CCW loop = tangent rotated -90 deg
    t = np.gradient(pts, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    nrm = np.column_stack([t[:, 1], -t[:, 0]])
    return pts, nrm, t, xs, upper


def section_levels(field, margin=0.12):
    """Levels carried by blocks near the section, finest first. With a
    refine_body_box the surface is covered by MORE THAN ONE level (the
    NACA nose case: level 12 at the nose, 11 elsewhere), so the caller
    must sweep them and merge."""
    with h5py.File(field, "r") as ff:
        blocks = ff["blocks"][...]
        nb = int(ff.attrs["block_nb_x"])
        nx = int(ff.attrs["nx"]); nz = int(ff.attrs["nz"])
        lx = float(ff.attrs["lx"]); lz = float(ff.attrs["lz"])
        m = np.asarray(ff.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        out = set()
        for ox, oy, oz, l in blocks:
            hx = lx/(nx*2**(int(l)*int(m[0]))); hz = lz/(nz*2**(int(l)*int(m[2])))
            if (XLE - margin < (ox + nb/2)*hx < XLE + CHORD + margin
                    and ZLE - margin < (oz + nb/2)*hz < ZLE + margin):
                out.add(int(l))
    return sorted(out, reverse=True)


def surface_levels(field, margin=0.12):
    """Levels whose LEAF blocks straddle the surface — the levels that
    actually carry wall cells, finest first. The other levels near the
    section are the 2:1 cascade and hold no wall. Normally one; with a
    refine_body_box it is two (the NACA nose case: 12 at the nose, 11
    over the rest of the chord)."""
    pts, nrm, tng, xoc, upper = naca0012_polyline()
    tree = cKDTree(pts)
    with h5py.File(field, "r") as ff:
        blocks = ff["blocks"][...]
        nb = int(ff.attrs["block_nb_x"])
        nx = int(ff.attrs["nx"]); nz = int(ff.attrs["nz"])
        lx = float(ff.attrs["lx"]); lz = float(ff.attrs["lz"])
        m = np.asarray(ff.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
    out = []
    for lev in sorted({int(l) for l in blocks[:, 3]}, reverse=True):
        s = blocks[:, 3] == lev
        hx = lx/(nx*2**(lev*int(m[0]))); hz = lz/(nz*2**(lev*int(m[2])))
        cx = (blocks[s, 0] + nb/2)*hx; cz = (blocks[s, 2] + nb/2)*hz
        near = ((np.abs(cx - (XLE + 0.5*CHORD)) < 0.5*CHORD + margin)
                & (np.abs(cz - ZLE) < margin))
        if not near.any():
            continue
        d, _ = tree.query(np.column_stack([cx[near], cz[near]]))
        if (d <= 0.5*nb*np.hypot(hx, hz)).any():      # centre within a block half-diagonal
            out.append(lev)
    return out


def build_lattice(field, coef, level=None, margin=0.12):
    """Span-averaged un, wn, pn and the IBM coefficients of every block at
    `level` (default: the finest) near the section, as dense (z, x) arrays
    over one window. Unrefined or absent cells are NaN.

    un and wn are the STAGGERED low-face values (un[i] at x_i - h/2,
    wn[k] at z_k - h/2, verified against the solid mask); pn is
    cell-centred. NOTE the geometry tiles are stored (x, y, z) while the
    field datasets are (z, y, x)."""
    with h5py.File(field, "r") as ff, h5py.File(coef, "r") as fc:
        blocks = ff["blocks"][...]
        if not np.array_equal(fc["blocks"][...], blocks):
            raise SystemExit(f"{coef}: blocks table differs from {field}")
        nb = int(ff.attrs["block_nb_x"])
        nx = int(ff.attrs["nx"]); nz = int(ff.attrs["nz"])
        lx = float(ff.attrs["lx"]); lz = float(ff.attrs["lz"])
        m = np.asarray(ff.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        lmax = int(blocks[:, 3].max()) if level is None else int(level)
        hx = lx/(nx*2**(lmax*int(m[0]))); hz = lz/(nz*2**(lmax*int(m[2])))
        sel = [b for b, (ox, oy, oz, l) in enumerate(blocks) if l == lmax
               and XLE - margin < (ox + nb/2)*hx < XLE + CHORD + margin
               and ZLE - margin < (oz + nb/2)*hz < ZLE + margin]
        if not sel:
            return None
        ox0 = min(blocks[b][0] for b in sel); ox1 = max(blocks[b][0] for b in sel) + nb
        oz0 = min(blocks[b][2] for b in sel); oz1 = max(blocks[b][2] for b in sel) + nb
        sh = (oz1 - oz0, ox1 - ox0)
        U, W, P, CU, CW = (np.full(sh, np.nan) for _ in range(5))
        for b in sel:
            ox, oy, oz, l = blocks[b]
            i0 = ox - ox0; k0 = oz - oz0
            U[k0:k0+nb, i0:i0+nb] = ff["un"][b].mean(axis=1)   # (z, y, x) -> (z, x)
            W[k0:k0+nb, i0:i0+nb] = ff["wn"][b].mean(axis=1)
            P[k0:k0+nb, i0:i0+nb] = ff["pn"][b].mean(axis=1)
            t = fc["coef_blocks"][b][1:nb+1, 1:nb+1, 1:nb+1, :].mean(axis=1)
            CU[k0:k0+nb, i0:i0+nb] = t[:, :, 0].T              # tiles are (x, y, z)
            CW[k0:k0+nb, i0:i0+nb] = t[:, :, 2].T
    return dict(U=U, W=W, P=P, CU=CU, CW=CW, ox0=ox0, oz0=oz0, hx=hx, hz=hz)


def far_field_p(field):
    """Mean level-0 pressure in the far-upstream reference box."""
    with h5py.File(field, "r") as f:
        blocks = f["blocks"][...]
        nb = int(f.attrs["block_nb_x"])
        h0x = float(f.attrs["lx"])/int(f.attrs["nx"])
        h0z = float(f.attrs["lz"])/int(f.attrs["nz"])
        px0, px1 = XLE + PBOX_DX[0], XLE + PBOX_DX[1]
        pz0, pz1 = ZLE + PBOX_DZ[0], ZLE + PBOX_DZ[1]
        tot, n = 0.0, 0
        for b, (ox, oy, oz, lev) in enumerate(blocks):
            if lev != 0:
                continue
            xc = (ox + 0.5 + np.arange(nb))*h0x
            zc = (oz + 0.5 + np.arange(nb))*h0z
            if xc[0] > px1 or xc[-1] < px0 or zc[0] > pz1 or zc[-1] < pz0:
                continue
            sx = (xc > px0) & (xc < px1); sz = (zc > pz0) & (zc < pz1)
            if sx.any() and sz.any():
                tot += float(f["pn"][b].mean(axis=1)[np.ix_(sz, sx)].sum())
                n += int(sx.sum()*sz.sum())
    if n == 0:
        raise SystemExit("no level-0 far-field box found for p_inf")
    return tot/n


def surface_geometry(P, valid, tree, pts, nrm, reach):
    """Signed normal distance (+ outside) and nearest polyline index for a
    set of points; NaN / -1 where invalid or out of reach."""
    P2 = P.reshape(-1, 2); v = valid.ravel()
    dn = np.full(v.size, np.nan); ci = np.full(v.size, -1)
    d, i = tree.query(P2[v], distance_upper_bound=reach)
    ok = np.isfinite(d) & (i < pts.shape[0])
    at = np.nonzero(v)[0][ok]
    dn[at] = np.einsum("ij,ij->i", P2[v][ok] - pts[i[ok]], nrm[i[ok]])
    ci[at] = i[ok]
    return dn, ci


def clean_mask(D):
    """Cells outside the penalization band entirely: all four in-plane
    faces carry coef = 0."""
    CU, CW = D["CU"], D["CW"]
    c = np.zeros_like(CU, dtype=bool)
    c[:-1, :-1] = ((CU[:-1, :-1] == 0) & (CU[:-1, 1:] == 0)
                   & (CW[:-1, :-1] == 0) & (CW[1:, :-1] == 0))
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--coef", required=True,
                    help="prepared IBM case file (the Cf estimator reads its "
                         "coefficient tiles to find the penalization band)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--stations", type=int, default=300,
                    help="equal arc-length stations around the section")
    ap.add_argument("--dmax-cp", type=float, default=12.0,
                    help="Cp extrapolation depth in fine cells")
    ap.add_argument("--levels", type=int, default=0,
                    help="how many of the finest levels carry the surface "
                         "(0 = auto-detect, which is what you want)")
    ap.add_argument("--re", type=float, default=4.0e5)
    ap.add_argument("--nose", type=float, nargs=2, default=[50.0, 48.0],
                    help="LE position (x z) in the chord-lift plane")
    ap.add_argument("--plot", default=None)
    a = ap.parse_args()
    global XLE, ZLE
    XLE, ZLE = a.nose
    nu = 1.0/a.re

    p_inf = far_field_p(a.h5)
    pts, nrm, tng, xoc, upper = naca0012_polyline()
    arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    tree = cKDTree(pts)
    nst = a.stations
    edges = np.linspace(0.0, arc[-1], nst + 1)
    mids = np.clip(np.searchsorted(arc, 0.5*(edges[:-1] + edges[1:])), 0, pts.shape[0] - 1)

    def station(ci):
        s = np.full(ci.shape, -1)
        g = ci >= 0
        s[g] = np.clip(np.searchsorted(edges, arc[ci[g]], side="right") - 1, 0, nst - 1)
        return s

    # Both estimators are SUMS over the points of a station, so a surface
    # spanning several levels is handled by accumulating each level's
    # contribution with its own cell size. (Needed since refine_body_box:
    # the NACA nose sits at level 12, the rest of the chord at 11.)
    cf_num = np.zeros(nst); cf_den = np.zeros(nst); cf_n = np.zeros(nst)
    S0 = np.zeros(nst); S1 = np.zeros(nst); S2 = np.zeros(nst)
    T0 = np.zeros(nst); T1 = np.zeros(nst)

    levels = surface_levels(a.h5)
    if a.levels:
        levels = section_levels(a.h5)[:a.levels]
    print(f"{a.h5}: surface carried by level(s) {levels}")
    h_fine = None
    # Cp depth is PHYSICAL, set by the coarsest surface level: the 12 h
    # default is a converged depth (sweep 2/4/8/12/16 h -> -1.737 ..
    # -1.7634), so letting it shrink with h on a refined patch would
    # read an unconverged Cp there (-1.751 instead of -1.775).
    with h5py.File(a.h5, "r") as _f:
        _m = np.asarray(_f.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        h_cp = float(_f.attrs["lx"])/(int(_f.attrs["nx"])*2**(min(levels)*int(_m[0])))
    depth_cp = a.dmax_cp*h_cp
    for lev in levels:
        D = build_lattice(a.h5, a.coef, lev)
        if D is None:
            continue
        h = D["hx"]
        if h_fine is None:
            h_fine = h
        U, W, P, CU, CW = D["U"], D["W"], D["P"], D["CU"], D["CW"]
        clean = clean_mask(D)
        kk, ii = np.mgrid[0:U.shape[0], 0:U.shape[1]]
        xc = (D["ox0"] + ii + 0.5)*D["hx"]; zc = (D["oz0"] + kk + 0.5)*D["hz"]
        reach = depth_cp + 6.0*h
        dnU, ciU = surface_geometry(np.stack([xc - 0.5*D["hx"], zc], -1),
                                    np.isfinite(U), tree, pts, nrm, reach)
        dnW, ciW = surface_geometry(np.stack([xc, zc - 0.5*D["hz"]], -1),
                                    np.isfinite(W), tree, pts, nrm, reach)
        dnP, ciP = surface_geometry(np.stack([xc, zc], -1),
                                    clean & np.isfinite(P), tree, pts, nrm, reach)
        uu, ww, pp = U.ravel(), W.ravel(), P.ravel()
        cu, cw = CU.ravel(), CW.ravel()
        tU = tng[ciU]*np.where(tng[ciU, 0:1] >= 0, 1.0, -1.0)
        tW = tng[ciW]*np.where(tng[ciW, 0:1] >= 0, 1.0, -1.0)
        sU, sW, sP = station(ciU), station(ciW), station(ciP)
        print(f"  level {lev}: {int(np.isfinite(U).sum())} cells near the section, "
              f"h = {h:.3e}")

        # ---- Cf: through-origin fit on the PENALIZATION BAND ----
        for s_, dn_, c_, q_, t_ in ((sU, dnU, cu, uu, tU[:, 0]),
                                    (sW, dnW, cw, ww, tW[:, 1])):
            m = (s_ >= 0) & (dn_ > 0) & (c_ > 0) & (c_ < SOLID_MIN) & (dn_ < 2*h)
            if not m.any():
                continue
            g = dn_[m]*t_[m]
            cf_num += np.bincount(s_[m], weights=q_[m]*g, minlength=nst)
            cf_den += np.bincount(s_[m], weights=g*g, minlength=nst)
            cf_n += np.bincount(s_[m], minlength=nst)
        # ---- Cp: linear wall extrapolation over the clean cells ----
        m = (sP >= 0) & (dnP > 0.05*h) & (dnP <= depth_cp)
        if m.any():
            d_ = dnP[m]; q_ = pp[m]; s_ = sP[m]
            S0 += np.bincount(s_, minlength=nst)
            S1 += np.bincount(s_, weights=d_, minlength=nst)
            S2 += np.bincount(s_, weights=d_*d_, minlength=nst)
            T0 += np.bincount(s_, weights=q_, minlength=nst)
            T1 += np.bincount(s_, weights=d_*q_, minlength=nst)

    if h_fine is None:
        raise SystemExit("no blocks near the section")
    h = h_fine
    out = {k: np.full(nst, np.nan) for k in ("cf", "cp", "xoc", "up", "ncf", "ncp")}
    out["xoc"] = xoc[mids]
    out["up"] = upper[mids].astype(float)
    ok = cf_den > 0.0
    out["cf"][ok] = 2.0*nu*cf_num[ok]/cf_den[ok]
    out["ncf"][ok] = cf_n[ok]
    det = S0*S2 - S1*S1                       # normal equations of p = p_w + b d
    ok = (S0 >= 3) & (np.abs(det) > 0.0)
    out["cp"][ok] = 2.0*((S2[ok]*T0[ok] - S1[ok]*T1[ok])/det[ok] - p_inf)
    out["ncp"][ok] = S0[ok]

    tag = re.sub(r"\.h5$", "", os.path.basename(a.h5))
    outp = a.out or f"cpcf_{tag}.npz"
    np.savez(outp, **out, h=h, p_inf=p_inf, re=a.re)
    good = np.isfinite(out["cf"]).sum(), np.isfinite(out["cp"]).sum()
    print(f"{outp}: {good[0]}/{nst} Cf stations, {good[1]}/{nst} Cp stations")

    for name in ("cp", "cf"):
        path = outp.replace(".npz", f"_{name}.dat")
        with open(path, "w") as f:
            f.write(f"# NACA0012 Re={a.re:.2e} {name} vs x/c "
                    f"(upper block, then lower)\n")
            for side, lab in ((1.0, "upper"), (0.0, "lower")):
                f.write(f"# {lab}\n")
                sel = (out["up"] == side) & np.isfinite(out[name])
                o = np.argsort(out["xoc"][sel])
                for xx, vv in zip(out["xoc"][sel][o], out[name][sel][o]):
                    f.write(f"{xx:.6f} {vv:+.6e}\n")
        print(f"wrote {path}")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
        for side, lab in ((1.0, "upper"), (0.0, "lower")):
            sel = out["up"] == side
            o = np.argsort(out["xoc"][sel])
            ax1.plot(out["xoc"][sel][o], out["cp"][sel][o], label=lab)
            ax2.plot(out["xoc"][sel][o], out["cf"][sel][o], label=lab)
        ax1.invert_yaxis(); ax1.set_xlabel("x/c"); ax1.set_ylabel(r"$C_p$")
        ax2.set_xlabel("x/c"); ax2.set_ylabel(r"$C_f$")
        for ax in (ax1, ax2):
            ax.grid(alpha=0.3); ax.legend()
        fig.tight_layout(); fig.savefig(a.plot, dpi=150)
        print(f"wrote {a.plot}")


if __name__ == "__main__":
    main()
