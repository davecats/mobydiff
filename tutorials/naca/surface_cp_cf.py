#!/usr/bin/env python3
"""Surface Cp and Cf on the immersed NACA 0012 from a block-table field
snapshot (span-y / refine_dims = xz orientation).

  surface_cp_cf.py <field.h5> [--out cpcf_aoaX.npz] [--stations 300]
                   [--dmax-cf 2.5] [--dmax-cp 4.0] [--plot cpcf.png]

Method (per surface station on the ANALYTIC section polyline):

  Cf — least-squares wall gradient USING ONLY FLUID POINTS and the
  EXACT immersed no-slip condition u_t(0) = 0: the tangential velocity
  of the span-averaged near-wall fluid cells is fitted THROUGH THE
  ORIGIN against the wall-normal distance d,
      g = sum(w u_t d) / sum(w d^2),   w = 1/d,
  over cells with d <= dmax_cf fine cells. The zero intercept IS the
  wall boundary condition, so no solid/forced cell value enters. The
  fit is then re-restricted to cells with d+ = d u_tau/nu <= 5 (viscous
  sublayer, one iteration with u_tau = sqrt(nu |g|)) so the linear law
  is actually valid. Cf = 2 nu g / U_inf^2, signed along the TE-ward
  tangent (attached flow -> positive; recirculation -> negative).

  Cp — least-squares WALL EXTRAPOLATION: p of the same near-wall fluid
  cells is fitted linearly in d (p = p_w + b d, d <= dmax_cp cells) and
  read at d = 0. This keeps the thin-BL dp/dn ~ 0 behaviour where it
  holds and still captures the finite normal gradient near the curved
  LE. p_inf is the mean level-0 pressure in a far-upstream box
  (x in [0.8, 1.8], z in [5.5, 6.5]); Cp = (p_w - p_inf)/(0.5 U^2).

Geometry: the analytic closed-TE NACA 0012 section (the make_airfoil_stl
formula) as a dense polyline; every FINEST-level span-averaged fluid
cell within reach of the surface is assigned to its nearest polyline
point (cKDTree) -> signed normal distance + arc station. Solid cells
are identified by the IBM interior signature |u| < 1e-20 (the
penalization damps buried velocities to ~1e-26) and never used.

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


def naca0012_polyline(n=4096):
    """Closed-TE NACA 0012 section as (points, outward normals, x/c,
    upper-side flag), ordered TE -> upper -> LE -> lower -> TE (CCW in
    the x-z plane)."""
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


def load_fine_plane(path):
    """Span-averaged u, w, p of every FINEST-level block, with cell
    centres in the chord-lift plane, plus the far-field p reference."""
    with h5py.File(path, "r") as f:
        blocks = f["blocks"][...]
        nb = int(f.attrs["block_nb_x"])
        nx = int(f.attrs["nx"]); nz = int(f.attrs["nz"])
        lx = float(f.attrs["lx"]); lz = float(f.attrs["lz"])
        mask = np.asarray(f.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        lmax = int(blocks[:, 3].max())
        hx = lx/(nx*2**(lmax*int(mask[0])))
        hz = lz/(nz*2**(lmax*int(mask[2])))
        U, W, P = f["un"], f["wn"], f["pn"]
        K = f["k"] if "k" in f else None

        xs, zs, us, ws, ps, ks = [], [], [], [], [], []
        p_far, n_far = 0.0, 0
        for bid, (ox, oy, oz, lev) in enumerate(blocks):
            if lev == 0:
                # far-field pressure reference box (level 0, upstream)
                h0x = lx/nx; h0z = lz/nz
                xc0 = (ox + 0.5 + np.arange(nb))*h0x
                zc0 = (oz + 0.5 + np.arange(nb))*h0z
                px0, px1 = XLE + PBOX_DX[0], XLE + PBOX_DX[1]
                pz0, pz1 = ZLE + PBOX_DZ[0], ZLE + PBOX_DZ[1]
                if xc0[0] > px1 or xc0[-1] < px0 or zc0[0] > pz1 or zc0[-1] < pz0:
                    continue
                pm = P[bid].mean(axis=1)      # span average -> (z, x)
                sel_x = (xc0 > px0) & (xc0 < px1)
                sel_z = (zc0 > pz0) & (zc0 < pz1)
                if sel_x.any() and sel_z.any():
                    p_far += float(pm[np.ix_(sel_z, sel_x)].sum())
                    n_far += int(sel_x.sum()*sel_z.sum())
                continue
            if lev != lmax:
                continue
            xc = (ox + 0.5 + np.arange(nb))*hx
            zc = (oz + 0.5 + np.arange(nb))*hz
            # quick reject: keep only blocks near the section
            if (xc[-1] < XLE - 0.1 or xc[0] > XLE + CHORD + 0.1 or
                    zc[-1] < ZLE - 0.2 or zc[0] > ZLE + 0.2):
                continue
            u = U[bid][...].mean(axis=1)      # (z, x) span-averaged
            w = W[bid][...].mean(axis=1)
            p = P[bid][...].mean(axis=1)
            kk = K[bid][...].mean(axis=1) if K is not None else np.zeros_like(u) + np.nan
            Z, X = np.meshgrid(zc, xc, indexing="ij")
            xs.append(X.ravel()); zs.append(Z.ravel())
            us.append(u.ravel()); ws.append(w.ravel()); ps.append(p.ravel())
            ks.append(kk.ravel())
    if n_far == 0:
        raise SystemExit("no level-0 far-field box found for p_inf")
    return (np.concatenate(xs), np.concatenate(zs), np.concatenate(us),
            np.concatenate(ws), np.concatenate(ps), np.concatenate(ks),
            p_far/n_far, hx)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--out", default=None)
    ap.add_argument("--stations", type=int, default=300,
                    help="arc stations around the section")
    ap.add_argument("--dmax-cf", type=float, default=2.5,
                    help="Cf fit depth in fine cells (before the d+ <= 5 cut)")
    ap.add_argument("--dmax-cp", type=float, default=4.0,
                    help="Cp extrapolation depth in fine cells")
    ap.add_argument("--re", type=float, default=4.0e5)
    ap.add_argument("--nose", type=float, nargs=2, default=[50.0, 48.0],
                    help="LE position (x z) in the chord-lift plane")
    ap.add_argument("--plot", default=None)
    a = ap.parse_args()
    global XLE, ZLE
    XLE, ZLE = a.nose

    x, z, u, w, p, kfield, p_inf, h = load_fine_plane(a.h5)
    nu = 1.0/a.re
    fluid = (np.abs(u) + np.abs(w)) > 1e-20
    print(f"{a.h5}: {x.size} finest-level cells near the section "
          f"({int(fluid.sum())} fluid), h = {h:.3e}, p_inf = {p_inf:.5f}")

    pts, nrm, tng, xoc, upper = naca0012_polyline()
    tree = cKDTree(pts)
    reach = max(a.dmax_cf, a.dmax_cp)*h + 2*h
    d_all, idx = tree.query(np.column_stack([x, z]),
                            distance_upper_bound=reach + 5*h)
    near = np.isfinite(d_all) & (idx < pts.shape[0]) & fluid
    ci = idx[near]
    dx = np.column_stack([x[near], z[near]]) - pts[ci]
    dn = np.einsum("ij,ij->i", dx, nrm[ci])          # signed normal distance
    ut = (u[near]*tng[ci, 0] + w[near]*tng[ci, 1])   # tangential velocity
    kn = kfield[near]
    # TE-ward sign: tangent with positive x component on both sides
    sgn = np.where(tng[ci, 0] >= 0.0, 1.0, -1.0)
    ut *= sgn
    pn = p[near]
    ok = dn > 0.05*h                                  # outside the surface
    ci, dn, ut, pn, kn = ci[ok], dn[ok], ut[ok], pn[ok], kn[ok]

    # station = arc bin over the polyline index
    nst = a.stations
    bins = np.linspace(0, pts.shape[0], nst + 1).astype(int)
    st_of = np.searchsorted(bins, ci, side="right") - 1

    out = {k: np.full(nst, np.nan) for k in ("cf", "cp", "xoc", "up", "ncf", "ncp", "d0")}
    for s in range(nst):
        sel = st_of == s
        mid = (bins[s] + min(bins[s+1], pts.shape[0] - 1))//2
        out["xoc"][s] = xoc[mid]
        out["up"][s] = 1.0 if upper[mid] else 0.0
        if not sel.any():
            continue
        d, v, q, kv = dn[sel], ut[sel], pn[sel], kn[sel]
        # ---- Cf: d0-gated HYBRID estimator ----
        # (1) the ORIGINAL anchored through-origin fit (u_t(0) = 0 at
        # the analytic polyline, d+ <= 5 iteration) — kept verbatim: it
        # is exact where the discrete wall coincides with the polyline.
        # (2) a robust free-intercept fit u_t = g (d - d0) over a
        # deeper range, which MEASURES the local effective-wall offset
        # d0 (staircase/penalization penetration: ~ -1.2 h at the nose
        # vs ~ -0.3 h mid-chord) and reads the slope from the effective
        # wall. Blend by |d0|: anchored where the wall is where the
        # polyline says (mid-chord), free-intercept where it is not
        # (the nose — anchoring at the wrong wall was the laminar-zone
        # Cf spike/jitter).
        m = d <= a.dmax_cf*h
        g_anc = np.nan
        for _ in range(3):
            if m.sum() < 2:
                break
            wgt = 1.0/d[m]
            g_anc = float(np.sum(wgt*v[m]*d[m])/np.sum(wgt*d[m]**2))
            utau = np.sqrt(nu*abs(g_anc)) if g_anc != 0.0 else 0.0
            m2 = m & (d*utau/nu <= 5.0) if utau > 0 else m
            if m2.sum() == m.sum() or m2.sum() < 2:
                break
            m = m2
        n_anc = m.sum()
        mf = (d <= 1.6*a.dmax_cf*h) & (d >= 0.75*h)
        g_free = np.nan
        d0 = 0.0
        for _ in range(4):
            if mf.sum() < 4:
                break
            wgt = 1.0/d[mf]
            Sw = wgt.sum(); Sd = np.sum(wgt*d[mf]); Sv = np.sum(wgt*v[mf])
            Sdd = np.sum(wgt*d[mf]**2); Sdv = np.sum(wgt*d[mf]*v[mf])
            det = Sw*Sdd - Sd*Sd
            if det <= 0.0:
                break
            g_free = (Sw*Sdv - Sd*Sv)/det
            c0 = (Sdd*Sv - Sd*Sdv)/det
            d0 = -c0/g_free if g_free != 0.0 else 0.0
            r = v[mf] - (g_free*d[mf] + c0)
            sig = np.std(r)
            if sig > 0 and np.any(np.abs(r) > 2.5*sig):
                keep = np.abs(r) <= 2.5*sig
                idx = np.nonzero(mf)[0]
                mf = np.zeros_like(mf); mf[idx[keep]] = True
                continue
            break
        if n_anc >= 2 and np.isfinite(g_anc):
            cf_val = 2.0*nu*g_anc
            # LAMINAR gate: the free-intercept branch is valid (and
            # needed) only where the profile is laminar — linear over
            # the whole fit depth, and the effective-wall offset is the
            # dominant error (the nose staircase). Turbulent stations
            # keep the ORIGINAL anchored estimator verbatim (at y+ > 5
            # samples the free fit aliases buffer-layer curvature into
            # d0). k < 1e-4 U^2 near the wall = laminar (pinned/decayed
            # zones are ~0; turbulent BL carries ~1e-2). No-RANS
            # snapshots have no k dataset -> NaN -> always anchored.
            kmed = np.nanmedian(kv[d <= 1.6*a.dmax_cf*h]) if kv.size else np.nan
            if np.isfinite(g_free) and np.isfinite(kmed) and kmed < 1.0e-4:
                cf_val = 2.0*nu*g_free
                out["d0"][s] = d0
            out["cf"][s] = cf_val
            out["ncf"][s] = n_anc
        # ---- Cp: linear wall extrapolation ----
        mc = d <= a.dmax_cp*h
        if mc.sum() >= 3:
            A = np.column_stack([np.ones(mc.sum()), d[mc]])
            coef, *_ = np.linalg.lstsq(A, q[mc], rcond=None)
            out["cp"][s] = 2.0*(coef[0] - p_inf)
            out["ncp"][s] = mc.sum()
        elif mc.any():
            out["cp"][s] = 2.0*(q[mc][np.argmin(d[mc])] - p_inf)
            out["ncp"][s] = mc.sum()

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
