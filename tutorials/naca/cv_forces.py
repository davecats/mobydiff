#!/usr/bin/env python3
"""Control-volume (momentum-balance) drag and lift from a block-table
snapshot — the AUTHORITATIVE force statistic for the B11 campaign: the
buried interior is REMOVED, so the penalization integral under-reads the
(pressure-dominated) lift by construction (validation/naca0012 README).

  cv_forces.py <field.h5> [--boxes 1.5 2.5 4.0] [--re 4e5]
               [--nose 50 48] [--span-y]

For each control box (margin m around the profile bbox, chord-lift
plane, span-averaged fields, unit span):

  F = - oint [ rho u (u.n) + (p - p_inf) n - tau.n ] dl
  tau = (mu + mu_t) (grad u + grad u^T)   (2D in-plane components)

evaluated on the box border by sampling the reassembled fields on the
border lines (bilinear from the painted window). d/dt terms vanish for
a converged state; box-to-box agreement is the convergence/consistency
check (the boxes should sit in well-resolved regions: with the B11
per-level boxes, margins 1.5-4c lie at Delta = c/384 .. c/96).
C = 2 F / (U^2 c).

p_inf is sampled far upstream on the border's own level to avoid any
global-constant ambiguity.
"""
import argparse

import h5py
import numpy as np


def load_plane(path, x0, x1, z0, z1, span_y=True):
    """Span-averaged u (chordwise), w (lift), p, nut painted on the finest
    lattice inside the window [x0,x1]x[z0,z1] (chord-lift plane)."""
    names = ["un", "wn", "pn", "nut"] if span_y else ["un", "vn", "pn", "nut"]
    with h5py.File(path, "r") as f:
        blocks = f["blocks"][...]
        nb = int(f.attrs["block_nb_x"])
        nx = int(f.attrs["nx"]); ny = int(f.attrs["ny"]); nz = int(f.attrs["nz"])
        lx = float(f.attrs["lx"]); lz = float(f.attrs["lz"])
        mask = np.asarray(f.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        lmax = int(blocks[:, 3].max())
        hx = lx/(nx*2**(lmax*int(mask[0])))
        hz = lz/(nz*2**(lmax*int(mask[2])))
        i0, i1 = int(x0/hx), int(np.ceil(x1/hx))
        j0, j1 = int(z0/hz), int(np.ceil(z1/hz))
        out = {n: np.full((j1-j0, i1-i0), np.nan) for n in names}
        data = {n: f[n] for n in names if n in f}
        for bid, (ox, oy, oz, lev) in enumerate(blocks):
            fx = 2**((lmax-int(lev))*int(mask[0]))
            fz = 2**((lmax-int(lev))*int(mask[2]))
            bx0, bz0 = ox*fx, oz*fz
            if bx0 >= i1 or bx0 + nb*fx <= i0 or bz0 >= j1 or bz0 + nb*fz <= j0:
                continue
            si0 = max(i0, bx0); si1 = min(i1, bx0 + nb*fx)
            sj0 = max(j0, bz0); sj1 = min(j1, bz0 + nb*fz)
            for n in out:
                if n not in data:
                    continue
                row = data[n][bid][...]
                if span_y:
                    m2 = row.mean(axis=1)          # (z, x)
                else:
                    m2 = row.mean(axis=0)
                rep = m2.repeat(fz, axis=0).repeat(fx, axis=1)
                out[n][sj0-j0:sj1-j0, si0-i0:si1-i0] = \
                    rep[sj0-bz0:sj1-bz0, si0-bx0:si1-bx0]
        xc = (np.arange(i0, i1) + 0.5)*hx
        zc = (np.arange(j0, j1) + 0.5)*hz
    if span_y:
        u, w = out["un"], out["wn"]
    else:
        u, w = out["un"], out["vn"]
    return xc, zc, u, w, out["pn"], out.get("nut"), hx, hz


def cv_force(path, nose, margin, re, span_y=True):
    x0 = nose[0] - margin
    x1 = nose[0] + 1.0 + margin
    z0 = nose[1] - margin
    z1 = nose[1] + margin
    pad = 0.15*margin
    xc, zc, u, w, p, nut, hx, hz = load_plane(path, x0-pad, x1+pad, z0-pad, z1+pad,
                                              span_y)
    nu = 1.0/re
    nueff = nu + (nut if nut is not None else 0.0)

    # gradients on the painted lattice (borders sit inside one level, so
    # central differences on the uniform paint are consistent there)
    dudx = np.gradient(u, xc, axis=1); dudz = np.gradient(u, zc, axis=0)
    dwdx = np.gradient(w, xc, axis=1); dwdz = np.gradient(w, zc, axis=0)

    def line(vals, x=None, z=None):
        """Sample a field row/column nearest to the requested border."""
        if x is not None:
            i = int(np.argmin(np.abs(xc - x)))
            return vals[:, i]
        j = int(np.argmin(np.abs(zc - z)))
        return vals[j, :]

    # p_inf: upstream border mean far from the wake
    p_inf = float(np.nanmean(line(p, x=x0)))

    Fx = Fz = 0.0
    zin = (zc >= z0) & (zc <= z1)
    xin = (xc >= x0) & (xc <= x1)
    # west (n = -x) and east (n = +x)
    for xb, sgn in ((x0, -1.0), (x1, +1.0)):
        uu, ww, pp = line(u, x=xb)[zin], line(w, x=xb)[zin], line(p, x=xb)[zin]
        ne = line(nueff, x=xb)[zin] if nut is not None else nu
        txx = 2.0*ne*line(dudx, x=xb)[zin]
        txz = ne*(line(dudz, x=xb)[zin] + line(dwdx, x=xb)[zin])
        Fx -= sgn*np.nansum((uu*uu + (pp - p_inf) - txx))*hz
        Fz -= sgn*np.nansum((uu*ww - txz))*hz
    # south (n = -z) and north (n = +z)
    for zb, sgn in ((z0, -1.0), (z1, +1.0)):
        uu, ww, pp = line(u, z=zb)[xin], line(w, z=zb)[xin], line(p, z=zb)[xin]
        ne = line(nueff, z=zb)[xin] if nut is not None else nu
        tzz = 2.0*ne*line(dwdz, z=zb)[xin]
        txz = ne*(line(dudz, z=zb)[xin] + line(dwdx, z=zb)[xin])
        Fx -= sgn*np.nansum((ww*uu - txz))*hx
        Fz -= sgn*np.nansum((ww*ww + (pp - p_inf) - tzz))*hx
    return 2.0*Fx, 2.0*Fz     # C_D, C_L (q = 0.5, c = 1, unit span)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5", nargs="+")
    ap.add_argument("--boxes", type=float, nargs="+", default=[1.5, 2.5, 4.0])
    ap.add_argument("--re", type=float, default=4.0e5)
    ap.add_argument("--nose", type=float, nargs=2, default=[50.0, 48.0])
    ap.add_argument("--aoa", type=float, default=0.0,
                    help="angle of attack in degrees: rotate the body-axis "
                         "(Fx, Fz) force into wind axes (drag along the "
                         "freestream, lift normal to it)")
    a = ap.parse_args()
    ca, sa = np.cos(np.radians(a.aoa)), np.sin(np.radians(a.aoa))
    for path in a.h5:
        print(f"== {path}")
        for m in a.boxes:
            fx, fz = cv_force(path, a.nose, m, a.re)
            cd = fx*ca + fz*sa      # wind axes
            cl = fz*ca - fx*sa
            print(f"  CV margin {m:4.1f} c: C_D = {cd:+.5f}   C_L = {cl:+.5f}")


if __name__ == "__main__":
    main()
