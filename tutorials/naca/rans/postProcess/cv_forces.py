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

FLUX-EXACT (finite-volume) evaluation: each border is SNAPPED to the
nearest cell-face coordinate of the coarsest level it crosses (a coarse
face is a face at every finer level — dyadic lines), and the fluxes are
integrated per leaf block at its NATIVE level using the staggered face
velocities directly (un on x-faces, wn/vn on lift-faces, span-averaged;
no repainting — memory stays per-block, the L11 OOM landmine of the old
painted implementation does not exist here). Face areas are the exact
level-dependent spacings, so a uniform field integrates to EXACTLY zero
around the closed box. p and tangential velocities are interpolated one
half-cell to the face (central where both cells are stored in the same
block, second-order one-sided at block edges); the viscous term uses
central/one-sided differences at the face's own spacing. On a border
that coincides with a 2:1 level line, the face DOFs come from the block
that STORES that face (the low-side owner's reconciled value — the
conservative copy). d/dt terms vanish for a converged state;
box-to-box agreement is the consistency check. C = 2 F / (U^2 c).

p_inf is the mean face pressure of the upstream (west) border.
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


class BlockTable:
    """Leaf-table geometry + lazy per-block dataset access in the
    chord-lift plane. axis c = chord (x), axis l = lift (z when span_y,
    y otherwise); the span axis is averaged out on read."""

    def __init__(self, path, span_y=True):
        self.f = h5py.File(path, "r")
        f = self.f
        self.blocks = f["blocks"][...]
        self.nb = int(f.attrs["block_nb_x"])
        mask = np.asarray(f.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        n = {"x": int(f.attrs["nx"]), "y": int(f.attrs["ny"]), "z": int(f.attrs["nz"])}
        L = {"x": float(f.attrs["lx"]), "y": float(f.attrs["ly"]), "z": float(f.attrs["lz"])}
        lift = "z" if span_y else "y"
        span = "y" if span_y else "z"
        self.h0 = {"c": L["x"]/n["x"], "l": L[lift]/n[lift]}
        self.mask = {"c": int(mask[0]), "l": int(mask[2] if span_y else mask[1])}
        # number of span-direction block layers per level: contributions
        # from multiple layers at the same (c, l) must be AVERAGED (the
        # span axis is homogeneous), not summed.
        self.n_span = n[span]
        self.mask_span = int(mask[1] if span_y else mask[2])
        # dataset row layout is (z, y, x); the span axis to average and the
        # lift axis inside the row:
        self.span_axis = 1 if span_y else 0
        self.names = {"un": "un", "wn": "wn" if span_y else "vn",
                      "pn": "pn", "nut": "nut"}
        # block origins in the (c, l) plane per axis in LEVEL cells
        oc = self.blocks[:, 0]
        ol = self.blocks[:, 2] if span_y else self.blocks[:, 1]
        self.lev = self.blocks[:, 3].astype(int)
        self.oc, self.ol = oc.astype(int), ol.astype(int)
        self.has_nut = "nut" in f
        self._cache = {}

    def h(self, ax, lev):
        return self.h0[ax]/2**(lev*self.mask[ax])

    def nspan(self, lev):
        return max(1, (self.n_span*2**(lev*self.mask_span))//self.nb)

    def row(self, name, bid):
        """Span-averaged (l, c) tile of a dataset for one block."""
        key = (name, bid)
        if key not in self._cache:
            d = self.f[self.names[name]][bid][...]
            self._cache[key] = d.mean(axis=self.span_axis)  # -> (l, c) or (c ...)
            if self.span_axis == 0:
                # (z,y,x) span=z averaged -> (y, x) = (l, c): already (l, c)
                pass
        return self._cache[key]


def snap_border(bt, ax, coord, t0, t1):
    """Snap `coord` (a border normal to axis `ax`) to the nearest face of
    the coarsest level among blocks the border crosses (a coarse face is
    a face at every finer level)."""
    for _ in range(2):
        lev_min = None
        for b in range(len(bt.blocks)):
            l = bt.lev[b]
            hc, hl = bt.h(ax, l), bt.h("l" if ax == "c" else "c", l)
            o = (bt.oc if ax == "c" else bt.ol)[b]
            ot = (bt.ol if ax == "c" else bt.oc)[b]
            if not (o*hc <= coord <= (o + bt.nb)*hc):
                continue
            if (ot + bt.nb)*hl < t0 or ot*hl > t1:
                continue
            lev_min = l if lev_min is None else min(lev_min, l)
        if lev_min is None:
            raise SystemExit(f"border {ax}={coord}: no blocks crossed")
        hs = bt.h(ax, lev_min)
        snapped = round(coord/hs)*hs
        if abs(snapped - coord) < 1e-12:
            return snapped
        coord = snapped
    return coord


def cv_force(path, nose, margin, re, span_y=True, verbose=False):
    bt = BlockTable(path, span_y)
    nu = 1.0/re
    nb = bt.nb
    c0, c1 = nose[0] - margin, nose[0] + 1.0 + margin
    l0, l1 = nose[1] - margin, nose[1] + margin
    c0 = snap_border(bt, "c", c0, l0, l1)
    c1 = snap_border(bt, "c", c1, l0, l1)
    l0 = snap_border(bt, "l", l0, c0, c1)
    l1 = snap_border(bt, "l", l1, c0, c1)
    if verbose:
        print(f"    snapped box: c [{c0:.6f}, {c1:.6f}]  l [{l0:.6f}, {l1:.6f}]")

    def faces_on_border(ax, coord, t0, t1):
        """Yield (bid, iloc, rows): blocks whose stored `ax` faces include
        the border, with the face's in-block index (1-based, 1..nb) and
        the tangential row indices inside [t0, t1]."""
        out = []
        for b in range(len(bt.blocks)):
            lev = bt.lev[b]
            hn = bt.h(ax, lev)
            ht = bt.h("l" if ax == "c" else "c", lev)
            o = (bt.oc if ax == "c" else bt.ol)[b]
            ot = (bt.ol if ax == "c" else bt.oc)[b]
            fi = coord/hn - o          # face index within the block, 0..nb
            i = int(round(fi))
            if abs(fi - i) > 1e-9 or not (0 <= i <= nb - 1):
                # stored faces are local 0..nb-1 (u(i) = west face of cell i);
                # the block whose WEST edge is the border stores it at i=0.
                continue
            tc = (ot + 0.5 + np.arange(nb))*ht
            rows = np.nonzero((tc > t0) & (tc < t1))[0]
            if rows.size:
                out.append((b, i + 1, rows, lev, hn, ht))
        return out

    def face_values(ax, coord, t0, t1):
        """Per tangential cell on the border: (dl, u_n, ut_face, p_face,
        nue_face, dun_dn, dun_dt, dut_dn). ax = normal axis."""
        un_name = "un" if ax == "c" else "wn"
        ut_name = "wn" if ax == "c" else "un"
        vals = []
        for b, iloc, rows, lev, hn, ht in faces_on_border(ax, coord, t0, t1):
            U = bt.row(un_name, b)     # (l, c)
            T = bt.row(ut_name, b)
            P = bt.row("pn", b)
            NUT = bt.row("nut", b) if bt.has_nut else None

            def A(M, t, n):
                """M[tangential t, normal n] respecting axis order (l, c)."""
                return M[t, n] if ax == "c" else M[n, t]

            for t in rows:
                i = iloc - 1           # 0-based normal index of face/cell east
                u_n = A(U, t, i)
                # neighbours along the normal (cells i-1 | i share the face)
                if i >= 1:
                    p_face = 0.5*(A(P, t, i-1) + A(P, t, i))
                    ut_face = 0.5*(A(T, t, i-1) + A(T, t, i))
                    nue = nu + (0.5*(A(NUT, t, i-1) + A(NUT, t, i)) if NUT is not None else 0.0)
                    dut_dn = (A(T, t, i) - A(T, t, i-1))/hn
                else:
                    p_face = 1.5*A(P, t, 0) - 0.5*A(P, t, 1)
                    ut_face = 1.5*A(T, t, 0) - 0.5*A(T, t, 1)
                    nue = nu + (max(1.5*A(NUT, t, 0) - 0.5*A(NUT, t, 1), 0.0) if NUT is not None else 0.0)
                    dut_dn = (A(T, t, 1) - A(T, t, 0))/hn
                # normal-velocity gradients at the face (native spacing)
                if 1 <= i <= nb - 2:
                    dun_dn = (A(U, t, i+1) - A(U, t, i-1))/(2.0*hn)
                elif i == 0:
                    dun_dn = (A(U, t, 1) - A(U, t, 0))/hn
                else:
                    dun_dn = (A(U, t, i) - A(U, t, i-1))/hn
                if 1 <= t <= nb - 2:
                    dun_dt = (A(U, t+1, i) - A(U, t-1, i))/(2.0*ht)
                elif t == 0:
                    dun_dt = (A(U, 1, i) - A(U, 0, i))/ht
                else:
                    dun_dt = (A(U, t, i) - A(U, t-1, i))/ht
                vals.append((ht/bt.nspan(lev), u_n, ut_face, p_face, nue,
                             dun_dn, dun_dt, dut_dn))
        return vals

    # p_inf: mean face pressure of the upstream border
    west = face_values("c", c0, l0, l1)
    p_inf = float(np.mean([v[3] for v in west]))

    Fc = Fl = 0.0
    for ax, coord, sgn, cached in (("c", c0, -1.0, west), ("c", c1, +1.0, None),
                                   ("l", l0, -1.0, None), ("l", l1, +1.0, None)):
        vals = cached if cached is not None else face_values(
            ax, coord, l0 if ax == "c" else c0, l1 if ax == "c" else c1)
        fn = ft = 0.0
        for dl, u_n, ut, p, nue, dun_dn, dun_dt, dut_dn in vals:
            tau_nn = 2.0*nue*dun_dn
            tau_nt = nue*(dun_dt + dut_dn)
            fn += (u_n*u_n + (p - p_inf) - tau_nn)*dl
            ft += (u_n*ut - tau_nt)*dl
        if ax == "c":
            Fc -= sgn*fn
            Fl -= sgn*ft
        else:
            Fl -= sgn*fn
            Fc -= sgn*ft
    bt.f.close()
    return 2.0*Fc, 2.0*Fl     # C_D, C_L (q = 0.5, c = 1, unit span)


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
