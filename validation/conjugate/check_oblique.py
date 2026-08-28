#!/usr/bin/env python3
"""C2 gate 1: the oblique plane interface, measured against Section 3.

    ./check_oblique.py flux  <case.h5> --theta T --kappa K --ratio R [...]
    ./check_oblique.py field <field.h5> --theta T --kappa K --ratio R [...]

THE MANUFACTURED SOLUTION. With n the unit normal of the interface plane, t
the in-plane tangential direction and xi = n.(x - x0) the signed distance,

    T = q_n xi/k(side) + A (t.x),        k = 1 (fluid, xi > 0) | kappa_s

is an exact conjugate solution for ANY kappa_s: it is linear in each
material, continuous across the plane, and its normal flux k d_n T = q_n is
continuous there. The CONTROLLED RATIO of Section 3 is

    |grad_t T|/|d_n T| (fluid side) = A/q_n = r,

which is what the whole increment turns on: the argument for dropping the
tangential term is that r ~ 1e-2 in a DNS thermal layer.

`flux` IS THE MEASUREMENT. Per cut face it forms, from the case file's own
phi (= +-dwall_blocks signed by coef_p_blocks) and the ANALYTIC field:

    q_n^base = (T_R - T_L)/(a R)                        the baseline's implied
    q_n^corr = (T_R - T_L - h_d s_t^num)/(a R)          normal flux, Eq. (qn)
    e_face   = h_d s_t/(T_R - T_L)                      the C2 indicator

with R the series resistance and s_t^num the solver's own discrete tangential
estimate. Because the exact field is piecewise linear, T_R - T_L = a q_n R +
h_d s_t EXACTLY, so the baseline's relative error reduces to Section 3's
closed form e/(1 - e) with NO discretisation of its own -- the check that the
prediction and the measurement agree is then a real check of both.

Nothing here runs the solver: the field is analytic and phi comes from the
case file, so this measures the SCHEME, not a transient. The solver's own
e_face (printed by [scalar] indicator_interval) is the cross-check that the
Fortran computes what this does.

`field` is the companion BVP statement: the r = infinity case (pure
tangential, q_n = 0) is the ONE ratio whose boundary data is constant on
every face -- grad T is then the same in both materials, so all six faces
carry the exact constant Neumann value -- and it is therefore the only ratio
at which the solved field can be compared with the analytic one. The
solution is defined up to a constant (pure Neumann), so the mean offset is
removed before comparing.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scalar"))
from scalar_tools import BlockGeometry                       # noqa: E402

SOLID_THRESHOLD = 1.0e20
MIN_COSINE = 5.0e-2          # scalar.f90 CONJ_MIN_COSINE
MIN_GRADPHI = 1.0e-2         # scalar.f90 CONJ_MIN_GRADPHI


class Plane:
    """The manufactured oblique-plane solution."""

    def __init__(self, theta_deg, x0, y0, kappa, q_n, amp):
        th = math.radians(theta_deg)
        self.n = np.array([-math.sin(th), math.cos(th), 0.0])
        self.t = np.array([math.cos(th), math.sin(th), 0.0])
        self.p0 = np.array([x0, y0, 0.0])
        self.kappa = kappa
        self.q_n = q_n
        self.amp = amp

    def xi(self, x, y, z):
        return (x - self.p0[0]) * self.n[0] + (y - self.p0[1]) * self.n[1]

    def temperature(self, x, y, z):
        xi = self.xi(x, y, z)
        k = np.where(xi < 0.0, self.kappa, 1.0)
        return self.q_n * xi / k + self.amp * ((x - self.p0[0]) * self.t[0]
                                               + (y - self.p0[1]) * self.t[1])

    def gradient(self, side_fluid):
        """The (constant) gradient in one material."""
        k = 1.0 if side_fluid else self.kappa
        return (self.q_n / k) * self.n + self.amp * self.t

    def s_t(self, d):
        """The exact tangential scalar of a d-direction arm. It is the same
        on both sides -- that continuity is what makes the correction
        evaluable at all -- so the fluid gradient serves."""
        g = self.gradient(True)
        return g[d] - self.n[d] * float(self.n @ g)


def load_phi(case):
    """phi and the cell centres of a SINGLE-LEVEL case file, assembled into
    one global array with a ghost layer: index m holds cell m-1, so 0 and
    n+1 are the ghosts. The ghost-inclusive tiles supply every ghost,
    corners included, straight from the geometry."""
    with h5py.File(case, "r") as h5:
        blocks = h5["blocks"][...]
        nb = int(h5.attrs["block_nb"])
        node = [h5["x_nodes"][...], h5["y_nodes"][...], h5["z_nodes"][...]]
        dwall = h5["dwall_blocks"][...]
        coef = h5["coef_p_blocks"][...]
    if int(blocks[:, 3].max()) != 0:
        raise SystemExit("single-level case files only")

    n = [len(a) - 1 for a in node]
    phi = np.zeros((n[2] + 2, n[1] + 2, n[0] + 2))
    for bid in range(blocks.shape[0]):
        ox, oy, oz = (int(v) for v in blocks[bid, :3])
        solid = np.abs(coef[bid]) > SOLID_THRESHOLD
        tile = np.where(solid, -np.maximum(dwall[bid], 1e-300), dwall[bid])
        # CASE-FILE tiles are stored (i, j, k) -- unlike the FIELD datasets,
        # which are (k, j, i). C1's weight checker could not tell the two
        # apart (its slab is x/z symmetric); an oblique plane can.
        phi[oz:oz + nb + 2, oy:oy + nb + 2, ox:ox + nb + 2] = tile.transpose(2, 1, 0)

    centres = []
    for d in range(3):
        ext = np.concatenate(([2 * node[d][0] - node[d][1]], node[d],
                              [2 * node[d][-1] - node[d][-2]]))
        centres.append(0.5 * (ext[:-1] + ext[1:]))      # 0..n+1, ghosts included
    return phi, centres


def face_measure(phi, centres, plane, d, interior_margin=0):
    """Every cut face normal to direction d, as flat arrays of the quantities
    Section 3 is written in. Array axes are (k, j, i), so direction
    d = 0, 1, 2 is axis 2, 1, 0."""
    axis = 2 - d
    tang = [t for t in range(3) if t != d]
    shape = phi.shape
    xc, yc, zc = centres
    gk, gj, gi = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]),
                             np.arange(shape[2]), indexing="ij")
    temp = plane.temperature(xc[gi], yc[gj], zc[gk])
    cd = [float(np.diff(c)[0]) for c in centres]     # uniform grids only
    hd = cd[d]
    invd = 1.0 / hd

    def sl(off, ax=None, shift=0):
        """The trimmed core, shifted by `off` along the face-normal axis and
        by `shift` along the array axis `ax`."""
        s = [slice(1, shape[m] - 1) for m in range(3)]
        s[axis] = slice(1 + off, shape[axis] - 1 + off)
        if ax is not None:
            base = s[ax]
            s[ax] = slice(base.start + shift, base.stop + shift)
        return tuple(s)

    # The face between the cell at index m-1 (L) and m (R), over the trimmed
    # core: the trim is what makes the +-1 tangential stencil addressable and
    # costs only the outermost ghost layer.
    lo, hi = sl(-1), sl(0)

    def tangential(field, direction):
        ax = 2 - direction
        return 0.5 * ((field[sl(-1, ax, 1)] - field[sl(-1, ax, -1)])
                      + (field[sl(0, ax, 1)] - field[sl(0, ax, -1)])) \
            / (2.0 * cd[direction])

    gtd = (temp[hi] - temp[lo]) * invd
    gt1 = tangential(temp, tang[0])
    gt2 = tangential(temp, tang[1])
    gpd = (phi[hi] - phi[lo]) * invd
    gp1 = tangential(phi, tang[0])
    gp2 = tangential(phi, tang[1])
    gn = np.sqrt(gpd**2 + gp1**2 + gp2**2)
    ok = gn >= MIN_GRADPHI
    safe = np.where(ok, gn, 1.0)
    dn = (gpd * gtd + gp1 * gt1 + gp2 * gt2) / safe
    nd = np.where(ok, gpd / safe, 0.0)
    st_num = np.where(ok, gtd - nd * dn, 0.0)

    pl, pr = phi[lo], phi[hi]
    cut = (pl < 0.0) != (pr < 0.0)
    if interior_margin > 0:
        keep = np.zeros_like(cut)
        m = interior_margin
        keep[m:-m, m:-m, m:-m] = True
        cut &= keep
    if not cut.any():
        return None

    pl, pr = pl[cut], pr[cut]
    # a = n.e_d is SIGNED (the note's n runs from the L cell to the R one),
    # and only the grazing guard takes its magnitude. Getting this wrong
    # flips the sign of the implied normal flux and reads as a 200 % error.
    a_cos = (pl - pr) * invd
    w = np.where(np.abs(a_cos) < MIN_COSINE, 0.5, pl / (pl - pr))
    kl = np.where(pl < 0.0, plane.kappa, 1.0)
    kr = np.where(pr < 0.0, plane.kappa, 1.0)
    res = hd * (w / kl + (1.0 - w) / kr)             # the series resistance
    return dict(d=d, n=int(cut.sum()), hd=hd, a=a_cos, w=w, res=res,
                kface=hd / res,
                kloc=np.where(pl + pr < 0.0, plane.kappa, 1.0),
                kl=kl, kr=kr, nd=nd[cut], gtd=gtd[cut],
                dt=(temp[hi] - temp[lo])[cut], st_num=st_num[cut],
                st_exact=plane.s_t(d))


def debias(m):
    """The SHIPPED s_t: the raw projection, corrected for the jump the arm
    difference straddles (scalar.f90 conjugate_tangential).

    For a piecewise-linear field with a kink mu in the normal derivative the
    face-centred gradient is G_L + mu[n/2 + (1/2 - w) n_d e_d]: the two
    TANGENTIAL differences weight the sides (1/2, 1/2), so the projection
    removes their share of the jump exactly -- that is the note's argument --
    but the ARM difference weights them (w, 1-w), and what survives is
        s_t^raw = s_t + mu n_d (1/2 - w)(1 - n_d^2).
    Estimating mu from the face's own normal flux makes n_d cancel and
    leaves a linear equation whose solution is below. 1 - c >= 1/2 always.
    """
    c = m["kface"] * (1.0 / m["kr"] - 1.0 / m["kl"]) \
        * (0.5 - m["w"]) * (1.0 - m["nd"] * m["nd"])
    return (m["st_num"] - c * m["gtd"]) / (1.0 - c)


def collect(case, plane, margin=0):
    phi, centres = load_phi(case)
    out = []
    for d in range(3):
        m = face_measure(phi, centres, plane, d, margin)
        if m is not None:
            out.append(m)
    if not out:
        raise SystemExit("no cut faces found")
    return out, phi, centres


def cmd_flux(a):
    plane = Plane(a.theta, a.x0, a.y0, a.kappa, a.q_n, a.amp)
    parts, phi, centres = collect(a.case, plane, a.margin)

    # How good the geometry itself is, on an OBLIQUE plane -- C1's weight
    # gate only ever saw a grid-aligned one. Measured in the INTERFACE BAND,
    # which is the only place phi is read: the STL is a closed body, so far
    # from the plane the nearest surface is one of its padding faces and
    # dwall is (correctly) the distance to that instead.
    xc, yc, zc = centres
    gz, gy, gx = np.meshgrid(np.arange(len(zc)), np.arange(len(yc)),
                             np.arange(len(xc)), indexing="ij")
    xi = plane.xi(xc[gx], yc[gy], zc[gz])
    band = np.abs(xi) < 3.0 * float(np.diff(xc)[0])
    phi_err = float(np.abs(phi - xi)[band].max())

    fr = []          # relative error of the interface-NORMAL flux
    fc = []
    ef = []          # the indicator, exact and discrete
    efn = []
    fe = []          # relative error of the FACE flux
    fec = []
    fd = []          # ...with the DE-BIASED s_t (see debias())
    scale = 0.0
    nfaces = 0
    # The note's q_n is the flux along n running from the L cell to the R
    # one; the manufactured field's q_n is along the interface normal
    # pointing into the FLUID. The two differ by a sign, and the difference
    # is face-independent because `a` carries the orientation.
    qn_ref = -a.q_n
    for m in parts:
        nfaces += m["n"]
        st = m["st_exact"]
        qn_base = m["dt"] / (m["a"] * m["res"])
        st_db = debias(m)
        qn_corr = (m["dt"] - m["hd"] * m["st_num"]) / (m["a"] * m["res"])
        qn_db = (m["dt"] - m["hd"] * st_db) / (m["a"] * m["res"])
        if qn_ref != 0.0:
            fr.append(qn_base / qn_ref - 1.0)
            fc.append(qn_corr / qn_ref - 1.0)
            fd.append(qn_db / qn_ref - 1.0)
        ef.append(m["hd"] * st / m["dt"])
        efn.append(m["hd"] * st_db / m["dt"])
        # The exact face flux is k times the directional derivative of the
        # material the face MIDPOINT is in -- no sign conventions involved.
        gsolid = plane.gradient(False)[m["d"]]
        gfluid = plane.gradient(True)[m["d"]]
        f_exact = m["kloc"] * np.where(m["kloc"] > 1.0, gsolid, gfluid)
        f_base = m["kface"] * m["dt"] / m["hd"]
        f_corr = f_base + st_db * (m["kloc"] - m["kface"])
        scale = max(scale, float(np.abs(f_exact).max()))
        fe.append(f_base - f_exact)
        fec.append(f_corr - f_exact)

    # The solver's indicator drops faces whose denominator is below 1e-6 of
    # the largest one -- a cut face carrying no normal flux has e_face =
    # infinity and means nothing -- so the checker uses the same floor, or
    # the two statistics would not be comparable.
    denom = np.concatenate([np.atleast_1d(np.abs(m["dt"] / m["hd"])) for m in parts])
    keep = denom > 1.0e-6 * denom.max()

    def stat(chunks, norm=1.0, mask=True):
        if not chunks:
            return float("nan"), float("nan")
        v = np.concatenate([np.atleast_1d(c) for c in chunks]) / norm
        if mask and v.size == keep.size:
            v = v[keep]
        return float(np.abs(v).max()), float(np.sqrt(np.mean(v * v)))

    ef_max, ef_rms = stat(ef)
    efn_max, efn_rms = stat(efn)
    qb_max, qb_rms = stat(fr)
    qc_max, qc_rms = stat(fc)
    qd_max, qd_rms = stat(fd)
    fb_max, fb_rms = stat(fe, scale)
    fc_max, fc_rms = stat(fec, scale)
    # Section 3's closed form, from the exact quantities alone.
    pred = np.concatenate([np.atleast_1d(e / (1.0 - e)) for e in ef])[keep]
    pred_max = float(np.abs(pred).max())

    ratio = float("inf") if a.q_n == 0.0 else a.amp / a.q_n
    print(f"   theta = {a.theta:g} deg  kappa_s = {a.kappa:g}  ratio = {ratio:g}"
          f"  h = {parts[0]['hd']:.6g}  cut faces = {nfaces}"
          f"  ({int((~keep).sum())} below the denominator floor)")
    print(f"   |phi - exact| (case file, oblique plane) = {phi_err:.3e}")
    print(f"   e_face                    max = {ef_max:.4e}   rms = {ef_rms:.4e}")
    print(f"   e_face (solver s_t)       max = {efn_max:.4e}   rms = {efn_rms:.4e}")
    print(f"   q_n rel err  BASELINE     max = {qb_max:.4e}   rms = {qb_rms:.4e}"
          f"   (Section 3 predicts {pred_max:.4e})")
    print(f"   q_n rel err  RAW s_t      max = {qc_max:.4e}   rms = {qc_rms:.4e}"
          "   (the note's construction 1, NOT shipped)")
    print(f"   q_n rel err  CORRECTED    max = {qd_max:.4e}   rms = {qd_rms:.4e}")
    print(f"   face flux rel err BASE    max = {fb_max:.4e}   rms = {fb_rms:.4e}")
    print(f"   face flux rel err CORR    max = {fc_max:.4e}   rms = {fc_rms:.4e}")
    if a.emit:
        with open(a.emit, "a") as fh:
            fh.write("%g %g %g %g %d %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n"
                     % (a.theta, a.kappa, ratio, parts[0]["hd"], nfaces,
                        ef_max, ef_rms, qb_max, qb_rms, qd_max, qd_rms,
                        fb_max, fb_rms, phi_err))
    ok = math.isfinite(ef_max)
    if a.tolerance is not None:
        ok = ok and fc_max <= a.tolerance
        print("   PASS" if ok else f"   FAIL (corrected face-flux tolerance {a.tolerance:g})")
    return 0 if ok else 1


def cmd_residual(a):
    """The discrete flux-divergence residual of each scheme on the EXACT
    field -- the LOCAL TRUNCATION ERROR, which is what the cell balance sees.

    The manufactured solution has zero divergence, so anything left is the
    scheme's. This is the mechanism behind the BVP result, and it is also the
    measurement that says what a fix has to do:

      * `k_face` (C1) and `k_loc` (C2 as shipped) both leave an O(1/h)
        residual at cut cells -- and `k_loc`, which makes the POINTWISE face
        flux exact, leaves the LARGER of the two;
      * `k_area`, the AREA-weighted mean over the face, applied at EVERY face
        the interface clips, makes it vanish. The reason is a two-line
        divergence-theorem argument: with the exact face AVERAGES the cell
        balance IS the exact surface integral, which is zero.
    """
    plane = Plane(a.theta, a.x0, a.y0, a.kappa, a.q_n, a.amp)
    phi, centres = load_phi(a.case)
    shape = phi.shape
    gk, gj, gi = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    temp = plane.temperature(centres[0][gi], centres[1][gj], centres[2][gk])
    cd = [float(np.diff(c)[0]) for c in centres]

    print(f"   theta = {a.theta:g}  kappa_s = {a.kappa:g}  h = {cd[0]:.6g}"
          f"  ratio = {float('inf') if a.q_n == 0 else a.amp / a.q_n:g}")
    for mode, label in (("base", "C1 baseline   k_face"),
                        ("mid",  "C2 shipped    k_loc "),
                        ("area", "PROPOSED      k_area")):
        total = np.zeros(shape)
        band = np.zeros(shape, dtype=bool)
        defined = np.ones(shape, dtype=bool)
        for d in range(3):
            axis = 2 - d
            F, cut = face_flux_field(phi, temp, centres, plane, d, mode)
            # the low-face flux of cell m, as a full-size array: the high face
            # of m is the low face of m+1. Rolling keeps every direction on
            # the SAME cell grid -- combining per-direction arrays of
            # different shapes silently drops directions.
            full = np.full(shape, np.nan)
            fcut = np.zeros(shape, dtype=bool)
            idx = tuple(slice(1, shape[m] - 1) for m in range(3))
            full[idx] = F
            fcut[idx] = cut
            hi = np.roll(full, -1, axis=axis)
            total = total + (hi - full) / cd[d]
            band = band | fcut | np.roll(fcut, -1, axis=axis)
            defined = defined & np.isfinite(full) & np.isfinite(hi)
        core = tuple(slice(2, shape[m] - 2) for m in range(3))
        r, b, ok = total[core], band[core], defined[core]
        cutm, unc = b & ok, (~b) & ok
        print(f"   {label}   cut-cell |div| max = {np.abs(r[cutm]).max():.4e}"
              f"   rms = {r[cutm].std():.4e}"
              f"   (uncut cells max {np.abs(r[unc]).max():.2e})")
    return 0


def face_area_fraction(phic, n1, n2, h1, h2):
    """Fraction of an h1 x h2 face whose phi > 0, for a plane through the
    face centre with in-plane normal components (n1, n2) -- the standard 2D
    plane-in-rectangle form (Scardovelli & Zaleski). Validated against brute
    force to the quadrature error (9e-4 at 400^2 samples).

    It needs only phi at the face centre (= the mean of the two cell values,
    exact for a plane) and the in-plane components of the face-centred
    grad phi -- both already formed for s_t. NO new data.
    """
    a1 = np.abs(n1) * h1
    a2 = np.abs(n2) * h2
    s = phic + 0.5 * (a1 + a2)                  # distance from the lowest corner
    amin, amax, tot = np.minimum(a1, a2), np.maximum(a1, a2), a1 + a2
    sa = np.where(amax > 0.0, amax, 1.0)
    pr = np.where(a1 * a2 > 0.0, a1 * a2, 1.0)
    f = np.where(s <= 0.0, 0.0,
        np.where(s >= tot, 1.0,
        np.where(s <= amin, s * s / (2.0 * pr),
        np.where(s <= amax, (2.0 * s - amin) / (2.0 * sa),
                 1.0 - (tot - s) ** 2 / (2.0 * pr)))))
    f = np.where(tot <= 0.0, np.where(phic > 0.0, 1.0, 0.0), f)   # face || interface
    return np.clip(f, 0.0, 1.0)


def face_flux_field(phi, temp, centres, plane, d, mode):
    """Every scheme's flux on EVERY low face along direction d.

    mode = "base"  C1: F = k_face (T_R - T_L)/h
           "mid"   C2 as shipped: + s_t (k_loc - k_face), k_loc = the material
                   at the face MIDPOINT, at marker-cut faces
           "area"  the proposal: + s_t (k_area - k_face) at EVERY face, with
                   k_area the AREA-weighted mean over the face -- which is
                   what the cell balance actually wants (see the README).
    """
    axis = 2 - d
    tang = [t for t in range(3) if t != d]
    shape = phi.shape
    cd = [float(np.diff(c)[0]) for c in centres]
    hd = cd[d]
    invd = 1.0 / hd

    def sl(off, ax=None, shift=0):
        s = [slice(1, shape[m] - 1) for m in range(3)]
        s[axis] = slice(1 + off, shape[axis] - 1 + off)
        if ax is not None:
            s[ax] = slice(s[ax].start + shift, s[ax].stop + shift)
        return tuple(s)

    lo, hi = sl(-1), sl(0)
    pl_, pr_ = phi[lo], phi[hi]
    cut = (pl_ < 0.0) != (pr_ < 0.0)
    kl = np.where(pl_ < 0.0, plane.kappa, 1.0)
    kr = np.where(pr_ < 0.0, plane.kappa, 1.0)
    gap = pl_ - pr_
    w = np.where(np.abs(gap) * invd < MIN_COSINE, 0.5,
                 pl_ / np.where(gap == 0.0, 1.0, gap))
    kface = np.where(cut, 1.0 / (w / kl + (1.0 - w) / kr), kl)
    gtd = (temp[hi] - temp[lo]) * invd
    F = kface * gtd
    if mode != "base":
        def tg(fld, direction):
            ax = 2 - direction
            return 0.5 * ((fld[sl(-1, ax, 1)] - fld[sl(-1, ax, -1)])
                          + (fld[sl(0, ax, 1)] - fld[sl(0, ax, -1)])) / (2.0 * cd[direction])
        gt1, gt2 = tg(temp, tang[0]), tg(temp, tang[1])
        gpd = (phi[hi] - phi[lo]) * invd
        gp1, gp2 = tg(phi, tang[0]), tg(phi, tang[1])
        gn = np.sqrt(gpd**2 + gp1**2 + gp2**2)
        ok = gn >= MIN_GRADPHI
        safe = np.where(ok, gn, 1.0)
        nd = np.where(ok, gpd / safe, 0.0)
        st = np.where(ok, gtd - nd * (gpd * gtd + gp1 * gt1 + gp2 * gt2) / safe, 0.0)
        c = kface * (1.0 / kr - 1.0 / kl) * (0.5 - w) * (1.0 - nd * nd)
        st = (st - c * gtd) / (1.0 - c)
        if mode == "mid":
            kloc = np.where(pl_ + pr_ < 0.0, plane.kappa, 1.0)
            F = F + np.where(cut, st * (kloc - kface), 0.0)
        else:
            karea = face_area_fraction(0.5 * (pl_ + pr_),
                                       np.where(ok, gp1 / safe, 0.0),
                                       np.where(ok, gp2 / safe, 0.0),
                                       cd[tang[0]], cd[tang[1]])
            karea = karea + (1.0 - karea) * plane.kappa
            F = F + st * (karea - kface)
    return F, cut


def cmd_field(a):
    """The r = infinity BVP: solved field vs the analytic linear one."""
    plane = Plane(a.theta, a.x0, a.y0, a.kappa, a.q_n, a.amp)
    with h5py.File(a.field, "r") as h5:
        geo = BlockGeometry(h5)
        th = h5[a.name][...]
        num, den, linf = 0.0, 0.0, 0.0
        off_n, off_d = 0.0, 0.0
        for bid in range(geo.n_blocks):
            x, y, z, dV = geo.mesh(bid)
            ref = plane.temperature(*np.broadcast_arrays(x, y, z))
            off_n += float(np.sum((th[bid] - ref) * dV))
            off_d += float(np.sum(dV))
        offset = off_n / off_d                       # pure Neumann: up to a constant
        for bid in range(geo.n_blocks):
            x, y, z, dV = geo.mesh(bid)
            ref = plane.temperature(*np.broadcast_arrays(x, y, z)) + offset
            err = th[bid] - ref
            num += float(np.sum(err * err * dV))
            den += float(np.sum(dV))
            linf = max(linf, float(np.abs(err).max()))
    l2 = math.sqrt(num / den)
    print(f"   theta = {a.theta:g}  kappa_s = {a.kappa:g}   L2 = {l2:.6e}   Linf = {linf:.6e}")
    if a.emit:
        with open(a.emit, "a") as fh:
            fh.write("%s %g %g %.8e %.8e\n" % (a.tag, a.theta, a.kappa, l2, linf))
    ok = math.isfinite(l2) and (a.tolerance is None or linf <= a.tolerance)
    if a.tolerance is not None:
        print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g})")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--theta", type=float, required=True)
        p.add_argument("--kappa", type=float, required=True)
        p.add_argument("--x0", type=float, default=0.5)
        p.add_argument("--y0", type=float, default=0.5)
        p.add_argument("--q-n", type=float, default=1.0, dest="q_n")
        p.add_argument("--amp", type=float, default=0.0)
        p.add_argument("--emit", default=None, help="append one row to this table")
        p.add_argument("--tolerance", type=float, default=None)

    p = sub.add_parser("flux")
    p.add_argument("case")
    p.add_argument("--margin", type=int, default=0,
                   help="drop faces within this many cells of the domain edge")
    common(p)
    p.set_defaults(func=cmd_flux)

    p = sub.add_parser("residual")
    p.add_argument("case")
    common(p)
    p.set_defaults(func=cmd_residual)

    p = sub.add_parser("field")
    p.add_argument("field")
    p.add_argument("--name", default="theta")
    p.add_argument("--tag", default="run")
    common(p)
    p.set_defaults(func=cmd_field)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
