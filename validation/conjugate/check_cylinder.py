#!/usr/bin/env python3
"""C2 gate 2: the cylindrical shell, and the O(kappa h/a) error in w.

    ./check_cylinder.py flux <case.h5> --radius A --kappa K --source S

THE MANUFACTURED SOLUTION. A solid cylinder of radius a carrying a uniform
volumetric source S, in a conducting fluid, has the exact steady solution

    T = T0 - S r^2/(4 kappa_s)                        r <= a   (solid)
    T = T0 - S a^2/(4 kappa_s) - (S a^2/2) ln(r/a)    r >= a   (fluid, log)

with the interface-normal flux q_n = -S a/2, the SAME at every point of the
interface. Both material laws are satisfied exactly and the flux balance
holds by construction, so this is a genuine conjugate solution -- for
D = 1/(Re Pr) = 1 and C_s = 1, which the ini sets.

WHY THIS GATE EXISTS. The field is RADIAL, so grad T is purely normal and
s_t = 0 identically: the C2 tangential term is inert here by construction.
What is left in the cut-face flux is the one approximation the level-set
weight makes -- w is exact for a PLANE, and for a curved interface phi_L and
phi_R measure to DIFFERENT nearest points, so w carries a relative error
O(kappa h/a) (LaTeX note Section 6.6). This gate isolates exactly that, and
it is also the check that the correction does NOT invent a flux where there
is no tangential gradient.

The measurement is per cut face on the ANALYTIC field with the case file's
own phi -- the machinery is check_oblique's, so both gates measure the same
way. No boundary condition is involved (and none would do: the log solution
is radial, a box boundary is not).
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np

from check_oblique import (load_phi, face_measure, debias,       # noqa: E402
                           face_flux_field, MIN_COSINE, MIN_GRADPHI)


class Shell:
    """The radial conjugate solution above; the interface of a 'field object'
    that check_oblique's face machinery consumes."""

    def __init__(self, cx, cy, radius, kappa, source, t0=0.0):
        self.c = (cx, cy)
        self.a = radius
        self.kappa = kappa
        self.source = source
        self.t0 = t0

    def radius(self, x, y):
        return np.sqrt((x - self.c[0])**2 + (y - self.c[1])**2)

    def temperature(self, x, y, z):
        r = np.maximum(self.radius(x, y), 1.0e-300)
        inner = self.t0 - self.source * r**2 / (4.0 * self.kappa)
        outer = (self.t0 - self.source * self.a**2 / (4.0 * self.kappa)
                 - 0.5 * self.source * self.a**2 * np.log(r / self.a))
        return np.where(r <= self.a, inner, outer)

    def s_t(self, d):
        return 0.0                       # radial: no tangential gradient at all

    @property
    def q_n(self):
        return -0.5 * self.source * self.a


class Dipole:
    """A conducting cylinder in a UNIFORM far-field gradient -- the curved
    counterpart of the oblique plane, and the right test for a face-flux
    scheme on a curved interface.

        T_out = -G cos(t) (r + beta a^2/r),   beta = (1 - kappa)/(1 + kappa)
        T_in  = -G cos(t) gamma r,            gamma = 2/(1 + kappa)

    Both pieces are HARMONIC (r^-1 cos t is harmonic in 2D), [T] = 0 and
    [k d_n T] = 0 hold at r = a, and there is NO source -- so the exact
    divergence is zero everywhere and the truncation-residual test applies
    directly. Unlike the radial log solution of `Shell`, this one carries a
    TANGENTIAL gradient, and its local ratio |grad_t T|/|d_n T| = |tan t|
    sweeps the whole range 0 -> infinity around the body: one run samples
    every ratio the decision turns on.

    The interface-normal flux is q_r(t) = -G cos(t) 2 kappa/(1 + kappa), and
    the nearest interface point of a face centre is its RADIAL projection, so
    the face centre's polar angle is the exact angle to compare at -- no O(h)
    error enters the reference.
    """

    def __init__(self, cx, cy, radius, kappa, grad=1.0):
        self.c = (cx, cy)
        self.a = radius
        self.kappa = kappa
        self.g = grad
        self.beta = (1.0 - kappa)/(1.0 + kappa)
        self.gamma = 2.0/(1.0 + kappa)

    def _xy(self, x, y):
        X = x - self.c[0]
        Y = y - self.c[1]
        return X, Y, np.maximum(X*X + Y*Y, 1.0e-300)

    def temperature(self, x, y, z):
        X, Y, r2 = self._xy(x, y)
        inner = -self.g*self.gamma*X
        outer = -self.g*(X + self.beta*self.a**2*X/r2)
        return np.where(r2 <= self.a**2, inner, outer)

    def theta(self, x, y):
        X, Y, _ = self._xy(x, y)
        return np.arctan2(Y, X)

    def q_r(self, theta):
        """k d_r T at r = a, outward -- continuous across the interface."""
        return -self.g*np.cos(theta)*2.0*self.kappa/(1.0 + self.kappa)

    def ratio(self, theta):
        """|grad_t T|/|d_n T| measured on the FLUID side -- the same quantity
        `r` that Plane controls directly, so that the plane and the body are
        binned on ONE definition.

        FIXED 2026-09-14: this returned |tan t|, which is the kappa_s = 1
        form. At r = a the tangential gradient is (1 + beta) and the fluid
        normal gradient (1 - beta), and (1+beta)/(1-beta) = 1/kappa_s, so

            |grad_t T|/|d_n T| = |tan t|/kappa_s.

        The old form overstated the ratio by the full contrast, i.e. by a
        DECADE at kappa_s = 10 -- which is what the README's position-resolved
        crossover was binned on. See the README's C2 correction note."""
        return np.abs(np.tan(theta))/self.kappa

    def grad_outer(self, theta):
        """grad T at r = a, from the fluid side (its TANGENTIAL part is the
        one s_t is after, and that part is continuous)."""
        gx = -self.g*(1.0 - self.beta*np.cos(2.0*theta))
        gy = self.g*self.beta*np.sin(2.0*theta)
        return gx, gy, np.zeros_like(gx)

    def s_t_exact(self, theta, d):
        """The exact tangential scalar of a d-direction arm at angle t."""
        g = self.grad_outer(theta)
        n = (np.cos(theta), np.sin(theta), np.zeros_like(theta))
        dn = sum(n[q]*g[q] for q in range(3))
        return g[d] - n[d]*dn

    def s_t(self, d):
        """face_measure's constant-s_t hook. A curved interface has no single
        value, and this measurement does not use it -- NaN so that any
        accidental consumer shows up instead of silently reading a wrong 0."""
        return float("nan")


class Multipole:
    """The cylinder in a mode-m harmonic far field -- the Dipole generalised,
    and the reason it exists is that the DIPOLE'S INTERIOR IS LINEAR.

        T_out = (r^m + beta a^2m r^-m) cos(m t),   beta = (1 - k)/(1 + k)
        T_in  = gamma r^m cos(m t),                gamma = 2/(1 + k)

    Both harmonic, [T] = 0 and [k d_r T] = 0 at r = a, no source -- so the
    exact divergence is zero and every measurement the Dipole supports
    carries over.

    WHY IT MATTERS FOR THE ONE-SIDED ESTIMATOR. At m = 1 the interior is
    gamma r cos t = gamma x, exactly LINEAR, so a same-side stencil inside
    the body is exact by construction and any solid-side score is flattered
    by the test problem rather than earned. At m = 2 the interior is
    gamma (x^2 - y^2) and NEITHER side is linear, so the two sides are
    scored on equal terms. Measure on m = 2 before believing a side.
    """

    def __init__(self, cx, cy, radius, kappa, mode=2, grad=1.0):
        self.c = (cx, cy)
        self.a = radius
        self.kappa = kappa
        self.m = int(mode)
        self.g = grad
        self.beta = (1.0 - kappa)/(1.0 + kappa)
        self.gamma = 2.0/(1.0 + kappa)

    def _polar(self, x, y):
        X, Y = x - self.c[0], y - self.c[1]
        return np.sqrt(np.maximum(X*X + Y*Y, 1.0e-300)), np.arctan2(Y, X)

    def temperature(self, x, y, z):
        r, t = self._polar(x, y)
        m, a = self.m, self.a
        rs = np.maximum(r, 1.0e-300)
        outer = (rs**m + self.beta*a**(2*m)*rs**(-m))*np.cos(m*t)
        inner = self.gamma*rs**m*np.cos(m*t)
        return self.g*np.where(r <= a, inner, outer)

    def theta(self, x, y):
        return self._polar(x, y)[1]

    def q_r(self, theta):
        """k d_r T at r = a, outward -- continuous across the interface."""
        m, a = self.m, self.a
        return self.g*m*a**(m - 1)*(1.0 - self.beta)*np.cos(m*theta)

    def ratio(self, theta):
        """|grad_t T|/|d_n T| on the FLUID side = |tan(m t)|/kappa_s."""
        return np.abs(np.tan(self.m*theta))/self.kappa

    def grad_outer(self, theta):
        """grad T at r = a from the fluid side."""
        m, a = self.m, self.a
        dr = self.g*m*a**(m - 1)*(1.0 - self.beta)*np.cos(m*theta)
        dt = -self.g*m*a**(m - 1)*(1.0 + self.beta)*np.sin(m*theta)
        c, s = np.cos(theta), np.sin(theta)
        return dr*c - dt*s, dr*s + dt*c, np.zeros_like(dr)

    def s_t_exact(self, theta, d):
        g = self.grad_outer(theta)
        n = (np.cos(theta), np.sin(theta), np.zeros_like(theta))
        dn = sum(n[q]*g[q] for q in range(3))
        return g[d] - n[d]*dn

    def s_t(self, d):
        return float("nan")


def cmd_dipole(a):
    """Cylinder in a uniform gradient: the two measurements that decide.

    (1) the local interface-flux error, BINNED BY THE LOCAL RATIO -- baseline
        against corrected, i.e. where on a real body the correction pays;
    (2) the cut-cell truncation residual for the three s_t multipliers, on a
        CURVED interface (the exact divergence is zero, so anything left is
        the scheme's; away from the body the background is the ordinary
        O(h^2) Laplacian truncation, which is what the cut cells must be
        compared against).
    """
    dip = Dipole(a.centre[0], a.centre[1], a.radius, a.kappa, a.grad)
    phi, centres = load_phi(a.case)
    shape = phi.shape
    gk, gj, gi = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    xx, yy, zz = centres[0][gi], centres[1][gj], centres[2][gk]
    temp = dip.temperature(xx, yy, zz)
    cd = [float(np.diff(c)[0]) for c in centres]

    # ---- (1) the interface flux, binned by the local ratio ----------------
    # The implied normal flux is (dT - h s_t)/(a R) whatever multiplier the
    # FLUX uses, so this half compares exactly two schemes: s_t = 0 and the
    # shipped s_t.
    rows = []
    for d in range(3):
        m = face_measure(phi, centres, dip, d)
        if m is None:
            continue
        axis = 2 - d
        sl = [slice(1, shape[q] - 1) for q in range(3)]
        sl[axis] = slice(1, shape[axis] - 1)
        # face centre = midpoint of the two cell centres
        lo = [slice(1, shape[q] - 1) for q in range(3)]
        lo[axis] = slice(0, shape[axis] - 2)
        xf = 0.5*(xx[tuple(sl)] + xx[tuple(lo)])
        yf = 0.5*(yy[tuple(sl)] + yy[tuple(lo)])
        cut = (phi[tuple(lo)] < 0.0) != (phi[tuple(sl)] < 0.0)
        th = dip.theta(xf[cut], yf[cut])
        qref = -dip.q_r(th)                      # the note's orientation
        base = m["dt"]/(m["a"]*m["res"])/qref - 1.0
        corr = (m["dt"] - m["hd"]*debias(m))/(m["a"]*m["res"])/qref - 1.0
        ste = dip.s_t_exact(th, d)
        rows.append((dip.ratio(th), base, corr, m["st_num"], debias(m), ste))
    rat = np.concatenate([r[0] for r in rows])
    base = np.concatenate([r[1] for r in rows])
    corr = np.concatenate([r[2] for r in rows])
    st_raw = np.concatenate([r[3] for r in rows])
    st_shp = np.concatenate([r[4] for r in rows])
    st_ex = np.concatenate([r[5] for r in rows])

    print(f"   radius = {a.radius:g}  kappa_s = {a.kappa:g}  h = {cd[0]:.6g}"
          f"   cut faces = {rat.size}")
    print("   interface-flux error by LOCAL RATIO |grad_t T|/|d_n T| = |tan t|:")
    print("     ratio band      faces   baseline rms   corrected rms   winner")
    edges = [(0.0, 0.1), (0.1, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, np.inf)]
    for l, h in edges:
        sel = (rat >= l) & (rat < h)
        if not sel.any():
            continue
        b = float(np.sqrt(np.mean(base[sel]**2)))
        c = float(np.sqrt(np.mean(corr[sel]**2)))
        print(f"     {l:5.2f} - {h:5.2f}   {int(sel.sum()):5d}   {b:12.4e}   {c:12.4e}"
              f"   {'corrected' if c < b else 'BASELINE'}")

    # ---- (1b) how good s_t itself is, on a CURVED interface ---------------
    # This is the quantity that sets the floor of any multiplier: the exact
    # s_t is known here, unlike on the shell where it is identically zero.
    scale = float(np.sqrt(np.mean(st_ex**2)))
    print(f"   s_t itself (rms over cut faces, exact rms = {scale:.4e}):")
    for lab, v in (("raw projection ", st_raw), ("shipped de-bias", st_shp)):
        err = float(np.sqrt(np.mean((v - st_ex)**2)))
        print(f"     {lab}   rms error = {err:.4e}   = {err/scale*100:6.2f} % of the signal")

    # ---- (2) the truncation residual, three multipliers -------------------
    print("   cut-cell truncation residual (exact divergence is zero):")
    for mode, label in (("base", "C1 baseline   k_face"),
                        ("mid",  "C2 shipped    k_loc "),
                        ("area", "PROPOSED      k_area")):
        total = np.zeros(shape)
        band = np.zeros(shape, dtype=bool)
        defined = np.ones(shape, dtype=bool)
        for d in range(3):
            axis = 2 - d
            F, cut = face_flux_field(phi, temp, centres, dip, d, mode)
            full = np.full(shape, np.nan)
            fcut = np.zeros(shape, dtype=bool)
            idx = tuple(slice(1, shape[q] - 1) for q in range(3))
            full[idx] = F
            fcut[idx] = cut
            hi = np.roll(full, -1, axis=axis)
            total = total + (hi - full)/cd[d]
            band = band | fcut | np.roll(fcut, -1, axis=axis)
            defined = defined & np.isfinite(full) & np.isfinite(hi)
        core = tuple(slice(2, shape[q] - 2) for q in range(3))
        r, b, ok = total[core], band[core], defined[core]
        cm, un = b & ok, (~b) & ok
        print(f"     {label}   max = {np.abs(r[cm]).max():.4e}   rms = {r[cm].std():.4e}"
              f"   (uncut background {np.abs(r[un]).max():.2e})")
        if a.emit:
            with open(a.emit, "a") as fh:
                fh.write("%g %g %.8e %s %.8e %.8e %.8e\n"
                         % (a.radius, a.kappa, cd[0], mode,
                            np.abs(r[cm]).max(), r[cm].std(), np.abs(r[un]).max()))
    return 0


def cmd_flux(a):
    shell = Shell(a.centre[0], a.centre[1], a.radius, a.kappa, a.source)
    phi, centres = load_phi(a.case)
    parts = [m for m in (face_measure(phi, centres, shell, d, a.margin)
                         for d in range(3)) if m is not None]
    if not parts:
        raise SystemExit("no cut faces found")

    # The note's q_n runs from the L cell of the arm to its R cell; the
    # shell's is along +r, i.e. into the fluid. The two differ by a sign,
    # face-independently, because `a` carries the orientation (the same
    # convention check_oblique makes).
    qn_ref = -shell.q_n
    qb, qc, qr, st, nfaces = [], [], [], [], 0
    for m in parts:
        nfaces += m["n"]
        st_db = debias(m)
        qb.append(m["dt"] / (m["a"] * m["res"]) / qn_ref - 1.0)
        qc.append((m["dt"] - m["hd"] * st_db) / (m["a"] * m["res"]) / qn_ref - 1.0)
        qr.append((m["dt"] - m["hd"] * m["st_num"]) / (m["a"] * m["res"]) / qn_ref - 1.0)
        st.append(m["hd"] * st_db / m["dt"])             # e_face, the shipped s_t

    def stat(chunks):
        v = np.concatenate([np.atleast_1d(c) for c in chunks])
        return float(np.abs(v).max()), float(np.sqrt(np.mean(v * v)))

    qb_max, qb_rms = stat(qb)
    qc_max, qc_rms = stat(qc)
    qr_max, qr_rms = stat(qr)
    ef_max, ef_rms = stat(st)
    h = parts[0]["hd"]
    print(f"   radius = {a.radius:g}  kappa_s = {a.kappa:g}  h = {h:.6g}"
          f"  cut faces = {nfaces}   kappa_curv h = {h / a.radius:.4f}")
    print(f"   q_n rel err  BASELINE     max = {qb_max:.4e}   rms = {qb_rms:.4e}")
    print(f"   q_n rel err  RAW s_t      max = {qr_max:.4e}   rms = {qr_rms:.4e}"
          "   (the note's construction 1, NOT shipped)")
    print(f"   q_n rel err  CORRECTED    max = {qc_max:.4e}   rms = {qc_rms:.4e}")
    print(f"   e_face (discrete, exact 0) max = {ef_max:.4e}   rms = {ef_rms:.4e}")
    if a.emit:
        with open(a.emit, "a") as fh:
            fh.write("%g %g %.8e %d %.8e %.8e %.8e %.8e %.8e\n"
                     % (a.radius, a.kappa, h, nfaces, qb_max, qb_rms,
                        qc_max, qc_rms, ef_max))
    ok = math.isfinite(qb_max) and (a.tolerance is None or qb_rms <= a.tolerance)
    if a.tolerance is not None:
        print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g})")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("flux")
    p.add_argument("case")
    p.add_argument("--radius", type=float, required=True)
    p.add_argument("--kappa", type=float, required=True)
    p.add_argument("--source", type=float, default=1.0)
    p.add_argument("--centre", type=float, nargs=2, default=(0.5, 0.5))
    p.add_argument("--margin", type=int, default=0)
    p.add_argument("--emit", default=None)
    p.add_argument("--tolerance", type=float, default=None)
    p.set_defaults(func=cmd_flux)

    p = sub.add_parser("dipole")
    p.add_argument("case")
    p.add_argument("--radius", type=float, required=True)
    p.add_argument("--kappa", type=float, required=True)
    p.add_argument("--grad", type=float, default=1.0)
    p.add_argument("--centre", type=float, nargs=2, default=(0.5, 0.5))
    p.add_argument("--emit", default=None)
    p.set_defaults(func=cmd_dipole)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
