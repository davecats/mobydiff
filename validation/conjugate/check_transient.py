#!/usr/bin/env python3
"""The TRANSIENT two-material slab (increment C3, gate 1).

C1 gated the steady state, where the capacity is irrelevant by construction.
This is the gate the fluid-fraction-weighted capacity exists for, and it is
built on the one transient of a two-material slab that is available in closed
form: a single DECAYING EIGENMODE.

    T(y, t) = exp(-mu t) X(y),
    X = A sin(k_s y)            (solid, y < y_w,  k_s = sqrt(mu C_s/kappa_s))
      = B sin(k_f (L - y))      (fluid, y > y_w,  k_f = sqrt(mu))

with homogeneous Dirichlet ends. [T] = 0 and [k dT/dy] = 0 at y_w give the
eigenvalue condition

    kappa_s k_s cot(k_s y_w) + k_f cot(k_f (L - y_w)) = 0 ,

whose lowest root is found here by bracketing and Brent (float64 to ~1e-15,
which is five orders below anything measured against it). Everything else is
then closed form -- the field at any time, and the interface flux

    q(t) = -kappa_s X_s'(y_w) exp(-mu t) ,

which is what the solver's interface-heat diagnostic must reproduce.

WHY A SINGLE MODE. The decay rate mu is a global functional of the capacity
distribution, so a cut cell carrying ONE material's whole heat capacity -- the
C1 pointwise value -- is an O(1) error on an O(h) band and shows up as a
FIRST-order error in mu, in the field and in the flux. The fraction-weighted
capacity removes it and restores second order. Both curves are measured
(run_gates_c3.sh drives the pointwise one with the archived C2 binary), which
is the whole content of the gate: not "the new capacity is defensible" but
"the old one costs an order, here it is".

    ./check_transient.py mode  --y-wall Y --kappa K --cap C [--ly L]
    ./check_transient.py seed  base.h5 out.h5 --y-wall Y --kappa K --cap C
    ./check_transient.py check field.h5 --y-wall Y --kappa K --cap C
                               [--heat heat.txt] [--emit table.dat] [--tag T]
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import sys

import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scalar"))
from scalar_tools import BlockGeometry                       # noqa: E402


class Mode:
    """The lowest eigenmode of the two-material slab, normalised to A = 1."""

    def __init__(self, y_wall, kappa, cap, ly=1.0):
        self.y_wall, self.kappa, self.cap, self.ly = y_wall, kappa, cap, ly
        self.alpha_s = kappa / cap
        self.mu = self._lowest_root()
        self.ks = math.sqrt(self.mu / self.alpha_s)
        self.kf = math.sqrt(self.mu)
        self.A = 1.0
        self.B = math.sin(self.ks * y_wall) / math.sin(self.kf * (ly - y_wall))

    def residual(self, mu):
        ks = math.sqrt(mu / self.alpha_s)
        kf = math.sqrt(mu)
        a = ks * self.y_wall
        b = kf * (self.ly - self.y_wall)
        return self.kappa * ks / math.tan(a) + kf / math.tan(b)

    def _lowest_root(self):
        """The first root of F. F falls from +infinity (mu -> 0) to
        -infinity at the first pole of either cotangent, so bracketing up to
        the first pole and bisecting is unconditionally safe -- no continuation
        in kappa_s, no initial guess to tune."""
        pole = min(math.pi / self.y_wall, math.pi / (self.ly - self.y_wall))
        # the first pole in mu: k_s y_w = pi  or  k_f (L - y_w) = pi
        mu_pole = min(self.alpha_s * (math.pi / self.y_wall) ** 2,
                      (math.pi / (self.ly - self.y_wall)) ** 2)
        lo, hi = 1e-12 * mu_pole, mu_pole * (1.0 - 1e-13)
        flo, fhi = self.residual(lo), self.residual(hi)
        if not (flo > 0.0 > fhi):
            raise SystemExit(f"eigenvalue not bracketed: F({lo:g}) = {flo:g}, "
                             f"F({hi:g}) = {fhi:g}, pole {pole:g}")
        for _ in range(200):                       # bisection: 200 halvings
            mid = 0.5 * (lo + hi)                  # is far past float64
            if self.residual(mid) > 0.0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def field(self, y, t=0.0):
        s = self.A * np.sin(self.ks * np.asarray(y))
        f = self.B * np.sin(self.kf * (self.ly - np.asarray(y)))
        return math.exp(-self.mu * t) * np.where(np.asarray(y) <= self.y_wall, s, f)

    def flux(self, t=0.0):
        """The PHYSICAL heat flux -k dT/dy at y_w, in the +y sense, i.e. from
        the solid into the fluid -- the sign convention of the solver's
        interface-heat diagnostic (positive = into the fluid). Continuous
        across the interface, which is the whole point of the problem."""
        return -(self.kappa * self.A * self.ks * math.cos(self.ks * self.y_wall)
                 * math.exp(-self.mu * t))


def describe(m):
    print(f"   mu = {m.mu:.16e}   k_s = {m.ks:.12f}   k_f = {m.kf:.12f}"
          f"   |F| = {abs(m.residual(m.mu)):.3e}")
    print(f"   A = {m.A:g}   B = {m.B:.12f}   q(0) = {m.flux(0.0):.16e}")


def cmd_mode(a):
    describe(Mode(a.y_wall, a.kappa, a.cap, a.ly))
    return 0


def cmd_seed(a):
    shutil.copyfile(a.base, a.out)
    with h5py.File(a.out, "r+") as h5:
        geo = BlockGeometry(h5)
        m = Mode(a.y_wall, a.kappa, a.cap, geo.leng[1])
        th = h5[a.name]
        for bid in range(geo.n_blocks):
            _, y, _, _ = geo.mesh(bid)
            th[bid] = m.field(np.broadcast_to(y, th[bid].shape))
    describe(m)
    print(f"{a.out}: eigenmode seeded")
    return 0


def read_heat(path, column=2):
    """The interface-heat file: step, time, then (staircase, graded, total)
    per scalar. Column 2 (1-based within the scalar's triple) is the total."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            v = [float(x) for x in line.split()]
            if len(v) >= 5:
                rows.append((v[1], v[2 + column]))
    return rows


def cmd_check(a):
    with h5py.File(a.field, "r") as h5:
        geo = BlockGeometry(h5)
        m = Mode(a.y_wall, a.kappa, a.cap, geo.leng[1])
        t = float(h5.attrs["t_current"])
        th = h5[a.name][...]
        num = den = 0.0
        linf = 0.0
        for bid in range(geo.n_blocks):
            _, y, _, dV = geo.mesh(bid)
            ref = m.field(np.broadcast_to(y, th[bid].shape), t)
            err = th[bid] - ref
            w = np.broadcast_to(dV, th[bid].shape)
            num += float(np.sum(err * err * w))
            den += float(np.sum(w))
            linf = max(linf, float(np.abs(err).max()))
        area = geo.leng[0] * geo.leng[2]
    l2 = math.sqrt(num / den)
    amp = abs(m.field(np.array([m.y_wall]), t)[0])
    print(f"   t = {t:.6f}   mu = {m.mu:.12f}   decay = {math.exp(-m.mu*t):.6e}")
    print(f"   field:  L2 = {l2:.6e}   Linf = {linf:.6e}   "
          f"relative to the surviving amplitude {linf/amp:.6e}")

    # The interface flux, from the SOLVER's own diagnostic. The heat file
    # reports the heat INTO THE FLUID over the whole interface, so dividing
    # by the interface area gives q -- with the sign of the +y flux, since
    # the solid is below.
    qrel = float("nan")
    murel = float("nan")
    if a.heat:
        rows = read_heat(a.heat)
        if not rows:
            raise SystemExit(f"{a.heat}: no data rows")
        th_t, heat = rows[-1]
        q_num = heat / area
        q_ref = m.flux(th_t)
        qrel = abs(q_num - q_ref) / abs(q_ref)
        print(f"   flux:   solver {q_num:.12e}   exact {q_ref:.12e}   "
              f"relative {qrel:.6e}   (t = {th_t:.6f})")

        # ...and the DECAY RATE the run actually realises, which is the
        # functional the capacity distribution controls: two flux samples an
        # interval apart give mu_num with no reference to the amplitude.
        if len(rows) >= 2:
            (t0, h0), (t1, h1) = rows[-2], rows[-1]
            if h0 * h1 > 0.0 and t1 > t0:
                mu_num = math.log(abs(h0 / h1)) / (t1 - t0)
                murel = abs(mu_num - m.mu) / m.mu
                print(f"   decay rate: solver {mu_num:.12f}   exact {m.mu:.12f}"
                      f"   relative {murel:.6e}")

    if a.emit:
        with open(a.emit, "a") as fh:
            fh.write("%s %.10e %.10e %.10e %.10e\n" % (a.tag, l2, linf, qrel, murel))
    ok = math.isfinite(l2) and (a.tolerance is None or linf <= a.tolerance)
    if a.tolerance is not None:
        print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g})")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--y-wall", type=float, required=True)
        p.add_argument("--kappa", type=float, required=True)
        p.add_argument("--cap", type=float, required=True)
        p.add_argument("--name", default="theta")

    p = sub.add_parser("mode")
    common(p)
    p.add_argument("--ly", type=float, default=1.0)
    p.set_defaults(func=cmd_mode)

    p = sub.add_parser("seed")
    p.add_argument("base")
    p.add_argument("out")
    common(p)
    p.set_defaults(func=cmd_seed)

    p = sub.add_parser("check")
    p.add_argument("field")
    common(p)
    p.add_argument("--heat", default=None)
    p.add_argument("--emit", default=None)
    p.add_argument("--tag", default="run")
    p.add_argument("--tolerance", type=float, default=None)
    p.set_defaults(func=cmd_check)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
