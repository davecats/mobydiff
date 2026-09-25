#!/usr/bin/env python3
"""F1 gate: [scalar.N] source_dir, the Kasagi source on a chosen direction.

Two checks, because they pin different things:

  closed-form   Uniform (u, v, w) in a triply-periodic box is an exact steady
                solution and a uniform scalar has no gradient, so the Kasagi
                source integrates exactly to theta = source * u_dir * t. The
                three directions give three DIFFERENT numbers, so a branch
                reading the wrong component lands on another one's answer
                rather than merely on a wrong magnitude.
  identical     In a channel driven along the x-z diagonal, u(y) = w(y) cell
                by cell, so the x and the z source are the same number and
                the two runs must agree BIT FOR BIT. The premise is checked
                too (max|u - w| must be 0), or the test proves nothing.

    ./check_source_dir.py closed-form ku_x.h5 ku_y.h5 ku_z.h5 \
        --components 1.0 2.0 3.0 --time 0.04
    ./check_source_dir.py identical kc_x.h5 kc_z.h5
"""

from __future__ import annotations

import argparse
import h5py
import numpy as np


def cmd_closed_form(a):
    ok = True
    for path, comp in zip(a.files, a.components):
        with h5py.File(path, "r") as h5:
            th = h5[a.name][...]
        exact = a.source * comp * a.time
        err = float(np.abs(th - exact).max())
        # The tolerance is a few ulp of the answer itself: the arithmetic is
        # an exact accumulation of identical increments, not an approximation.
        tol = max(a.tolerance, 8.0 * np.spacing(abs(exact)))
        good = err <= tol
        ok &= good
        print(f"   {path}: theta = {th.min():.16f} .. {th.max():.16f}"
              f"   exact {exact:.16f}   max|err| {err:.3e}  {'ok' if good else 'FAIL'}")
    # The three answers must also be DISTINCT, or the test is blind.
    vals = [a.source * c * a.time for c in a.components]
    if len(set(vals)) != len(vals):
        print("   the reference values are not distinct -- the test is blind")
        ok = False
    print("   PASS" if ok else "   FAIL")
    return 0 if ok else 1


def cmd_identical(a):
    with h5py.File(a.files[0], "r") as A, h5py.File(a.files[1], "r") as B:
        names = [n for n in ("un", "vn", "wn", "pn", a.name) if n in A and n in B]
        worst = {n: float(np.abs(A[n][...] - B[n][...]).max()) for n in names}
        scale = {n: float(np.abs(A[n][...]).max()) for n in names}
        u, w = A["un"][...], A["wn"][...]
        premise = float(np.abs(u - w).max())
        umax = float(np.abs(u).max())
    for n in names:
        print(f"   {n:6s} max|diff| {worst[n]:.3e}   scale {scale[n]:.3e}")
    print(f"   premise: max|u - w| = {premise:.3e} on max|u| = {umax:.3e}"
          f"   ({'u == w' if premise == 0.0 else 'NOT EQUAL -- the test proves nothing'})")
    ok = premise == 0.0 and umax > 0.0 and all(v == 0.0 for v in worst.values())
    print("   PASS (bit-identical)" if ok else "   FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("closed-form")
    p.add_argument("files", nargs="+")
    p.add_argument("--components", type=float, nargs="+", required=True)
    p.add_argument("--time", type=float, required=True)
    p.add_argument("--source", type=float, default=1.0)
    p.add_argument("--name", default="theta")
    p.add_argument("--tolerance", type=float, default=0.0)
    p.set_defaults(func=cmd_closed_form)

    p = sub.add_parser("identical")
    p.add_argument("files", nargs=2)
    p.add_argument("--name", default="theta")
    p.set_defaults(func=cmd_identical)

    a = ap.parse_args()
    if a.cmd == "closed-form" and len(a.files) != len(a.components):
        raise SystemExit("closed-form: one --components value per file")
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
