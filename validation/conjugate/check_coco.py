#!/usr/bin/env python3
"""THE COCO CEILING: how much of the cut-cell truncation error is captured by
a LOCAL analytic model, flat versus locally curved?

    ./check_coco.py --kappa 10 --grids 64 128 256 --scheme base

A COCO/deferred correction subtracts -L_h(T_e), with T_e a rescaled local
analytic solution of the dominant local operator. Since L(T_e) = 0 exactly,
L_h(T_e) IS a truncation error, and the correction removes as much of the
scheme's own truncation error as T_e resembles the true local behaviour. This
measures that resemblance directly -- the CEILING of the method, before any
question of how the amplitudes would be estimated from a discrete field.

THE TWO LOCAL MODELS, both exact two-material conjugate solutions anchored at
the cell's own nearest interface point (a, th0), inside = solid:

  flat    the 180/180-degree wedge pair, i.e. what the SHIPPED scheme's
          series resistance already is (design note Sec. 5: "the COCO
          construction with 180/180-degree wedges IS the series-resistance
          scheme"). Piecewise linear, kink on the TANGENT PLANE:
              T = (q_n/kappa) xi + a_t tau

  curved  the same two modes on a circular arc of the LOCAL curvature:
              T = (q_n a/kappa) ln(r/a)              radial mode
                + A [r or gamma r or (r + beta a^2/r)] sin(th - th0)
          with A = a_t/(1 + beta) = a_t (1 + kappa_s)/2 -- note that the
          2/(1 + kappa_s) scaling which defeated every DISCRETE estimate of
          s_t (README stage 1b) is here an ANALYTIC property of the basis.
          `flat` is its a -> infinity limit (verified: the difference falls
          exactly as 1/a).

Both amplitudes are taken from the EXACT solution. That is deliberate: this
is a ceiling, not a scheme. If the ceiling is low the idea is dead whatever
the estimator.

L_h is the scheme's own operator, and it is available EXACTLY as the matrix
check_convergence.py assembles by colouring, so L_h(T_e) at cell C is one
sparse row dotted with T_e on that row's stencil.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from check_convergence import build
from check_cylinder import Multipole


def local_models(xs, ys, th0, a, ks, qn, at, cx=0.5, cy=0.5):
    """(flat, curved) evaluated at the stencil points of one cell."""
    b1 = (1.0 - ks) / (1.0 + ks)
    g1 = 2.0 / (1.0 + ks)
    nx, ny = np.cos(th0), np.sin(th0)
    gx, gy = cx + a * nx, cy + a * ny
    xi = (xs - gx) * nx + (ys - gy) * ny
    ta = -(xs - gx) * ny + (ys - gy) * nx
    flat = (qn / np.where(xi <= 0.0, ks, 1.0)) * xi + at * ta

    X, Y = xs - cx, ys - cy
    r = np.maximum(np.hypot(X, Y), 1.0e-300)
    th = np.arctan2(Y, X)
    kap = np.where(r <= a, ks, 1.0)
    A = at / (1.0 + b1)
    curved = ((qn * a / kap) * np.log(r / a)
              + np.where(r <= a, A * g1 * r, A * (r + b1 * a * a / r))
              * np.sin(th - th0))
    return flat, curved


def run(case, model, scheme, a_rad, ks):
    A, b, ex, phi, U, cd, _ = build(case, model, scheme)
    nj, ni = ex.shape
    Ac = A.tocsr()
    resid = Ac @ ex.ravel() + b.ravel()

    # cell-centre coordinates of the unknown grid
    import h5py  # noqa: F401  (load_phi already did the io; reuse its centres)
    from check_oblique import load_phi
    _, centres = load_phi(case)
    xc = centres[0][2:2 + ni]
    yc = centres[1][2:2 + nj]
    XX, YY = np.meshgrid(xc, yc)

    p = phi[2][U[1], U[2]]
    band = np.zeros_like(p, dtype=bool)
    for ax in (0, 1):
        for sh in (-1, 1):
            band |= (np.roll(p, sh, axis=ax) < 0.0) != (p < 0.0)
    rows = np.flatnonzero(band.ravel())

    m = model.m
    g, beta = model.g, model.beta
    e_flat = np.zeros(rows.size)
    e_curv = np.zeros(rows.size)
    ref = resid[rows]
    for t, C in enumerate(rows):
        lo, hi = Ac.indptr[C], Ac.indptr[C + 1]
        cols, w = Ac.indices[lo:hi], Ac.data[lo:hi]
        xs, ys = XX.ravel()[cols], YY.ravel()[cols]
        th0 = np.arctan2(YY.ravel()[C] - 0.5, XX.ravel()[C] - 0.5)
        # the exact two mode amplitudes at that interface point
        qn = g * m * a_rad ** (m - 1) * (1.0 - beta) * np.cos(m * th0)
        at = -g * m * a_rad ** (m - 1) * (1.0 + beta) * np.sin(m * th0)
        fl, cu = local_models(xs, ys, th0, a_rad, ks, qn, at)
        e_flat[t] = w @ fl
        e_curv[t] = w @ cu
    rms = lambda v: float(np.sqrt(np.mean(v * v)))
    return dict(n=rows.size, ref=rms(ref),
                flat=rms(ref - e_flat), curved=rms(ref - e_curv),
                cf=1.0 - rms(ref - e_flat) / rms(ref),
                cc=1.0 - rms(ref - e_curv) / rms(ref))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="cyl")
    ap.add_argument("--grids", type=int, nargs="+", default=[64, 128, 256])
    ap.add_argument("--kappa", type=float, default=10.0)
    ap.add_argument("--radius", type=float, default=0.25)
    ap.add_argument("--mode", type=int, default=2)
    ap.add_argument("--scheme", default="base")
    a = ap.parse_args()
    model = Multipole(0.5, 0.5, a.radius, a.kappa, a.mode)
    print(f"== COCO ceiling: scheme '{a.scheme}', kappa_s = {a.kappa:g},"
          f" harmonic mode {a.mode}")
    print("      n   cells   |L_h(T_exact)|   residual after FLAT   after CURVED"
          "      captured")
    out = {}
    for n in a.grids:
        r = run(f"{a.tag}_{n}.h5", model, a.scheme, a.radius, a.kappa)
        out[n] = r
        print(f"   {n:4d}   {r['n']:5d}   {r['ref']:14.4e}   {r['flat']:19.4e}"
              f"   {r['curved']:12.4e}   flat {r['cf']*100:5.1f}%"
              f"  curved {r['cc']*100:5.1f}%")
    if len(a.grids) >= 2:
        h = np.array([1.0 / n for n in a.grids])
        for k, lab in (("ref", "uncorrected"), ("flat", "after FLAT  "),
                       ("curved", "after CURVED")):
            v = np.array([out[n][k] for n in a.grids])
            p = np.polyfit(np.log(h), np.log(np.maximum(v, 1e-300)), 1)[0]
            print(f"   {lab} residual: observed order {p:5.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
