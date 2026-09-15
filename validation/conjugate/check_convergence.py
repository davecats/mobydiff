#!/usr/bin/env python3
"""A PROPER convergence test on a CURVED interface: solve, refine, measure.

Everything in the C2/stage-0/1/1b campaign measured a local TRUNCATION
residual. That is not a convergence statement: the residual lives on a
codimension-one set, and a conservative scheme can damp it. This solves the
discrete conjugate BVP and measures the SOLUTION error under refinement.

    ./check_convergence.py --kappa 10 --mode 2 --grids 64 128 256

THE PROBLEM. The Multipole is harmonic in BOTH materials with no source, so
the exact solution satisfies div(k grad T) = 0 and the discrete problem is
L_h T = 0 with Dirichlet data taken from the exact solution on the two
outermost cell layers (the divergence stencil reaches one cell in, and the
interface is far from the boundary -- radius 0.25 in a unit box).

WHY IT CAN BE SOLVED EXACTLY RATHER THAN ITERATIVELY. Every scheme here is
LINEAR in T with purely geometric coefficients (k_face, the area fraction, the
closure's A and B, the extension's own averaging pattern all depend on phi
alone, and phi is static). So the operator has a matrix, and it is recovered
EXACTLY by graph colouring -- P^2 applications with P = 2R+1 -- and solved
directly. No Krylov tolerance enters the measured error.

The geometry is z-invariant and the z-direction flux of a z-invariant field is
identically zero (checked: 1.4e-13), so the unknowns are two-dimensional.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from check_oblique import load_phi, face_flux_field
from check_cylinder import Dipole, Multipole
from check_oblique import Plane

RADIUS = 6                      # colouring half-width; verified a posteriori


def build(case, model, mode):
    """The operator's matrix A and the boundary vector b, with L(T) = A u + b
    for T = (exact outside the unknowns) + (u inside)."""
    phi, centres = load_phi(case)
    shape = phi.shape
    cd = [float(np.diff(c)[0]) for c in centres]
    gk, gj, gi = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    exact = model.temperature(centres[0][gi], centres[1][gj], centres[2][gk])
    U = (slice(None), slice(2, shape[1] - 2), slice(2, shape[2] - 2))
    nj, ni = shape[1] - 4, shape[2] - 4

    def divergence(T):
        total = np.zeros(shape)
        idx = tuple(slice(1, shape[m] - 1) for m in range(3))
        for d in range(3):
            axis = 2 - d
            F, _ = face_flux_field(phi, T, centres, model, d, mode)
            full = np.zeros(shape)
            full[idx] = F
            total = total + (np.roll(full, -1, axis=axis) - full) / cd[d]
        return total[U][2]                      # any z-slice: the field is invariant

    def apply2d(u2d, with_exact):
        T = exact.copy() if with_exact else np.zeros(shape)
        T[U] = u2d[None, :, :]
        return divergence(T)

    b = apply2d(np.zeros((nj, ni)), True)       # exact data, zero unknowns
    zero = apply2d(np.zeros((nj, ni)), False)   # must be exactly 0 (linearity)
    assert np.abs(zero).max() == 0.0, "operator is not linear/homogeneous"

    # --- the matrix, by colouring -------------------------------------------
    P = 2 * RADIUS + 1
    jj, ii = np.meshgrid(np.arange(nj), np.arange(ni), indexing="ij")
    rows, cols, vals = [], [], []
    for a in range(P):
        for bcol in range(P):
            e = np.zeros((nj, ni))
            e[(jj % P == a) & (ii % P == bcol)] = 1.0
            y = apply2d(e, False)
            nz = np.abs(y) > 0.0
            if not nz.any():
                continue
            jm, im = jj[nz], ii[nz]
            # the unique same-colour cell within the radius, in closed form
            js = jm + ((a - jm + RADIUS) % P) - RADIUS
            iss = im + ((bcol - im + RADIUS) % P) - RADIUS
            good = (js >= 0) & (js < nj) & (iss >= 0) & (iss < ni)
            if not good.all():
                # a contribution from outside the same-colour window means the
                # stencil is wider than RADIUS -- the recovery would be wrong
                raise SystemExit(f"stencil exceeds RADIUS={RADIUS}; raise it")
            rows.append(jm * ni + im)
            cols.append(js * ni + iss)
            vals.append(y[nz])
    A = sp.csr_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(nj * ni, nj * ni))
    return A, b, exact[U][2], phi, U, cd, apply2d


def run(case, model, mode, verify):
    A, b, ex, phi, U, cd, apply2d = build(case, model, mode)
    nj, ni = ex.shape
    # THE MATRIX MUST BE THE OPERATOR, and this is checked rather than
    # argued: the colouring recovery picks the unique same-colour cell within
    # RADIUS, so a stencil WIDER than RADIUS silently attributes an entry to
    # the wrong column while staying in range. Only a direct comparison
    # catches that.
    rng = np.random.default_rng(0)
    v = rng.standard_normal((nj, ni))
    ref = apply2d(v, False)
    rel = np.abs(A @ v.ravel() - ref.ravel()).max() / max(np.abs(ref).max(), 1e-300)
    if rel > 1e-12:
        raise SystemExit(f"matrix != operator (rel {rel:.3e}) -- raise RADIUS")
    if verify:
        print(f"      matrix-vs-operator: {rel:.3e}")
    u = spla.spsolve(A.tocsc(), -b.ravel()).reshape(nj, ni)
    e = u - ex
    # the interface band: cells with a sign change among their 4 in-plane
    # neighbours, i.e. the cells the cut faces feed
    p = phi[2][U[1], U[2]]
    band = np.zeros_like(e, dtype=bool)
    for ax in (0, 1):
        for sh in (-1, 1):
            band |= (np.roll(p, sh, axis=ax) < 0.0) != (p < 0.0)
    out = {}
    for name, m in (("all", np.ones_like(band)), ("fluid", p > 0), ("solid", p < 0),
                    ("band", band), ("away", ~band)):
        out[name] = (float(np.sqrt(np.mean(e[m]**2))), float(np.abs(e[m]).max()))
    out["scale"] = float(np.sqrt(np.mean(ex**2)))
    out["nnz"] = A.nnz / (nj * ni)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="cyl")
    ap.add_argument("--grids", type=int, nargs="+", default=[64, 128, 256])
    ap.add_argument("--kappa", type=float, default=10.0)
    ap.add_argument("--radius", type=float, default=0.25)
    ap.add_argument("--mode", type=int, default=2)
    ap.add_argument("--geometry", choices=["cylinder", "plane"], default="cylinder")
    ap.add_argument("--theta", type=float, default=30.0)
    ap.add_argument("--x0", type=float, default=0.5)
    ap.add_argument("--q-n", type=float, default=1.0, dest="q_n")
    ap.add_argument("--amp", type=float, default=0.01)
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--summary", action="store_true",
                    help="one machine-readable line per scheme")
    ap.add_argument("--schemes", nargs="+",
                    default=["base", "area", "area1sext"])
    a = ap.parse_args()

    if a.geometry == "plane":
        # The manufactured oblique plane of check_oblique: linear in each
        # material, so harmonic in both and div(k grad T) = 0 exactly -- the
        # same BVP statement the cylinder satisfies. This is the geometry the
        # as-built note's "oblique face, r > 0" row is about.
        model = Plane(a.theta, a.x0, 0.5117, a.kappa, a.q_n, a.amp)
    else:
        model = (Dipole(0.5, 0.5, a.radius, a.kappa) if a.mode == 1 else
                 Multipole(0.5, 0.5, a.radius, a.kappa, a.mode))
    print(f"== SOLUTION convergence, curved interface: kappa_s = {a.kappa:g},"
          f" harmonic mode {a.mode}" if a.geometry == "cylinder" else
          f"== SOLUTION convergence, OBLIQUE PLANE: theta = {a.theta:g} deg,"
          f" kappa_s = {a.kappa:g}, r = {a.amp/a.q_n if a.q_n else float(chr(105)+chr(110)+chr(102)):g}")
    res = {s: {} for s in a.schemes}
    for n in a.grids:
        for s in a.schemes:
            res[s][n] = run(f"{a.tag}_{n}.h5", model, s, a.verify)
        sc = res[a.schemes[0]][n]["scale"]
        print(f"   n = {n:4d}   |T| = {sc:.4e}")
        for s in a.schemes:
            r = res[s][n]
            print(f"     {s:12s} L2 {r['all'][0]:.4e}  Linf {r['all'][1]:.4e}"
                  f"   band L2 {r['band'][0]:.4e}   away L2 {r['away'][0]:.4e}"
                  f"   (nnz/row {r['nnz']:.1f})")
    if len(a.grids) >= 2:
        h = np.array([1.0 / n for n in a.grids])
        print("   observed order (least squares):")
        for s in a.schemes:
            line = f"     {s:12s}"
            for key in ("all", "band", "away"):
                for q, lab in ((0, "L2"), (1, "Li")):
                    v = np.array([res[s][n][key][q] for n in a.grids])
                    p = np.polyfit(np.log(h), np.log(np.maximum(v, 1e-300)), 1)[0]
                    line += f"   {key}-{lab} {p:5.2f}"
            print(line)
        if a.summary:
            for s2 in a.schemes:
                v = np.array([res[s2][n]["all"][0] for n in a.grids])
                vi = np.array([res[s2][n]["all"][1] for n in a.grids])
                p = np.polyfit(np.log(h), np.log(v), 1)[0]
                pi = np.polyfit(np.log(h), np.log(vi), 1)[0]
                print(f"SUMMARY {a.kappa:g} {a.mode} {s2} {p:.4f} {pi:.4f} "
                      f"{v[-1]:.6e} {vi[-1]:.6e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
