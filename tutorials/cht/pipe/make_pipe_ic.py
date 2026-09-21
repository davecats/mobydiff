#!/usr/bin/env python3
"""Initial condition for the conjugate pipe: a turbulent-like mean profile
plus a divergence-free, large-scale perturbation.

WHY THIS EXISTS, because it cost a run to learn: the solver's own cold start
([flow] initial = uniform + initial_noise) puts a PLUG profile in the pipe,
and a plug has all of its shear inside the wall cell. White noise on top of it
has nothing to feed on -- transient growth needs the mean shear -- so it
dissipates, and pipe flow, being linearly stable, then relaxes to the laminar
state the forcing supports. MEASURED on the first attempt: after t = 25 the
fluctuations were down to 1e-3 and <w> on the axis was 1.47 and rising toward
the laminar 3.1 (u_b = 1.55), with u'/u_b ~ 1e-3. That is not slow transition,
it is relaminarisation.

WHAT THIS WRITES instead:

  mean      w = u_c (1 - r/R)^(1/n) inside the pipe, zero outside, with u_c
            set so that the bulk velocity is exactly --u-bulk. A 1/7 power law
            is already close to the turbulent profile at Re_b = 5300, so the
            run starts near the attractor rather than having to diffuse its
            way there over the viscous time R^2/nu = 1325 D/u_b.

  perturbation   u' = curl(g(r) A(x)), which is divergence-free EXACTLY (any
            curl is), so the first projection has almost nothing to remove and
            the perturbation survives the first step at the amplitude it was
            given. A is a sum of a few random Fourier modes with wavelengths
            of order the diameter -- LARGE scales, the ones that grow -- and
            g(r) = (1 - (r/R)^2)^2 confines the perturbation to the fluid and
            makes it vanish at the wall, with the product differentiated
            analytically (curl(gA) = g curl A + grad g x A), so the envelope
            does not spoil the divergence.

Each velocity component is evaluated at ITS OWN staggered location (q(i) is
the LOW face of cell i), the same convention the solver and pipe_stats.py use.

    ./make_pipe_ic.py p_coarse_dev_2500.h5 IC_pipe.h5 --amp 0.15
"""

from __future__ import annotations

import argparse
import shutil
import sys

import h5py
import numpy as np


def block_geometry(h5):
    blocks = h5["blocks"][...]
    nb = (int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]),
          int(h5.attrs["block_nb_z"]))
    if int(blocks[:, 3].max()) != 0:
        raise SystemExit("make_pipe_ic: single-level files only")
    return blocks, nb


def mean_profile(r, R, u_bulk, n):
    """u_c (1 - r/R)^(1/n), with u_c fixed by the bulk velocity."""
    s = np.linspace(0.0, 1.0, 20001)
    shape = (1.0 - s) ** (1.0 / n)
    u_c = u_bulk / (2.0 * np.trapezoid(shape * s, s))
    out = np.zeros_like(r)
    m = r < R
    out[m] = u_c * (1.0 - r[m] / R) ** (1.0 / n)
    return out, u_c


class Potential:
    """A = sum_n a_n sin(k_n.x + phi_n), with the curl taken analytically."""

    def __init__(self, nmodes, lz, R, rng):
        # wavelengths: a few diameters along the axis, ~R across it
        self.k = np.empty((nmodes, 3))
        self.a = np.empty((nmodes, 3))
        self.p = rng.uniform(0.0, 2.0 * np.pi, nmodes)
        for n in range(nmodes):
            kz = 2.0 * np.pi / lz * rng.integers(1, 7)
            kx = np.pi / R * rng.integers(1, 4) * rng.choice([-1.0, 1.0])
            ky = np.pi / R * rng.integers(1, 4) * rng.choice([-1.0, 1.0])
            self.k[n] = (kx, ky, kz)
            self.a[n] = rng.normal(size=3)

    def A_and_grad(self, x, y, z):
        """A (3 arrays) and dA_i/dx_j (3x3 arrays) at the given points."""
        A = [np.zeros_like(x) for _ in range(3)]
        dA = [[np.zeros_like(x) for _ in range(3)] for _ in range(3)]
        for n in range(len(self.p)):
            kx, ky, kz = self.k[n]
            ph = kx * x + ky * y + kz * z + self.p[n]
            s, c = np.sin(ph), np.cos(ph)
            for i in range(3):
                A[i] += self.a[n, i] * s
                for j, kj in enumerate((kx, ky, kz)):
                    dA[i][j] += self.a[n, i] * kj * c
        return A, dA


def curl_component(comp, x, y, z, pot, R, cx, cy):
    """Component `comp` of curl(g A) at the given points, g = (1-(r/R)^2)^2."""
    A, dA = pot.A_and_grad(x, y, z)
    px, py = x - cx, y - cy
    r2 = px * px + py * py
    t = np.maximum(1.0 - r2 / (R * R), 0.0)
    g = t * t
    # dg/dx = 2 t * (-2 px / R^2), same for y; dg/dz = 0
    gx = -4.0 * t * px / (R * R)
    gy = -4.0 * t * py / (R * R)
    gz = np.zeros_like(g)
    dg = (gx, gy, gz)
    # curl(gA)_i = g (curl A)_i + (grad g x A)_i
    i, j, k = comp, (comp + 1) % 3, (comp + 2) % 3
    curlA = dA[k][j] - dA[j][k]
    cross = dg[j] * A[k] - dg[k] * A[j]
    return g * curlA + cross


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("template", help="a snapshot of the case (layout + metadata)")
    ap.add_argument("out")
    ap.add_argument("--centre", type=float, nargs=2, default=(0.65, 0.65))
    ap.add_argument("--radius", type=float, default=0.5)
    ap.add_argument("--u-bulk", type=float, default=1.0)
    ap.add_argument("--power", type=float, default=7.0, help="1/n power law")
    ap.add_argument("--amp", type=float, default=0.15,
                    help="perturbation rms, as a fraction of u_bulk")
    ap.add_argument("--modes", type=int, default=24)
    ap.add_argument("--seed", type=int, default=20260915)
    a = ap.parse_args()

    shutil.copyfile(a.template, a.out)
    cx, cy = a.centre
    R = a.radius
    rng = np.random.default_rng(a.seed)

    with h5py.File(a.out, "r+") as h5:
        blocks, nb = block_geometry(h5)
        node = {d: h5[d][...] for d in "xyz"}
        cent = {d: 0.5 * (node[d][:-1] + node[d][1:]) for d in "xyz"}
        lz = float(h5.attrs["lz"])
        pot = Potential(a.modes, lz, R, rng)

        # Pass 1: the perturbation, to measure its rms before scaling.
        sums = np.zeros(2)
        pert = {}
        for name, comp in (("un", 0), ("vn", 1), ("wn", 2)):
            out = np.empty((blocks.shape[0], nb[2], nb[1], nb[0]))
            for bid, (ox, oy, oz, _l) in enumerate(blocks):
                # component `comp` sits on the LOW face in its own direction
                gx = (node["x"][ox:ox + nb[0]] if comp == 0
                      else cent["x"][ox:ox + nb[0]])
                gy = (node["y"][oy:oy + nb[1]] if comp == 1
                      else cent["y"][oy:oy + nb[1]])
                gz = (node["z"][oz:oz + nb[2]] if comp == 2
                      else cent["z"][oz:oz + nb[2]])
                Z, Y, X = np.meshgrid(gz, gy, gx, indexing="ij")
                q = curl_component(comp, X, Y, Z, pot, R, cx, cy)
                out[bid] = q
                inside = ((X - cx) ** 2 + (Y - cy) ** 2) < R * R
                sums += (float((q[inside] ** 2).sum()), float(inside.sum()))
            pert[name] = out
        rms = np.sqrt(sums[0] / max(sums[1], 1.0))
        scale = a.amp * a.u_bulk / rms
        print(f"   perturbation rms {rms:.4g} -> scaled by {scale:.4g}"
              f" for {a.amp:g} u_bulk")

        # Pass 2: write mean + scaled perturbation.
        u_c = None
        for name, comp in (("un", 0), ("vn", 1), ("wn", 2)):
            dset = h5[name]
            for bid, (ox, oy, oz, _l) in enumerate(blocks):
                gx = (node["x"][ox:ox + nb[0]] if comp == 0
                      else cent["x"][ox:ox + nb[0]])
                gy = (node["y"][oy:oy + nb[1]] if comp == 1
                      else cent["y"][oy:oy + nb[1]])
                gz = (node["z"][oz:oz + nb[2]] if comp == 2
                      else cent["z"][oz:oz + nb[2]])
                Z, Y, X = np.meshgrid(gz, gy, gx, indexing="ij")
                q = scale * pert[name][bid]
                if comp == 2:
                    r = np.hypot(X - cx, Y - cy)
                    w, u_c = mean_profile(r, R, a.u_bulk, a.power)
                    q = q + w
                # nothing moves in the solid: the envelope already vanishes at
                # the wall, but the mean profile is cut there exactly.
                dset[bid] = q
        h5["pn"][...] = 0.0
        h5.attrs["t_current"] = 0.0
        h5.attrs["step"] = 0
    print(f"   {a.out}: 1/{a.power:g} power law, u_c = {u_c:.4f},"
          f" u_bulk = {a.u_bulk:g}, perturbation {a.amp:g} u_bulk,"
          f" {a.modes} modes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
