#!/usr/bin/env python3
"""Controlled initial conditions for the AoA = 0 startup tests
(2026-07-16): both are restart files carrying the run's exact
12088-leaf block layout (template = a completed aoa0 snapshot).

  ic_zero.h5    u = v = w = p = 0 everywhere: the flow enters only
                through the inlet boundary conditions (no impulsive
                interior, no body-interior velocity to annihilate).
  ic_masked.h5  uniform freestream (1, 0, 0) but ZEROED inside the body
                and within a 2-fine-cell layer around it (analytic 0012
                signed distance <= 2*Delta_fine): the impulsive interior
                velocity that the penalization would otherwise slam to
                zero -- the suspected ringing source -- is removed while
                the far field starts established.

k/omega/nut are DELETED so RANS reinitializes cleanly from tu/nut_ratio
(a zeroed omega would poison the 1/omega terms); step/t_current reset
to 0.

Usage: make_ic_tests.py <template_snapshot.h5>
"""
import shutil
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

XLE, ZLE, CHORD = 4.5, 6.0, 1.0


def polyline(n=4096):
    tt = 0.12
    beta = np.linspace(0.0, np.pi, n//2)
    x = 0.5*(1.0 - np.cos(beta))
    yt = 5.0*tt*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2
                 + 0.2843*x**3 - 0.1036*x**4)
    xs = np.concatenate([x[::-1], x[1:-1]])
    zs = np.concatenate([yt[::-1], -yt[1:-1]])
    pts = np.column_stack([XLE + CHORD*xs, ZLE + CHORD*zs])
    t = np.gradient(pts, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    nrm = np.column_stack([t[:, 1], -t[:, 0]])   # outward (CCW loop)
    return pts, nrm


def reset_common(h5):
    for name in ("k", "omega", "nut", "gamma", "rethetat", "fd"):
        if name in h5:
            del h5[name]
    h5.attrs["step"] = np.int32(0)
    h5.attrs["t_current"] = 0.0


def main():
    template = sys.argv[1]

    # ---- ic_zero.h5 -------------------------------------------------
    shutil.copyfile(template, "ic_zero.h5")
    with h5py.File("ic_zero.h5", "r+") as h5:
        for name in ("un", "vn", "wn", "pn"):
            h5[name][...] = 0.0
        reset_common(h5)
    print("wrote ic_zero.h5 (fields = 0 everywhere)")

    # ---- ic_masked.h5 -----------------------------------------------
    shutil.copyfile(template, "ic_masked.h5")
    with h5py.File("ic_masked.h5", "r+") as h5:
        blocks = h5["blocks"][...]
        nb = int(h5.attrs["block_nb_x"])
        nx = int(h5.attrs["nx"]); ny = int(h5.attrs["ny"]); nz = int(h5.attrs["nz"])
        lx = float(h5.attrs["lx"]); lz = float(h5.attrs["lz"])
        mask = np.asarray(h5.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        lmax = int(blocks[:, 3].max())
        hfx = lx/(nx*2**(lmax*int(mask[0])))     # finest spacing
        layer = 2.0*hfx                          # the 2-cell layer

        pts, nrm = polyline()
        tree = cKDTree(pts)

        U, V, W, P = h5["un"], h5["vn"], h5["wn"], h5["pn"]
        nzero = 0
        for bid, (ox, oy, oz, lev) in enumerate(blocks):
            hx = lx/(nx*2**(int(lev)*int(mask[0])))
            hz = lz/(nz*2**(int(lev)*int(mask[2])))
            xc = (ox + 0.5 + np.arange(nb))*hx
            zc = (oz + 0.5 + np.arange(nb))*hz
            # freestream everywhere first
            U[bid] = 1.0
            V[bid] = 0.0
            W[bid] = 0.0
            P[bid] = 0.0
            # bbox gate: keep pure freestream outside the section bbox + margin
            if (xc[-1] < XLE - 0.05 or xc[0] > XLE + CHORD + 0.05 or
                    zc[-1] < ZLE - 0.12 or zc[0] > ZLE + 0.12):
                continue
            Z, X = np.meshgrid(zc, xc, indexing="ij")
            q = np.column_stack([X.ravel(), Z.ravel()])
            _, idx = tree.query(q)
            dn = np.einsum("ij,ij->i", q - pts[idx], nrm[idx]).reshape(nb, nb)
            m = dn <= layer                       # inside or within 2 cells
            if not m.any():
                continue
            u = U[bid][...]
            u[np.broadcast_to(m[:, None, :], (nb, nb, nb))] = 0.0
            U[bid] = u
            nzero += int(m.sum())*nb
        reset_common(h5)
    print(f"wrote ic_masked.h5 (freestream, {nzero} cells zeroed in "
          f"body + {layer:.2e} layer)")


if __name__ == "__main__":
    main()
