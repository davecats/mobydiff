#!/usr/bin/env python3
"""Metrics for the A0 freestream gates (docs/next_session_airfoil.md).

  oblique <h5> --u0 U --v0 V   gate (b): max |u-u0|, |v-v0|, |w|, interior div
  pois <io_h5> <ref_h5> [--drift <earlier_io_h5>]
                               gate (c): profile vs the periodic reference,
                               pressure linearity + outlet pin, optional drift
  vortex <h5...>               gate (d): perturbation energy per snapshot;
                               reflected fraction = E(last)/E(first)

Fields are single-level block-table files; rows are reassembled globally.
The domain-boundary HIGH staggered face (the outlet face) lives in solver
halos and is not written, so the divergence check covers cells with all six
stored faces (i < nx-1 etc.); constants + interior exactness pin the rest.
"""
import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from compare_fields import load_field  # noqa: E402


def load(path):
    with h5py.File(path, "r") as f:
        fields = {n: load_field(f, n) for n in ("un", "vn", "wn", "pn")}
        nodes = (f["x"][...], f["y"][...], f["z"][...])
    return fields, nodes


def interior_div(fields, nodes):
    """max |div| over cells whose six faces are all stored (z periodic)."""
    u, v, w = fields["un"], fields["vn"], fields["wn"]
    dx = np.diff(nodes[0]); dy = np.diff(nodes[1]); dz = np.diff(nodes[2])
    nz, ny, nx = u.shape
    up = np.roll(u, -1, axis=2); vp = np.roll(v, -1, axis=1); wp = np.roll(w, -1, axis=0)
    div = ((up - u)/dx[np.newaxis, np.newaxis, :]
           + (vp - v)/dy[np.newaxis, :, np.newaxis]
           + (wp - w)/dz[:, np.newaxis, np.newaxis])
    return float(np.max(np.abs(div[:, :ny - 1, :nx - 1])))


def cmd_oblique(a):
    fields, nodes = load(a.h5)
    du = float(np.max(np.abs(fields["un"] - a.u0)))
    dv = float(np.max(np.abs(fields["vn"] - a.v0)))
    dw = float(np.max(np.abs(fields["wn"])))
    dd = interior_div(fields, nodes)
    print(f"max|u-u0| = {du:.3e}   max|v-v0| = {dv:.3e}   max|w| = {dw:.3e}")
    print(f"interior max|div| = {dd:.3e}")
    ok = du == 0.0 and dv == 0.0 and dw == 0.0 and dd < 1e-12
    print("oblique gate:", "PASS (exact)" if ok else "FAIL")
    return 0 if ok else 1


def profile(fields, xslice):
    """u(y) averaged over z at one x index."""
    return fields["un"][:, :, xslice].mean(axis=0)


def cmd_pois(a):
    io, nodes = load(a.io)
    ref, _ = load(a.ref)
    ny = io["un"].shape[1]
    nx = io["un"].shape[2]
    pr = ref["un"].mean(axis=(0, 2))          # periodic: x-invariant
    peak = float(np.max(pr))
    worst = 0.0
    for frac in (0.5, 0.9):
        xi = int(frac*nx)
        d = float(np.max(np.abs(profile(io, xi) - pr)))/peak
        print(f"profile dev vs reference at x/lx={frac}: {d:.3e}")
        worst = max(worst, d)
    # pressure: z,y-averaged p(x) should be linear, ~0 at the outlet end
    px = io["pn"].mean(axis=(0, 1))
    xc = 0.5*(nodes[0][:-1] + nodes[0][1:])
    fit = np.polyfit(xc, px, 1)
    resid = float(np.max(np.abs(px - np.polyval(fit, xc))))
    print(f"p(x): slope = {fit[0]:.4e} (theory {-8.0/100.0:.4e}), "
          f"nonlinearity = {resid:.2e}, last-cell p = {px[-1]:.3e}")
    status = 0 if worst < 2e-2 else 1
    if a.drift:
        prev, _ = load(a.drift)
        d = float(np.max(np.abs(io["un"] - prev["un"])))
        print(f"drift max|u(t2)-u(t1)| = {d:.3e}")
        status |= 0 if d < 1e-8 else 1
    print("pois gate:", "PASS" if status == 0 else "FAIL")
    return status


def cmd_vortex(a):
    e0 = None
    for path in a.h5:
        fields, nodes = load(path)
        with h5py.File(path, "r") as f:
            t = float(f.attrs.get("t", np.nan))
        du = fields["un"] - 1.0
        e = float(np.sum(du*du) + np.sum(fields["vn"]**2) + np.sum(fields["wn"]**2))
        if e0 is None:
            e0 = e
        print(f"{os.path.basename(path):32s} t={t:8.4f}  E_pert={e:.6e}  E/E0={e/e0:.3e}")
    print(f"reflected fraction (last/first) = {e/e0:.3e}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    o = sub.add_parser("oblique")
    o.add_argument("h5")
    o.add_argument("--u0", type=float, required=True)
    o.add_argument("--v0", type=float, required=True)
    o.set_defaults(func=cmd_oblique)
    p = sub.add_parser("pois")
    p.add_argument("io")
    p.add_argument("ref")
    p.add_argument("--drift", default=None)
    p.set_defaults(func=cmd_pois)
    v = sub.add_parser("vortex")
    v.add_argument("h5", nargs="+")
    v.set_defaults(func=cmd_vortex)
    a = ap.parse_args()
    sys.exit(a.func(a))


if __name__ == "__main__":
    main()
