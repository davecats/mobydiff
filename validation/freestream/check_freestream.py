#!/usr/bin/env python3
"""Metrics for the A0 freestream gates (docs/next_session_airfoil.md).

  oblique <h5> --u0 U --v0 V   gate (b): max |u-u0|, |v-v0|, |w|, interior div
  pois <io_h5> <ref_h5> [--drift <earlier_io_h5>]
                               gate (c): profile vs the periodic reference,
                               pressure linearity + outlet pin, optional drift
  vortex <h5...>               gate (d): perturbation energy per snapshot;
                               reflected fraction = E(last)/E(first)
  mirror <high_h5> <low_h5>    outlet gate: a run with the outlet at x_max
                               against its mirror image (outlet at x_min),
                               u(i) <-> -u(nx-i), v, p (i) <-> (nx-1-i)
  uniform <h5> --u0 U --v0 V   outlet gate: deviation from uniform flow and
                               the stored pressure range, on any block layout
                               (the zero-initialised refined legs)

pois and vortex take --mirror for a run whose outlet is at x_min (flow in
-x): the fields are mapped into the high-outlet frame first, so the numbers
printed are directly those of the x_max case.

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


def to_high_frame(fields, inlet_u=0.0):
    """Map a run with the outlet at x_min (flow in -x) into the frame of its
    mirror image: x -> lx - x, u -> -u. Stored u faces are 0..nx-1 (face 0 on
    the low boundary), so face i of the mirrored field is face nx - i of the
    run; its face 0 is the run's INLET face nx, which the file does not hold
    (a pinned value: inlet_u, in the mirrored frame)."""
    out = {}
    u = fields["un"]
    m = np.empty_like(u)
    m[:, :, 1:] = -u[:, :, :0:-1]
    m[:, :, 0] = inlet_u
    out["un"] = m
    for n in ("vn", "wn", "pn"):
        out[n] = fields[n][:, :, ::-1]
    return out


def load(path, mirror=False, inlet_u=0.0):
    with h5py.File(path, "r") as f:
        fields = {n: load_field(f, n) for n in ("un", "vn", "wn", "pn")}
        nodes = (f["x"][...], f["y"][...], f["z"][...])
    if mirror:
        fields = to_high_frame(fields, inlet_u)
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
    io, nodes = load(a.io, a.mirror)
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
        prev, _ = load(a.drift, a.mirror)
        d = float(np.max(np.abs(io["un"] - prev["un"])))
        print(f"drift max|u(t2)-u(t1)| = {d:.3e}")
        status |= 0 if d < 1e-8 else 1
    print("pois gate:", "PASS" if status == 0 else "FAIL")
    return status


def cmd_vortex(a):
    e0 = None
    for path in a.h5:
        fields, nodes = load(path, a.mirror, inlet_u=1.0)
        with h5py.File(path, "r") as f:
            t = float(f.attrs.get("t_current", np.nan))
        du = fields["un"] - 1.0
        e = float(np.sum(du*du) + np.sum(fields["vn"]**2) + np.sum(fields["wn"]**2))
        if e0 is None:
            e0 = e
        print(f"{os.path.basename(path):32s} t={t:8.4f}  E_pert={e:.6e}  E/E0={e/e0:.3e}")
    print(f"reflected fraction (last/first) = {e/e0:.3e}")
    return 0


def cmd_mirror(a):
    hi, _ = load(a.high)
    lo, _ = load(a.low, mirror=True)
    # face 0 of either frame is an inlet face the other file does not hold
    du = float(np.max(np.abs(hi["un"][:, :, 1:] - lo["un"][:, :, 1:])))
    dv = float(np.max(np.abs(hi["vn"] - lo["vn"])))
    dw = float(np.max(np.abs(hi["wn"] - lo["wn"])))
    dp = float(np.max(np.abs(hi["pn"] - lo["pn"])))
    print(f"high vs mirrored low: max|du| = {du:.3e}  max|dv| = {dv:.3e}  "
          f"max|dw| = {dw:.3e}  max|dp| = {dp:.3e}")
    ok = max(du, dv, dw, dp) <= a.tol
    print(f"mirror gate (tol {a.tol:.1e}):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_uniform(a):
    fields, _ = load(a.h5)
    du = float(np.max(np.abs(fields["un"] - a.u0)))
    dv = float(np.max(np.abs(fields["vn"] - a.v0)))
    dw = float(np.max(np.abs(fields["wn"])))
    pmin, pmax = float(fields["pn"].min()), float(fields["pn"].max())
    print(f"max|u-u0| = {du:.3e}   max|v-v0| = {dv:.3e}   max|w| = {dw:.3e}   "
          f"p in [{pmin:.3e}, {pmax:.3e}]")
    ok = max(du, dv, dw, abs(pmin), abs(pmax)) <= a.tol
    print(f"uniform gate (tol {a.tol:.1e}):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


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
    p.add_argument("--mirror", action="store_true")
    p.set_defaults(func=cmd_pois)
    v = sub.add_parser("vortex")
    v.add_argument("h5", nargs="+")
    v.add_argument("--mirror", action="store_true")
    v.set_defaults(func=cmd_vortex)
    m = sub.add_parser("mirror")
    m.add_argument("high")
    m.add_argument("low")
    m.add_argument("--tol", type=float, default=1e-12)
    m.set_defaults(func=cmd_mirror)
    u = sub.add_parser("uniform")
    u.add_argument("h5")
    u.add_argument("--u0", type=float, required=True)
    u.add_argument("--v0", type=float, required=True)
    u.add_argument("--tol", type=float, default=1e-9)
    u.set_defaults(func=cmd_uniform)
    a = ap.parse_args()
    sys.exit(a.func(a))


if __name__ == "__main__":
    main()
