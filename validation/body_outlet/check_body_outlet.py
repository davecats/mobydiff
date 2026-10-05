#!/usr/bin/env python3
"""Checks for the body-at-outlet gates (docs/next_session_body_at_outlet.md).

    ./check_body_outlet.py profile  pois_body_20000.h5 pois_body_periodic_20000.h5 pois_body.case.h5
    ./check_body_outlet.py solid    pois_body_20000.h5 pois_body.case.h5
    ./check_body_outlet.py mirror   pois_body_20000.h5 pois_body_low_20000.h5
    ./check_body_outlet.py permute  pois_body_20000.h5 pois_body_y_20000.h5 y   (or z)

profile  the inflow/outflow profile against the periodic twin (same immersed
         walls) near the outlet and ON the outlet face; the pressure along the
         gap centre: linear, pinned to the outlet value, last-cell p = G dx/2.
solid    the velocity at every solid-centred staggered location (coef >=
         1e20): interior and on the outlet plane (the un_xmax rows).
mirror   the outlet on the low face: u(i) = -u(nx + 2 - i) etc. against the
         x_max case.
permute  the stream along y (walls normal to x) or z (walls normal to y)
         against the x case under the coordinate permutation.
same     two snapshots equal at max_abs 0 incl. the outlet planes (rank,
         device, restart and zero-coefficient-twin gates).

Fields are reassembled with tools/compare_fields.py (arrays as [z, y, x]);
the solid marker comes from the case file's coef_blocks (leaf, i, j, k, var,
ghost-inclusive). Thresholds are set from the measurements recorded in
README.md; every number is printed.
"""

import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from compare_fields import load_field, block_geometry  # noqa: E402

SOLID = 1.0e20


def load(fn):
    h = h5py.File(fn, "r")
    u, v, w, p = (load_field(h, k) for k in ("un", "vn", "wn", "pn"))
    return h, u, v, w, p


def outlet_plane(h, name="un_xmax"):
    """The high-face plane rows reassembled as [z, y] (x) / [z, x] (y) / [y, x] (z)."""
    blocks, nb, lmax, mzyx, shape = block_geometry(h)
    nz, ny, nx = shape
    if name not in h:
        return None
    rows = h[name][...]
    if name == "un_xmax":
        out = np.zeros((nz, ny)); n1, n2 = nb[1], nb[2]
    elif name == "vn_ymax":
        out = np.zeros((nz, nx)); n1, n2 = nb[0], nb[2]
    else:
        out = np.zeros((ny, nx)); n1, n2 = nb[0], nb[1]
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        if name == "un_xmax" and ox + nb[0] == nx:
            out[oz:oz + nb[2], oy:oy + nb[1]] = rows[bid].reshape(n2, n1)
        if name == "vn_ymax" and oy + nb[1] == ny:
            out[oz:oz + nb[2], ox:ox + nb[0]] = rows[bid].reshape(n2, n1)
        if name == "wn_zmax" and oz + nb[2] == nz:
            out[oy:oy + nb[1], ox:ox + nb[0]] = rows[bid].reshape(n2, n1)
    return out


def solid_masks(case_fn, h):
    """Boolean [z, y, x] masks of solid-centred u, v, w locations (interior
    indices 1..nb of every surviving leaf) and the x_max plane rows (u at
    nb+1) -- from coef_blocks (leaf, i, j, k, var), ghost-inclusive."""
    c = h5py.File(case_fn, "r")
    coef = c["coef_blocks"][...]
    cblocks = c["blocks"][...]
    blocks, nb, lmax, mzyx, shape = block_geometry(h)
    nz, ny, nx = shape
    masks = [np.zeros(shape, dtype=bool) for _ in range(3)]
    plane = np.zeros((nz, ny), dtype=bool)      # u at nb+1 of the x_max blocks
    plane_lo = np.zeros((nz, ny), dtype=bool)   # u at 1 of the x_min blocks
    # match case-file leaves to field leaves on (origin, level)
    key = {tuple(int(t) for t in row): i for i, row in enumerate(cblocks)}
    for (ox, oy, oz, lev) in blocks:
        n = key[(int(ox), int(oy), int(oz), int(lev))]
        for var in range(3):
            tile = coef[n, 1:nb[0] + 1, 1:nb[1] + 1, 1:nb[2] + 1, var]   # (i, j, k)
            masks[var][oz:oz + nb[2], oy:oy + nb[1], ox:ox + nb[0]] = np.transpose(tile, (2, 1, 0)) > SOLID
        if ox + nb[0] == nx:
            tile = coef[n, nb[0] + 1, 1:nb[1] + 1, 1:nb[2] + 1, 0]       # (j, k) at i = nb+1
            plane[oz:oz + nb[2], oy:oy + nb[1]] = tile.T > SOLID
        if ox == 0:
            tile = coef[n, 1, 1:nb[1] + 1, 1:nb[2] + 1, 0]               # (j, k) at i = 1
            plane_lo[oz:oz + nb[2], oy:oy + nb[1]] = tile.T > SOLID
    return masks, (plane, plane_lo)


def cmd_profile(a):
    h, u, v, w, p = load(a[0])
    hp, up, vp, wp, pp = load(a[1])
    masks, _ = solid_masks(a[2], h)
    nz, ny, nx = u.shape
    x = h["x"][...]; y = h["y"][...]
    dx = x[1] - x[0]; dy = y[1] - y[0]
    yc = 0.5 * (y[:-1] + y[1:])
    # the periodic twin: x-averaged developed profile (the field is x-invariant)
    uper = up.mean(axis=(0, 2))
    uio_90 = u[:, :, int(0.9 * nx)].mean(axis=0)
    uio_75 = u[:, :, int(0.75 * nx)].mean(axis=0)
    uout = outlet_plane(h).mean(axis=0)
    fluid = ~masks[0].any(axis=(0, 2))
    def dev(a_, b_):
        return np.abs(a_ - b_)[fluid].max()
    bulk_io = uout[fluid].sum() * dy
    bulk_per = uper[fluid].sum() * dy
    print(f"bulk flux (fluid rows): io outlet {bulk_io:.6e}  periodic {bulk_per:.6e}  ratio {bulk_io / bulk_per:.6f}")
    print(f"max|u_io(x/lx=0.9) - u_per|   {dev(uio_90, uper):.3e}   (peak {uper.max():.4f})")
    print(f"max|u_io(outlet face) - u_per| {dev(uout, uper):.3e}")
    print(f"max|u_io(outlet) - u_io(0.75)| {dev(uout, uio_75):.3e}   (developed-ness)")
    print(f"max|u_io(outlet) - u_io(0.9)|  {dev(uout, uio_90):.3e}")
    # pressure along the gap centre
    jc = ny // 2
    pc = p[:, jc, :].mean(axis=0)
    xc = 0.5 * (x[:-1] + x[1:])
    sel = slice(int(0.25 * nx), int(0.75 * nx))
    G = -np.polyfit(xc[sel], pc[sel], 1)[0]
    resid = np.abs(pc[sel] - np.polyval(np.polyfit(xc[sel], pc[sel], 1), xc[sel])).max()
    print(f"pressure gradient G (fit x/lx 0.25..0.75) {G:.6e}  (forcing of the twin 1.2); linearity residual {resid:.2e}")
    print(f"last-cell p {pc[-1]:.6e}   G dx/2 {G * dx / 2:.6e}   ratio {pc[-1] / (G * dx / 2):.5f}")
    print(f"p(x) at the gap centre, last 4 cells: {pc[-4:]}")
    print(f"max|v| in the fluid {np.abs(v[~masks[1]]).max():.3e}   max|w| {np.abs(w[~masks[2]]).max():.3e}")


def cmd_solid(a):
    h, u, v, w, p = load(a[0])
    masks, (plane, plane_lo) = solid_masks(a[1], h)
    for name, f, m in (("u", u, masks[0]), ("v", v, masks[1]), ("w", w, masks[2])):
        if m.any():
            print(f"max|{name}| at solid-centred {name} locations: {np.abs(f[m]).max():.3e}   ({m.sum()} of them)")
        else:
            print(f"no solid-centred interior {name} location")
    # the outlet plane: x_max (un_xmax rows) or, for the mirrored case, x_min (index 0 of un)
    uout = outlet_plane(h)
    which = "x_max"
    if uout is None:
        uout = u[:, :, 0]; plane = plane_lo; which = "x_min"
    if plane.any():
        print(f"max|u| on the {which} plane, solid rows: {np.abs(uout[plane]).max():.3e}   ({plane.sum()} faces)")
    else:
        print(f"no solid face on the {which} plane")
    print(f"max|u| on the {which} plane, fluid rows: {np.abs(uout[~plane]).max():.3e}")
    # the neighbour plane (the last interior u face): the "n solid, f fluid" limit shows here
    un_ = u[:, :, -1] if which == "x_max" else u[:, :, 1]
    nsol = masks[0][:, :, -1] if which == "x_max" else masks[0][:, :, 1]
    if nsol.any():
        print(f"max|u| on the neighbour plane, solid rows: {np.abs(un_[nsol]).max():.3e}   ({nsol.sum()} faces); "
              f"the outlet face in those rows: max|u| {np.abs(uout[nsol]).max():.3e}")
    if masks[0].any():
        print(f"max|p| in solid-centred cells (u-marker): {np.abs(p[masks[0]]).max():.3e}; in the fluid {np.abs(p[~masks[0]]).max():.3e}")
    # mass balance
    y = h["y"][...]; z = h["z"][...]
    dy = y[1] - y[0]; dz = z[1] - z[0]
    # the inlet flux: the x_min face (index 0 of un), or for the mirrored
    # case the last interior u face (the x_max inlet face itself is in no
    # dataset; the projection makes the two equal to its residual)
    qin = (u[:, :, 0] if which == "x_max" else u[:, :, -1]).sum() * dy * dz
    qout = uout.sum() * dy * dz
    print(f"flux in {qin:.10e}  out {qout:.10e}  imbalance (out-in)/in {(qout - qin) / qin:.3e}")


def cmd_mirror(a):
    h, u, v, w, p = load(a[0])
    hl, ul, vl, wl, pl = load(a[1])
    # the mirrored case: x -> lx - x. u-face i <-> nx - i (interior 1..nx-1); cells i <-> nx-1-i
    d_u = np.abs(ul[:, :, 1:] + u[:, :, 1:][:, :, ::-1]).max()
    d_v = np.abs(vl - v[:, :, ::-1]).max()
    d_w = np.abs(wl - w[:, :, ::-1]).max()
    d_p = np.abs(pl - p[:, :, ::-1]).max()
    uout = outlet_plane(h)
    # the low case's outlet face is index 0 of un; the high case's is un_xmax
    d_f = np.abs(ul[:, :, 0] + uout).max()
    print(f"mirror: u {d_u:.3e}  v {d_v:.3e}  w {d_w:.3e}  p {d_p:.3e}  outlet face {d_f:.3e}")


def cmd_permute(a):
    h, u, v, w, p = load(a[0])
    hq, uq, vq, wq, pq = load(a[1])
    ax = a[2]
    # x case arrays [z, y, x]; stream x, walls normal y, periodic z.
    if ax == "y":
        # stream y, walls normal x, periodic z: (x, y, z) -> (y, x, z): arrays [z, x', y'] = swap axes 1,2
        d = [np.abs(vq - np.swapaxes(u, 1, 2)).max(), np.abs(uq - np.swapaxes(v, 1, 2)).max(),
             np.abs(wq - np.swapaxes(w, 1, 2)).max(), np.abs(pq - np.swapaxes(p, 1, 2)).max()]
        # x case plane [z, y]; y case plane vn_ymax [z, x'] with x' = y
        f = np.abs(outlet_plane(hq, "vn_ymax") - outlet_plane(h)).max()
    else:
        # stream z, walls normal y, periodic x: (x, y, z) -> (z, y, x): arrays swap axes 0,2
        d = [np.abs(wq - np.swapaxes(u, 0, 2)).max(), np.abs(vq - np.swapaxes(v, 0, 2)).max(),
             np.abs(uq - np.swapaxes(w, 0, 2)).max(), np.abs(pq - np.swapaxes(p, 0, 2)).max()]
        # x case plane [z, y]; z case plane wn_zmax [y, x'] with x' = z
        f = np.abs(outlet_plane(hq, "wn_zmax") - outlet_plane(h).T).max()
    print(f"permute {ax}: stream {d[0]:.3e}  normal {d[1]:.3e}  span {d[2]:.3e}  p {d[3]:.3e}  outlet face {f:.3e}")


def cmd_same(a):
    """max_abs over un/vn/wn/pn AND the high-outlet planes of two snapshots;
    exit 1 unless every one is exactly 0 (the bit-exactness form)."""
    h, u, v, w, p = load(a[0])
    hq, uq, vq, wq, pq = load(a[1])
    worst = 0.0
    for name, f, g in (("un", u, uq), ("vn", v, vq), ("wn", w, wq), ("pn", p, pq)):
        d = np.abs(f - g).max(); worst = max(worst, d)
        print(f"   {name} max_abs={d:.16e}")
    for name in ("un_xmax", "vn_ymax", "wn_zmax"):
        if name in h or name in hq:
            a_ = outlet_plane(h, name); b_ = outlet_plane(hq, name)
            d = np.abs(a_ - b_).max() if a_ is not None and b_ is not None else np.inf
            worst = max(worst, d)
            print(f"   {name} max_abs={d:.16e}")
    print("   PASS (max_abs 0)" if worst == 0.0 else "   FAIL")
    sys.exit(0 if worst == 0.0 else 1)


if __name__ == "__main__":
    cmd = sys.argv[1]
    {"profile": cmd_profile, "solid": cmd_solid, "mirror": cmd_mirror,
     "permute": cmd_permute, "same": cmd_same}[cmd](sys.argv[2:])
