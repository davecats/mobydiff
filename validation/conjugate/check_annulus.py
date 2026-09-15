#!/usr/bin/env python3
"""F5 gate: the ANNULAR pipe body, checked out of the prepared case file.

The pipe's immersed body is everything OUTSIDE r_inner -- the fluid is the
hole -- and make_geometry_stl.py closes that solid with an outer surface that
is pure padding. Three things then have to hold, and this checks all three
against the analytic geometry, using no solver output at all:

  class   every cell is solid exactly when its centre is outside the FACETED
          inner polygon. The facets ARE the geometry (the cylinder command's
          docstring), so the reference is the m-gon, not the ideal circle.
  phi     |dwall| equals the analytic distance to that polygon. This is the
          one number the conjugate scheme reads: its sign picks the material
          and its magnitude sets the cut-face level-set weight w.
  near    inside the domain the PIPE WALL is the nearest surface everywhere,
          i.e. the padding never wins. docs/next_session_pipe_cht.md F2 says
          to verify this rather than assume it, because phi in a shell is the
          distance to the nearest of the two surfaces -- were the padding
          nearer anywhere, phi there would describe a fictitious wall.

Case-file tiles are ghost-inclusive, shape (nBlocks, nb+2, nb+2, nb+2) in
(i, j, k) order -- NOT the (k, j, i) of the field datasets; see
check_oblique.py load_phi, which is where that trap was first paid for. A
pipe can tell them apart (it is invariant in z but not in x or y), and does:
with the axes swapped the phi error below is 0.18 rather than round-off.
The ghost cells are included in every check: a cut arm may sit on a block
boundary, which is what the ghost layer exists for.

    ./check_annulus.py pipe64.h5 --centre 0.65 0.65 --r-inner 0.5 \
        --facets 16384 --box-half 1.25 --domain-half 0.65
"""

from __future__ import annotations

import argparse
import h5py
import numpy as np

SOLID_THRESHOLD = 1e20


def polygon_distance(px, py, R, m):
    """Signed distance to a regular m-gon of CIRCUMRADIUS R centred at the
    origin, positive outside. Exact, by folding into one sector: the edge is
    the segment x = R cos(pi/m), |y| <= R sin(pi/m) in the frame whose x axis
    runs through the edge midpoint."""
    r = np.hypot(px, py)
    sector = 2.0 * np.pi / m
    # Angle measured from the nearest edge MIDPOINT, i.e. folded into
    # [-pi/m, pi/m]. Vertices sit at angles 0, sector, ... (make_geometry_stl),
    # so edge midpoints are offset by half a sector.
    th = np.arctan2(py, px) - 0.5 * sector
    th = th - sector * np.round(th / sector)
    a = R * np.cos(0.5 * sector)          # apothem
    hw = R * np.sin(0.5 * sector)         # half edge length
    x, y = r * np.cos(th), r * np.abs(np.sin(th))
    # Nearest point on the edge segment: its interior if |y| <= hw, else the
    # vertex. The sign is that of (x - a) only in the first case; outside the
    # slab the point is beyond the vertex and therefore outside the polygon
    # whenever x >= a, which the explicit inside test below settles.
    d_edge = np.abs(x - a)
    d_vert = np.hypot(x - a, y - hw)
    d = np.where(y <= hw, d_edge, d_vert)
    inside = (y <= hw) & (x <= a) | (y > hw) & (x * hw + y * a <= a * hw + hw * a - 1e-300) & (x < a)
    # The robust inside test: a convex polygon point is inside iff it is on
    # the inner side of the nearest edge line, which in this frame is x <= a.
    inside = x <= a
    return np.where(inside, -d, d)


def box_distance(px, py, h):
    """Distance to the boundary of the square [-h, h]^2, for points inside
    it (the padding is outside the domain, so that is the only case)."""
    return np.minimum(h - np.abs(px), h - np.abs(py))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--centre", type=float, nargs=2, required=True)
    ap.add_argument("--r-inner", type=float, required=True)
    ap.add_argument("--facets", type=int, required=True)
    ap.add_argument("--box-half", type=float, default=None)
    ap.add_argument("--r-outer", type=float, default=None)
    ap.add_argument("--domain-half", type=float, required=True)
    ap.add_argument("--axis", choices=("x", "y", "z"), default="z")
    ap.add_argument("--tol-phi", type=float, default=1e-9)
    a = ap.parse_args()

    with h5py.File(a.case, "r") as h5:
        blocks = h5["blocks"][...]
        nb = int(h5.attrs["block_nb"])
        node = {d: h5[f"{d}_nodes"][...] for d in "xyz"}
        dwall = h5["dwall_blocks"][...]
        coef = h5["coef_p_blocks"][...]
    if int(blocks[:, 3].max()) != 0:
        raise SystemExit("check_annulus: single-level case files only")

    # Ghost-inclusive cell centres per direction, as in check_conjugate.py.
    cent = {}
    for d in "xyz":
        n = node[d]
        ext = np.concatenate(([2 * n[0] - n[1]], n, [2 * n[-1] - n[-2]]))
        cent[d] = 0.5 * (ext[:-1] + ext[1:])

    iw = "xyz".index(a.axis)
    ic = [d for d in range(3) if d != iw]
    cx, cy = a.centre

    worst_phi, n_cells, n_flip = 0.0, 0, 0
    worst_near, n_dom = -1e30, 0
    for bid in range(blocks.shape[0]):
        org = [int(blocks[bid, d]) for d in range(3)]
        # tile axes are (i, j, k) = (x, y, z)
        c = [cent["xyz"[d]][org[d]:org[d] + nb + 2] for d in range(3)]
        grid = list(np.meshgrid(c[0], c[1], c[2], indexing="ij"))
        pu, pv = grid[ic[0]] - cx, grid[ic[1]] - cy

        phi_ref = polygon_distance(pu, pv, a.r_inner, a.facets)
        solid = np.abs(coef[bid]) > SOLID_THRESHOLD
        phi = np.where(solid, -np.maximum(dwall[bid], 1e-300), dwall[bid])

        # class: the body is OUTSIDE the polygon, so solid <=> phi_ref > 0.
        # Cells whose centre lands within the chord error of the surface are
        # genuinely ambiguous (the STL is a polygon, the reference the same
        # polygon, but the ray parity and the distance are separate codes);
        # count them rather than let them fail the test.
        amb = np.abs(phi_ref) < 1e-12
        n_flip += int((solid != (phi_ref > 0))[~amb].sum())

        # phi: sign and magnitude together, i.e. phi itself.
        worst_phi = max(worst_phi, float(np.abs(phi - (-phi_ref)).max()))
        n_cells += phi.size

        # near: only cells INSIDE the domain, where the check has meaning.
        ind = (np.abs(pu) <= a.domain_half) & (np.abs(pv) <= a.domain_half)
        if ind.any():
            if a.box_half is not None:
                d_pad = box_distance(pu[ind], pv[ind], a.box_half)
            else:
                d_pad = a.r_outer - np.hypot(pu[ind], pv[ind])
            # margin > 0 means the pipe wall wins everywhere.
            worst_near = max(worst_near, float((np.abs(phi_ref[ind]) - d_pad).max()))
            n_dom += int(ind.sum())

    print(f"   cells = {n_cells} (ghost-inclusive), in-domain = {n_dom}")
    print(f"   class: solid <=> outside the {a.facets}-gon   flips = {n_flip}")
    print(f"   phi:   max|phi - phi_exact| = {worst_phi:.3e}  (tol {a.tol_phi:g})")
    print(f"   near:  max(d_wall - d_padding) = {worst_near:.4f}"
          f"   ({'pipe wall nearest everywhere' if worst_near < 0 else 'PADDING WINS'})")
    ok = n_flip == 0 and worst_phi <= a.tol_phi and worst_near < 0.0
    print("   PASS" if ok else "   FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
