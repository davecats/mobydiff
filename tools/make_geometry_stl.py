#!/usr/bin/env python3
"""ASCII STLs for the conjugate gates: a TILTED half-space, a CIRCULAR
cylinder, and an ANNULAR shell (the pipe).

    ./make_geometry_stl.py plane    out.stl --theta 30 --x0 0.5 --y0 0.5 ...
    ./make_geometry_stl.py cylinder out.stl --centre 0.5 0.5 --radius 0.2 ...
    ./make_geometry_stl.py annulus  out.stl --axis z --r-inner 0.5 --box-half 1.1 ...

ASCII for the same reason as `make_slab_stl.py`: moby_prepare parses ASCII
vertices straight to float64, so the surface is exactly the analytic one and
the reference needs no quantisation dance.

Both bodies are CLOSED solids that extend beyond the domain, so that inside
the domain the nearest surface point is always on the surface of interest --
the tilted plane, or the cylinder's skin -- and never on the padding.

KEEP THE BODIES TIGHT. The BVH point-triangle distance forms d^2 as a
difference of terms of the vertices' magnitude, so padding costs precision
quadratically (measured in C1: +-4 loses ~64x more than +-0.5). The defaults
here are the smallest that still keep the padding faces out of range, and the
gates report |phi - exact| so the floor is visible rather than assumed.
"""

from __future__ import annotations

import argparse
import math


def write_stl(path, triangles, name="body"):
    """triangles = iterable of (normal, (v0, v1, v2))."""
    with open(path, "w") as fh:
        fh.write(f"solid {name}\n")
        n_tri = 0
        for n, tri in triangles:
            fh.write("  facet normal %.17g %.17g %.17g\n" % tuple(n))
            fh.write("    outer loop\n")
            for p in tri:
                fh.write("      vertex %.17g %.17g %.17g\n" % tuple(p))
            fh.write("    endloop\n")
            fh.write("  endfacet\n")
            n_tri += 1
        fh.write(f"endsolid {name}\n")
    return n_tri


def box_faces(v):
    """The 12 triangles of a hexahedron given its 8 corners in the standard
    order (low face 0..3 counter-clockwise seen from outside, high face 4..7).
    Normals are computed from the vertices, so a ROTATED box needs no extra
    care."""
    quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (3, 7, 6, 2), (0, 4, 7, 3), (1, 2, 6, 5)]
    for a, b, c, d in quads:
        for tri in ((a, b, c), (a, c, d)):
            p, q, r = (v[t] for t in tri)
            e1 = [q[m] - p[m] for m in range(3)]
            e2 = [r[m] - p[m] for m in range(3)]
            nx = e1[1] * e2[2] - e1[2] * e2[1]
            ny = e1[2] * e2[0] - e1[0] * e2[2]
            nz = e1[0] * e2[1] - e1[1] * e2[0]
            s = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
            yield (nx / s, ny / s, nz / s), (p, q, r)


def cmd_plane(a):
    """A half-space below the plane through (x0, y0) with normal
    n = (-sin t, cos t, 0), realised as a rotated box whose TOP face is that
    plane and whose other five faces sit outside the domain."""
    th = math.radians(a.theta)
    nrm = (-math.sin(th), math.cos(th), 0.0)
    tan = (math.cos(th), math.sin(th), 0.0)

    def pt(s, d, z):
        """s along the plane, d along -n (into the solid), z as given."""
        return (a.x0 + s * tan[0] - d * nrm[0],
                a.y0 + s * tan[1] - d * nrm[1], z)

    corners = [pt(-a.span, a.depth, a.z0), pt(a.span, a.depth, a.z0),
               pt(a.span, 0.0, a.z0), pt(-a.span, 0.0, a.z0),
               pt(-a.span, a.depth, a.z1), pt(a.span, a.depth, a.z1),
               pt(a.span, 0.0, a.z1), pt(-a.span, 0.0, a.z1)]
    n = write_stl(a.out, box_faces(corners), "plane")
    print(f"{a.out}: half-space below the plane through ({a.x0!r}, {a.y0!r}) "
          f"at {a.theta:g} deg, {n} triangles")


def cmd_cylinder(a):
    """A circular cylinder with its axis along z, faceted with `--facets`
    segments. The caps sit outside the domain, so the skin is the only
    surface the distance can see -- and the facets are the geometry: the
    gate's reference must use the SAME polygon, not the ideal circle, or it
    measures the faceting instead of the scheme."""
    cx, cy = a.centre
    m = a.facets
    ang = [2.0 * math.pi * i / m for i in range(m)]
    ring0 = [(cx + a.radius * math.cos(t), cy + a.radius * math.sin(t), a.z0) for t in ang]
    ring1 = [(cx + a.radius * math.cos(t), cy + a.radius * math.sin(t), a.z1) for t in ang]
    tris = []
    for i in range(m):
        j = (i + 1) % m
        # ang[i] + pi/m, NOT the mean of ang[i] and ang[j]: the latter is 180
        # degrees out at the wrap-around facet, where ang[j] has come back to
        # 0. Cosmetic for the solver -- geometry_stl.f90 ignores the stored
        # normal, parity being orientation-independent -- but wrong is wrong.
        nx = math.cos(ang[i] + math.pi / m)
        ny = math.sin(ang[i] + math.pi / m)
        tris.append(((nx, ny, 0.0), (ring0[i], ring0[j], ring1[j])))
        tris.append(((nx, ny, 0.0), (ring0[i], ring1[j], ring1[i])))
    for i in range(1, m - 1):            # the two caps, fanned
        tris.append(((0.0, 0.0, -1.0), (ring0[0], ring0[i + 1], ring0[i])))
        tris.append(((0.0, 0.0, 1.0), (ring1[0], ring1[i], ring1[i + 1])))
    n = write_stl(a.out, tris, "cylinder")
    print(f"{a.out}: cylinder r = {a.radius!r} at ({cx!r}, {cy!r}), "
          f"{m} facets, {n} triangles")


# Cyclic frames: (u, v, w) with w along the axis, all right-handed, so the
# winding built below survives the mapping and the normals come out outward.
# `centre` is taken in ASCENDING coordinate order, which for the y axis is the
# other way round from the cyclic pair -- swapped here, once.
AXIS_FRAME = {"x": (1, 2, 0), "y": (2, 0, 1), "z": (0, 1, 2)}


def cmd_annulus(a):
    """A closed ANNULAR shell: the solid between `--r-inner` and an outer
    surface, with its axis along `--axis` and its ends beyond the domain.

    This is the pipe's immersed body. The FLUID is the hole -- the body is
    everything outside r_inner -- so the outer surface is padding whose only
    job is to close the solid, and it must stay far enough out that no cell
    centre INSIDE the domain is nearer to it than to the pipe wall. Otherwise
    phi = +-dwall, which in a shell is the distance to the NEAREST of the two
    surfaces, would report a fictitious wall (the conjugate scheme reads phi's
    sign to pick the material and its magnitude only at CUT faces, so the
    damage is confined -- but --domain-half checks it rather than assuming).

    A BOX outer surface (--box-half) is tighter than a cylindrical one
    (--r-outer) for a square cross-section: the binding case is the domain
    corner, and the box clears it at a smaller vertex magnitude. Keep it tight
    -- the module docstring's precision argument applies here too.

    The facets ARE the geometry: an inscribed m-gon, chord error
    r (1 - cos(pi/m)), printed below. Any reference must use the same polygon.
    """
    if (a.r_outer is None) == (a.box_half is None):
        raise SystemExit("annulus: give exactly one of --r-outer / --box-half")
    iu, iv, iw = AXIS_FRAME[a.axis]
    cu, cv = (a.centre[1], a.centre[0]) if a.axis == "y" else a.centre

    def P(u, v, w):
        """Place a canonical-frame point into world coordinates."""
        p = [0.0, 0.0, 0.0]
        p[iu], p[iv], p[iw] = u, v, w
        return tuple(p)

    def N(nu, nv, nw):
        return P(nu, nv, nw)

    m = a.facets
    ang = [2.0 * math.pi * i / m for i in range(m)]
    ri = a.r_inner
    inner = [[P(cu + ri * math.cos(t), cv + ri * math.sin(t), w) for t in ang]
             for w in (a.a0, a.a1)]
    if a.r_outer is not None:
        ro = a.r_outer
        outer = [[P(cu + ro * math.cos(t), cv + ro * math.sin(t), w) for t in ang]
                 for w in (a.a0, a.a1)]
    else:
        # The outer surface is a square; the end caps below still connect
        # ring to ring, so the square is sampled at the SAME m angles by
        # projecting each ray onto the box. Corners land exactly on a facet
        # boundary only if m % 4 == 0, which is required.
        if m % 4:
            raise SystemExit("annulus: --box-half needs --facets divisible by 4")
        h = a.box_half

        def on_box(t):
            c, s_ = math.cos(t), math.sin(t)
            return h / max(abs(c), abs(s_)) * c, h / max(abs(c), abs(s_)) * s_

        outer = [[P(cu + q[0], cv + q[1], w) for q in map(on_box, ang)]
                 for w in (a.a0, a.a1)]

    i0, i1 = inner
    o0, o1 = outer
    tris = []
    for i in range(m):
        j = (i + 1) % m
        # Inner skin: normals point INTO the hole, -r. (Mid-angle as
        # ang[i] + pi/m -- see cmd_cylinder on the wrap-around facet.)
        nu = -math.cos(ang[i] + math.pi / m)
        nv = -math.sin(ang[i] + math.pi / m)
        tris.append((N(nu, nv, 0.0), (i0[i], i1[j], i0[j])))
        tris.append((N(nu, nv, 0.0), (i0[i], i1[i], i1[j])))
        # Outer skin: normals point out. Computed from the vertices, because
        # for the box form the facet normal is not the mid-angle direction.
        p, q, r = o0[i], o0[j], o1[j]
        du, dv = q[iu] - p[iu], q[iv] - p[iv]
        sn = math.hypot(du, dv) or 1.0
        tris.append((N(dv / sn, -du / sn, 0.0), (o0[i], o0[j], o1[j])))
        tris.append((N(dv / sn, -du / sn, 0.0), (o0[i], o1[j], o1[i])))
        # The two annular end caps.
        tris.append((N(0.0, 0.0, -1.0), (o0[i], i0[i], i0[j])))
        tris.append((N(0.0, 0.0, -1.0), (o0[i], i0[j], o0[j])))
        tris.append((N(0.0, 0.0, 1.0), (o1[i], i1[j], i1[i])))
        tris.append((N(0.0, 0.0, 1.0), (o1[i], o1[j], i1[j])))
    n = write_stl(a.out, tris, "annulus")

    chord = ri * (1.0 - math.cos(math.pi / m))
    outer_desc = (f"r_outer = {a.r_outer!r}" if a.r_outer is not None
                  else f"box half-width {a.box_half!r}")
    print(f"{a.out}: annulus about {a.axis}, r_inner = {ri!r}, {outer_desc}, "
          f"{a.axis} in [{a.a0!r}, {a.a1!r}], {m} facets, {n} triangles")
    print(f"  inner chord error {chord:.3g}")
    if a.domain_half is not None:
        # The binding point is the cross-section corner: furthest from the
        # pipe wall, nearest to the outer surface.
        d = a.domain_half
        corner_r = math.hypot(d, d)
        to_inner = corner_r - ri
        to_outer = (a.r_outer - corner_r) if a.r_outer is not None else (a.box_half - d)
        ok = to_outer > to_inner
        print(f"  nearest-surface check at the domain corner (|x| = |y| = {d!r}): "
              f"to the pipe wall {to_inner:.4f}, to the outer surface {to_outer:.4f}"
              f"  -> {'OK' if ok else 'THE OUTER SURFACE WINS -- pad it further out'}")
        if not ok:
            raise SystemExit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("plane")
    p.add_argument("out")
    p.add_argument("--theta", type=float, required=True, help="tilt from the x axis, degrees")
    p.add_argument("--x0", type=float, default=0.5)
    p.add_argument("--y0", type=float, default=0.5)
    p.add_argument("--span", type=float, default=2.0, help="half-extent along the plane")
    p.add_argument("--depth", type=float, default=2.0, help="extent into the solid")
    p.add_argument("--z0", type=float, default=-0.5)
    p.add_argument("--z1", type=float, default=0.5)
    p.set_defaults(func=cmd_plane)

    p = sub.add_parser("cylinder")
    p.add_argument("out")
    p.add_argument("--centre", type=float, nargs=2, default=(0.5, 0.5))
    p.add_argument("--radius", type=float, required=True)
    p.add_argument("--facets", type=int, default=2048)
    p.add_argument("--z0", type=float, default=-0.5)
    p.add_argument("--z1", type=float, default=0.5)
    p.set_defaults(func=cmd_cylinder)

    p = sub.add_parser("annulus")
    p.add_argument("out")
    p.add_argument("--axis", choices=("x", "y", "z"), default="z")
    p.add_argument("--centre", type=float, nargs=2, default=(0.0, 0.0),
                   help="cross-plane centre, in ascending coordinate order")
    p.add_argument("--r-inner", type=float, required=True)
    p.add_argument("--r-outer", type=float, default=None,
                   help="cylindrical outer surface (exclusive with --box-half)")
    p.add_argument("--box-half", type=float, default=None,
                   help="square outer surface of this half-width; tighter than "
                        "--r-outer for a square cross-section")
    p.add_argument("--facets", type=int, default=16384)
    p.add_argument("--a0", type=float, default=-0.5, help="axial start (pad past the domain)")
    p.add_argument("--a1", type=float, default=0.5, help="axial end (pad past the domain)")
    p.add_argument("--domain-half", type=float, default=None,
                   help="cross-section half-width; checks that the pipe wall is "
                        "the nearest surface everywhere inside it")
    p.set_defaults(func=cmd_annulus)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    main()
