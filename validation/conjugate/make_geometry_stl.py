#!/usr/bin/env python3
"""ASCII STLs for the C2 gates: a TILTED half-space and a CIRCULAR cylinder.

    ./make_geometry_stl.py plane    out.stl --theta 30 --x0 0.5 --y0 0.5 ...
    ./make_geometry_stl.py cylinder out.stl --centre 0.5 0.5 --radius 0.2 ...

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
        nx = math.cos(0.5 * (ang[i] + ang[j]))
        ny = math.sin(0.5 * (ang[i] + ang[j]))
        tris.append(((nx, ny, 0.0), (ring0[i], ring0[j], ring1[j])))
        tris.append(((nx, ny, 0.0), (ring0[i], ring1[j], ring1[i])))
    for i in range(1, m - 1):            # the two caps, fanned
        tris.append(((0.0, 0.0, -1.0), (ring0[0], ring0[i + 1], ring0[i])))
        tris.append(((0.0, 0.0, 1.0), (ring1[0], ring1[i], ring1[i + 1])))
    n = write_stl(a.out, tris, "cylinder")
    print(f"{a.out}: cylinder r = {a.radius!r} at ({cx!r}, {cy!r}), "
          f"{m} facets, {n} triangles")


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

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    main()
