#!/usr/bin/env python3
"""Axis-aligned box as an ASCII STL -- the conjugate gates' flat solid slab.

ASCII is deliberate. `moby_prepare`'s STL reader parses ASCII vertices
straight to float64 (P1b), while a binary STL carries float32 vertices, so an
ASCII slab puts the interface at EXACTLY the requested plane and the analytic
reference needs no quantisation dance (the T1 lesson in CLAUDE.md, where the
float32 rounding of the les_ibm wall planes was part of the as-built geometry).

The slab is padded well beyond the domain in the periodic directions on
purpose: geometry_stl.f90 images a periodic dimension only when the mesh is
NARROWER than a cell, so a full-span padded slab is its own periodic
continuation and its side skin is never seen as a wall. Only the top face
y = y_wall is a real interface.

    ./make_slab_stl.py out.stl --y-top 0.234 --pad 4.0
"""

from __future__ import annotations

import argparse


def box_triangles(lo, hi):
    """The 12 triangles of an axis-aligned box, outward normals."""
    (x0, y0, z0), (x1, y1, z1) = lo, hi
    v = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
         (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    quads = [((0, 3, 2, 1), (0.0, 0.0, -1.0)),      # z = z0
             ((4, 5, 6, 7), (0.0, 0.0, 1.0)),       # z = z1
             ((0, 1, 5, 4), (0.0, -1.0, 0.0)),      # y = y0
             ((3, 7, 6, 2), (0.0, 1.0, 0.0)),       # y = y1  (the interface)
             ((0, 4, 7, 3), (-1.0, 0.0, 0.0)),      # x = x0
             ((1, 2, 6, 5), (1.0, 0.0, 0.0))]       # x = x1
    for (a, b, c, d), n in quads:
        yield n, (v[a], v[b], v[c])
        yield n, (v[a], v[c], v[d])


def write_stl(path, lo, hi, name="slab"):
    with open(path, "w") as fh:
        fh.write(f"solid {name}\n")
        for n, tri in box_triangles(lo, hi):
            fh.write("  facet normal %.17g %.17g %.17g\n" % n)
            fh.write("    outer loop\n")
            for p in tri:
                fh.write("      vertex %.17g %.17g %.17g\n" % p)
            fh.write("    endloop\n")
            fh.write("  endfacet\n")
        fh.write(f"endsolid {name}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--y-top", type=float, required=True,
                    help="the interface plane; the slab fills everything below it")
    # KEEP THE BOX SMALL. The BVH point-triangle distance forms d^2 as a
    # difference of terms of the vertices' magnitude, so a slab padded to
    # +-4 loses ~64x more in that cancellation than one padded to +-0.5 --
    # measured: |dwall - exact| 3.2e-13 vs ~1e-14 at a cut 0.003 from the
    # plane. It is a geometric-precision floor, not a scheme error, but it
    # sets the tolerance of every gate downstream, so do not inflate it.
    ap.add_argument("--y-bottom", type=float, default=-1.0,
                    help="deep enough that the bottom face is never the nearest")
    ap.add_argument("--pad", type=float, default=0.5,
                    help="extent beyond the domain in x and z")
    a = ap.parse_args()

    y0 = a.y_bottom
    write_stl(a.out, (-a.pad, y0, -a.pad), (a.pad, a.y_top, a.pad))
    print(f"{a.out}: slab y in [{y0!r}, {a.y_top!r}], padded +-{a.pad} in x,z")


if __name__ == "__main__":
    main()
