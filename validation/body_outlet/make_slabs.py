#!/usr/bin/env python3
"""Two plane wall slabs as ONE ASCII STL, for the body-at-outlet gates
(docs/next_session_body_at_outlet.md): the solid below `lo` and above `hi`
along the wall-normal axis, both slabs padded beyond the domain in every
direction so that the body crosses the inlet and the outlet planes (and
their ghost rows: an STL that ends at the plane would leave the tangential
ghosts fluid -- B0, question 6).

    ./make_slabs.py out.stl --normal y --lo 0.259375 --hi 1.259375 \
        --box 4 1.5 0.125 [--pad 2] [--x0 -2 --x1 6] ...

`--box` is the domain (lx ly lz); the slabs span [-pad, L + pad] in the two
tangential directions and `pad` deep into the solid. `--lo-x0/--lo-x1`
(and `--hi-x0/--hi-x1`) limit the LOWER (UPPER) slab's streamwise extent for
the body-ends-before-the-plane / starts-at-the-plane limits; the default is
the full padded span. ASCII for the reason tools/make_geometry_stl.py
gives: moby_prepare parses ASCII vertices to float64 directly.
"""

import argparse
import sys

sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), "..", "..", "tools"))
from make_geometry_stl import write_stl, box_faces  # noqa: E402

AX = {"x": 0, "y": 1, "z": 2}


def box_corners(lo, hi):
    """The 8 corners of the axis-aligned box [lo, hi] in box_faces' order."""
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    return [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
            (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("out")
    p.add_argument("--normal", choices=("x", "y", "z"), default="y")
    p.add_argument("--stream", choices=("x", "y", "z"), default="x",
                   help="the streamwise axis (for --lo-x0/--lo-x1 limits)")
    p.add_argument("--lo", type=float, required=True, help="top of the lower slab")
    p.add_argument("--hi", type=float, required=True, help="bottom of the upper slab")
    p.add_argument("--box", type=float, nargs=3, required=True, help="lx ly lz")
    p.add_argument("--pad", type=float, default=2.0)
    p.add_argument("--lo-x0", type=float, default=None)
    p.add_argument("--lo-x1", type=float, default=None)
    p.add_argument("--hi-x0", type=float, default=None)
    p.add_argument("--hi-x1", type=float, default=None)
    a = p.parse_args()

    n = AX[a.normal]
    s = AX[a.stream]
    tris = []
    for which, x0lim, x1lim in (("lo", a.lo_x0, a.lo_x1), ("hi", a.hi_x0, a.hi_x1)):
        lo = [-a.pad] * 3
        hi = [a.box[d] + a.pad for d in range(3)]
        if which == "lo":
            hi[n] = a.lo
            lo[n] = a.lo - a.pad
        else:
            lo[n] = a.hi
            hi[n] = a.hi + a.pad
        if x0lim is not None:
            lo[s] = x0lim
        if x1lim is not None:
            hi[s] = x1lim
        tris.extend(box_faces(box_corners(lo, hi)))
    ntri = write_stl(a.out, tris, "slabs")
    print(f"{a.out}: two slabs normal to {a.normal} below {a.lo!r} and above {a.hi!r}, {ntri} triangles")


if __name__ == "__main__":
    main()
