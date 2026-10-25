#!/usr/bin/env python3
"""Slice, measure and rotate the sailplane STL (ASCII, millimetres).

Step 1 of docs/next_session_sailplane_junction.md: the lifting-line model
needs the planform (chord, leading edge, incidence, dihedral per station),
XFOIL needs the section shapes, and the calibration box needs the geometry
rotated to the angle of attack of the run (there is no [ibm] stl_rotate key;
the rotation is done here, once, before moby_prepare).

    stl_sections.py planform  FILE.stl                 # table over the span
    stl_sections.py section   FILE.stl Y OUT.dat        # XFOIL loop at y = Y m
    stl_sections.py rotate    FILE.stl ALPHA OUT.stl   [--pivot X Z]
                                        # nose-up ALPHA deg about the y axis

Axes as in the STL: x streamwise (nose at 0), y spanwise, z up.  Everything
is converted to METRES on reading (stl_scale 0.001 in the ini), the rotated
STL is written back in millimetres.  The wing is separated from the fuselage
and fin by its x range (x < 6 m); the fin starts at x = 7.7 m.
"""
import sys
import numpy as np

WING_XMAX = 6.0  # m: everything forward of this at |y| > 0.45 is wing


def read_ascii_stl(path):
    v = []
    with open(path) as f:
        for line in f:
            if "vertex" in line:
                v.append([float(t) for t in line.split()[1:4]])
    return np.array(v).reshape(-1, 3, 3) * 1e-3


def slice_plane(tri, axis, c0):
    """Intersection segments' endpoints of all triangles with axis = c0."""
    c0 = c0 + 1e-7  # off the vertices of a mesh symmetric about y = 0
    lo, hi = tri[:, :, axis].min(1), tri[:, :, axis].max(1)
    pts = []
    for t in tri[(lo < c0) & (hi > c0)]:
        for i in range(3):
            a, b = t[i], t[(i + 1) % 3]
            if (a[axis] - c0) * (b[axis] - c0) < 0:
                s = (c0 - a[axis]) / (b[axis] - a[axis])
                pts.append(a + s * (b - a))
    return np.array(pts)


def wing_section(tri, y0):
    p = slice_plane(tri, 1, y0)
    return p[p[:, 0] < WING_XMAX] if len(p) else p


def planform(tri):
    print("# y[m]  xle[m]  chord[m]  zle[m]  incidence[deg]  t/c")
    ys = np.arange(0.5, 9.75, 0.25)
    chords = []
    for y0 in ys:
        p = wing_section(tri, y0)
        if len(p) < 10:
            chords.append(np.nan)
            continue
        x, z = p[:, 0], p[:, 2]
        ile, ite = np.argmin(x), np.argmax(x)
        c = x[ite] - x[ile]
        chords.append(c)
        inc = -np.degrees(np.arctan2(z[ite] - z[ile], x[ite] - x[ile]))
        xn = (x - x[ile]) / c
        thick = max((z[s].max() - z[s].min()) / c
                    for s in ((xn >= i / 50) & (xn < (i + 1) / 50) for i in range(50))
                    if s.sum() > 1)
        print("%5.2f  %.3f  %.3f  %.3f  %+.2f  %.3f" % (y0, x[ile], c, z[ile], inc, thick))
    chords = np.array(chords)
    m = ~np.isnan(chords)
    b = 2 * tri[:, :, 1].max()
    s_half = np.trapz(chords[m], ys[m]) + chords[m][0] * ys[m][0]  # carry-over to y = 0
    print("# span %.3f m, S %.2f m2 (root chord carried to y = 0), AR %.1f, mean chord %.3f m"
          % (b, 2 * s_half, b * b / (2 * s_half), 2 * s_half / b))
    print("# fuselage: x  |y|max[m]  z range[m]")
    for x0 in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0]:
        p = slice_plane(tri, 0, x0)
        body = p[(np.abs(p[:, 1]) < 0.45) & (p[:, 2] > -0.3) & (p[:, 2] < 1.0)]
        print("# %.1f  %.3f  [%.3f, %.3f]" % (x0, np.abs(body[:, 1]).max(), body[:, 2].min(), body[:, 2].max()))


def section(tri, y0, out):
    """Chain the slice segments into one closed loop, chord-normalised, and write
    it in the XFOIL convention (from the trailing edge over the upper surface
    to the leading edge and back along the lower surface)."""
    p = wing_section(tri, y0)
    x, z = p[:, 0], p[:, 2]
    ile = np.argmin(x)
    c = x.max() - x[ile]
    q = np.round(np.c_[(x - x[ile]) / c, (z - z[ile]) / c], 7)
    # slice_plane appends the two endpoints of each triangle's cut segment in
    # order: q[2k] -- q[2k+1] is one segment; chain them by shared endpoints
    seg = {}
    for k in range(0, len(q) - 1, 2):
        a, b = tuple(q[k]), tuple(q[k + 1])
        seg.setdefault(a, []).append(b)
        seg.setdefault(b, []).append(a)
    start = max(seg, key=lambda t: t[0])  # the trailing edge
    loop, prev, cur = [start], None, start
    while True:
        nxt = [t for t in seg[cur] if t != prev and t != cur]
        if not nxt or nxt[0] == start:
            break
        prev, cur = cur, nxt[0]
        loop.append(cur)
    pts = np.array(loop)
    if pts[1, 1] < pts[-1, 1]:  # walk the upper surface first
        pts = np.r_[pts[:1], pts[1:][::-1]]
    np.savetxt(out, pts, fmt="%.6f", header="sailplane section y = %.3f m, chord %.4f m" % (y0, c))
    print("wrote %s: %d points, chord %.4f m" % (out, len(pts), c))


def rotate(path, alpha_deg, out, pivot):
    """Nose-up rotation about the y axis through (pivot_x, pivot_z), in metres."""
    a = np.radians(alpha_deg)
    ca, sa = np.cos(a), np.sin(a)
    px, pz = pivot
    with open(path) as f, open(out, "w") as g:
        for line in f:
            w = line.split()
            if w and w[0] in ("vertex", "normal") or (len(w) == 5 and w[0] == "facet"):
                i = 1 if w[0] == "vertex" else 2
                x, y, z = (float(t) for t in w[i:i + 3])
                if w[0] == "vertex":  # metres about the pivot, back to millimetres
                    x, z = x * 1e-3 - px, z * 1e-3 - pz
                    # pitch about +y: a point ahead of the pivot (x < 0) moves up (z > 0)
                    xr, zr = ca * x + sa * z, -sa * x + ca * z
                    x, z = (xr + px) * 1e3, (zr + pz) * 1e3
                else:
                    x, z = ca * x + sa * z, -sa * x + ca * z
                g.write("%s%s %.6e %.6e %.6e\n" % (line[:line.find(w[0])], " ".join(w[:i]), x, y, z))
            else:
                g.write(line)
    print("wrote %s: rotated %+.3f deg nose-up about (x, z) = (%.3f, %.3f) m" % (out, alpha_deg, px, pz))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "planform":
        planform(read_ascii_stl(sys.argv[2]))
    elif cmd == "section":
        section(read_ascii_stl(sys.argv[2]), float(sys.argv[3]), sys.argv[4])
    elif cmd == "rotate":
        piv = (2.84, 0.47)  # root quarter chord, wing plane
        if "--pivot" in sys.argv:
            k = sys.argv.index("--pivot")
            piv = (float(sys.argv[k + 1]), float(sys.argv[k + 2]))
        rotate(sys.argv[2], float(sys.argv[3]), sys.argv[4], piv)
    else:
        sys.exit(__doc__)
