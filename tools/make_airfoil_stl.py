#!/usr/bin/env python3
"""Extruded cylinder / NACA 4-digit airfoil STL for the quasi-2D IBM cases.

The section polygon is extruded along z past both domain faces (the periodic
images close the geometry; same convention as the les_ibm wall slabs) and
exported as STL (float32 vertices -- part of the as-built geometry). Inside
the STL = solid, mobygeom's default.

  python3 make_airfoil_stl.py cylinder --xc 8 --yc 8 --d 1 --lz 0.5 --out cyl.stl
  python3 make_airfoil_stl.py naca --code 0012 --xc 8 --yc 8 --chord 1 \
      --lz 0.5 --out n0012.stl
  python3 make_airfoil_stl.py selig --file sd7003.dat --xc 8 --yc 8 --chord 1 \
      --lz 0.5 --out sd7003.stl

Run with the geometry venv (/home/davide/ibmc/bin/python: trimesh + shapely).
"""
import argparse

import numpy as np
import trimesh


def cylinder_section(a):
    t = np.linspace(0.0, 2.0*np.pi, a.n, endpoint=False)
    x = a.xc + 0.5*a.d*np.cos(t)
    y = a.yc + 0.5*a.d*np.sin(t)
    return np.column_stack([x, y])


def naca4_section(a):
    """NACA 4-digit closed section, chord from (xc,yc) along +x."""
    m = int(a.code[0])/100.0
    p = int(a.code[1])/10.0
    tt = int(a.code[2:])/100.0
    # cosine-clustered chordwise stations
    beta = np.linspace(0.0, np.pi, a.n//2)
    x = 0.5*(1.0 - np.cos(beta))
    yt = 5.0*tt*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2
                 + 0.2843*x**3 - 0.1036*x**4)   # closed TE coefficient
    if m > 0.0:
        yc = np.where(x < p, m/p**2*(2.0*p*x - x**2),
                      m/(1.0 - p)**2*((1.0 - 2.0*p) + 2.0*p*x - x**2))
        dyc = np.where(x < p, 2.0*m/p**2*(p - x), 2.0*m/(1.0 - p)**2*(p - x))
    else:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
    th = np.arctan(dyc)
    xu, yu = x - yt*np.sin(th), yc + yt*np.cos(th)
    xl, yl = x + yt*np.sin(th), yc - yt*np.cos(th)
    # upper TE->LE, lower LE->TE, drop duplicated LE/TE points
    xs = np.concatenate([xu[::-1], xl[1:-1]])
    ys = np.concatenate([yu[::-1], yl[1:-1]])
    return np.column_stack([a.xc + a.chord*xs, a.yc + a.chord*ys])


def selig_section(a):
    """Section from a Selig-format coordinate file (one name line, then
    normalized x y pairs TE -> upper -> LE -> lower -> TE). SD7003 and the
    rest of the UIUC database. Optional --resample fits a periodic cubic
    spline through the published points (LE clustering preserved via the
    chordwise parametrization) for IBM grids finer than the file's segments;
    default keeps the coordinates verbatim."""
    rows = []
    with open(a.file) as f:
        lines = f.read().splitlines()
    for line in lines[1:]:                       # first line is the name
        parts = line.split()
        if len(parts) != 2:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    pts = np.array(rows)
    if pts.shape[0] < 10:
        raise SystemExit(f"{a.file}: no Selig coordinate block found")
    if np.any(pts[:, 0] > 1.5):
        raise SystemExit(f"{a.file}: looks like Lednicer format (point counts "
                         "on line 2); convert to Selig ordering first")
    # drop an exactly repeated closing TE point; shapely closes the ring
    if np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]
    if a.resample:
        from scipy.interpolate import splev, splprep
        tck, _ = splprep([pts[:, 0], pts[:, 1]], s=0.0, per=True)
        u = np.linspace(0.0, 1.0, a.resample, endpoint=False)
        pts = np.column_stack(splev(u, tck))
    return np.column_stack([a.xc + a.chord*pts[:, 0], a.yc + a.chord*pts[:, 1]])


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="shape", required=True)
    for name in ("cylinder", "naca", "selig"):
        p = sub.add_parser(name)
        p.add_argument("--xc", type=float, required=True)
        p.add_argument("--yc", type=float, required=True,
                       help="section-centre coordinate in the LIFT direction "
                            "(y for --span z, z for --span y)")
        p.add_argument("--lz", type=float, required=True,
                       help="domain extent along the SPAN axis (Lz for --span z, "
                            "Ly for --span y)")
        p.add_argument("--span", choices=("z", "y"), default="z",
                       help="extrusion (span) axis: z (default) or y (the "
                            "xz-quadtree airfoil orientation: chord x, lift z)")
        p.add_argument("--pad", type=float, default=0.25, help="span overhang past both faces")
        p.add_argument("--n", type=int, default=720, help="section points")
        p.add_argument("--out", required=True)
        if name == "cylinder":
            p.add_argument("--d", type=float, required=True, help="diameter")
        elif name == "naca":
            p.add_argument("--code", default="0012", help="NACA 4-digit code")
            p.add_argument("--chord", type=float, default=1.0)
        else:
            p.add_argument("--file", required=True, help="Selig-format coordinate file")
            p.add_argument("--chord", type=float, default=1.0)
            p.add_argument("--resample", type=int, default=0,
                           help="periodic-spline resample to this many points (0 = verbatim)")
    a = ap.parse_args()

    from shapely.geometry import Polygon
    sections = {"cylinder": cylinder_section, "naca": naca4_section, "selig": selig_section}
    pts = sections[a.shape](a)
    mesh = trimesh.creation.extrude_polygon(Polygon(pts), height=a.lz + 2.0*a.pad)
    mesh.apply_translation([0.0, 0.0, -a.pad])
    if a.span == "y":
        # Section plane becomes x-z (chord x, LIFT z), span along y: swap
        # the y and z axes. The swap mirrors the mesh, so re-invert the
        # face winding to keep inside = solid.
        mesh.vertices = mesh.vertices[:, [0, 2, 1]]
        mesh.invert()
    mesh.export(a.out)
    print(f"{a.out}: {len(mesh.vertices)} vertices, watertight={mesh.is_watertight}, "
          f"span {a.span} in [{-a.pad}, {a.lz + a.pad}]")


if __name__ == "__main__":
    main()
