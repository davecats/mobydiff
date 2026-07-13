#!/usr/bin/env python3
"""Per-level dwall/yeff cross-check for the 3-level cylinder ransgeom dump.

Independent reference: exact distance to the AS-BUILT extruded prism --
the float32 STL section polygon (T1 lesson: the quantized vertices are the
geometry). The prism is extruded past both domain z faces, so for every
domain point the nearest surface point is either on the lateral surface
(= the 2D distance to the section ring at the same z) or, for points whose
xy lies inside the ring, on a cap (= |z - z_cap|); for outside-xy points
the cap rim distance hypot(d2, dz) >= d2 never wins. All x/y domain faces
are declared patches in dwall.ini, so no domain wall is min'ed in and the
dump is the pure body distance.

Usage: check_dwall_cylinder.py <ransgeom.h5> <cylinder.stl> [--tolerance 1e-9]
"""
import argparse
import sys

import h5py
import numpy as np


def read_binary_stl_vertices(path: str) -> np.ndarray:
    with open(path, "rb") as f:
        f.seek(80)
        n = int(np.frombuffer(f.read(4), dtype="<u4")[0])
        rec = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)
    tri = rec[:, 12:48].copy().view("<f4").reshape(n, 3, 3)
    return tri.astype(np.float64)


def section_ring(tri: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Ring of unique section-xy vertices ordered by angle (the section is a
    convex 720-gon, so the angular sort around the centroid is the ring),
    plus the extrusion z span."""
    v = tri.reshape(-1, 3)
    z_lo, z_hi = float(v[:, 2].min()), float(v[:, 2].max())
    xy = np.unique(v[:, 0:2], axis=0)
    c = xy.mean(axis=0)
    order = np.argsort(np.arctan2(xy[:, 1] - c[1], xy[:, 0] - c[0]))
    return xy[order], z_lo, z_hi


def ring_distance_inside(pts_xy: np.ndarray, ring: np.ndarray,
                         chunk: int = 4096) -> tuple[np.ndarray, np.ndarray]:
    """Exact 2D distance to the closed ring polyline + point-in-polygon."""
    a = ring
    b = np.roll(ring, -1, axis=0)
    ab = b - a
    ab2 = np.einsum("ij,ij->i", ab, ab)
    d = np.empty(pts_xy.shape[0])
    inside = np.empty(pts_xy.shape[0], dtype=bool)
    for s in range(0, pts_xy.shape[0], chunk):
        p = pts_xy[s:s + chunk]
        ap = p[:, None, :] - a[None, :, :]
        t = np.clip(np.einsum("mij,ij->mi", ap, ab) / ab2, 0.0, 1.0)
        dv = ap - t[:, :, None] * ab[None, :, :]
        d[s:s + chunk] = np.sqrt(np.einsum("mij,mij->mi", dv, dv).min(axis=1))
        # crossing-number test against the same edges
        ya, yb = a[None, :, 1], b[None, :, 1]
        cond = (ya <= p[:, 1:2]) != (yb <= p[:, 1:2])
        xi = a[None, :, 0] + (p[:, 1:2] - ya) / (yb - ya + (~cond)) * ab[None, :, 0]
        inside[s:s + chunk] = (np.sum(cond & (xi > p[:, 0:1]), axis=1) % 2) == 1
    return d, inside


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("stl")
    ap.add_argument("--tolerance", type=float, default=1.0e-9)
    a = ap.parse_args()

    with h5py.File(a.dump, "r") as h5:
        blocks = h5["blocks"][...]
        dwall = h5["dwall"][...]
        yeff = h5["yeff"][...]
        xc, yc, zc = h5["xc"][...], h5["yc"][...], h5["zc"][...]

    ring, z_lo, z_hi = section_ring(read_binary_stl_vertices(a.stl))
    print(f"STL section: {ring.shape[0]} vertices, z span [{z_lo}, {z_hi}]")

    levels = blocks[:, 3]
    n_lev = int(levels.max()) + 1
    print(f"{a.dump}: {blocks.shape[0]} blocks, levels {np.bincount(levels).tolist()}")

    max_err_dwall = np.zeros(n_lev)
    max_err_yeff = np.zeros(n_lev)
    for lev in range(n_lev):
        sel = np.where(levels == lev)[0]
        if sel.size == 0:
            continue
        # all cell centres of this level, in the dump's (z, y, x) order;
        # the ring distance is cached on unique xy (cells repeat across z
        # and across blocks) -- z enters only through the cap distance.
        pts_list = []
        for b in sel:
            Zb, Yb, Xb = np.meshgrid(zc[b], yc[b], xc[b], indexing="ij")
            pts_list.append(np.column_stack([Xb.ravel(), Yb.ravel(), Zb.ravel()]))
        pts = np.concatenate(pts_list)
        xy, inv = np.unique(pts[:, 0:2], axis=0, return_inverse=True)
        d2, inside = ring_distance_inside(xy, ring)
        d2p, insp = d2[inv], inside[inv]
        dz_cap = np.minimum(np.abs(pts[:, 2] - z_lo), np.abs(pts[:, 2] - z_hi))
        ref = np.where(insp, np.minimum(d2p, dz_cap), d2p)
        # yeff floor from the block's own centre spacing (the gate grid is
        # uniform per direction, so centre spacing == cell size exactly)
        b0 = sel[0]
        h = min(float(np.diff(xc[b0]).min()), float(np.diff(yc[b0]).min()),
                float(np.diff(zc[b0]).min()))
        ref_yeff = np.maximum(ref, 0.5 * h)

        got = dwall[sel].reshape(sel.size, -1).ravel()
        goty = yeff[sel].reshape(sel.size, -1).ravel()
        max_err_dwall[lev] = np.abs(got - ref).max()
        max_err_yeff[lev] = np.abs(goty - ref_yeff).max()

    for lev in range(n_lev):
        print(f"level {lev}: max|dwall-ref| = {max_err_dwall[lev]:.3e}, "
              f"max|yeff-ref| = {max_err_yeff[lev]:.3e}")
    ok = (max_err_dwall.max() <= a.tolerance and max_err_yeff.max() <= a.tolerance
          and np.count_nonzero(np.bincount(levels)) == n_lev)
    print("dwall gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
