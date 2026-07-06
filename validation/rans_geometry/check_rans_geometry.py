#!/usr/bin/env python3
"""Gate checker for the RANS T1 geometry state (docs/next_session_iddes.md).

Reads a <prefix>_ransgeom.h5 dump (rans.f90: blocks table, interior
dwall/yeff/wallcell in the block-table layout, per-block cell-centre
coordinates) and compares against an independent reference:

  --mode flat   the les_ibm plane-wall channel: exact point-to-box surface
                distance to the two wall slabs (make_walls_stl.py geometry),
                min'ed with the y domain walls; wallcell predicted exactly
                from the slab planes.
  --mode wavy   the analytic wavy wall: high-precision minimization of the
                distance to the sine curve (scipy), min'ed with the y domain
                walls; wallcell sanity checks only (curved surface).

y_eff reference = max(dwall_ref, half the smallest local spacing at the
block's level). Uniform grids only (all T1 gate cases are uniform).
"""
import argparse
import sys

import h5py
import numpy as np
from scipy.optimize import minimize_scalar

FLAT = dict(
    lx=12.566370614359172, ly=2.5, lz=6.283185307179586,
    nx=64, ny=80, nz=64,
    pad=0.25, y_lo=8.3*2.5/80, y_hi=72.3*2.5/80,
)
WAVY = dict(
    lx=1.0, ly=1.0, lz=0.25,
    nx=64, ny=64, nz=16,
    amp=2.5e-2, n_wave=1, offset=1.0e-2,
)


def box_surface_distance(pts: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Unsigned distance from points (N,3) to the surface of an axis box."""
    d = np.maximum(lo - pts, pts - hi)          # per-axis signed exterior distance
    outside = np.linalg.norm(np.maximum(d, 0.0), axis=1)
    inside = -np.max(d, axis=1)                 # > 0 strictly inside
    return np.where(outside > 0.0, outside, inside)


def flat_body_distance(pts: np.ndarray) -> np.ndarray:
    # STL stores float32 vertices, so the as-built wall planes are the
    # float32-rounded box bounds (the ~1e-7 quantization is part of the
    # committed geometry, not a solver error).
    f32 = lambda v: np.float64(np.float32(v))
    p = FLAT
    lo_box = (f32(np.array([-p["pad"], -p["pad"], -p["pad"]])),
              f32(np.array([p["lx"] + p["pad"], p["y_lo"], p["lz"] + p["pad"]])))
    hi_box = (f32(np.array([-p["pad"], p["y_hi"], -p["pad"]])),
              f32(np.array([p["lx"] + p["pad"], p["ly"] + p["pad"], p["lz"] + p["pad"]])))
    return np.minimum(box_surface_distance(pts, *lo_box),
                      box_surface_distance(pts, *hi_box))


def wavy_height(s: np.ndarray | float) -> np.ndarray | float:
    p = WAVY
    return p["amp"]*0.5*(1.0 + np.sin(2.0*np.pi*p["n_wave"]*s/p["lx"])) + p["offset"]


def wavy_curve_distance(x: float, y: float) -> float:
    """Distance to the sine curve y = wavy_height(s), scanning one period
    centred on x then refining with a bounded scalar minimization."""
    period = WAVY["lx"]/WAVY["n_wave"]
    s = np.linspace(x - 0.5*period, x + 0.5*period, 2049)
    d2 = (x - s)**2 + (y - wavy_height(s))**2
    k = int(np.argmin(d2))
    a, b = s[max(k - 1, 0)], s[min(k + 1, s.size - 1)]
    res = minimize_scalar(lambda t: (x - t)**2 + (y - wavy_height(t))**2,
                          bounds=(a, b), method="bounded",
                          options={"xatol": 1.0e-13})
    return float(np.sqrt(min(res.fun, d2[k])))


def load_dump(path: str):
    with h5py.File(path, "r") as h5:
        nbx = int(h5.attrs["block_nb_x"])
        blocks = h5["blocks"][...]
        dwall = h5["dwall"][...]      # (n, nbz, nby, nbx)
        yeff = h5["yeff"][...]
        wallcell = h5["wallcell"][...]
        xc, yc, zc = h5["xc"][...], h5["yc"][...], h5["zc"][...]
    return nbx, blocks, dwall, yeff, wallcell, xc, yc, zc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--mode", choices=("flat", "wavy"), required=True)
    ap.add_argument("--tolerance", type=float, default=1.0e-9)
    args = ap.parse_args()

    p = FLAT if args.mode == "flat" else WAVY
    nb, blocks, dwall, yeff, wallcell, xc, yc, zc = load_dump(args.dump)
    n_blocks = blocks.shape[0]
    levels = blocks[:, 3]
    print(f"{args.dump}: {n_blocks} blocks, levels {np.bincount(levels)}")

    # Per-point reference on each block's own cell centres (dump order z,y,x).
    max_err_dwall = np.zeros(int(levels.max()) + 1)
    max_err_yeff = np.zeros_like(max_err_dwall)
    bad_wallcell = 0

    # The wavy minimization is per unique (x, y); cache it.
    wavy_cache: dict[tuple[float, float], float] = {}

    for b in range(n_blocks):
        lev = int(levels[b])
        Z, Y, X = np.meshgrid(zc[b], yc[b], xc[b], indexing="ij")
        pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

        if args.mode == "flat":
            body = flat_body_distance(pts)
        else:
            body = np.empty(pts.shape[0])
            for i, (x, y) in enumerate(zip(pts[:, 0], pts[:, 1])):
                key = (round(x, 14), round(y, 14))
                if key not in wavy_cache:
                    wavy_cache[key] = wavy_curve_distance(x, y)
                body[i] = wavy_cache[key]

        # y domain walls (both gate cases: y non-periodic, no-slip).
        dom = np.minimum(np.abs(pts[:, 1]), np.abs(p["ly"] - pts[:, 1]))
        ref_dwall = np.minimum(body, dom)
        spacing = np.array([p["lx"]/p["nx"], p["ly"]/p["ny"], p["lz"]/p["nz"]])/2**lev
        ref_yeff = np.maximum(ref_dwall, 0.5*spacing.min())

        shape = dwall[b].shape
        max_err_dwall[lev] = max(max_err_dwall[lev],
                                 np.abs(dwall[b] - ref_dwall.reshape(shape)).max())
        max_err_yeff[lev] = max(max_err_yeff[lev],
                                np.abs(yeff[b] - ref_yeff.reshape(shape)).max())

        if args.mode == "flat":
            # Exact wall-cell prediction: u/w faces sit at the cell-centre y
            # (4 faces), v faces at the two y nodes; a face is solid iff its
            # y is outside the fluid gap.
            h = p["ly"]/p["ny"]/2**lev
            solid = lambda y: (y < p["y_lo"]) | (y > p["y_hi"])
            n_solid = (4*solid(Y).astype(int) + solid(Y - 0.5*h).astype(int)
                       + solid(Y + 0.5*h).astype(int))
            ref_marker = np.where(n_solid == 0, 0, np.where(n_solid == 6, 2, 1))
            bad_wallcell += int((wallcell[b] != ref_marker).sum())

    n_wall = int((wallcell == 1).sum())
    n_solid = int((wallcell == 2).sum())
    print(f"wall cells: {n_wall}, solid cells: {n_solid}")
    for lev in range(max_err_dwall.size):
        print(f"level {lev}: max|dwall-ref| = {max_err_dwall[lev]:.3e}, "
              f"max|yeff-ref| = {max_err_yeff[lev]:.3e}")

    ok = max_err_dwall.max() <= args.tolerance and max_err_yeff.max() <= args.tolerance
    if args.mode == "flat":
        print(f"wallcell mismatches: {bad_wallcell}")
        ok = ok and bad_wallcell == 0
    else:
        # Curved surface: sanity only — a wall band exists and cells well
        # above the surface (max height 0.06) stay fluid.
        high_nonfluid = 0
        for b in range(n_blocks):
            mask = yc[b] > 0.1
            high_nonfluid += int((wallcell[b][:, mask, :] != 0).sum())
        ok = ok and n_wall > 0 and high_nonfluid == 0
        print(f"wall band present: {n_wall > 0}; non-fluid cells above y=0.1: {high_nonfluid}")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
