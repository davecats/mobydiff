#!/usr/bin/env python3
"""Build the finer wall-normal (blayer) node line for tbl_finey and save it.

    make_finey_grid.py            # -> new_y_256.npy

Ports src/modules/init.f90 build_blayer_line / natural_wall_coordinate EXACTLY
(validated bit-for-bit against the solver's ny=176 line, max|err| 1e-13), so the
line here equals what the solver rebuilds from [grid] ny + [case.boundaryLayer]
dyw_plus / jb / resolved_height on restart. Design: keep the WALL spacing fixed
(same physical dy_wall = same dy+ ~ 0.23 -> the diffusive/Peclet dt limit, set by
the smallest cell, is unchanged so dt stays 0.02) and raise ny so the largest
in-BL spacing drops from dy+ ~ 6 to ~ 4 (isotropic with dx+ ~ dz+ ~ 4). Matching
dy_wall requires raising dyw_plus 0.15 -> 0.270 as ny 176 -> 256 (the knob is a
clustering strength, not the wall spacing itself).
"""
import numpy as np

LY, H, JB, RMAX = 100.0, 36.0, 45.0, 1.2      # ly, resolved_height, blend index, max outer ratio
NY, DYW = 224, 0.221                           # target: dy_wall matches the ny=176/dyw=0.15 line
#   (ny=224 -> dy+_max within delta99 ~4.3 typ / 4.9 worst, 176 M cells fits the 5090's 32 GB;
#    ny=256/dyw=0.270 gives ~3.7 but 201 M ~= 32 GB would OOM a single 5090.)


def ypsi(j, jb, dyw):
    j = np.asarray(j, float)
    blend = (j / jb) ** 2
    outer = (0.75 * j) ** (4.0 / 3.0)          # 0.75 = alpha*c_eta*0.75 = 1.25*0.8*0.75
    return np.where(j <= 0, 0.0, (dyw * j + outer * blend) / (1.0 + blend))


def solve_geom_ratio(dy, m, span):
    if m < 1 or dy <= 0 or span <= 0 or dy * m >= span:
        return 1.0, False
    hi = 3.0
    if dy * (hi ** m - 1) / (hi - 1) < span:
        return 1.0, False
    lo = 1.0 + 1e-9
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        hi, lo = (mid, lo) if dy * (mid ** m - 1) / (mid - 1) > span else (hi, mid)
    return 0.5 * (lo + hi), True


def build_blayer_line(n, dyw, length=LY, jb=JB, h=H):
    for n_in in range(n - 2, max(8, n // 5) - 1, -1):
        n_out = n - n_in
        yb = ypsi(np.arange(n_in + 1), jb, dyw)
        dy_edge = h - h * yb[n_in - 1] / yb[n_in]
        if dy_edge <= 0:
            continue
        r, ok = solve_geom_ratio(dy_edge, n_out, length - h)
        if ok and r <= RMAX:
            node = np.zeros(n + 1)
            node[:n_in + 1] = h * yb / yb[n_in]
            node[0], node[n_in] = 0.0, h
            for i in range(1, n_out + 1):
                node[n_in + i] = h + dy_edge * (r ** i - 1) / (r - 1)
            node[n] = length
            return node, n_in, r
    raise RuntimeError("no valid blayer split")


if __name__ == "__main__":
    node, n_in, r = build_blayer_line(NY, DYW)
    np.save("new_y.npy", node)
    dy = np.diff(node)
    print(f"ny={NY} dyw={DYW}: n_in={n_in} r={r:.4f}")
    print(f"  dy_wall = {dy[0]:.6e}  (dy+ ~ {dy[0]*450*0.05:.3f})")
    print(f"  dy at resolved edge (y={H}) = {dy[n_in-1]:.4f}")
    print(f"  wrote new_y.npy ({len(node)} nodes, y_max={node[-1]:.1f})")
