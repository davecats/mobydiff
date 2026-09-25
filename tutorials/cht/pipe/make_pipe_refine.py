#!/usr/bin/env python3
"""Interpolate a converged coarse-grid field onto the fine grid.

WHY: the fine grid's thermal legs cost about 9x the coarse ones per unit time
-- 1.6x the cells per direction AND a Peclet-limited dt ~ h^2 -- so repeating
the campaign's transients there would take days for fields we have already
converged once. The physics does not change with the grid; only the
discretisation error does, which is the whole point of the two-grid study. So
the fine run starts from the coarse ANSWER and only has to re-adjust to its
own truncation error, which happens on the fluid's own time scale, not on the
shell's diffusion time.

Trilinear interpolation, each velocity component taken at ITS OWN staggered
location (`q(i)` is the LOW face of cell i) and each scalar at the cell
centre. The result is NOT discretely divergence-free -- the first projection
removes that, and the correction is small because the coarse field was. The
scalars are interpolated straight through the immersed interface, which
smears it over one coarse cell; that too is an initial condition, not an
answer, and the conjugate coefficient re-establishes the jump within a few
steps.

The TEMPLATE is a fine-grid snapshot (run the solver for one step on the fine
case) -- it supplies the block table, the node lines and every attribute, so
this script only rewrites dataset contents.

    ./make_pipe_refine.py p_coarse2_stat_635000.h5 fine_template.h5 IC_fine.h5
"""

from __future__ import annotations

import argparse
import shutil
import sys

import h5py
import numpy as np


def block_geometry(h5):
    blocks = h5["blocks"][...]
    nb = (int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]),
          int(h5.attrs["block_nb_z"]))
    if int(blocks[:, 3].max()) != 0:
        raise SystemExit("make_pipe_refine: single-level files only")
    shape = (int(h5.attrs["nz"]), int(h5.attrs["ny"]), int(h5.attrs["nx"]))
    return blocks, nb, shape


def load_global(h5, name):
    blocks, nb, shape = block_geometry(h5)
    arr = np.empty(shape)
    d = h5[name]
    for bid, (ox, oy, oz, _l) in enumerate(blocks):
        arr[oz:oz + nb[2], oy:oy + nb[1], ox:ox + nb[0]] = d[bid]
    return arr


def axis_weights(src, dst, periodic, length):
    """Index pair and weight for interpolating from the `src` sample line to
    the `dst` positions. Outside the source range the value is clamped, except
    on a periodic axis where it wraps -- the jacket, where clamping applies,
    holds a constant field anyway."""
    n = len(src)
    h = src[1] - src[0]
    if periodic:
        f = (dst - src[0]) / h
        i0 = np.floor(f).astype(int)
        w = f - i0
        return i0 % n, (i0 + 1) % n, w
    f = np.clip((dst - src[0]) / h, 0.0, n - 1.0000000001)
    i0 = np.floor(f).astype(int)
    return i0, np.minimum(i0 + 1, n - 1), f - i0


def interp3(arr, wz, wy, wx):
    """Separable trilinear interpolation; each w is (i0, i1, w)."""
    (kz0, kz1, tz), (ky0, ky1, ty), (kx0, kx1, tx) = wz, wy, wx
    out = np.zeros((len(kz0), len(ky0), len(kx0)))
    for kz, cz in ((kz0, 1 - tz), (kz1, tz)):
        for ky, cy in ((ky0, 1 - ty), (ky1, ty)):
            for kx, cx in ((kx0, 1 - tx), (kx1, tx)):
                out += (cz[:, None, None] * cy[None, :, None]
                        * cx[None, None, :]) * arr[np.ix_(kz, ky, kx)]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", help="converged coarse snapshot")
    ap.add_argument("template", help="fine-grid snapshot (one step) to reuse")
    ap.add_argument("out")
    a = ap.parse_args()

    shutil.copyfile(a.template, a.out)
    with h5py.File(a.source, "r") as src, h5py.File(a.out, "r+") as dst:
        snode = {d: src[d][...] for d in "xyz"}
        scent = {d: 0.5 * (snode[d][:-1] + snode[d][1:]) for d in "xyz"}
        dnode = {d: dst[d][...] for d in "xyz"}
        dcent = {d: 0.5 * (dnode[d][:-1] + dnode[d][1:]) for d in "xyz"}
        blocks, nb, _ = block_geometry(dst)
        names = [k for k, v in src.items()
                 if isinstance(v, h5py.Dataset) and v.ndim == 4 and k != "blocks"
                 and k in dst]
        # staggered axis per variable: the LOW face in its own direction
        stag = {"un": (0,), "vn": (1,), "wn": (2,)}
        print(f"   {len(names)} datasets: {' '.join(names)}")
        for nm in names:
            if nm == "vfrac":                    # geometry, already correct
                continue
            field = load_global(src, nm)
            comp = stag.get(nm, ())
            sx = snode["x"][:-1] if 0 in comp else scent["x"]
            sy = snode["y"][:-1] if 1 in comp else scent["y"]
            sz = snode["z"][:-1] if 2 in comp else scent["z"]
            out = dst[nm]
            for bid, (ox, oy, oz, _l) in enumerate(blocks):
                dx = (dnode["x"][ox:ox + nb[0]] if 0 in comp
                      else dcent["x"][ox:ox + nb[0]])
                dy = (dnode["y"][oy:oy + nb[1]] if 1 in comp
                      else dcent["y"][oy:oy + nb[1]])
                dz = (dnode["z"][oz:oz + nb[2]] if 2 in comp
                      else dcent["z"][oz:oz + nb[2]])
                wx = axis_weights(sx, dx, False, 0.0)
                wy = axis_weights(sy, dy, False, 0.0)
                wz = axis_weights(sz, dz, True, 0.0)
                out[bid] = interp3(field, wz, wy, wx)
            print(f"   {nm}: interpolated", flush=True)
        with h5py.File(a.source, "r") as s2:
            dst.attrs["t_current"] = 0.0
            dst.attrs["step"] = 0
    print(f"   wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
