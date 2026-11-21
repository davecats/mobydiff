#!/usr/bin/env python3
"""Tile a boundary-layer restart field N times in the periodic z direction.

    tile_restart.py SRC.h5 DST.h5 [--ntile 2] [--noise 1e-3] [--seed 1]

SRC is a block-table restart (4D datasets (nblocks, nz, ny, nx) + `blocks`, the
origin (i, j, k, level) of every block) or a legacy global-3D one. DST is
legacy global-3D (no `blocks` dataset): the reader takes the global path, which
works on any rank count and block layout (make_global_restart.py).

The copy is exact: z is periodic, so un, vn, wn, pn repeated ntile times are
the same divergence-free field on the wider box (a z-face of wn at index k is
the same face in every copy). The attributes nz, lz and the z node line are
rewritten; every other attribute is copied (the solver checks only the grid
keys against the ini, docs/configuration.md "[restart]"). The outlet planes
(un_xmax, vn_ymax: solver state of the run that wrote SRC) are NOT carried:
an IC tool writes the plane or deletes it (CLAUDE.md), and the run predicts
them again from its first step.

--noise adds a uniform random perturbation of that amplitude (U_inf units) to
un, vn, wn of every copy but the first, away from the inlet, the outlet and
the wall rows (x index 100 .. nx-100, y index 1 .. 0.8 ny), so the copies do
not stay identical: the dynamics keep an exact Lz/ntile symmetry otherwise
until the trip-born flow has replaced the layer. The first projection removes
the perturbation's divergence.
"""
import argparse

import h5py
import numpy as np

VARS = ("un", "vn", "wn", "pn")


def read_global(src, name):
    """The global (nz, ny, nx) array of one variable."""
    d = src[name]
    if d.ndim == 3:
        return d[...]
    nz, ny, nx = (int(src.attrs[k]) for k in ("nz", "ny", "nx"))
    out = np.empty((nz, ny, nx), dtype=d.dtype)
    bnz, bny, bnx = d.shape[1:]
    for b, (i0, j0, k0, lev) in enumerate(src["blocks"][...]):
        assert lev == 0, "level-0 (uniform) block tables only"
        out[k0:k0 + bnz, j0:j0 + bny, i0:i0 + bnx] = d[b]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--ntile", type=int, default=2)
    ap.add_argument("--noise", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    with h5py.File(a.src, "r") as src, h5py.File(a.dst, "w") as dst:
        nz, ny, nx = (int(src.attrs[k]) for k in ("nz", "ny", "nx"))
        lz = float(src.attrs["lz"])
        for k, v in src.attrs.items():
            dst.attrs[k] = v
        dst.attrs["nz"] = nz * a.ntile
        dst.attrs["lz"] = lz * a.ntile
        dst.create_dataset("x", data=src["x"][...])
        dst.create_dataset("y", data=src["y"][...])
        z = src["z"][...]
        assert len(z) == nz + 1 and abs(z[-1] - z[0] - lz) < 1e-9 * lz, "z is not the node line"
        dst.create_dataset("z", data=np.concatenate([z] + [z[1:] + c * lz for c in range(1, a.ntile)]))
        i0, i1, j1 = 100, nx - 100, int(0.8 * ny)
        for name in VARS:
            g = read_global(src, name)
            t = np.concatenate([g] * a.ntile, axis=0)
            del g
            if a.noise > 0 and name != "pn":
                for c in range(1, a.ntile):
                    sl = t[c * nz:(c + 1) * nz, 1:j1, i0:i1]
                    sl += a.noise * (2.0 * rng.random(sl.shape) - 1.0)
            dst.create_dataset(name, data=t)
            print(f"  {name}: ({nz}, {ny}, {nx}) -> {t.shape}", flush=True)
            del t
    print(f"wrote {a.dst}: legacy global-3D, nz {nz} -> {nz * a.ntile}, lz {lz} -> {lz * a.ntile}")


if __name__ == "__main__":
    main()
