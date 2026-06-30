#!/usr/bin/env python3
"""Gates 1-2 for the LES<->IBM coupling: measure the time-mean nut field at the
IBM wall from the stats-leg snapshots.

  Gate 1 (HARD): every cell flagged SOLID by the IBM mask (any of its 6 staggered
                 faces has |coef| > 1e20) has nut == 0 in EVERY snapshot.
  Gate 2:        the time-mean nut(y) profile -> 0 INTO the wall and shows no
                 spurious spike at the band cells (the first fluid cells adjacent
                 to the solid). The physical behaviour is nut->0 at the wall.

The mask uses the SAME rule as src/modules/les.f90 (solid_cell = any of coef(VAR_U,i),
coef(VAR_U,i+1), coef(VAR_V,j), coef(VAR_V,j+1), coef(VAR_W,k), coef(VAR_W,k+1)
exceeds 1e20). For the single-level run the coefficients are the global `coef`
dataset of ibm_coeff.h5; nut snapshots are the block-table layout.

Usage:  python3 measure_nut.py [--run runs/a_wale/stats] [--coef ibm_coeff.h5]
"""
from __future__ import annotations
import argparse
import glob
import os

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SOLID_TH = 1e20
RETAU = 180.0
NB = 8


def load_coef_solid_mask(coef_path, ny):
    """Per global cell (x,y,z): True if the les.f90 ibm_aware rule marks it solid.
    coef dataset shape (nx+2, ny+2, nz+2, 3) with a one-cell ghost layer."""
    with h5py.File(coef_path, "r") as f:
        c = f["coef"][...]                       # (nx+2, ny+2, nz+2, 3) ghost-incl
    sol = np.abs(c) > SOLID_TH                    # (.,.,.,var)
    nx, nyt, nz = c.shape[0] - 2, c.shape[1] - 2, c.shape[2] - 2
    # interior cells 1..n ; faces: u(i),u(i+1); v(j),v(j+1); w(k),w(k+1)
    u, v, w = sol[..., 0], sol[..., 1], sol[..., 2]
    i = slice(1, nx + 1)
    j = slice(1, nyt + 1)
    k = slice(1, nz + 1)
    ip = slice(2, nx + 2)
    jp = slice(2, nyt + 2)
    kp = slice(2, nz + 2)
    solid = (u[i, j, k] | u[ip, j, k] | v[i, j, k] | v[i, jp, k] |
             w[i, j, k] | w[i, j, kp])           # (nx, ny, nz)
    return solid                                  # axes (x, y, z)


def reassemble_nut(snap, ny_base):
    """time-mean nut(x,y,z) reassembled to the base lattice (single level)."""
    with h5py.File(snap, "r") as f:
        bl = f["blocks"][...]
        n = f["nut"][...]                         # (n_blocks, NBz, NBy, NBx)? see below
    nb = n.shape[1]
    nx = (bl[:, 0].max() // nb + 1) * nb
    nz = (bl[:, 2].max() // nb + 1) * nb
    out = np.zeros((nx, ny_base, nz))
    for bid, (ox, oy, oz, lev) in enumerate(bl):
        # block dataset axes mirror the writer: (z, y, x) within the block
        out[ox:ox + nb, oy:oy + nb, oz:oz + nb] = np.transpose(n[bid], (2, 1, 0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=os.path.join(HERE, "runs", "a_wale", "stats"),
                    help="stats-leg dir with channel_ibm_*.h5 snapshots")
    ap.add_argument("--coef", default=os.path.join(HERE, "ibm_coeff.h5"))
    ap.add_argument("--prefix", default="channel_ibm")
    a = ap.parse_args()

    snaps = sorted(glob.glob(os.path.join(a.run, f"{a.prefix}_*.h5")),
                   key=lambda p: int("".join(filter(str.isdigit, os.path.basename(p)))))
    snaps = [s for s in snaps if "stats" not in os.path.basename(s)]
    if not snaps:
        raise SystemExit(f"no snapshots in {a.run}")
    with h5py.File(snaps[0], "r") as f:
        yt = f["y"][...]
    ny = len(yt) - 1
    ycen = 0.5 * (yt[:-1] + yt[1:])

    solid = load_coef_solid_mask(a.coef, ny)     # (x, y, z)

    # Gate 1: nut == 0 in every solid cell, every snapshot
    worst = 0.0
    nut_sum = np.zeros((solid.shape[0], ny, solid.shape[2]))
    for s in snaps:
        nut = reassemble_nut(s, ny)
        worst = max(worst, float(np.abs(nut[solid]).max()))
        nut_sum += nut
    nut_mean = nut_sum / len(snaps)
    print(f"== {len(snaps)} snapshots from {os.path.relpath(a.run, HERE)}")
    print(f"GATE 1 (solid nut==0): max|nut| over all solid cells = {worst:.3e}  "
          f"{'PASS' if worst == 0.0 else 'FAIL'}")

    # Gate 2: time+plane-mean nut(y); flag any band-cell spike
    solid_y = solid.all(axis=(0, 2))             # whole-plane solid rows
    nut_y = nut_mean.mean(axis=(0, 2)) / (1.0 / RETAU)   # nut/nu_mol
    print(f"GATE 2 (no spurious wall spike): nut/nu_mol(y), wall band marked [B]")
    # band cells = fluid rows adjacent to a solid row
    fluid = ~solid_y
    adj = np.zeros(ny, bool)
    adj[1:] |= solid_y[:-1] & fluid[1:]
    adj[:-1] |= solid_y[1:] & fluid[:-1]
    core = fluid & ~adj
    core_max = nut_y[core].max() if core.any() else float("nan")
    for j in range(ny):
        if solid_y[j] and not (j > 0 and fluid[j - 1]) and not (j < ny - 1 and fluid[j + 1]):
            continue                              # skip deep solid rows
        tag = "[B]" if adj[j] else ("[S]" if solid_y[j] else "")
        if abs(ycen[j] - 0.26) < 0.4 or abs(ycen[j] - 2.26) < 0.4:
            print(f"   y={ycen[j]:6.4f} y+_wall~{abs(ycen[j]-0.259375)*RETAU:6.2f} "
                  f"nut/nu={nut_y[j]:8.4f} {tag}")
    band_max = nut_y[adj].max() if adj.any() else float("nan")
    print(f"   band-cell max nut/nu = {band_max:.4f};  core max nut/nu = {core_max:.4f};  "
          f"ratio band/core = {band_max / core_max:.2f}")
    print("   (a ratio >> 1 = spurious band spike -> band-aware damping needed; "
          "~<1 = physical nut->0 into the wall)")


if __name__ == "__main__":
    main()
