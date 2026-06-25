#!/usr/bin/env python3
"""Per-y-row fluctuation-rms profile of a mobydiff channel block-table field,
across a 2:1 wall-band interface. Reference-free: for every cell-row (level,j)
gather all (x,z) cells across the blocks sharing it, remove the (x,z) mean, and
RMS the fluctuation. The spurious interface band shows up as a localized spike in
u'/v'/w'/p' rms at exactly the coarse cell-row bordering the refined wall band.

Usage: python3 tools/channel_band_profile.py FIELD.h5 [--ymax 0.9]
"""
from __future__ import annotations
import argparse, h5py, numpy as np


def level_line(base, lev):
    """Midpoint-subdivided node line at level `lev` from the base node line."""
    line = base.copy()
    for _ in range(lev):
        mid = 0.5 * (line[:-1] + line[1:])
        new = np.empty(2 * len(line) - 1)
        new[0::2] = line
        new[1::2] = mid
        line = new
    return line


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("field")
    ap.add_argument("--ymax", type=float, default=0.95, help="print rows with y-centre < ymax")
    a = ap.parse_args()
    f = h5py.File(a.field, "r")
    nb = int(f.attrs["block_nb_x"])
    yb = f["y"][...]
    blocks = f["blocks"][...]
    D = {v: f[{"u": "un", "v": "vn", "w": "wn", "p": "pn"}[v]][...] for v in "uvwp"}
    f.close()

    # group cell values by (level, global-row j); store list of (k,i) planes
    rows = {}
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        for jj in range(nb):
            key = (int(lev), int(oy) + jj)
            rows.setdefault(key, []).append({v: D[v][bid, :, jj, :] for v in "uvwp"})

    lines = {L: level_line(yb, L) for L in set(k[0] for k in rows)}
    recs = []
    for (lev, gj), planes in rows.items():
        yc = 0.5 * (lines[lev][gj] + lines[lev][gj + 1])
        rms = {}
        for v in "uvwp":
            allcells = np.concatenate([pl[v].ravel() for pl in planes])
            rms[v] = np.sqrt(np.mean((allcells - allcells.mean()) ** 2))
        recs.append((yc, lev, gj, rms))
    recs.sort()

    print(f"{a.field}")
    print(f"  {'y':>7} {'lev':>3} {'gj':>4}   {'u_rms':>9} {'v_rms':>9} {'w_rms':>9} {'p_rms':>9}")
    prev_lev = None
    for yc, lev, gj, rms in recs:
        if yc > a.ymax:
            continue
        mark = "  <== LEVEL CHANGE (interface)" if (prev_lev is not None and lev != prev_lev) else ""
        print(f"  {yc:7.4f} {lev:3d} {gj:4d}   "
              + " ".join(f"{rms[v]:9.4e}" for v in "uvwp") + mark)
        prev_lev = lev


if __name__ == "__main__":
    main()
