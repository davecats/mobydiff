#!/usr/bin/env python3
"""Edge/corner interface banding metric for an embedded refined patch.

With a single refined patch floating in the coarse interior there is no
homogeneous direction to average for fluctuations, and the 2:1 interface band is
LOW-wavenumber (a reference-free roughness/Laplacian proxy misses it -- see the
interface-validation-suite memory). So this metric is REFERENCE-BASED: the
patch run's coarse cells OUTSIDE the patch start bit-identical to a base-grid
control (the patch IC's level-0 leaves are the same base interpolation as
`make_channel_restart --mode base`). After the same number of steps the
cell-by-cell difference, classified by each coarse leaf's adjacency to the
patch, isolates the interface artifact:

  interior : coarse leaf with no refined neighbour (background difference level,
             from the global pressure coupling; should be small/decay with
             distance from the patch)
  face     : a 6-axis face-neighbour is refined  (the patch's 6 faces)
  edge     : an edge-neighbour is refined         (12 edges)
  corner   : only a corner-neighbour is refined   (8 corners)

GATE: edge/corner classes not worse than the face class, and the face class not
much above the interior background (no band/spurious energy at the edges/corners
beyond the known flat-face residual). Precedence face>edge>corner: the strongest
contact a leaf has defines its class.

NOTE: statistically meaningful only on developed (time-averaged) fields. On an
early snapshot it is a first, caveated reading -- the flow has barely advected,
so the coarse region away from the patch is still close to the control and the
difference is dominated by the patch's local interface footprint.

Usage:
  python3 tools/patch_interface_diff.py PATCH.h5 BASE.h5 [--periodic-y]
"""
from __future__ import annotations
import argparse
import h5py
import numpy as np

VARS = {"u": "un", "v": "vn", "w": "wn", "p": "pn"}


def read_blocks(path):
    f = h5py.File(path, "r")
    nb = int(f.attrs["block_nb_x"])
    blocks = f["blocks"][...]
    data = {v: f[VARS[v]][...] for v in VARS}
    f.close()
    return nb, blocks, data


def refined_base_blocks(blocks, nb):
    """Set of base-block coords (cx,cy,cz) that are refined (have level>=1
    children). origin is in level-l cells: base coord = origin//nb//2**l."""
    fine = set()
    for ox, oy, oz, lev in blocks:
        if lev >= 1:
            f = nb * (2 ** lev)
            fine.add((int(ox) // f, int(oy) // f, int(oz) // f))
    return fine


def coarse_leaf_index(blocks, nb):
    """Map base-block coord -> bid for every level-0 leaf."""
    idx = {}
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        if lev == 0:
            idx[(int(ox) // nb, int(oy) // nb, int(oz) // nb)] = bid
    return idx


def classify(c, fine, gnbt, periodic):
    """Return 'face'|'edge'|'corner'|'interior' for base-block coord c by the
    strongest contact it has with a refined block among its 26 neighbours."""
    def wrap(cn):
        out = list(cn)
        for d in range(3):
            if periodic[d]:
                out[d] %= gnbt[d]
            elif out[d] < 0 or out[d] >= gnbt[d]:
                return None
        return tuple(out)

    best = 0  # 0 none, 1 corner, 2 edge, 3 face
    for dz in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                man = abs(dx) + abs(dy) + abs(dz)
                if man == 0:
                    continue
                cn = wrap((c[0] + dx, c[1] + dy, c[2] + dz))
                if cn is None or cn not in fine:
                    continue
                kind = {1: 3, 2: 2, 3: 1}[man]  # face/edge/corner
                best = max(best, kind)
    return {0: "interior", 1: "corner", 2: "edge", 3: "face"}[best]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("patch", help="refined-patch block-table field")
    ap.add_argument("base", help="base-grid control field (same steps, no patch)")
    ap.add_argument("--periodic-y", action="store_true",
                    help="wrap y when classifying (default: walls in y)")
    a = ap.parse_args()

    nb, pblocks, pdata = read_blocks(a.patch)
    nbb, bblocks, bdata = read_blocks(a.base)
    assert nb == nbb, "block sizes differ"

    # global base-block lattice
    gnbt = tuple(int(b[:, d].max()) // nb + 1 for d, b in
                 ((0, bblocks), (1, bblocks), (2, bblocks)))
    periodic = (True, a.periodic_y, True)

    fine = refined_base_blocks(pblocks, nb)
    pidx = coarse_leaf_index(pblocks, nb)
    bidx = coarse_leaf_index(bblocks, nb)

    # accumulate squared difference per class per variable
    classes = ("interior", "face", "edge", "corner")
    acc = {cl: {v: [0.0, 0] for v in VARS} for cl in classes}
    nleaf = {cl: 0 for cl in classes}
    for c, pbid in pidx.items():
        if c not in bidx:
            continue
        cl = classify(c, fine, gnbt, periodic)
        nleaf[cl] += 1
        bbid = bidx[c]
        for v in VARS:
            d = (pdata[v][pbid] - bdata[v][bbid]).ravel()
            acc[cl][v][0] += float(np.sum(d * d))
            acc[cl][v][1] += d.size

    print(f"patch:  {a.patch}")
    print(f"base :  {a.base}")
    print(f"refined base blocks: {len(fine)}   coarse leaves classified: "
          f"{sum(nleaf.values())}")
    print(f"  {'class':>8} {'nleaf':>6}   "
          + " ".join(f"{v+'_rmsΔ':>11}" for v in VARS))
    for cl in classes:
        rms = {}
        for v in VARS:
            s, n = acc[cl][v]
            rms[v] = np.sqrt(s / n) if n else 0.0
        print(f"  {cl:>8} {nleaf[cl]:6d}   "
              + " ".join(f"{rms[v]:11.4e}" for v in VARS))

    # ratios vs interior background and vs face
    print("\n  ratios (rms diff relative to FACE class; gate: edge,corner <= ~face):")
    base_rms = {v: np.sqrt(acc['face'][v][0] / max(acc['face'][v][1], 1)) for v in VARS}
    for cl in ("face", "edge", "corner"):
        r = {}
        for v in VARS:
            s, n = acc[cl][v]
            x = np.sqrt(s / n) if n else 0.0
            r[v] = x / base_rms[v] if base_rms[v] > 0 else 0.0
        print(f"  {cl:>8}        "
              + " ".join(f"{r[v]:11.3f}" for v in VARS))


if __name__ == "__main__":
    main()
