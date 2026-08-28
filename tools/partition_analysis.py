#!/usr/bin/env python3
"""Offline analysis of how a block-to-rank partitioning cuts the 2:1 interface.

P1 of docs/next_session_2to1_penalty.md. The measured defect: cross-level peer
points are IDENTICAL at 2 and 4 ranks (923700 / 923692), because `refine_dims =
xz` puts the y tile in the high Morton bits, so every fine block precedes every
coarse one and a linear split of the leaf list always cuts the level change.
Every interface transfer that could be a device-local copy becomes an MPI
message.

This script answers, without running the solver, how much of that is avoidable
and what it would cost in load balance. It reads a leaf table -- either the
`blocks` dataset of any case/snapshot file, or the stdout of
`build_cpu/leaftable_test` -- reconstructs the 26-neighbour block graph in
finest-lattice index space, and scores candidate partitionings.

  ./partition_analysis.py <leaves.txt | file.h5> --nb 64 44 48 --dims xz \\
      --grid 2048 176 96 --periodic 0 0 1 --ranks 2 4 8 16

The traffic proxy is the SHARED AREA in finest cells, which is what an exchange
entry actually carries; pair counts alone would weight a corner like a face.
"""

import argparse
import sys
from itertools import product


def read_leaves(path):
    """-> (ids, origins[level-l cells], levels). Accepts leaftable_test stdout
    or an HDF5 file with a `blocks` dataset (origin x,y,z + level per row)."""
    if path.endswith(".h5"):
        import h5py
        import numpy as np
        with h5py.File(path, "r") as f:
            b = np.asarray(f["blocks"])
        return [(int(r[0]), int(r[1]), int(r[2]), int(r[3])) for r in b]
    rows = []
    for line in open(path):
        f = line.split()
        if len(f) == 5 and f[0].isdigit():
            rows.append(tuple(int(v) for v in f[1:]))
    if not rows:
        sys.exit(f"no leaf rows found in {path}")
    return rows


def extents(rows, nb, mask, nlev):
    """Finest-lattice [lo, hi) per block. A level-l cell spans 2**((L-l)*mask)
    finest cells in a refined direction and 1 in an unrefined one."""
    out = []
    for (ox, oy, oz, lv) in rows:
        lo, hi = [], []
        for d, o in enumerate((ox, oy, oz)):
            s = 2 ** ((nlev - 1 - lv) * mask[d])
            lo.append(o * s)
            hi.append((o + nb[d]) * s)
        out.append((lo, hi, lv))
    return out


def neighbour_pairs(ext, dom, periodic, tile):
    """26-neighbour pairs with their shared area in finest cells.

    Two blocks are neighbours when, in every direction, they either OVERLAP or
    TOUCH (hi_a == lo_b, or across a periodic seam), and touch in at least one.
    The exchanged surface is the product of the overlap extents in the
    directions that do not touch -- a face pair gets an area, an edge pair a
    length, a corner a single point.

    Bucketed by level-0 tile: refinement only ever subdivides a tile, so a
    leaf's neighbours all lie in its own or an adjacent tile. Exact, and it
    turns an O(n^2) sweep into O(n) -- needed at 25k leaves.
    """
    ntile = [max(1, -(-dom[d] // tile[d])) for d in range(3)]
    buckets = {}
    for i, (lo, _hi, _lv) in enumerate(ext):
        key = tuple(lo[d] // tile[d] for d in range(3))
        buckets.setdefault(key, []).append(i)

    def adjacent(a0, a1, b0, b1, d):
        if (a1 == b0) or (b1 == a0):
            return True
        return periodic[d] and ((a1 % dom[d] == b0 % dom[d])
                                or (b1 % dom[d] == a0 % dom[d]))

    pairs, seen = [], set()
    for key, members in buckets.items():
        cand = []
        for off in product((-1, 0, 1), repeat=3):
            k = tuple((key[d] + off[d]) % ntile[d] if periodic[d]
                      else key[d] + off[d] for d in range(3))
            cand.extend(buckets.get(k, ()))
        for i in members:
            loA, hiA, lvA = ext[i]
            for j in cand:
                if j <= i:
                    continue
                if (i, j) in seen:
                    continue
                loB, hiB, lvB = ext[j]
                touch, span, ok = 0, [], True
                for d in range(3):
                    a0, a1, b0, b1 = loA[d], hiA[d], loB[d], hiB[d]
                    ov = min(a1, b1) - max(a0, b0)
                    if ov > 0:
                        span.append(ov)
                        continue
                    if adjacent(a0, a1, b0, b1, d):
                        touch += 1
                    else:
                        ok = False
                        break
                if ok and touch >= 1:
                    seen.add((i, j))
                    area = 1
                    for s in span:
                        area *= s
                    pairs.append((i, j, area, lvA != lvB, touch))
    return pairs


def morton_key(cx, cz):
    k = 0
    for b in range(20):
        k |= ((cx >> b) & 1) << (2 * b)
        k |= ((cz >> b) & 1) << (2 * b + 1)
    return k


def assign_current(rows, nranks):
    """The shipped rule: leaf ids split linearly (blocks.f90 zorder_owner)."""
    n = len(rows)
    q, r = divmod(n, nranks)
    split = r * (q + 1)
    out = []
    for i in range(n):
        out.append(i // (q + 1) if i < split else r + (i - split) // q)
    return out


def assign_columns(rows, ext, nb, mask, nlev, nranks, morton=True):
    """Column-wise: every leaf in the same level-0 (x,z) column -- all y tiles,
    both levels -- goes to one rank. The 2:1 interface here is a horizontal
    plane, so cutting VERTICALLY keeps every cross-level pair inside a rank."""
    colw = [nb[d] * 2 ** ((nlev - 1) * mask[d]) for d in range(3)]
    cols = {}
    for i, (lo, hi, lv) in enumerate(ext):
        c = (lo[0] // colw[0], lo[2] // colw[2])
        cols.setdefault(c, []).append(i)
    keys = sorted(cols, key=lambda c: morton_key(*c) if morton else c)
    # Hand columns out along the curve, advancing a rank once it holds its
    # share. Compare the RUNNING TOTAL against the cumulative target -- against
    # a per-rank load it would advance far too early.
    owner = [0] * len(ext)
    total, placed, r = len(ext), 0, 0
    for c in keys:
        if r < nranks - 1 and placed + len(cols[c]) / 2 > total * (r + 1) / nranks:
            r += 1
        for i in cols[c]:
            owner[i] = r
        placed += len(cols[c])
    return owner


def assign_ykey_low(ext, nb, mask, nlev, nranks):
    """The cheap fix: keep the LINEAR leaf split exactly as it is, and move the
    y tile from the HIGH bits of the xz Morton key to the LOW bits. The leaf
    order then runs down each (x,z) column before moving on, so the existing
    closed-form zorder_owner split becomes column-wise for free -- no
    partitioner, no ownership table, no load-balance heuristic.

    Scored here as "reorder the key, then split leaves linearly", which is what
    the implementation would actually do; a rank boundary may therefore fall
    INSIDE a column, unlike the idealised column-wise assignment."""
    colw = [nb[d] * 2 ** ((nlev - 1) * mask[d]) for d in range(3)]
    order = sorted(range(len(ext)), key=lambda i: (
        morton_key(ext[i][0][0] // colw[0], ext[i][0][2] // colw[2]),
        ext[i][0][1], ext[i][0][0], ext[i][0][2]))
    n = len(ext)
    q, r = divmod(n, nranks)
    split = r * (q + 1)
    owner = [0] * n
    for pos, i in enumerate(order):
        owner[i] = pos // (q + 1) if pos < split else r + (pos - split) // q
    return owner


def score(owner, pairs, nranks, label):
    tot = {True: 0, False: 0}
    cut = {True: 0, False: 0}
    for i, j, area, cross, _t in pairs:
        tot[cross] += area
        if owner[i] != owner[j]:
            cut[cross] += area
    load = [0] * nranks
    for o in owner:
        load[o] += 1
    imb = max(load) / (sum(load) / nranks)
    print(f"  {label:<26} cross-level cut {cut[True]:>9,} / {tot[True]:>9,}"
          f" ({100*cut[True]/max(1,tot[True]):5.1f} %)   same-level cut"
          f" {cut[False]:>9,}   total peer area {cut[True]+cut[False]:>9,}"
          f"   load imb {imb:.3f}")
    return cut[True] + cut[False]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("leaves")
    p.add_argument("--nb", nargs=3, type=int, required=True)
    p.add_argument("--grid", nargs=3, type=int, required=True)
    p.add_argument("--dims", default="xyz")
    p.add_argument("--periodic", nargs=3, type=int, default=[1, 1, 1])
    p.add_argument("--ranks", nargs="+", type=int, default=[2, 4, 8, 16])
    a = p.parse_args()

    mask = [1, 1, 1] if a.dims == "xyz" else [1, 0, 1]
    rows = read_leaves(a.leaves)
    nlev = max(r[3] for r in rows) + 1
    ext = extents(rows, a.nb, mask, nlev)
    dom = [a.grid[d] * 2 ** ((nlev - 1) * mask[d]) for d in range(3)]
    tile = [a.nb[d] * 2 ** ((nlev - 1) * mask[d]) for d in range(3)]
    pairs = neighbour_pairs(ext, dom, [bool(x) for x in a.periodic], tile)

    nfine = sum(1 for r in rows if r[3] > 0)
    xarea = sum(ar for _i, _j, ar, cr, _t in pairs if cr)
    print(f"{len(rows)} leaves ({len(rows)-nfine} coarse, {nfine} fine), "
          f"{nlev} levels, {len(pairs)} neighbour pairs")
    print(f"cross-level shared area: {xarea:,} finest cells "
          f"({100*xarea/sum(ar for _i,_j,ar,_c,_t in pairs):.1f} % of all "
          f"neighbour area)\n")

    for R in a.ranks:
        print(f"ranks = {R}")
        score(assign_current(rows, R), pairs, R, "current (Morton linear)")
        score(assign_columns(rows, ext, a.nb, mask, nlev, R), pairs, R,
              "column-wise (Morton x,z)")
        score(assign_ykey_low(ext, a.nb, mask, nlev, R), pairs, R,
              "y-tile in LOW key bits")
        print()


if __name__ == "__main__":
    main()
