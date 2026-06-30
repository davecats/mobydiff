#!/usr/bin/env python3
"""Time-averaged edge/corner banding metric for an embedded refined patch
(gate 5 / DEVELOPED flow). Companion to patch_interface_diff.py (which is the
instantaneous-snapshot version).

Given a set of 3D field snapshots over the averaging window for the refined-patch
run AND for a base-grid control (same channel, no patch), it computes per coarse
cell the time-mean and the fluctuation rms (sqrt(<q^2> - <q>^2)), then aggregates
two diagnostics by each coarse leaf's adjacency to the patch
(interior/face/edge/corner -- see patch_interface_diff.py):

  MEAN footprint   rms over class cells of (<q>_patch - <q>_base).
                   Patch and base share the base lattice outside the patch, so
                   this is a matched cell-by-cell difference. Away from the patch
                   the time-mean is statistically the same -> interior ~ 0; a
                   nonzero face/edge/corner value is the mean interface transport
                   footprint (the -<u'v'>-type effect). GATE: edge/corner not
                   worse than face; small.

  BAND ratio       <fluct_rms>_patch / <fluct_rms>_base over class cells, matched
                   cell-by-cell (controls for the strong y-dependence of the
                   stresses). interior ~ 1.0 (no patch nearby). A 2:1 interface
                   BAND shows as a ratio > 1 at face/edge/corner. GATE: ratios
                   ~1; edge/corner not worse than face; consistent with the flat-
                   face ~5% small-scale v'/w' loss (ratio ~0.95 if anything).

Usage:
  python3 tools/patch_interface_stats.py \
      --patch 'runs/gate5/patch/stats/channel_field_*.h5' \
      --base  'runs/gate5/base/stats/base_field_*.h5' [--periodic-y]
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from patch_interface_diff import (  # noqa: E402
    VARS, read_blocks, refined_base_blocks, coarse_leaf_index, classify)

import h5py  # noqa: E402


# var -> HDF5 dataset name. "nut" (LES eddy viscosity) is included only when the
# snapshot files carry it (auto-detected), so this tool still works on non-LES runs.
BASE_VARMAP = {"u": "un", "v": "vn", "w": "wn", "p": "pn"}


def varmap_for(path):
    with h5py.File(path, "r") as f:
        has_nut = "nut" in f
    vm = dict(BASE_VARMAP)
    if has_nut:
        vm["nut"] = "nut"
    return vm


def accumulate(paths, varmap):
    """Time mean and variance per cell over the snapshot set, on the run's own
    leaf layout. Returns nb, blocks, mean{v}, fluct_rms{v} (each (nleaf,nb,nb,nb))."""
    paths = sorted(paths)
    if not paths:
        sys.exit("no snapshots matched")
    nb, blocks, _ = read_blocks(paths[0])
    s1 = {v: None for v in varmap}
    s2 = {v: None for v in varmap}
    n = 0
    for p in paths:
        f = h5py.File(p, "r")
        for v in varmap:
            q = f[varmap[v]][...].astype(np.float64)
            s1[v] = q.copy() if s1[v] is None else s1[v] + q
            s2[v] = q * q if s2[v] is None else s2[v] + q * q
        f.close()
        n += 1
    mean = {v: s1[v] / n for v in varmap}
    fluct = {v: np.sqrt(np.maximum(s2[v] / n - mean[v] ** 2, 0.0)) for v in varmap}
    print(f"  {os.path.dirname(paths[0])}: {n} snapshots")
    return nb, blocks, mean, fluct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patch", required=True, help="glob of refined-patch snapshots")
    ap.add_argument("--base", required=True, help="glob of base-control snapshots")
    ap.add_argument("--periodic-y", action="store_true")
    a = ap.parse_args()

    ppaths, bpaths = sorted(glob.glob(a.patch)), sorted(glob.glob(a.base))
    if not ppaths or not bpaths:
        sys.exit("no snapshots matched --patch/--base")
    # nut is compared only when BOTH runs carry it (e.g. both LES); otherwise the
    # metric reduces to u,v,w,p exactly as before.
    varmap = {v: ds for v, ds in varmap_for(ppaths[0]).items() if v in varmap_for(bpaths[0])}
    VARS = list(varmap)
    print(f"variables: {VARS}")

    print("accumulating patch run...")
    nb, pblocks, pmean, pfluct = accumulate(ppaths, varmap)
    print("accumulating base run...")
    nbb, bblocks, bmean, bfluct = accumulate(bpaths, varmap)
    assert nb == nbb, "block sizes differ"

    gnbt = tuple(int(bblocks[:, d].max()) // nb + 1 for d in range(3))
    periodic = (True, a.periodic_y, True)
    fine = refined_base_blocks(pblocks, nb)
    pidx = coarse_leaf_index(pblocks, nb)
    bidx = coarse_leaf_index(bblocks, nb)

    classes = ("interior", "face", "edge", "corner")
    classmap = {c: classify(c, fine, gnbt, periodic) for c in pidx}
    nleaf = {cl: 0 for cl in classes}
    for c in pidx:
        nleaf[classmap[c]] += 1

    def rms(acc):
        s, k = acc
        return np.sqrt(s / k) if k else 0.0

    # --- BAND ratio: cross-run, matched cell-by-cell (controls for y(stresses)).
    fp = {cl: {v: [0.0, 0] for v in VARS} for cl in classes}
    fb = {cl: {v: [0.0, 0] for v in VARS} for cl in classes}
    for c, pbid in pidx.items():
        if c not in bidx:
            continue
        cl = classmap[c]
        bbid = bidx[c]
        for v in VARS:
            qp = pfluct[v][pbid].ravel(); qb = bfluct[v][bbid].ravel()
            fp[cl][v][0] += float(np.sum(qp ** 2)); fp[cl][v][1] += qp.size
            fb[cl][v][0] += float(np.sum(qb ** 2)); fb[cl][v][1] += qb.size

    # --- MEAN footprint: SINGLE-RUN deviation of each coarse cell's time-mean
    # from its x,z-HOMOGENEOUS time-mean at the same y-row (interior cells only).
    # No cross-run decorrelation: each run is referenced to itself. The base run
    # (no patch) gives the method's noise floor; the patch run's face/edge/corner
    # excess above that floor is the real mean interface footprint.
    def homog_footprint(blocks, mean, idx):
        # x,z-homogeneous mean per (v, global y-row) from interior coarse cells
        hsum = {v: {} for v in VARS}
        hcnt = {}
        for c, bid in idx.items():
            if classmap.get(c) != "interior":
                continue
            oy = int(blocks[bid][1])
            for jj in range(nb):
                gj = oy + jj
                hcnt[gj] = hcnt.get(gj, 0) + nb * nb
                for v in VARS:
                    hsum[v][gj] = hsum[v].get(gj, 0.0) + float(mean[v][bid][:, jj, :].sum())
        hom = {v: {gj: hsum[v][gj] / hcnt[gj] for gj in hcnt} for v in VARS}
        acc = {cl: {v: [0.0, 0] for v in VARS} for cl in classes}
        for c, bid in idx.items():
            cl = classmap.get(c)
            if cl is None:
                continue
            oy = int(blocks[bid][1])
            for jj in range(nb):
                gj = oy + jj
                for v in VARS:
                    if gj not in hom[v]:
                        continue
                    dev = (mean[v][bid][:, jj, :] - hom[v][gj]).ravel()
                    acc[cl][v][0] += float(np.sum(dev ** 2)); acc[cl][v][1] += dev.size
        return acc

    pf = homog_footprint(pblocks, pmean, pidx)
    bf = homog_footprint(bblocks, bmean, bidx)

    print(f"\nrefined base blocks: {len(fine)}   coarse leaves: {sum(nleaf.values())}")
    print("  class    nleaf:  " + "  ".join(f"{cl} {nleaf[cl]}" for cl in classes))

    print("\nMEAN footprint  rms(<q>_cell - <q>_row-homog) per class "
          "(PATCH run; BASE = noise floor):")
    print(f"  {'class':>8} {'run':>5}   " + " ".join(f"{v:>10}" for v in VARS))
    for cl in classes:
        print(f"  {cl:>8} {'patch':>5}   " + " ".join(f"{rms(pf[cl][v]):10.3e}" for v in VARS))
        print(f"  {cl:>8} {'base':>5}   " + " ".join(f"{rms(bf[cl][v]):10.3e}" for v in VARS))

    print("\nBAND ratio  <fluct_rms>_patch / <fluct_rms>_base per class "
          "(gate: ~1; edge,corner <= face):")
    print(f"  {'class':>8}          " + " ".join(f"{v:>10}" for v in VARS))
    for cl in classes:
        row = []
        for v in VARS:
            rp, rb = rms(fp[cl][v]), rms(fb[cl][v])
            row.append(rp / rb if rb > 0 else 0.0)
        print(f"  {cl:>8}          " + " ".join(f"{x:10.3f}" for x in row))


if __name__ == "__main__":
    main()
