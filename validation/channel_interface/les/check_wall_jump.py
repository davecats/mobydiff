#!/usr/bin/env python3
"""Eddy viscosity of the wall rows next to a level jump (wall_jump.ini).

The flow is x- and z-invariant, so within one level the eddy viscosity of a
wall row has one value; the columns next to the 2:1 interface show what the
SGS kernel makes of the halo it reads there. Reports, per level and wall, the
relative deviation from the level's median of the column on the LOW-x side of
a block face at the jump (its halo index nb+1 is the edge halo of the wall
ghost row) and gates it:

    fine side    <= 2e-2   (prolongation of the coarse ghost row)
    coarse side  <= 8e-2   (restriction of the ONE fine ghost row: first
                            order, see the README)

Before the exchange wrote that halo (main up to 154e48f) it kept its initial
value 0 and the same numbers were 0.26 and 0.13.
"""
import sys
import h5py
import numpy as np

LIMIT = {0: 8e-2, 1: 2e-2}
NAME = {0: "coarse", 1: "fine"}

with h5py.File(sys.argv[1], "r") as f:
    blocks = f["blocks"][...]                      # origin x, y, z (level cells), level
    nb = round(f["nut"].shape[-1])
    nut = f["nut"][...].reshape(len(blocks), nb, nb, nb)   # block, k, j, i

status = 0
for level in (0, 1):
    of_level = np.where(blocks[:, 3] == level)[0]
    y_top = blocks[of_level, 1].max()
    for wall, y0, row in (("lower", 0, 0), ("upper", y_top, nb - 1)):
        column = {}                                # global x index -> values of the row
        for b in of_level[blocks[of_level, 1] == y0]:
            for i in range(nb):
                column.setdefault(blocks[b, 0] + i, []).append(nut[b, :, row, i])
        gx = np.array(sorted(column))
        lo = np.array([np.min(np.concatenate(column[g])) for g in gx])
        hi = np.array([np.max(np.concatenate(column[g])) for g in gx])
        median = np.median(hi)
        dev = np.maximum(abs(hi/median - 1.0), abs(lo/median - 1.0))
        # The column whose +x neighbour is at the other level: the last
        # column of the fine band, or the coarse column just below it.
        jump = gx.max() if level == 1 else max(g for g in gx if g + 1 not in column)
        if level == 0 and jump == gx.max():        # periodic wrap: the band is interior
            jump = max(g for g in gx[:-1] if g + 1 not in column)
        at_jump = dev[gx == jump][0]
        ok = at_jump <= LIMIT[level]
        status |= not ok
        print(f"{'PASS' if ok else 'FAIL'} {NAME[level]:6s} side, {wall} wall: nut {median:.4e}, "
              f"column next to the jump (x index {jump}) off by {at_jump:.2e} (limit {LIMIT[level]:.0e}); "
              f"largest elsewhere {np.max(dev[gx != jump]):.2e}")
sys.exit(status)
