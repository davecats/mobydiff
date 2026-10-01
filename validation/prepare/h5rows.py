"""Row matching between two block-table files.

A block-layout dataset holds one row per leaf, in the leaf order of the case
file its writer ran from -- and that order is no longer unique: a case file
built since 2026-10-01 carries the minimum-surface key-bit order in xyz mode,
an older or mobygeom-written one the legacy interleave (blocks.f90 leaf_key).
The leaves themselves are what must agree; rows are matched on
(origin, level) through each file's own `blocks` table.
"""
import numpy as np


def match_rows(cand_blocks, ref_blocks):
    """perm with cand_blocks[perm] == ref_blocks row for row, or None when
    the two tables do not hold the same leaves. The identity when the row
    orders already agree."""
    cb, rb = np.asarray(cand_blocks), np.asarray(ref_blocks)
    if cb.shape != rb.shape:
        return None
    if np.array_equal(cb, rb):
        return np.arange(len(rb))
    where = {tuple(int(v) for v in row): i for i, row in enumerate(cb)}
    try:
        perm = np.array([where[tuple(int(v) for v in row)] for row in rb])
    except KeyError:
        return None
    return perm if len(where) == len(rb) else None


def read_rows(dset, rows):
    """dset[rows] for an arbitrary (not increasing) row list -- h5py wants
    increasing indices."""
    order = np.argsort(rows, kind="stable")
    out = dset[rows[order]]
    back = np.empty_like(order)
    back[order] = np.arange(len(order))
    return out[back]
