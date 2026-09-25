#!/usr/bin/env python3
"""Shared geometry helpers for the passive-scalar gates.

A solver field file is in the block-table layout: per-variable datasets of
shape (nBlocksGlobal, nbz, nby, nbx) plus a `blocks` table (ox, oy, oz, level
in level-l cells) and the BASE-level node lines x/y/z. Level-l node lines are
l rounds of midpoint subdivision of the base line in every REFINED direction
([blocks] refine_dims, attribute `refine_dims`, absent = xyz octree) -- the
same construction blocks.f90 uses, so cell centres and widths here are the
solver's to the last bit.
"""

from __future__ import annotations

import numpy as np
import h5py


def subdivide(line: np.ndarray) -> np.ndarray:
    """Midpoint subdivision of a node line (blocks.f90 subdivide_node_line)."""
    fine = np.empty(2 * (line.size - 1) + 1, dtype=np.float64)
    fine[0::2] = line
    fine[1::2] = 0.5 * (line[:-1] + line[1:])
    return fine


class BlockGeometry:
    """Per-level node lines and the block table of one field file."""

    def __init__(self, h5: h5py.File):
        self.blocks = h5["blocks"][...]
        self.nb = (int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]),
                   int(h5.attrs["block_nb_z"]))
        self.mask = np.asarray(h5.attrs.get("refine_dims", [1, 1, 1]), dtype=int)
        self.leng = (float(h5.attrs["lx"]), float(h5.attrs["ly"]), float(h5.attrs["lz"]))
        lmax = int(self.blocks[:, 3].max())
        base = [h5["x"][...], h5["y"][...], h5["z"][...]]
        self.lines = [[base[d].copy()] for d in range(3)]
        for d in range(3):
            for _ in range(lmax):
                nxt = subdivide(self.lines[d][-1]) if self.mask[d] else self.lines[d][-1]
                self.lines[d].append(nxt)

    def block_axes(self, bid: int):
        """(centres, widths) triple for block bid, one array per direction."""
        ox, oy, oz, lev = (int(v) for v in self.blocks[bid])
        out = []
        for d, o in enumerate((ox, oy, oz)):
            line = self.lines[d][lev]
            lo = line[o:o + self.nb[d]]
            hi = line[o + 1:o + self.nb[d] + 1]
            out.append((0.5 * (lo + hi), hi - lo))
        return out

    def mesh(self, bid: int):
        """(x, y, z, dV) broadcast to the (nbz, nby, nbx) dataset order."""
        (xc, dx), (yc, dy), (zc, dz) = self.block_axes(bid)
        x = xc[None, None, :]
        y = yc[None, :, None]
        z = zc[:, None, None]
        dV = dz[:, None, None] * dy[None, :, None] * dx[None, None, :]
        return x, y, z, dV

    @property
    def n_blocks(self) -> int:
        return self.blocks.shape[0]


def integrate(h5: h5py.File, name: str) -> float:
    """Volume integral of a cell-centred dataset over every stored cell."""
    geo = BlockGeometry(h5)
    data = h5[name]
    total = 0.0
    for bid in range(geo.n_blocks):
        _, _, _, dV = geo.mesh(bid)
        total += float(np.sum(data[bid] * dV))
    return total


def volume(h5: h5py.File) -> float:
    geo = BlockGeometry(h5)
    total = 0.0
    for bid in range(geo.n_blocks):
        _, _, _, dV = geo.mesh(bid)
        total += float(np.sum(dV))
    return total


def field_error(h5: h5py.File, name: str, exact):
    """(L2, Linf) of dataset `name` against exact(x, y, z), volume-weighted L2."""
    geo = BlockGeometry(h5)
    data = h5[name]
    num = 0.0
    vol = 0.0
    linf = 0.0
    for bid in range(geo.n_blocks):
        x, y, z, dV = geo.mesh(bid)
        err = data[bid] - exact(x, y, z)
        num += float(np.sum(err * err * dV))
        vol += float(np.sum(dV))
        linf = max(linf, float(np.max(np.abs(err))))
    return np.sqrt(num / vol), linf
