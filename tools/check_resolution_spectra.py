#!/usr/bin/env python3
"""Diagnose grid resolution from 1D u-spectra of a field snapshot.

For a turbulent-channel field, plot the streamwise (E_uu vs kx) and spanwise
(E_uu vs kz) one-dimensional energy spectra at one or more y+ heights, averaged
over the homogeneous directions (and over several snapshots if given). A
well-resolved field drops several decades with a steepening dissipative roll-off
near the grid cut-off; a *high-k pile-up* (the spectrum flattening or turning up
at the largest wavenumbers) is the unambiguous signature of under-resolution.
The kx plot diagnoses dx, the kz plot diagnoses dz.

Handles the block-table and legacy global-3D field layouts (uses the same
plane extraction as channel_interface_validation.py). For block files only the
requested level is assembled; a uniform/legacy file is level 0.

Usage:
  python3 tools/check_resolution_spectra.py FIELD1.h5 [FIELD2.h5 ...] \
      [--yplus 15 90 175] [--retau 180] [--level 0] [--out spectra.png]
"""

from __future__ import annotations

import argparse
import glob
import os

import h5py
import numpy as np


class Snapshot:
    """x-z plane extraction at a y-cell row (block-table or legacy layout)."""

    def __init__(self, path):
        self.h5 = h5py.File(path, "r")
        self.block = self.h5["un"].ndim == 4
        self.y = self.h5["y"][...]
        self.lx = float(self.h5.attrs["lx"])
        self.lz = float(self.h5.attrs["lz"])
        if self.block:
            self.blocks = self.h5["blocks"][...]   # (ox, oy, oz, level), level-l cells
            # Per-direction block sizes (cubic for the lattice; the whole rank box
            # for a single-block RANKBOX/uniform file).
            self.nbx = int(self.h5.attrs["block_nb_x"])
            self.nby = int(self.h5.attrs["block_nb_y"])
            self.nbz = int(self.h5.attrs["block_nb_z"])

    def yline(self, level):
        """Node line at a refinement level (midpoint subdivision of level 0)."""
        line = self.y
        for _ in range(level):
            fine = np.empty(2*(line.size - 1) + 1)
            fine[0::2] = line
            fine[1::2] = 0.5*(line[:-1] + line[1:])
            line = fine
        return line

    def plane(self, level, jrow):
        """Global (z, x) plane of u at level-l y-cell row jrow (1-based)."""
        data = self.h5["un"]
        if not self.block:
            return data[:, jrow-1, :]
        nbx, nby, nbz = self.nbx, self.nby, self.nbz
        nxl = int(self.h5.attrs["nx"])*2**level
        nzl = int(self.h5.attrs["nz"])*2**level
        out = np.full((nzl, nxl), np.nan)
        jc = jrow - 1
        for bid, (ox, oy, oz, lev) in enumerate(self.blocks):
            if lev != level or not (oy <= jc < oy + nby):
                continue
            out[oz:oz+nbz, ox:ox+nbx] = data[bid][:, jc - oy, :]
        return out

    def close(self):
        self.h5.close()


def row_of_yplus(line, yplus, retau):
    """1-based cell row whose centre is nearest the target y+ (lower half)."""
    cent = 0.5*(line[:-1] + line[1:])
    return int(np.argmin(np.abs(cent*retau - yplus))) + 1


def spectra(files, level, jrow):
    """Snapshot- and homogeneous-direction-averaged 1D u-spectra at a y row."""
    ex = ez = None
    n = 0
    lx = lz = None
    for path in files:
        s = Snapshot(path)
        lx, lz = s.lx, s.lz
        p = s.plane(level, jrow)
        s.close()
        if np.isnan(p).any():
            continue
        nz, nx = p.shape
        u = p - p.mean()
        ax = np.mean(np.abs(np.fft.rfft(u, axis=1))**2, axis=0)/nx**2
        az = np.mean(np.abs(np.fft.rfft(u, axis=0))**2, axis=1)/nz**2
        ex = ax if ex is None else ex + ax
        ez = az if ez is None else ez + az
        n += 1
    if n == 0:
        return None
    kx = 2*np.pi/lx*np.arange(ex.size)
    kz = 2*np.pi/lz*np.arange(ez.size)
    return (kx, ex/n), (kz, ez/n)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fields", nargs="+", help="field HDF5 file(s); averaged together")
    parser.add_argument("--yplus", type=float, nargs="+", default=[15.0, 90.0, 175.0],
                        help="y+ heights to probe (default: 15 90 175)")
    parser.add_argument("--retau", type=float, default=180.0)
    parser.add_argument("--level", type=int, default=0,
                        help="refinement level to assemble (block files; default 0)")
    parser.add_argument("--out", default="resolution_spectra.png")
    args = parser.parse_args()

    files = []
    for f in args.fields:
        files += sorted(glob.glob(f)) if any(c in f for c in "*?[") else [f]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    line = Snapshot(files[0]).yline(args.level)
    ny = line.size - 1
    fig, axes = plt.subplots(len(args.yplus), 2, figsize=(11, 3.4*len(args.yplus)),
                             squeeze=False)
    for r, yp in enumerate(args.yplus):
        jrow = row_of_yplus(line, yp, args.retau)
        ycent = 0.5*(line[jrow-1] + line[jrow])
        sp = spectra(files, args.level, jrow)
        for c, (comp, lab) in enumerate((("kx", "E_uu(kx)"), ("kz", "E_uu(kz)"))):
            ax = axes[r][c]
            if sp:
                k, e = sp[c]
                ax.loglog(k[1:], e[1:], "o-", ms=3, label="spectrum")
                # -5/3 guide anchored to the low-k end
                kk = k[1:k.size//3 + 1]
                ax.loglog(kk, e[1]*(kk/k[1])**(-5.0/3.0), "k--", lw=0.8, label="-5/3")
            ax.set_xlabel(comp); ax.set_ylabel("E_uu")
            ax.set_title(f"y+ ~ {yp:.0f} (row {jrow}/{ny}, y={ycent:.3f})  {lab}")
            ax.legend(fontsize=8)
    fig.suptitle(f"u-spectra: {os.path.basename(files[0])}"
                 + (f"  (+{len(files)-1} more)" if len(files) > 1 else "")
                 + f"   level {args.level}")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")
    print("under-resolution shows as a high-k pile-up (spectrum flat/rising near "
          "the right edge); kx diagnoses dx, kz diagnoses dz.")


if __name__ == "__main__":
    main()
