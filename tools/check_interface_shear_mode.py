#!/usr/bin/env python3
"""Interface shear-mode stability gate (docs/interface_review.md vi).

Reports max|vn| (interface-normal velocity) versus step for the
tutorials/interface_shear_mode case. The structured checkerboard seed on the
fine-side interface-normal velocity is a zero-coarse-average mode invisible to
the coarse pressure; if the projection coupling at the 2:1 interface is
non-contractive it GROWS, otherwise it CONTRACTS.

Usage:
    check_interface_shear_mode.py <rundir> [--prefix shear]

Verdict: PASS if max|vn| contracts monotonically (final < first), FAIL if it
grows. Run projection-only (MOBY_PROJECTION_ONLY=1) to isolate the linear
projection instability from the predictor.
"""
from __future__ import annotations
import argparse, glob, re, sys
import h5py
import numpy as np


def reassemble(h5, fld):
    blocks = h5["blocks"][...]
    nb = (int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]), int(h5.attrs["block_nb_z"]))
    lmax = int(blocks[:, 3].max())
    gx, gy, gz = (int(h5.attrs[a]) for a in ("nx", "ny", "nz"))
    A = np.full((gz * 2**lmax, gy * 2**lmax, gx * 2**lmax), np.nan)
    data = h5[fld][...]
    for b, (ox, oy, oz, lv) in enumerate(blocks):
        f = 2**(lmax - lv)
        blk = np.repeat(np.repeat(np.repeat(data[b], f, 0), f, 1), f, 2)
        A[oz*f:oz*f+nb[2]*f, oy*f:oy*f+nb[1]*f, ox*f:ox*f+nb[0]*f] = blk
    return A


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("rundir")
    p.add_argument("--prefix", default="shear")
    args = p.parse_args()

    files = sorted(glob.glob(f"{args.rundir}/{args.prefix}_*.h5"),
                   key=lambda q: int(re.findall(r"_(\d+)\.h5", q)[0]))
    if not files:
        print("no fields found"); return 1

    print(f"{'step':>6} {'max|vn|':>13} {'max|vn|@iface':>15}")
    series = []
    for f in files:
        s = int(re.findall(r"_(\d+)\.h5", f)[0])
        with h5py.File(f, "r") as h5:
            vn = reassemble(h5, "vn")
            ny = vn.shape[1]
            # interface y-planes: the refined band top (y=ny/2) and the wrap (y=0)
            iface = max(np.nanmax(np.abs(vn[:, ny//2 - 1:ny//2 + 1, :])),
                        np.nanmax(np.abs(vn[:, 0:1, :])), np.nanmax(np.abs(vn[:, ny-1:ny, :])))
            mx = np.nanmax(np.abs(vn))
        print(f"{s:6d} {mx:13.4e} {iface:15.4e}")
        series.append((s, mx))

    first, last = series[0][1], series[-1][1]
    if not np.isfinite(last):
        print("\nVERDICT: FAIL (diverged to NaN/Inf)"); return 1
    growth = last / first if first > 0 else float("inf")
    if last < first:
        print(f"\nVERDICT: PASS - max|vn| contracted {first:.3e} -> {last:.3e} (x{growth:.2e})")
        return 0
    print(f"\nVERDICT: FAIL - max|vn| grew {first:.3e} -> {last:.3e} (x{growth:.2e})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
