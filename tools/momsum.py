#!/usr/bin/env python3
"""Momentum-conservation check at the 2:1 interface (MOBY_TERMDUMP).

The divergence-form advection -div(q_i q) and the viscous flux div(nu grad q_i)
conserve total momentum: summed over a periodic domain every interior face flux
telescopes (flux_p of a cell = flux_m of its neighbour), so

    Sum_cells vol * term_i        (net momentum source for component i)

is round-off on a single-level grid. The ONLY faces that do not cancel are the
2:1 interface faces, where the coarse face flux differs from the summed fine
sub-face fluxes -- the un-refluxed momentum (Berger-Colella). So this GLOBAL sum
is the exact momentum-conservation error, and (since the uniform grid gives
round-off) it is entirely the interface imbalance -- the momentum analogue of the
mass gate Sum(vol*div) (tools/divsum.py).

LOCAL vs GLOBAL: per-cell vol*term is the *physical* flux divergence (nonzero
everywhere -- advection transports momentum), so there is no per-cell error to
report; conservation is the integral. To localize to ONE interface, run a case
with a single coarse-fine interface: its global sum IS that interface's local
imbalance. Run `uniform` as the zero reference and `slab`/`patch` for the
interface error.

Pass the MOBY_TERMDUMP `<prefix>_adv_<step>.h5` (run with MOBY_TERMDUMP=1/2/3 for
u/v/w) and optionally the matching `_dif`. Each component's sum must be
~round-off for momentum to be conserved.

Usage: momsum.py FILE_adv.h5 [FILE_dif.h5] ...
"""
import sys
import os
import h5py
import numpy as np


def analyze(fn):
    h5 = h5py.File(fn, "r")
    lx = float(h5.attrs["lx"]); nx = int(h5.attrs["nx"])
    blocks = h5["blocks"][...]
    t = h5["un"][...] + h5["vn"][...] + h5["wn"][...]   # total flux-div per cell (nblk,k,j,i)
    h5.close()
    s = 0.0; scale = 0.0
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        vol = (lx / (nx * 2 ** lev)) ** 3
        s += vol * t[bid].sum()
        scale += vol * np.abs(t[bid]).sum()
    return s, scale


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    for fn in sys.argv[1:]:
        s, scale = analyze(fn)
        kind = "adv" if "_adv" in fn else ("dif" if "_dif" in fn else "term")
        rel = abs(s) / scale if scale > 0 else 0.0
        print(f"{os.path.basename(fn)}  [{kind}]")
        print(f"  Sum(vol*{kind}) = {s:+.6e}   (rel {rel:.2e})   ~0 = momentum conserved")


if __name__ == "__main__":
    main()
