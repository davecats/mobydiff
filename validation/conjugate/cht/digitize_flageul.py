#!/usr/bin/env python3
"""Digitise the CONJUGATE temperature variance from Flageul et al. 2015 fig. 5.

    ./digitize_flageul.py [--pdf ../../../literature/flageul.pdf] [--dpi 600]

Writes `flageul_fig5a_conjug.dat` (y+, <T'^2>/T_tau^2), the black solid
"Conjug" curve of figure 5 -- their G = G_2 = 1 case, i.e. solid conductivity
AND diffusivity equal to the fluid's, which is exactly our k1.

WHY DIGITISE AT ALL. The paper tabulates nothing, so the comparison had been
resting on four values read off the plot by eye. Two of them were meaningfully
wrong: the peak is 6.21, not the 6.3 assumed, and -- the one that mattered --
the "wall" value read where the curve meets the axis is 1.1, but their first
data point sits at y+ = 0.49 and the curve is still climbing, so the value at
OUR first cell (y+ = 0.75) is 1.27. Comparing our y+ = 0.75 cell against their
1.1 mixed up two different heights.

HOW, AND HOW IT IS CHECKED. The figure is vector art, but pdftocairo drops the
figure-5 paths (they do not survive the SVG conversion), so this renders the
page as pixels instead and separates the curves by COLOUR: the other three
series are saturated red/green/blue, the Conjug line is the only achromatic
dark ink inside the axes. Frames, tick marks and the legend rectangle are
masked out; the curve is then one connected vertical run per column.

Axes are calibrated on the TICK MARKS, not on the frame: the major ticks are
found as the long inward marks on an axis the curves do not touch (the top
frame for x, the left frame for y), and their uniform spacing is verified
before use.

THE VALIDATION IS FREE, AND IS THE POINT: figure 5 plots the same curve TWICE,
left on a log x-axis and right on a linear one. Both panels are digitised
independently -- different axis type, different calibration, different legend
mask -- and must agree. They do, to 0.019 max and 0.007 rms in <T'^2>, with
the peak landing at 6.208 in both. A calibration slip would not survive that.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile

import numpy as np

PAGE = 11                      # the page carrying figure 5
# panel frames, in pixels at 600 dpi (found by the long dark rows/columns)
LEFT = dict(XL=1328, XR=2567, YT=4276, YB=5182, legend=(4840, 5175, 1980, 2362))
RIGHT = dict(XL=2857, XR=4102, YT=4276, YB=5182, legend=(4568, 4899, 3575, 3955))
# y calibration is shared: ticks 129.4 px apart, y = 6 at px 4404.5, y = 1 at 5053.5
Y_PX, Y_VAL = (4404.5, 5053.5), (6.0, 1.0)


def render(pdf, dpi, out):
    subprocess.run(["pdftoppm", "-f", str(PAGE), "-l", str(PAGE), "-r", str(dpi),
                    "-png", pdf, out], check=True)
    return f"{out}-{PAGE}.png"


def black_mask(im, P):
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    chroma = (np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b))
    blk = (im.mean(2) < 110) & (chroma < 40)       # dark AND achromatic
    y0, y1, x0, x1 = P["legend"]
    blk[y0:y1, x0:x1] = False
    M = 16                                          # tick marks live here
    blk[:P["YT"] + M, :] = False
    blk[P["YB"] - M:, :] = False
    blk[:, :P["XL"] + M] = False
    blk[:, P["XR"] - M:] = False
    return blk


def ticks(dark, P, axis):
    """Major ticks as the long inward marks on the frame the curves avoid."""
    out, keep = [], []
    if axis == "x":
        for x in range(P["XL"] + 1, P["XR"]):
            col = dark[P["YT"] + 1:P["YT"] + 60, x]
            n = 0
            while n < col.size and col[n]:
                n += 1
            if n >= 10:
                keep.append(x)
    else:
        for y in range(P["YT"] + 1, P["YB"]):
            row = dark[y, P["XL"] + 1:P["XL"] + 60]
            n = 0
            while n < row.size and row[n]:
                n += 1
            if n >= 10:
                keep.append(y)
    for i in keep:
        if out and i - out[-1][-1] <= 3:
            out[-1].append(i)
        else:
            out.append([i])
    return [float(np.mean(g)) for g in out]


def trace(im, P, x_of_px):
    blk = black_mask(im, P)
    a = (Y_VAL[0] - Y_VAL[1]) / (Y_PX[0] - Y_PX[1])
    xs, vs = [], []
    for x in range(P["XL"] + 16, P["XR"] - 16):
        ys = np.flatnonzero(blk[:, x])
        if ys.size == 0:
            continue
        grp = max(np.split(ys, np.flatnonzero(np.diff(ys) > 3) + 1), key=len)
        xs.append(x_of_px(x))
        vs.append(Y_VAL[0] + (grp.mean() - Y_PX[0]) * a)
    return np.array(xs), np.array(vs)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--pdf", default=os.path.join(here, "..", "..", "..",
                                                  "literature", "flageul.pdf"))
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--out", default=os.path.join(here, "flageul_fig5a_conjug.dat"))
    a = ap.parse_args()

    from PIL import Image
    with tempfile.TemporaryDirectory() as td:
        png = render(a.pdf, a.dpi, os.path.join(td, "pg"))
        im = np.array(Image.open(png).convert("RGB")).astype(int)
    dark = im.mean(2) < 128

    # ---- left panel: log x. Decades from the three major ticks -------------
    tx = ticks(dark, LEFT, "x")
    assert len(tx) == 3, f"expected 3 x decades in fig 5a, got {tx}"
    d = np.diff(tx)
    assert np.allclose(d, d.mean(), rtol=0.02), f"decade spacing not uniform: {d}"
    dec = d.mean()
    xs, vs = trace(im, LEFT, lambda px: 10.0 ** ((px - tx[0]) / dec))

    # ---- right panel: linear x, the SAME curve -> an independent check -----
    tr = ticks(dark, RIGHT, "x")
    tr = [t for t in tr if RIGHT["XL"] + 20 < t < RIGHT["XR"] - 20]   # drop frames
    dr = np.diff(tr)
    assert np.allclose(dr, dr.mean(), rtol=0.02), f"tick spacing not uniform: {dr}"
    step = dr.mean()                       # px per 20 wall units (labels 20..140)
    xs2, vs2 = trace(im, RIGHT, lambda px: 20.0 + (px - tr[0]) * 20.0 / step)

    g = np.logspace(np.log10(max(xs.min(), xs2.min())),
                    np.log10(min(xs.max(), xs2.max())), 300)
    dev = np.interp(g, xs2, vs2) - np.interp(g, xs, vs)
    print(f"  panel 5a: {xs.size} columns, y+ {xs.min():.2f}..{xs.max():.0f}")
    print(f"  panel 5b: {xs2.size} columns, y+ {xs2.min():.1f}..{xs2.max():.0f}")
    print(f"  CROSS-CHECK 5a vs 5b: max {np.abs(dev).max():.3f}, "
          f"rms {np.sqrt(np.mean(dev ** 2)):.3f}  (peak {vs.max():.3f} vs {vs2.max():.3f})")
    assert np.abs(dev).max() < 0.06, "the two panels disagree -- calibration is wrong"

    np.savetxt(a.out, np.c_[xs, vs], fmt="%.6g", header=(
        "Flageul, Benhamadouche, Lamballais & Laurence, IJHFF 55 (2015) 34-44,\n"
        "figure 5a, the CONJUGATE case (black solid 'Conjug', G = G_2 = 1).\n"
        "Re_tau = 149, Pr = 0.71.  Digitised by digitize_flageul.py; validated\n"
        "against the same curve in panel 5b (linear x) to 0.02 in <T'^2>.\n"
        "y+   <T'^2>/T_tau^2"))
    print(f"  wrote {a.out}")
    print(f"  peak {vs.max():.3f} at y+ {xs[int(np.argmax(vs))]:.1f}; "
          f"value at y+ 0.75 = {np.interp(0.75, xs, vs):.3f}")


if __name__ == "__main__":
    main()
