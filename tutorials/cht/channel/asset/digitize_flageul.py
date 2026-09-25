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


# The four series are separable by INK COLOUR alone. Conjug is the only
# achromatic dark ink inside the axes; the others are saturated primaries
# (measured in the render: exactly (0,128,0), (0,0,255), (255,0,0)).
SERIES = {
    # tol = the 5a-vs-5b agreement each series is required to meet. A 2-px
    # LINE is located to a fraction of a pixel; a SYMBOL centroid is much
    # coarser, and its bias differs between the panels because the same ~20 px
    # marker spans very different y+ widths on a log and on a linear axis. The
    # measured disagreement is therefore spread over the steep regions rather
    # than being one outlier -- checked, not assumed -- so the brackets are
    # good to ~0.2-0.4 in <T'^2> where the conjugate line is good to 0.015.
    "conjug": dict(kind="line",   label="Conjug (black solid)", tol=0.05, floor=0.0),
    "isoq":   dict(kind="marker", label="isoQ (green x)",       tol=0.45, floor=0.0),
    # isoT is UNRELIABLE below ~0.3: a '+' spans about +-0.08 in <T'^2>, so
    # where the curve is that small the symbol straddles the axis and its lower
    # half is CLIPPED IN THEIR FIGURE. That ink does not exist in the image and
    # no estimator recovers it; the digitised values there go non-monotone
    # (0.094, 0.079, 0.122 at y+ 0.5, 0.75, 1.0) instead of following the y^2
    # an ideal Dirichlet wall must. Use the analytic limit instead: isoT -> 0
    # at the wall, exactly, by definition of the boundary condition.
    "isot":   dict(kind="marker", label="isoT (blue +)",        tol=0.45, floor=0.3),
}


def ink_mask(im, P, which):
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    if which == "conjug":
        chroma = (np.maximum(np.maximum(r, g), b)
                  - np.minimum(np.minimum(r, g), b))
        m = (im.mean(2) < 110) & (chroma < 40)      # dark AND achromatic
    elif which == "isoq":
        m = (g > r + 40) & (g > b + 40)
    elif which == "isot":
        m = (b > r + 60) & (b > g + 60)
    else:
        raise ValueError(which)
    y0, y1, x0, x1 = P["legend"]                    # the legend carries samples
    m[y0:y1, x0:x1] = False
    # The margin exists to drop the frame and its TICK MARKS, which are black.
    # A colour mask already excludes them, so a coloured series needs only to
    # be clipped to the panel -- and MUST NOT get the wide margin: it would eat
    # the lower half of every symbol sitting near the axis and bias the
    # centroid UPWARD exactly where the isoT curve is smallest. (Symptom: the
    # digitised isoT went flat, 0.149 -> 0.163, where an ideal Dirichlet wall
    # must vary like y^2.)
    M = 16 if which == "conjug" else 2
    m[:P["YT"] + M, :] = False
    m[P["YB"] - M:, :] = False
    m[:, :P["XL"] + M] = False
    m[:, P["XR"] - M:] = False
    return m


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


def trace(im, P, x_of_px, which="conjug"):
    """Column-wise read of one series.

    A LINE is one connected vertical run, so take the longest (stray ink from a
    crossing series is a shorter run). MARKERS (x and +) are symmetric about
    their data point, so the column-wise CENTROID of the ink sits on the curve
    -- and it keeps doing so where neighbouring markers overlap, which they do
    densely near the peak, whereas a connected-run rule would merge them and
    read the union's extent instead.
    """
    msk = ink_mask(im, P, which)
    a = (Y_VAL[0] - Y_VAL[1]) / (Y_PX[0] - Y_PX[1])
    kind = SERIES[which]["kind"]
    xs, vs = [], []
    for x in range(P["XL"] + 16, P["XR"] - 16):
        ys = np.flatnonzero(msk[:, x])
        if ys.size == 0:
            continue
        if kind == "line":
            ys = max(np.split(ys, np.flatnonzero(np.diff(ys) > 3) + 1), key=len)
        elif ys.size < 3:
            continue                                # antialias speck, not a marker
        xs.append(x_of_px(x))
        vs.append(Y_VAL[0] + (ys.mean() - Y_PX[0]) * a)
    return np.array(xs), np.array(vs)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--pdf", default=os.path.join(here, "..", "..", "..",
                                                  "literature", "flageul.pdf"))
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--outdir", default=here)
    a = ap.parse_args()

    from PIL import Image
    with tempfile.TemporaryDirectory() as td:
        png = render(a.pdf, a.dpi, os.path.join(td, "pg"))
        im = np.array(Image.open(png).convert("RGB")).astype(int)
    dark = im.mean(2) < 128

    # ---- axis calibration, once, from the tick marks -----------------------
    tx = ticks(dark, LEFT, "x")
    assert len(tx) == 3, f"expected 3 x decades in fig 5a, got {tx}"
    d = np.diff(tx)
    assert np.allclose(d, d.mean(), rtol=0.02), f"decade spacing not uniform: {d}"
    dec = d.mean()
    tr = [t for t in ticks(dark, RIGHT, "x")
          if RIGHT["XL"] + 20 < t < RIGHT["XR"] - 20]          # drop the frames
    dr = np.diff(tr)
    assert np.allclose(dr, dr.mean(), rtol=0.02), f"tick spacing not uniform: {dr}"
    step = dr.mean()                     # px per 20 wall units (labels 20..140)

    for name, spec in SERIES.items():
        xs, vs = trace(im, LEFT, lambda px: 10.0 ** ((px - tx[0]) / dec), name)
        xs2, vs2 = trace(im, RIGHT,
                         lambda px: 20.0 + (px - tr[0]) * 20.0 / step, name)
        # THE CHECK: the same series, digitised off two panels whose axes have
        # nothing in common but the y scale.
        g = np.logspace(np.log10(max(xs.min(), xs2.min())),
                        np.log10(min(xs.max(), xs2.max())), 300)
        dev = np.interp(g, xs2, vs2) - np.interp(g, xs, vs)
        mx, rms = float(np.abs(dev).max()), float(np.sqrt(np.mean(dev ** 2)))
        out = os.path.join(a.outdir, f"flageul_fig5a_{name}.dat")
        np.savetxt(out, np.c_[xs, vs], fmt="%.6g", header=(
            "Flageul, Benhamadouche, Lamballais & Laurence, IJHFF 55 (2015) 34-44,\n"
            f"figure 5a, the {spec['label']} case.  Re_tau = 149, Pr = 0.71.\n"
            "Digitised by digitize_flageul.py and cross-validated against the same\n"
            f"series in panel 5b: max {mx:.3f}, rms {rms:.3f} in <T'^2>.\n"
            + (f"UNRELIABLE below <T'^2> = {spec['floor']}: the symbol is clipped by\n"
               "their own axis there. The exact limit is 0 at the wall.\n"
               if spec["floor"] else "")
            + "y+   <T'^2>/T_tau^2"))
        print(f"  {spec['label']:<22} {xs.size:5d} cols  y+ {xs.min():6.2f}..{xs.max():5.0f}"
              f"   peak {vs.max():6.3f} at y+ {xs[int(np.argmax(vs))]:5.1f}"
              f"   at y+0.75 {np.interp(0.75, xs, vs):6.3f}"
              f"   | 5a-vs-5b max {mx:.3f} rms {rms:.3f}")
        assert mx < spec["tol"], (f"{name}: the two panels disagree by {mx:.3f} "
                                  f"(> {spec['tol']}) -- calibration is wrong")
    # keep the historical filename for the conjugate case
    src = os.path.join(a.outdir, "flageul_fig5a_conjug.dat")
    print(f"  wrote {src} and the two brackets")


if __name__ == "__main__":
    main()
