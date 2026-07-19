#!/usr/bin/env python3
"""rms of the u wiggle (deviation from z-mean-smooth) at y~23 per snapshot."""
import sys, h5py, numpy as np
sys.path.insert(0, "../../tools")
from compare_fields import load_field
for fn in sys.argv[1:]:
    with h5py.File(fn, "r") as f:
        u = load_field(f, "un").mean(axis=0)
        t = float(f.attrs.get("t_current", np.nan)) if "t_current" in f.attrs else np.nan
    row = u[60, :] - 1.0
    # 2-dx content: difference against 3-point smoothed row
    sm = np.convolve(row, [0.25, 0.5, 0.25], mode="same")
    hi = row - sm
    print(f"{fn:24s} t={t:8.2f} rms(total)={np.sqrt((row**2).mean()):.3e} "
          f"rms(2dx)={np.sqrt((hi**2).mean()):.3e}")
