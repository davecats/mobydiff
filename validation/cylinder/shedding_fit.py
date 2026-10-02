#!/usr/bin/env python3
"""Shedding frequency of a cylinder force trace (README, "Which Strouhal
number is this domain's").

    ./shedding_fit.py forces_a.txt [forces_b.txt ...]

St is the maximum of the Hann-windowed Fourier amplitude of C_L over a
CONTINUOUS frequency in [0.14, 0.20], over the last 50 time units and over
the 50 before them (the two agree once the shedding has saturated). It is
the fundamental whatever else the trace carries, and it resolves 1e-4 where
the FFT bin of a 50-unit window is 0.02. The strongest line in [0.21, 0.60]
is printed too: on the scheme before 2026-10-01 the lift carried one, at a
frequency that depended on the domain.
"""
import sys
import numpy as np

def line(tt, cc, lo, hi):
    w = np.hanning(len(cc))
    fs = np.linspace(lo, hi, 6001)
    amp = np.array([abs(np.sum(w*cc*np.exp(-2j*np.pi*f*tt))) for f in fs])*2/np.sum(w)
    k = int(np.argmax(amp))
    return fs[k], amp[k]

for p in sys.argv[1:]:
    win = 50.0
    d = np.loadtxt(p, skiprows=1)
    t, cl, cd = d[:, 1], d[:, 2], d[:, 3]
    m = t >= t[-1] - win
    tt, cc = t[m], cl[m] - cl[m].mean()
    f1, a1 = line(tt, cc, 0.14, 0.20)
    f2, a2 = line(tt, cc, 0.21, 0.60)
    m2 = (t >= t[-1] - 2*win) & (t < t[-1] - win)
    f1b, a1b = line(t[m2], cl[m2] - cl[m2].mean(), 0.14, 0.20) if m2.sum() > 100 else (np.nan, np.nan)
    print(f"{p.split('/')[-1][:30]:30s} St {f1:.4f} amp {a1:.3f} | previous window St {f1b:.4f} amp {a1b:.3f} | "
          f"second line f {f2:.3f} amp {a2:.3f} | mean C_D {cd[m].mean():.4f}  C_L range {cl[m].min():+.2f}..{cl[m].max():+.2f}")
