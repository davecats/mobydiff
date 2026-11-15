#!/usr/bin/env python3
"""Lift history of the c/2048 calibration box from its strip read-outs.

Reads every ``l7_<step>_strips.csv`` written by ``strip_loads.py`` (one per
snapshot, fetched to ``~/horeka_stage/junction``), prints the trend table of
the calib README (area-weighted c_l, c_dp, c_m, root strip, outboard strips
y > 1.5 c, the mid-span separation point), fits the relaxation
``a - b exp(-(t - t0)/tau)`` to the area-weighted c_l and draws
``viz/l7_cl_history.png``.  The time of each snapshot comes from the chunk
logs' trace lines (the generic case has no runtime force output) and is
tabulated here.

    python3 cl_history.py [--dir ~/horeka_stage/junction] [--out viz/l7_cl_history.png]
"""
import argparse, glob, os, re
import numpy as np
from scipy.optimize import curve_fit

# step -> t of the snapshot (trace lines of the chunk logs; 154547 = t_final)
T = {61925: 16.05, 70000: 16.37, 80000: 16.77, 90000: 17.17, 100000: 17.57,
     110000: 18.005, 120000: 18.445, 127000: 18.757,
     130000: 18.892, 140000: 19.343, 150000: 19.795, 154547: 20.0}
T0 = 16.37            # fit origin (the first snapshot of the continuation chain)
REF = {'c/512 wall-function box': 0.382, '2D c/1024 resolved twin': 0.479,
       'XFOIL turbulent': 0.79}
PALETTE = ['#2a78d6', '#eb6834', '#8a63d2', '#1baf7a']


def relax(t, a, b, tau):
    return a - b * np.exp(-(t - T0) / tau)


def read_strips(path):
    rows = np.loadtxt(path, delimiter=',', comments='#')
    y0, A, cl, cdp, cm = rows[:, 1], rows[:, 4], rows[:, 5], rows[:, 6], rows[:, 8]
    w = A / A.sum()
    out = w @ cl, w @ cdp, w @ cm, cl[0]
    ob = y0 > 1.5
    mid = np.argmin(np.abs(0.5 * (rows[:, 1] + rows[:, 2]) - 2.0))
    return out + (A[ob] @ cl[ob] / A[ob].sum(), rows[mid, 9])   # x_sep1: the separation line (x_sep2 is the trailing edge)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default=os.path.expanduser('~/horeka_stage/junction'))
    ap.add_argument('--out', default=os.path.join(os.path.dirname(__file__), 'viz', 'l7_cl_history.png'))
    args = ap.parse_args()
    files = {int(re.search(r'l7_(\d+)_strips', f).group(1)): f
             for f in glob.glob(os.path.join(args.dir, 'l7_*_strips.csv'))}
    steps = sorted(s for s in files if s in T)
    missing = sorted(s for s in files if s not in T)
    if missing:
        print('no time tabulated for steps', missing, '-- add them to T from the trace lines')
    t = np.array([T[s] for s in steps])
    data = np.array([read_strips(files[s]) for s in steps])
    cl, cdp, cm, root, ob, xsep = data.T

    print('| step | t | c_l | dc_l/dt | c_dp | c_m | root strip | outboard y > 1.5 c | separation (mid-span) |')
    print('|---|---|---|---|---|---|---|---|---|')
    for i, s in enumerate(steps):
        slope = '' if i == 0 else '%+.3f' % ((cl[i] - cl[i - 1]) / (t[i] - t[i - 1]))
        print('| %d | %.2f | %.3f | %s | %.4f | %.3f | %.3f | %.3f | %.2f c |'
              % (s, t[i], cl[i], slope, cdp[i], cm[i], root[i], ob[i], xsep[i]))

    fits = {}
    for name, sel in (('all points from t %.2f' % T0, t >= T0 - 1e-9), ('last four', np.arange(len(t)) >= len(t) - 4)):
        if sel.sum() < 4:
            continue
        p, cov = curve_fit(relax, t[sel], cl[sel], p0=(0.5, 0.15, 1.5))
        fits[name] = (p, np.sqrt(np.diag(cov)))
        print('fit %-24s: asymptote %.3f +- %.3f, tau %.2f +- %.2f c/U, c_l(20) %.3f, c_l(21) %.3f'
              % (name, p[0], np.sqrt(cov[0, 0]), p[2], np.sqrt(cov[2, 2]), relax(20.0, *p), relax(21.0, *p)))
    sel = t >= T0 - 1e-9
    p, cov = curve_fit(relax, t[sel], ob[sel], p0=(0.5, 0.15, 1.5))
    print('fit outboard strips      : asymptote %.3f +- %.3f, tau %.2f c/U' % (p[0], np.sqrt(cov[0, 0]), p[2]))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    tt = np.linspace(T0, 21.0, 200)
    for (name, (p, _)), col, ls in zip(fits.items(), PALETTE[1:3], ('-', '--')):
        ax.plot(tt, relax(tt, *p), ls, color=col, lw=1.6,
                label='fit %s: %.3f, tau %.2f' % (name, p[0], p[2]))
    ax.plot(t, cl, 'o', color=PALETTE[0], ms=5, label='area-weighted c_l (strips)')
    ax.plot(t, ob, 's', color=PALETTE[3], ms=4, mfc='none', label='outboard strips y > 1.5 c')
    for (name, v), dy in zip(REF.items(), (-0.012, -0.012, -0.012)):
        if v < 0.6:
            ax.axhline(v, color='#888', lw=0.8, ls=':')
            ax.text(16.1, v + dy, name + ' %.3f' % v, fontsize=7, color='#555', va='top')
    ax.axvline(20.0, color='#bbb', lw=0.8)
    ax.text(20.02, 0.59, 't_final 20', fontsize=7, color='#555', va='top')
    ax.set_xlim(16.0, 21.0); ax.set_ylim(0.34, 0.60)
    ax.set_xlabel('t  (c/U)'); ax.set_ylabel('c_l of the box wing')
    ax.set_title('c/2048 calibration box: lift relaxation (XFOIL turbulent 0.79 off scale)', fontsize=9)
    ax.grid(alpha=0.25); ax.legend(loc='upper left', fontsize=7)
    fig.tight_layout(); fig.savefig(args.out)
    print('wrote', args.out)


if __name__ == '__main__':
    main()
