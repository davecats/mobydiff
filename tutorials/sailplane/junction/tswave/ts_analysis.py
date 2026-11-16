#!/usr/bin/env python3
"""Read-out of the TS-wave test: the immersed-plate Blasius base flow against
Blasius, and the harmonic disturbance's streamwise amplitude A(x) against the
spatial Orr-Sommerfeld growth (orr_sommerfeld.py, parallel flow).

    python3 ts_analysis.py --base tsbase_N.h5 --forced tsw_*.h5 --ini tswave.ini [--out viz/]

The forced snapshots are K phases of ONE forcing period (field_interval =
period/K steps); the disturbance u' = u_k - u_base is projected on the forcing
frequency, A(x, y) = |2/K sum_k u'_k exp(-i omega t_k)|, and its maximum over
y is the amplitude whose log-derivative is the measured -alpha_i(x).  Lengths
are inlet momentum thicknesses (the ini's units), Re = Re_theta,in.
"""
import argparse, configparser, glob, os, sys
import numpy as np
import h5py
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'tools'))
from compare_fields import load_field
import orr_sommerfeld as osm

BETA = 0.664112           # 2 f''(0): theta = beta * sqrt(nu x / U) ... Re_theta = beta sqrt(Re_x)
DISP = 1.720788           # delta*/sqrt(nu x/U)


def read_ini(path):
    cp = configparser.ConfigParser(inline_comment_prefixes=(';',), strict=False)
    cp.read(path)
    g = lambda s, k, d=None: cp.get(s, k, fallback=d)
    return dict(re=float(g('flow', 're')), y0=float(g('boundary', 'blasius_y0', 0.0)),
                omega=float(g('force', 'trip_omega', 0.0)), x0=float(g('force', 'trip_x0', 0.0)))


def base_flow(path, y0):
    with h5py.File(path, 'r') as f:
        u = load_field(f, 'un')[:, :, :].mean(axis=0)      # (ny, nx) spanwise mean, u at x faces
        x = f['x'][...]; y = f['y'][...]
        t = float(f.attrs['t_current'])
    yc = 0.5 * (y[1:] + y[:-1])                             # u rows sit at cell centres in y
    xf = x[1:]                                              # u columns at the faces x_{i+1}
    return u, xf, yc, t


def integral_thicknesses(u, yc, y0):
    """delta*, theta, H per column, integrating above the plate to the edge."""
    ds, th = np.zeros(u.shape[1]), np.zeros(u.shape[1])
    dy = np.gradient(yc)
    for i in range(u.shape[1]):
        ui = u[:, i]; ue = ui[yc > y0].max()
        sel = yc > y0
        f1 = np.where(sel, 1 - ui / ue, 0.0); f2 = np.where(sel, ui / ue * (1 - ui / ue), 0.0)
        ds[i] = (f1 * dy).sum(); th[i] = (f2 * dy).sum()
    return ds, th, ds / th


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', required=True); ap.add_argument('--forced', nargs='*', default=[])
    ap.add_argument('--ini', required=True); ap.add_argument('--out', default='viz')
    ap.add_argument('--xlim', nargs=2, type=float, default=(100.0, None),
                    help='x range of the read-out (default: 100 to 85 %% of lx)')
    args = ap.parse_args()
    p = read_ini(args.ini)
    re_t = p['re']; y0 = p['y0']
    xv = re_t / BETA ** 2                                   # virtual origin upstream of the inlet
    ub, xf, yc, tb = base_flow(args.base, y0)
    lx = xf[-1]
    x1 = args.xlim[1] or 0.85 * lx
    sel = (xf >= args.xlim[0]) & (xf <= x1)
    ds, th, H = integral_thicknesses(ub, yc, y0)
    rex = re_t * (xf + xv)
    th_bl = BETA * np.sqrt(xf + xv) / np.sqrt(re_t) * 1.0   # theta = beta sqrt(nu (x+xv)/U), nu = 1/Re_theta
    ds_bl = DISP * np.sqrt((xf + xv) / re_t)
    os.makedirs(args.out, exist_ok=True)
    print('BASE FLOW (t %.1f): worst over x in [%.0f, %.0f]: theta %.2f %%, delta* %.2f %%, H %.3f..%.3f (Blasius 2.591)'
          % (tb, args.xlim[0], x1, 100 * np.abs(th[sel] / th_bl[sel] - 1).max(),
             100 * np.abs(ds[sel] / ds_bl[sel] - 1).max(), H[sel].min(), H[sel].max()))
    np.savetxt(os.path.join(args.out, 'base_thicknesses.csv'),
               np.c_[xf, rex, ds, ds_bl, th, th_bl, H], delimiter=',',
               header='x,Re_x,delta*,delta*_Blasius,theta,theta_Blasius,H', comments='')
    if not args.forced:
        return

    omega = p['omega']
    ts, up = [], []
    for fp in sorted(args.forced):
        with h5py.File(fp, 'r') as f:
            u = load_field(f, 'un').mean(axis=0); ts.append(float(f.attrs['t_current']))
        up.append(u - ub)
    ts = np.array(ts); up = np.array(up)                   # (K, ny, nx)
    # the phases must be equispaced over one period: a run whose restart step
    # is not a multiple of field_interval also writes its final step, off
    # the lattice -- drop such a file rather than bias the projection
    d = np.diff(ts)
    while len(ts) > 2 and abs(d[-1] - np.median(d)) > 0.01 * np.median(d):
        print('dropping the off-lattice phase file at t %.1f (spacing %.2f against %.2f)' % (ts[-1], d[-1], np.median(d)))
        ts, up, d = ts[:-1], up[:-1], d[:-1]
    K = len(ts)
    print('phases: K %d, spacing %.2f, covering %.3f of the period 2 pi / omega' % (K, d.mean(), K * d.mean() * omega / (2 * np.pi)))
    ph = np.exp(-1j * omega * ts)[:, None, None]
    A = np.abs(2.0 / K * (up * ph).sum(axis=0))             # (ny, nx) amplitude at omega
    Amax = A.max(axis=0)                                    # over y
    lnA = np.log(np.maximum(Amax, 1e-300))
    # OS at the local Re_delta* on every 40th column (alpha per delta*,
    # converted to the theta units of x); interpolated to every column after
    re_ds = DISP * np.sqrt(rex)
    out = []
    guess = None
    for i in np.where(sel)[0][::40]:
        al = osm.os_spatial(re_ds[i], omega * ds_bl[i], N=120, guess=guess); guess = al
        out.append((xf[i], re_ds[i], -al.imag / ds_bl[i], al.real / ds_bl[i], Amax[i]))
    out = np.array(out)
    ai = np.interp(xf, out[:, 0], out[:, 2]); ar = np.interp(xf, out[:, 0], out[:, 3])
    # the measured growth rate: the least-squares slope of ln A over a window
    # of +- one local TS wavelength (2 pi / alpha_r, ~ 90 theta_in), which
    # averages the wavelength-scale ripple of the y-maximum
    gr = np.full_like(lnA, np.nan)
    for i in np.where(sel)[0]:
        lam = 2 * np.pi / ar[i]
        w = (xf >= xf[i] - lam) & (xf <= xf[i] + lam)
        gr[i] = np.polyfit(xf[w], lnA[w], 1)[0]
    # the reference station: the amplitude minimum behind the forcing's near
    # field (the non-modal part decays within ~60 theta_in of the ribbon)
    near = (xf > p['x0'] + 20) & (xf < p['x0'] + 200)
    i0 = np.where(near)[0][np.argmin(Amax[near])]
    N_meas = lnA - lnA[i0]
    N_os = np.concatenate([[0], np.cumsum(0.5 * (ai[1:] + ai[:-1]) * np.diff(xf))]); N_os -= N_os[i0]
    np.savetxt(os.path.join(args.out, 'amax.csv'), np.c_[xf, re_ds, Amax, gr, ai, N_meas, N_os][sel],
               delimiter=',', header='x,Re_dstar,A_max,dlnA_dx_fit,minus_alpha_i_OS,N_measured,N_OS',
               comments='')
    print('reference station x %.0f (Re_delta* %.0f), A %.3e; forcing amplitude %s'
          % (xf[i0], re_ds[i0], Amax[i0], 'trip_amp'))
    print('Re_delta*    x   -alpha_i(OS)  dlnA/dx(fit)   N_OS  N_meas')
    for r in range(800, 1401, 50):
        if r > re_ds[sel].max() - 20: break                 # the read-out range ends here
        i = np.argmin(np.abs(re_ds - r))
        print('%8.0f %6.0f   %+.5f     %+.5f   %6.3f %6.3f' % (r, xf[i], ai[i], gr[i], N_os[i], N_meas[i]))
    far = sel & (xf > xf[i0] + 200)
    if ai[far].min() < 0:                                   # branch II inside the range
        i2 = np.where(far)[0][np.argmax(Amax[far])]
        ib = np.where(far)[0][np.argmin(np.abs(ai[far]))]
        print('branch II: measured A_max peak at x %.0f (Re_delta* %.0f), OS neutral at x %.0f (Re_delta* %.0f)'
              % (xf[i2], re_ds[i2], xf[ib], re_ds[ib]))
        print('N at the peak: measured %.3f, OS %.3f at the same x; OS max %.3f' % (N_meas[i2], N_os[i2], N_os[far].max()))
    r1 = min(1300.0, re_ds[sel].max() - 20)
    band = sel & (re_ds > 850) & (re_ds < r1)
    print('growth rate, Re_delta* 850..%.0f: measured / OS mean ratio %.3f' % (r1, gr[band].sum() / ai[band].sum()))
    # the amplitude profile |u'|(y) at Re_delta* 1000 against the OS eigenfunction
    i1 = np.argmin(np.abs(re_ds - 1000.0))
    al1, yos, uos = osm.os_spatial(re_ds[i1], omega * ds_bl[i1], N=120, eigenfunction=True)
    prof = A[:, i1] / A[:, i1].max(); yprof = (yc - y0) / ds_bl[i1]
    np.savetxt(os.path.join(args.out, 'eigenfunction_re1000.csv'), np.c_[yprof, prof], delimiter=',',
               header='y_over_dstar,u_amplitude_normalised', comments='')
    uos_i = np.interp(yprof, yos, uos)
    print('eigenfunction at Re_delta* %.0f: max |measured - OS| %.3f of the peak, inner peak at y/delta* %.2f (OS %.2f)'
          % (re_ds[i1], np.abs(prof - uos_i)[(yprof > 0) & (yprof < 8)].max(), yprof[np.argmax(prof)], yos[np.argmax(uos)]))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4), dpi=140)
    ax[2].plot(prof[yc > y0], yprof[yc > y0], 'o', color='#2a78d6', ms=3, label='measured |u\'|, x %.0f' % xf[i1])
    ax[2].plot(uos, yos, '-', color='#eb6834', lw=1.6, label='Orr-Sommerfeld, Re_delta* %.0f' % re_ds[i1])
    ax[2].set_ylim(0, 8); ax[2].set_xlabel('|u\'| / max'); ax[2].set_ylabel('(y - y_plate) / delta*')
    ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
    s0 = sel & (xf >= xf[i0])
    ax[0].plot(xf[s0], N_meas[s0], color='#2a78d6', lw=1.8, label='measured ln(A/A_ref), u\' max over y')
    ax[0].plot(xf[s0], N_os[s0], '--', color='#eb6834', lw=1.6, label='Orr-Sommerfeld (parallel), integrated')
    ax[0].set_xlabel('x (theta_in)'); ax[0].set_ylabel('N from x %.0f' % xf[i0]); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    ax[1].plot(re_ds[s0], gr[s0], color='#2a78d6', lw=1.8, label='measured d ln A / dx (fit over +- one wavelength)')
    ax[1].plot(re_ds[s0], ai[s0], '--', color='#eb6834', lw=1.6, label='-alpha_i, Orr-Sommerfeld')
    ax[1].axhline(0, color='#888', lw=.8); ax[1].set_xlabel('Re_delta*'); ax[1].set_ylabel('growth rate (1/theta_in)')
    ax[1].set_ylim(-0.002, 0.003); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    fig.suptitle('TS wave through the immersed plate, omega %.4f, K = %d phases' % (omega, K), fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(args.out, 'ts_growth.png')); print('wrote', args.out)


if __name__ == '__main__':
    main()
