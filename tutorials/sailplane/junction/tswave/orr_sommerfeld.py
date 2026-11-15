#!/usr/bin/env python3
"""Spatial Orr-Sommerfeld eigenvalues of the Blasius boundary layer.

The reference of the step-3 TS-wave test (docs/next_session_sailplane_junction.md
section 5.1): for a given Reynolds number Re = U delta*/nu and a real
frequency omega (delta*/U units) the complex streamwise wavenumber alpha of
the least-damped Tollmien-Schlichting mode, alpha_i < 0 = growing, from the
Orr-Sommerfeld equation under the parallel-flow assumption:

    (U - c)(phi'' - alpha^2 phi) - U'' phi = (phi'''' - 2 alpha^2 phi'' + alpha^4 phi)/(i alpha Re),
    c = omega/alpha,  phi = phi' = 0 at the wall,  phi -> 0 aloft.

Chebyshev collocation on y = y_max (1 - cos)/2... mapped algebraically so that
half the points sit inside y < L/2 (the wave lives in the layer), the quartic
polynomial eigenvalue problem in alpha linearised to a companion matrix
(scipy.linalg.eig), the physical branch picked as the eigenvalue nearest a
given guess or, without one, the one with alpha_r ~ omega/0.35 (the TS phase
speed) and the smallest |alpha_i|.

Checks (run the file): Jordinson (1970) Re 998, omega 0.1122 -> alpha =
0.30866 - 0.00573 i  (the classic spatial value).

    python3 orr_sommerfeld.py                      # the check
    python3 orr_sommerfeld.py --re 1000 --F 1.1e-4 # alpha at one point
    python3 orr_sommerfeld.py --F 1.1e-4 --re-range 500 1800 50   # alpha_i(Re) along the plate
"""
import argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eig


def blasius_profile(eta_max=12.0, n=4001):
    """f, f', f'' on a uniform eta grid; f''(0) by shooting."""
    def rhs(eta, y):
        return [y[1], y[2], -0.5 * y[0] * y[2]]
    def shoot(fpp0):
        s = solve_ivp(rhs, (0, eta_max), [0, 0, fpp0], rtol=1e-11, atol=1e-13)
        return s.y[1, -1] - 1.0
    a, b = 0.3, 0.4
    for _ in range(60):
        m = 0.5 * (a + b)
        if shoot(a) * shoot(m) <= 0:
            b = m
        else:
            a = m
    fpp0 = 0.5 * (a + b)
    eta = np.linspace(0, eta_max, n)
    s = solve_ivp(rhs, (0, eta_max), [0, 0, fpp0], t_eval=eta, rtol=1e-11, atol=1e-13)
    return eta, s.y[0], s.y[1], s.y[2], fpp0


def cheb(N):
    """Chebyshev differentiation matrix and points on [-1, 1] (Trefethen)."""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def os_spatial(re, omega, N=120, ymax=40.0, yhalf=4.0, guess=None, eigenfunction=False):
    """Least-damped spatial TS eigenvalue alpha (delta* units) at Re, omega.
    With eigenfunction=True also the (y, |u'|) shape, u' = phi' normalised to
    max 1, for the comparison with the measured amplitude profile."""
    D, xi = cheb(N)
    # algebraic map xi in [-1,1] -> y in [0, ymax], half the points below yhalf
    a = yhalf * ymax / (ymax - 2 * yhalf)
    b = 1 + 2 * a / ymax
    y = a * (1 + xi) / (b - xi)
    dydxi = a * (1 + b) / (b - xi) ** 2
    D1 = np.diag(1 / dydxi) @ D
    D2 = D1 @ D1
    D4 = D2 @ D2
    # base flow on y (delta* units): eta = y * eta_max_disp... delta*/theta scaling
    eta_t, f, fp, fpp, fpp0 = blasius_profile()
    disp = eta_t[-1] - f[-1]                       # delta* in eta units (1.7208)
    eta = y * disp
    U = np.interp(eta, eta_t, fp, right=1.0)
    Upp = np.interp(eta, eta_t, -0.5 * f * fpp, right=0.0) * disp ** 2
    I = np.eye(N + 1)
    U = np.diag(U); Upp = np.diag(Upp)
    # OS as a polynomial in alpha:  sum_k alpha^k A_k phi = 0
    # i alpha Re [(U - omega/alpha)(D2 - a^2) - Upp] phi = (D4 - 2 a^2 D2 + a^4) phi
    # multiply out: i Re [alpha U D2 - alpha^3 U - omega D2 + omega alpha^2 - alpha Upp] - D4 + 2 alpha^2 D2 - alpha^4 = 0
    A0 = -1j * re * omega * D2 - D4
    A1 = 1j * re * (U @ D2 - Upp)
    A2 = 1j * re * omega * I + 2 * D2
    A3 = -1j * re * U
    A4 = -I
    # boundary conditions: phi = phi' = 0 at wall (last index, xi = -1) and top (index 0)
    n = N + 1
    rows_bc = [0, 1, n - 2, n - 1]
    for A in (A0, A1, A2, A3, A4):
        A[rows_bc, :] = 0
    A0[0, :] = I[0, :]; A0[1, :] = D1[0, :]
    A0[n - 1, :] = I[n - 1, :]; A0[n - 2, :] = D1[n - 1, :]
    # companion linearisation of the quartic pencil
    Z = np.zeros((n, n)); 
    M = np.block([[Z, I, Z, Z], [Z, Z, I, Z], [Z, Z, Z, I], [-A0, -A1, -A2, -A3]])
    B = np.block([[I, Z, Z, Z], [Z, I, Z, Z], [Z, Z, I, Z], [Z, Z, Z, A4]])
    w, V = eig(M, B)
    ok = np.isfinite(w) & (w.real > 0) & (np.abs(w.imag) < 0.2) & (w.real < 2)
    w, V = w[ok], V[:, ok]
    if guess is not None:
        k = np.argmin(np.abs(w - guess))
    else:
        # TS branch: phase speed 0.2-0.5, least damped
        cand = np.where((omega / w.real > 0.2) & (omega / w.real < 0.5))[0]
        if len(cand) == 0:
            return np.nan
        k = cand[np.argmin(w[cand].imag)] if w[cand].imag.min() < 0 else cand[np.argmin(np.abs(w[cand].imag))]
    if not eigenfunction:
        return w[k]
    phi = V[:n, k]                      # the first companion block is phi itself
    up = np.abs(D1 @ phi); up /= up.max()
    return w[k], y[::-1], up[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--re', type=float); ap.add_argument('--omega', type=float)
    ap.add_argument('--F', type=float, help='omega/Re (the dimensional frequency parameter)')
    ap.add_argument('--re-range', nargs=3, type=float, metavar=('RE0', 'RE1', 'DRE'))
    ap.add_argument('--N', type=int, default=160)
    args = ap.parse_args()
    if args.re_range:
        re0, re1, dre = args.re_range
        print('# Re_dstar  alpha_r  alpha_i   (F = %.3e; alpha in delta* units; alpha_i < 0 grows)' % args.F)
        guess = None
        for re in np.arange(re0, re1 + 0.5 * dre, dre):
            al = os_spatial(re, args.F * re, N=args.N, guess=guess)
            guess = al
            print('%7.0f  %.5f  %+.5f' % (re, al.real, al.imag))
        return
    if args.re:
        om = args.omega if args.omega else args.F * args.re
        al = os_spatial(args.re, om, N=args.N)
        print('Re %.1f omega %.5f: alpha = %.5f %+.5f i' % (args.re, om, al.real, al.imag))
        return
    for N in (120, 160, 200):
        al = os_spatial(998.0, 0.1122, N=N)
        print('Jordinson check N=%d: alpha = %.5f %+.5f i  (reference 0.30866 -0.00573 i)' % (N, al.real, al.imag))


if __name__ == '__main__':
    main()
