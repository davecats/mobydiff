#!/usr/bin/env python3
"""Checkers for the passive-scalar (S1) gates.

Subcommands
  uniform    max|s - c| over every STORED cell (no reassembly, so coarse and
             fine cells are each checked once) + the leaf level histogram.
  conserve   volume integral of a scalar in two snapshots (global
             conservation: d/dt int s dV = 0 in a periodic box).
  wave       error against the analytic advection-diffusion solution
             s = A exp(-D |k|^2 t) sin(k.(x - u t)) of a uniform flow.
  parabola   error against the steady conduction solution with a constant
             volumetric source between Dirichlet walls.
  profile    y profile (x,z average) of a scalar, optionally against the
             transient conduction series for the Pr sweep.
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np

from scalar_tools import BlockGeometry, integrate, volume, field_error


def cmd_uniform(a) -> int:
    with h5py.File(a.h5, "r") as f:
        levels = f["blocks"][:, 3]
        hist = np.bincount(levels, minlength=a.levels).tolist()
        dev = {}
        for spec in a.scalar:
            name, value = spec.split("=")
            dev[name] = float(np.max(np.abs(f[name][...] - float(value))))
    print(f"{a.h5}: {levels.size} leaves, level histogram {hist}")
    for name, d in dev.items():
        print(f"max|{name} - const| = {d:.3e}")
    ok = all(d == 0.0 for d in dev.values()) and len([h for h in hist if h > 0]) >= a.levels
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_conserve(a) -> int:
    totals, times, vols = [], [], []
    for path in (a.first, a.second):
        with h5py.File(path, "r") as f:
            totals.append(integrate(f, a.scalar))
            vols.append(volume(f))
            times.append(float(f.attrs["t_current"]))
            absmax = float(np.max(np.abs(f[a.scalar][...])))
    scale = abs(absmax) * vols[-1]
    drift = (totals[1] - totals[0]) / scale
    print(f"volume {vols[0]:.12e} (both files {vols[0] == vols[1]})")
    print(f"int s dV: {totals[0]:.16e} (t = {times[0]:g})"
          f" -> {totals[1]:.16e} (t = {times[1]:g})")
    print(f"relative drift (per max|s| V) = {drift:.3e}")
    ok = abs(drift) <= a.tolerance
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_wave(a) -> int:
    t0 = 0.0
    if a.ic is not None:
        with h5py.File(a.ic, "r") as f0:
            t0 = float(f0.attrs["t_current"])
    with h5py.File(a.h5, "r") as f:
        geo = BlockGeometry(f)
        # Elapsed time SINCE the manufactured IC was imposed (the seed run
        # leaves a nonzero t_current in the IC file).
        t = float(f.attrs["t_current"]) - t0
        k = np.array([2.0 * np.pi * n / L for n, L in zip(a.wave, geo.leng)])
        u = np.array(a.velocity, dtype=float)
        decay = np.exp(-a.diffusivity * float(k @ k) * t)
        phase = float(k @ u) * t

        def exact(x, y, z):
            return a.amp * decay * np.sin(k[0] * x + k[1] * y + k[2] * z - phase)

        l2, linf = field_error(f, a.scalar, exact)
        amp_t = a.amp * decay
    print(f"{a.h5}: t = {t:g}, analytic amplitude {amp_t:.8e}")
    print(f"L2 = {l2:.6e}   Linf = {linf:.6e}   L2/amp = {l2/amp_t:.6e}")
    return 0


def cmd_parabola(a) -> int:
    with h5py.File(a.h5, "r") as f:
        L = float(f.attrs["ly"])

        def exact(x, y, z):
            return a.s0 + (a.s1 - a.s0) * y / L \
                + a.source / (2.0 * a.diffusivity) * y * (L - y)

        l2, linf = field_error(f, a.scalar, exact)
    print(f"{a.h5}: L2 = {l2:.6e}   Linf = {linf:.6e}")
    return 0


def y_profile(f, name):
    """(y, s) profile of a cell-centred scalar, x/z-averaged, single level."""
    geo = BlockGeometry(f)
    data = f[name]
    acc = {}
    for bid in range(geo.n_blocks):
        (xc, dx), (yc, dy), (zc, dz) = geo.block_axes(bid)
        block = data[bid]
        for jj, yv in enumerate(yc):
            key = round(float(yv), 12)
            s = float(np.mean(block[:, jj, :]))
            n = block.shape[0] * block.shape[2]
            tot, cnt = acc.get(key, (0.0, 0))
            acc[key] = (tot + s * n, cnt + n)
    ys = np.array(sorted(acc))
    ss = np.array([acc[y][0] / acc[y][1] for y in ys])
    return ys, ss


def conduction_series(y, t, D, L, nterm=400):
    """s(y,t) for s(0)=s(L)=1, s(y,0)=0 (transient conduction between walls)."""
    s = np.ones_like(y)
    for n in range(1, 2 * nterm, 2):
        s -= (4.0 / (n * np.pi)) * np.sin(n * np.pi * y / L) \
            * np.exp(-D * (n * np.pi / L) ** 2 * t)
    return s


def cmd_profile(a) -> int:
    with h5py.File(a.h5, "r") as f:
        y, s = y_profile(f, a.scalar)
        t = float(f.attrs["t_current"])
        L = float(f.attrs["ly"])
    if a.diffusivity is None:
        for yy, ss in zip(y, s):
            print(f"{yy:.8f} {ss:.10e}")
        return 0
    ref = conduction_series(y, t, a.diffusivity, L)
    err = np.max(np.abs(s - ref))
    mid = np.argmin(np.abs(y - 0.5 * L))
    # Wall heat flux (Nusselt-like): D ds/dy at the wall, discrete one-sided
    # against the analytic series derivative.
    q_num = a.diffusivity * (1.0 - s[0]) / y[0]
    q_ref = a.diffusivity * (1.0 - ref[0]) / y[0]
    print(f"{a.h5}: t = {t:g}  D = {a.diffusivity:g}  ny = {y.size}")
    print(f"  centre s = {s[mid]:.6f} (analytic {ref[mid]:.6f}), "
          f"max|s - series| = {err:.3e}")
    print(f"  wall flux D(s_w - s_1)/y_1 = {q_num:.6e} (analytic {q_ref:.6e}, "
          f"ratio {q_num/q_ref:.6f})")
    ok = err <= a.tolerance
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("uniform")
    p.add_argument("h5")
    p.add_argument("--scalar", action="append", required=True,
                   help="name=value, repeatable")
    p.add_argument("--levels", type=int, default=1)
    p.set_defaults(func=cmd_uniform)

    p = sub.add_parser("conserve")
    p.add_argument("first")
    p.add_argument("second")
    p.add_argument("--scalar", default="s1")
    p.add_argument("--tolerance", type=float, default=1e-13)
    p.set_defaults(func=cmd_conserve)

    p = sub.add_parser("wave")
    p.add_argument("h5")
    p.add_argument("--scalar", default="s1")
    p.add_argument("--wave", type=int, nargs=3, default=[1, 0, 0])
    p.add_argument("--amp", type=float, default=1.0)
    p.add_argument("--velocity", type=float, nargs=3, default=[0.0, 0.0, 0.0])
    p.add_argument("--diffusivity", type=float, required=True)
    p.add_argument("--ic", default=None,
                   help="IC file: its t_current is the time origin of the wave")
    p.set_defaults(func=cmd_wave)

    p = sub.add_parser("parabola")
    p.add_argument("h5")
    p.add_argument("--scalar", default="s1")
    p.add_argument("--source", type=float, required=True)
    p.add_argument("--diffusivity", type=float, required=True)
    p.add_argument("--s0", type=float, default=0.0)
    p.add_argument("--s1", type=float, default=0.0)
    p.set_defaults(func=cmd_parabola)

    p = sub.add_parser("profile")
    p.add_argument("h5")
    p.add_argument("--scalar", default="s1")
    p.add_argument("--diffusivity", type=float, default=None)
    p.add_argument("--tolerance", type=float, default=1e-2)
    p.set_defaults(func=cmd_profile)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
