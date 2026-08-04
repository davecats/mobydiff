#!/usr/bin/env python3
"""Read the passive-scalar statistics the solver accumulates (increment S4).

`[scalar] stats_sample_interval` / `stats_write_interval` make the solver
write `[scalar] stats_file` in one of two layouts (`[scalar] stats_layout`):

  profile   wall-normal rows, x-z averaged, one file per refinement level
            (`name.h5`, `name_l1.h5`, ...) -- the channel form, written with
            the channel_stats HDF5 writer;
  plane     rows of the global (x,y) plane, z averaged -- the boundary-layer
            form, written with the bl_stats writer.

Either way the row tables carry SCALAR_NSTAT = 7 columns per scalar, in the
order scalar_stats.f90 defines them:

    <s>  <s^2>  <u_c s>  <v s>|lo  <J>|lo  <v s>|hi  <J>|hi

with `lo`/`hi` the cell's low/high y face, J = v s - D ds/dy the TOTAL
wall-normal flux built with the transport kernel's own face diffusivity
(molecular + nut/Pr_t). A wall row's J is therefore the exact discrete wall
flux, which is what makes theta_tau and the Nusselt number below exact rather
than reconstructed.

Subcommands
  profile  <stats.h5>   mean profile, rms, turbulent flux, wall flux ->
                        theta_tau / Nusselt (the channel form)
  plane    <stats.h5>   the same at chosen x stations (the boundary-layer
                        form): local wall flux, Nu_x, and the profile dump
  heat     <heat.txt>   the immersed body's heat release -> Nusselt, from the
                        runtime file `[scalar] heat_interval` writes
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np

NSTAT = 7
S, SS, US, CLO, JLO, CHI, JHI = range(7)
COLUMNS = ["s", "ss", "us", "conv_lo", "flux_lo", "conv_hi", "flux_hi"]


def read_stats(path: str) -> dict:
    """One statistics file, either layout, as plain arrays."""
    with h5py.File(path, "r") as f:
        nstat = int(f.attrs["nstat"])
        out = dict(nstat=nstat, nscalar=nstat // NSTAT,
                   step=int(f.attrs["step"]), t=float(f.attrs["t_current"]),
                   re=float(f.attrs["re"]),
                   profile=f["profile"][...], raw_sum=f["raw_sum"][...],
                   count=f["count"][...])
        if nstat % NSTAT:
            raise SystemExit(f"{path}: nstat {nstat} is not a multiple of {NSTAT}")
        if "nwall" in f.attrs:                      # profile layout
            out.update(layout="profile", y=f["coord"][...])
        else:                                       # plane layout
            out.update(layout="plane", x=f["xcoord"][...], y=f["ycoord"][...],
                       nx=int(f.attrs["nx"]), ny=int(f.attrs["ny"]))
    return out


def column(st: dict, index: int, which: int) -> np.ndarray:
    """Column `which` of scalar `index` (1-based) over all rows."""
    if not 1 <= index <= st["nscalar"]:
        raise SystemExit(f"scalar index {index} outside 1..{st['nscalar']}")
    return st["profile"][:, NSTAT * (index - 1) + which]


def rms(st: dict, index: int) -> np.ndarray:
    return np.sqrt(np.maximum(column(st, index, SS) - column(st, index, S) ** 2, 0.0))


def cmd_profile(a) -> int:
    st = read_stats(a.file)
    if st["layout"] != "profile":
        raise SystemExit("this file has the plane layout -- use the `plane` subcommand")
    y = st["y"]
    mean = column(st, a.index, S)
    fluc = rms(st, a.index)
    conv = column(st, a.index, CLO)
    flux = column(st, a.index, JLO)
    fhi = column(st, a.index, JHI)

    print(f"{a.file}: step {st['step']}, t = {st['t']:.6f}, "
          f"{st['nscalar']} scalar(s), {y.size} rows, Re = {st['re']:g}")

    # Both walls: the low face of the first row and the high face of the last.
    jlo_wall, jhi_wall = float(flux[0]), float(fhi[-1])
    theta_tau = 0.5 * abs(jlo_wall + jhi_wall)
    print(f"wall flux: y_min {jlo_wall:+.6e}   y_max {jhi_wall:+.6e}"
          f"   -> theta_tau = {theta_tau:.6f}")

    # Total-flux constancy: in a statistically steady antisymmetric channel
    # J(y) is exactly constant, so this measures convergence AND the
    # transport operator's conservation in one number.
    if theta_tau > 0.0:
        dev = float(np.max(np.abs(flux - jlo_wall))) / theta_tau
        print(f"total flux J(y): {flux.min():+.6e} .. {flux.max():+.6e},"
              f"  max|J - J_wall|/theta_tau = {dev:.4f}")

    if a.walls is not None:
        dT = abs(a.walls[0] - a.walls[1])
        height = a.height if a.height else float(y[-1] + y[0])
        dmol = 1.0 / (st["re"] * a.pr)
        nu = abs(jlo_wall) * height / (dmol * dT) if dT > 0.0 else float("nan")
        print(f"Nusselt (q_w H / (D dT), H = {height:g}, dT = {dT:g}, "
              f"D = {dmol:.6e}) = {nu:.4f}")

    if a.rows:
        n = min(a.rows, y.size)
        print(f"{'y':>12} {'<s>':>14} {'s_rms':>12} {'<v s>|lo':>13} {'J|lo':>13}")
        for j in list(range(n)) + ([] if y.size <= 2 * n else [-1]):
            print(f"{y[j]:12.6f} {mean[j]:14.6e} {fluc[j]:12.6e} "
                  f"{conv[j]:13.6e} {flux[j]:13.6e}")

    if a.dump:
        np.savetxt(a.dump, np.column_stack([y, mean, fluc, conv, flux, fhi]),
                   header="y <s> s_rms <v s>|lo J|lo J|hi")
        print(f"profile written to {a.dump}")
    return 0


def cmd_plane(a) -> int:
    st = read_stats(a.file)
    if st["layout"] != "plane":
        raise SystemExit("this file has the profile layout -- use the `profile` subcommand")
    nx, ny = st["nx"], st["ny"]
    x, y = st["x"], st["y"]
    mean = column(st, a.index, S).reshape(nx, ny)
    fluc = rms(st, a.index).reshape(nx, ny)
    flux = column(st, a.index, JLO).reshape(nx, ny)
    dmol = 1.0 / (st["re"] * a.pr)

    print(f"{a.file}: step {st['step']}, t = {st['t']:.6f}, "
          f"{st['nscalar']} scalar(s), {nx} x {ny} plane, Re = {st['re']:g}")
    stations = a.stations if a.stations else [x[nx // 4], x[nx // 2], x[3 * nx // 4]]
    print(f"{'x':>12} {'q_wall':>14} {'<s>(y_1)':>13} {'s_rms(y_1)':>13} {'Nu_x':>12}")
    for xs in stations:
        i = int(np.argmin(np.abs(x - xs)))
        qw = float(flux[i, 0])
        dT = abs(a.walls[0] - a.walls[1]) if a.walls else 1.0
        nux = abs(qw) * x[i] / (dmol * dT)
        print(f"{x[i]:12.6f} {qw:14.6e} {mean[i,0]:13.6e} {fluc[i,0]:13.6e} {nux:12.4f}")

    if a.dump:
        i = int(np.argmin(np.abs(x - (a.stations[0] if a.stations else x[nx // 2]))))
        np.savetxt(a.dump, np.column_stack([y, mean[i], fluc[i], flux[i]]),
                   header=f"x = {x[i]:g}: y <s> s_rms J|lo")
        print(f"profile written to {a.dump}")
    return 0


def cmd_heat(a) -> int:
    """The body heat release [scalar] heat_interval writes, as a Nusselt
    number. The columns are the cancellation-free pair (staircase interface
    flux, graded-cell penalization) -- see scalar_stats.f90 and the S3 FINDING
    in docs/next_session_scalar.md for why the penalization integral ALONE is
    structurally incomplete for a Dirichlet body."""
    data = np.loadtxt(a.file, ndmin=2)
    if data.size == 0:
        raise SystemExit(f"{a.file} holds no samples")
    dth = 1.0 / (a.re * a.pr)
    col = 2 + 3 * (a.index - 1)
    print(f"{'step':>10} {'t':>12} {'staircase/Lz':>15} {'graded/Lz':>13} "
          f"{'Q/Lz':>13} {'Nu':>10}")
    for row in data[-a.last:]:
        stair, graded, total = row[col] / a.lz, row[col + 1] / a.lz, row[col + 2] / a.lz
        nu = total / (np.pi * dth * a.value)
        print(f"{int(row[0]):10d} {row[1]:12.4f} {stair:15.6e} {graded:13.6e} "
              f"{total:13.6e} {nu:10.4f}")
    if a.band:
        nu = data[-1, col + 2] / a.lz / (np.pi * dth * a.value)
        ok = a.band[0] <= nu <= a.band[1]
        print(f"literature band [{a.band[0]}, {a.band[1]}]: " + ("PASS" if ok else "FAIL"))
        return 0 if ok else 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("profile")
    p.add_argument("file")
    p.add_argument("--index", type=int, default=1, help="scalar index (1-based)")
    p.add_argument("--pr", type=float, default=1.0)
    p.add_argument("--walls", type=float, nargs=2, default=None,
                   help="wall values: enables the Nusselt number")
    p.add_argument("--height", type=float, default=None,
                   help="reference length for Nu (default: the row extent)")
    p.add_argument("--rows", type=int, default=0, help="print this many near-wall rows")
    p.add_argument("--dump", default=None)
    p.set_defaults(fn=cmd_profile)

    p = sub.add_parser("plane")
    p.add_argument("file")
    p.add_argument("--index", type=int, default=1)
    p.add_argument("--pr", type=float, default=1.0)
    p.add_argument("--walls", type=float, nargs=2, default=None)
    p.add_argument("--stations", type=float, nargs="+", default=None)
    p.add_argument("--dump", default=None)
    p.set_defaults(fn=cmd_plane)

    p = sub.add_parser("heat")
    p.add_argument("file")
    p.add_argument("--index", type=int, default=1)
    p.add_argument("--re", type=float, required=True)
    p.add_argument("--pr", type=float, required=True)
    p.add_argument("--lz", type=float, default=1.0, help="span (Q is a volume integral)")
    p.add_argument("--value", type=float, default=1.0, help="body-fluid scalar difference")
    p.add_argument("--last", type=int, default=5, help="print the last N samples")
    p.add_argument("--band", type=float, nargs=2, default=None)
    p.set_defaults(fn=cmd_heat)

    a = ap.parse_args()
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
