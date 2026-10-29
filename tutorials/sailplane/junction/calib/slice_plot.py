#!/usr/bin/env python3
"""Flow visualisation of a block-table snapshot: planar slices sampled on
a regular grid with strip_loads.Field (trilinear; the staggered velocities
are shifted by half a cell onto cell centres). Four views of the
calibration box by default (chord units, body = zero velocity, masked):
  y = 0.6 m and 1.0 m sections (|u| + streamlines),
  plan view z = 10.11 just above the wing (u),
  cross-plane x = 7.6 through the junction (u + in-plane vectors).
Usage: slice_plot.py SNAP.h5 --case CASE.h5 --out DIR [--h 0.002]"""
import argparse, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from strip_loads import Field

ap = argparse.ArgumentParser()
ap.add_argument("field"); ap.add_argument("--case", required=True); ap.add_argument("--out", default="viz")
ap.add_argument("--h", type=float, default=0.002, help="sampling spacing (chords)")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
F = Field(a.field, a.case); t = float(F.f.attrs["t_current"]); tag = os.path.basename(a.field).replace(".h5", "")
hmin = F.h(F.lmax)

def grid(xr, yr, zr, h):
    axes = [np.arange(r[0], r[1] + 1e-9, h) if r[1] > r[0] else np.array([r[0]]) for r in (xr, yr, zr)]
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    return axes, np.c_[X.ravel(), Y.ravel(), Z.ravel()]

def fields(pts):
    row, lev = F.locate(pts)
    u = F.sample("un", pts + [0.5 * hmin[0], 0, 0], row, lev)   # low-face value -> centre of the cell behind it
    v = F.sample("vn", pts + [0, 0.5 * hmin[1], 0], row, lev)
    w = F.sample("wn", pts + [0, 0, 0.5 * hmin[2]], row, lev)
    p = F.sample("pn", pts, row, lev)
    solid = ((np.abs(u) + np.abs(v) + np.abs(w)) < 1e-9) | np.isnan(u)   # penalized cells, and the removed buried blocks (no leaf)
    u, v, w, p = (np.where(np.isnan(x), 0.0, x) for x in (u, v, w, p))
    return u, v, w, p, solid

def section(y, name):
    (xs, _, zs), pts = grid((6.8, 8.8), (y, y), (9.55, 10.45), a.h)
    u, v, w, p, solid = fields(pts); sh = (len(xs), len(zs))
    U = np.ma.array(u.reshape(sh), mask=solid.reshape(sh)); W = np.ma.array(w.reshape(sh), mask=solid.reshape(sh))
    mag = np.sqrt(U**2 + W**2)
    fig, ax = plt.subplots(figsize=(13, 6.2))
    pc = ax.pcolormesh(xs, zs, mag.T, cmap="viridis", vmin=0, vmax=1.6, shading="auto")
    ax.streamplot(xs, zs, U.T.filled(np.nan), W.T.filled(np.nan), density=2.2, color="w", linewidth=0.5, arrowsize=0.6)
    ax.contourf(xs, zs, solid.reshape(sh).T.astype(float), levels=[0.5, 1.5], colors=["#444444"])
    ax.set_aspect("equal"); ax.set_xlabel("x / c"); ax.set_ylabel("z / c")
    ax.set_title(f"{tag}, t = {t:.2f}: section y = {y:.3f} c ({name}), |u| and streamlines")
    fig.colorbar(pc, ax=ax, label="|u| / U", shrink=0.8); fig.tight_layout()
    f = os.path.join(a.out, f"{tag}_section_y{y:.3f}.png"); fig.savefig(f, dpi=130); plt.close(fig); print("wrote", f)

def plan(z):
    (xs, ys, _), pts = grid((6.6, 9.2), (0.0, 3.5), (z, z), a.h * 2)
    u, v, w, p, solid = fields(pts); sh = (len(xs), len(ys))
    U = np.ma.array(u.reshape(sh), mask=solid.reshape(sh))
    fig, ax = plt.subplots(figsize=(13, 7))
    pc = ax.pcolormesh(xs, ys, U.T, cmap="RdBu_r", vmin=0.0, vmax=1.6, shading="auto")
    ax.contourf(xs, ys, solid.reshape(sh).T.astype(float), levels=[0.5, 1.5], colors=["#444444"])
    ax.set_aspect("equal"); ax.set_xlabel("x / c"); ax.set_ylabel("y / c (symmetry plane at 0, cut at 3.5)")
    ax.set_title(f"{tag}, t = {t:.2f}: plan view z = {z:.2f} c (just above the upper surface), u / U")
    fig.colorbar(pc, ax=ax, label="u / U", shrink=0.8); fig.tight_layout()
    f = os.path.join(a.out, f"{tag}_plan_z{z:.2f}.png"); fig.savefig(f, dpi=130); plt.close(fig); print("wrote", f)

def cross(x):
    (_, ys, zs), pts = grid((x, x), (0.0, 1.6), (9.3, 10.7), a.h)
    u, v, w, p, solid = fields(pts); sh = (len(ys), len(zs))
    U = np.ma.array(u.reshape(sh), mask=solid.reshape(sh)); V = v.reshape(sh); W = w.reshape(sh); S = solid.reshape(sh)
    fig, ax = plt.subplots(figsize=(11, 9))
    pc = ax.pcolormesh(ys, zs, U.T, cmap="RdBu_r", vmin=0.0, vmax=1.6, shading="auto")
    k = 12; Yq, Zq = np.meshgrid(ys[::k], zs[::k], indexing="ij"); Vq = np.where(S, np.nan, V)[::k, ::k]; Wq = np.where(S, np.nan, W)[::k, ::k]
    ax.quiver(Yq, Zq, Vq, Wq, color="k", scale=6, width=0.0025)
    ax.contourf(ys, zs, S.T.astype(float), levels=[0.5, 1.5], colors=["#444444"])
    ax.set_aspect("equal"); ax.set_xlabel("y / c"); ax.set_ylabel("z / c")
    ax.set_title(f"{tag}, t = {t:.2f}: cross-plane x = {x:.1f} c through the junction, u / U and (v, w)")
    fig.colorbar(pc, ax=ax, label="u / U", shrink=0.8); fig.tight_layout()
    f = os.path.join(a.out, f"{tag}_cross_x{x:.1f}.png"); fig.savefig(f, dpi=130); plt.close(fig); print("wrote", f)

section(0.722, "0.6 m, junction station"); section(1.203, "1.0 m"); plan(10.11); cross(7.6)
