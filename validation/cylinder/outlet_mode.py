#!/usr/bin/env python3
"""The outlet pressure mode on the Re 100 cylinder: what one time unit does to
the stored pressure, per projection (run_outlet_mode.sh).

  outlet_mode.py <dir with om_<leg>_<k>_*.h5 and forces_om_<leg>_<k>.txt>

Per leg: the pressure at T (rms about its mean, max), where it sits (outlet
band, lateral bands, body box, the rest), its 2-dx content (p minus the mean
of its four x,y neighbours), the divergence the last projection left, the
control-volume C_L range over the leg, and the velocity difference to the
most converged leg of the same step size, with the place of its maximum.
"""
import glob
import os
import re
import sys

import h5py
import numpy as np

d = sys.argv[1] if len(sys.argv) > 1 else "."
ORDER = ["cheb12", "cheb60", "cheb240", "jac60", "jac240", "rb60", "rb240"]


def assemble(path):
    """-> dict of global (nz, ny, nx) arrays, t, (lx, ly, lz)."""
    with h5py.File(path) as h:
        blocks = h["blocks"][...]
        assert (blocks[:, 3] == 0).all(), "single-level case expected"
        nx, ny, nz = (int(h.attrs[k]) for k in ("nx", "ny", "nz"))
        out = {}
        for name in ("un", "vn", "wn", "pn"):
            a = h[name][...]
            nb = a.shape[1:]                       # (nbz, nby, nbx)
            g = np.empty((nz, ny, nx))
            for b, (ox, oy, oz, _l) in enumerate(blocks):
                g[oz:oz+nb[0], oy:oy+nb[1], ox:ox+nb[2]] = a[b]
            out[name] = g
        return out, float(h.attrs["t_current"]), tuple(float(h.attrs[k]) for k in ("lx", "ly", "lz"))


def final(pfx):
    f = sorted(glob.glob(os.path.join(d, pfx + "_[0-9]*.h5")),
               key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
    return f[-1] if f else None


legs, longs = {}, []
for path in glob.glob(os.path.join(d, "om_*_[0-9].ini")):
    m = re.match(r"om_(\w+)_(\d)\.ini", os.path.basename(path))
    if not m or not final(f"om_{m.group(1)}_{m.group(2)}"):
        continue
    if re.search(r"T\d+$", m.group(1)):
        longs.append(f"om_{m.group(1)}_{m.group(2)}")       # a multi-time-unit leg
    else:
        legs[(m.group(1), int(m.group(2)))] = f"om_{m.group(1)}_{m.group(2)}"
if not legs and not longs:
    sys.exit("no legs found in " + d)

# Time evolution of a long leg: one row per snapshot. The mode is the part of
# p that is ANTISYMMETRIC in y, so it is measured directly: A = half the
# difference of the mean pressure of the upper and lower quarter of the
# domain (x < lx - 1), next to the y-rms of p in the first and the last cell
# column, the outflow imbalance (mean u through the upper minus the lower
# half of the outlet) and the C_L range over the preceding time unit.
for pfx in sorted(longs):
    snaps = sorted(glob.glob(os.path.join(d, pfx + "_[0-9]*.h5")),
                   key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
    f = np.loadtxt(os.path.join(d, f"forces_{pfx}.txt"), skiprows=1)
    print(f"\n== {pfx}: time evolution")
    print(f"{'t':>8s} | {'p rms':>7s} {'A(antisym)':>10s} {'std in':>7s} {'std out':>7s} |"
          f" {'u_out up-lo':>11s} | {'C_L range, last t.u.':>20s} | {'C_D mean':>8s}")
    tprev = None
    for sp in snaps:
        q, t, (lx, ly, lz) = assemble(sp)
        p = (q["pn"] - q["pn"].mean()).mean(axis=0)               # (ny, nx)
        ny, nx = p.shape
        nxin = int(nx*(lx - 1.0)/lx)
        amp = 0.5*(p[3*ny//4:, :nxin].mean() - p[:ny//4, :nxin].mean())
        uo = q["un"].mean(axis=0)[:, -1]
        sel = (f[:, 1] <= t + 1e-9) & (f[:, 1] > (t - 1.0 if tprev is None else tprev) + 1e-9)
        cl = f[sel, 2] if sel.any() else np.array([np.nan])
        cd = f[sel, 3] if sel.any() else np.array([np.nan])
        print(f"{t:8.3f} | {np.sqrt((p**2).mean()):7.3f} {amp:+10.3f} {p[:, 0].std():7.3f} {p[:, -1].std():7.3f} |"
              f" {uo[ny//2:].mean() - uo[:ny//2].mean():+11.4f} | {cl.min():+9.3f} {cl.max():+9.3f}  | {cd.mean():8.4f}")
        tprev = t

for k in sorted({key[1] for key in legs}):
    names = [n for n in ORDER if (n, k) in legs]
    data = {n: assemble(final(legs[(n, k)])) for n in names}
    ref = next(n for n in ("rb240", "jac240", "cheb240", "rb60", "jac60", "cheb60") if n in data)
    q0, t0, (lx, ly, lz) = data[ref]
    nz, ny, nx = q0["pn"].shape
    dx, dy, dz = lx/nx, ly/ny, lz/nz
    xc = (np.arange(nx) + 0.5)*dx
    yc = (np.arange(ny) + 0.5)*dy
    X, Y = np.meshgrid(xc, yc)                                   # (ny, nx)
    regions = {
        "outlet": X > lx - 1.0,
        "lateral": (np.abs(Y - ly/2) > ly/2 - 1.0) & (X <= lx - 1.0),
        "body": (X > 4.0) & (X < 8.0) & (Y > 6.5) & (Y < 9.5),
    }
    regions["rest"] = ~(regions["outlet"] | regions["lateral"] | regions["body"])
    print(f"\n== dt = {5e-3/2**k:.3e} ({200 << k} steps), t = {t0:.4f}; velocity reference: {ref}")
    print(f"{'leg':8s} | {'p rms':>7s} {'max|p|':>7s} | " + " ".join(f"{r:>8s}" for r in regions)
          + f" | {'2dx rms':>8s} | {'div rms':>8s} {'div max':>8s} | {'C_L range':>17s} |"
          f" {'du rms':>8s} {'du max':>8s}  at (x, y)")
    for n in names:
        q, t, _ = data[n]
        assert abs(t - t0) < 1e-9, f"{n} ends at {t}"
        p = q["pn"] - q["pn"].mean()
        p2 = p.mean(axis=0)                                       # z-mean plane
        reg = " ".join(f"{np.sqrt((p2[m]**2).mean()):8.3f}" for m in regions.values())
        # 2-dx content: p minus the mean of its four in-plane neighbours
        chk = p2[1:-1, 1:-1] - 0.25*(p2[:-2, 1:-1] + p2[2:, 1:-1] + p2[1:-1, :-2] + p2[1:-1, 2:])
        # divergence of the stored velocity (u(i) on the LOW face of cell i);
        # the last x column and y row need a face the file does not hold
        u, v, w = q["un"], q["vn"], q["wn"]
        div = ((u[:, :-1, 1:] - u[:, :-1, :-1])/dx + (v[:, 1:, :-1] - v[:, :-1, :-1])/dy
               + (np.roll(w, -1, axis=0) - w)[:, :-1, :-1]/dz)
        cl = np.loadtxt(os.path.join(d, f"forces_{legs[(n, k)]}.txt"), skiprows=1)[:, 2]
        du = np.sqrt(sum((q[c] - q0[c])**2 for c in ("un", "vn", "wn"))).mean(axis=0)
        j, i = np.unravel_index(np.argmax(du), du.shape)
        print(f"{n:8s} | {np.sqrt((p**2).mean()):7.3f} {np.abs(p).max():7.3f} | {reg} | "
              f"{np.sqrt((chk**2).mean()):8.2e} | {np.sqrt((div**2).mean()):8.2e} {np.abs(div).max():8.2e} | "
              f"{cl.min():+8.3f} {cl.max():+8.3f} | {np.sqrt((du**2).mean()):8.2e} {du.max():8.2e}"
              f"  ({xc[i]:5.2f}, {yc[j]:5.2f})")
