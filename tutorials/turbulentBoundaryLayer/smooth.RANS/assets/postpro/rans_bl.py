#!/usr/bin/env python3
"""2D RANS boundary layer against the DNS of ../smooth.

    python3 rans_bl.py FIELD.h5 [FIELD2.h5 ...] [--plot out.png] [--labels a b]

FIELD.h5 is a moby_solve snapshot of a smooth.RANS run (block layout, any
block size). The fields are z-averaged (the RANS is spanwise uniform; nz = 4
only because a block needs 4 cells) and brought to the cell centres. The
integral quantities use EXACTLY the definitions of the DNS export
(../smooth/assets/matlab/export_mat.py): U_e = U at the top row, delta*,
theta by the midpoint rule, delta_99 by linear interpolation, tau_w from the
second-order one-sided gradient through U(0) = 0, Re_x from the virtual
Blasius origin. The DNS reference is ../smooth/assets/matlab/tbl_dns_smooth.mat
(written by export_mat.py; kept locally, not in git).
"""
import argparse
import os

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DNS_MAT = os.path.join(HERE, "..", "..", "..", "smooth", "assets", "matlab", "tbl_dns_smooth.mat")
BLASIUS_DSTAR = 1.7207876573


def load_field(h5, name):
    """Reassemble a single-level block-layout dataset to (nz, ny, nx)."""
    blocks = h5["blocks"][...]
    nb = [int(h5.attrs[f"block_nb_{d}"]) for d in "xyz"]
    shape = tuple(int(h5.attrs[a]) for a in ("nz", "ny", "nx"))
    assert blocks[:, 3].max() == 0, "single-level layouts only"
    arr = np.empty(shape)
    data = h5[name]
    for bid, (ox, oy, oz, _) in enumerate(blocks):
        arr[oz:oz + nb[2], oy:oy + nb[1], ox:ox + nb[0]] = data[bid]
    return arr


def load(path):
    """Cell-centred, z-averaged fields of one snapshot: dict of (nx, ny)
    arrays (x first, as the DNS export before its transpose)."""
    with h5py.File(path, "r") as f:
        xn, yn = f["x"][...], f["y"][...]
        re = float(f.attrs["re"])
        q = {n: load_field(f, n).mean(axis=0) for n in
             ("un", "vn", "pn", "k", "omega", "nut", "gamma", "rethetat") if n in f}
        step, t = int(f.attrs["step"]), float(f.attrs["t_current"])
    # un(i) / vn(j) sit on the LOW face of cell i / j. The high face of the
    # last column / row is a boundary plane stored apart (un_xmax, vn_ymax);
    # the last cell copies its low face instead (the outlet zone and the top
    # row, neither used by the comparisons).
    u, v = q.pop("un"), q.pop("vn")             # (ny, nx)
    U = 0.5 * (u + np.c_[u[:, 1:], u[:, -1:]])
    V = 0.5 * (v + np.r_[v[1:, :], v[-1:, :]])
    xc, yc = 0.5 * (xn[1:] + xn[:-1]), 0.5 * (yn[1:] + yn[:-1])
    out = dict(x=xc, y=yc, re=re, step=step, t=t, U=U.T, V=V.T)
    for n, a in q.items():
        out[n] = a.T
    return out


def integrals(x, y, U, re):
    """The DNS export's 1D quantities from a (nx, ny) mean-velocity array."""
    nu = 1.0 / re
    nx = len(x)
    yf = np.r_[0.0, 0.5 * (y[1:] + y[:-1]), 2 * y[-1] - 0.5 * (y[-1] + y[-2])]
    dy = np.diff(yf)
    Ue = U[:, -1]
    r = U / Ue[:, None]
    dstar = np.sum((1 - r) * dy, axis=1)
    theta = np.sum(r * (1 - r) * dy, axis=1)
    j99 = np.argmax(r >= 0.99, axis=1)
    jm = j99 - 1
    i = np.arange(nx)
    d99 = y[jm] + (0.99 - r[i, jm]) * (y[j99] - y[jm]) / (r[i, j99] - r[i, jm])
    y0, y1 = y[0], y[1]
    tauw = nu * (U[:, 0] * y1**2 - U[:, 1] * y0**2) / (y0 * y1 * (y1 - y0))
    utau = np.sqrt(np.abs(tauw))
    rex = (re / BLASIUS_DSTAR) ** 2 + x * re
    return dict(x=x, Ue=Ue, dstar=dstar, theta=theta, delta99=d99, H=dstar / theta,
                tauw=tauw, utau=utau, cf=2 * tauw / Ue**2, Re_x=rex,
                Re_theta=Ue * theta / nu, Re_tau=utau * d99 / nu)


def dns():
    import scipy.io as sio
    m = sio.loadmat(DNS_MAT, squeeze_me=True)
    return {k: m[k] for k in ("x", "y", "U", "uv", "k", "nut_over_nu", "cf", "Re_theta", "Re_x",
                              "Re_tau", "H", "theta", "delta99", "dstar", "utau")}


def table(name, r, d):
    print(f"== {name}")
    print("    x  | Re_th RANS  DNS  |  cf RANS   DNS    err% |  H RANS  DNS  | Re_tau RANS DNS")
    for xs in (30, 50, 70, 100, 150, 200, 300, 400, 500, 550):
        i = np.argmin(abs(r["x"] - xs))
        j = np.argmin(abs(d["x"] - xs))
        print(f"  {xs:4.0f} | {r['Re_theta'][i]:6.0f} {d['Re_theta'][j]:6.0f} | "
              f"{r['cf'][i]:.5f} {d['cf'][j]:.5f} {100 * (r['cf'][i] / d['cf'][j] - 1):+6.1f} | "
              f"{r['H'][i]:.3f} {d['H'][j]:.3f} | {r['Re_tau'][i]:5.0f} {d['Re_tau'][j]:5.0f}")
    # c_f at matched Re_theta (the transition-independent comparison)
    print("  cf at matched Re_theta:", end="")
    for rt in (400, 500, 600, 700, 800):
        if rt > r["Re_theta"].max():
            continue
        c1 = np.interp(rt, r["Re_theta"][r["x"] > 100], r["cf"][r["x"] > 100])
        c2 = np.interp(rt, d["Re_theta"][d["x"] > 100], d["cf"][d["x"] > 100])
        print(f"  {rt}: {100 * (c1 / c2 - 1):+.1f}%", end="")
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fields", nargs="+")
    ap.add_argument("--labels", nargs="*")
    ap.add_argument("--plot")
    ap.add_argument("--profiles", help="second figure: wall-unit profiles at matched Re_theta")
    ap.add_argument("--maps", help="third figure: gamma (or log10 omega), k, nu_t/nu maps of the FIRST field")
    a = ap.parse_args()
    d = dns()
    runs = []
    for p in a.fields:
        f = load(p)
        r = integrals(f["x"], f["y"], f["U"], f["re"])
        runs.append((os.path.basename(p), f, r))
        table(f"{p} (step {f['step']}, t {f['t']:.0f})", r, d)
    if not a.plot:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = a.labels or [n for n, _, _ in runs]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    ax = ax.ravel()
    ax[0].plot(d["x"], d["cf"], "k", lw=2, label="DNS")
    ax[1].plot(d["Re_theta"], d["cf"], "k", lw=2, label="DNS")
    ax[2].plot(d["Re_x"], d["Re_theta"], "k", lw=2, label="DNS")
    ax[3].plot(d["Re_theta"], d["H"], "k", lw=2, label="DNS")
    # reference curves: Blasius and Coles-Fernholz (kappa 0.384, C 4.127)
    rth = np.linspace(150, 1000, 200)
    rex = np.linspace(d["Re_x"].min(), d["Re_x"].max(), 200)
    ax[0].plot((rex - d["Re_x"][0]) / 450 + d["x"][0], 0.664 / np.sqrt(rex), "k:", lw=1, label="Blasius")
    ax[1].plot(rth, 0.441 / rth, "k:", lw=1, label="Blasius")
    ax[1].plot(rth, 2 * (np.log(rth) / 0.384 + 4.127) ** -2, "k--", lw=1, label="Coles-Fernholz")
    ax[2].plot(rex, 0.664 * np.sqrt(rex), "k:", lw=1, label="Blasius")
    for lab, f, r in zip(labels, [x[1] for x in runs], [x[2] for x in runs]):
        keep = r["x"] <= 550                     # the DNS export's range (outlet zone dropped)
        ax[0].plot(r["x"][keep], r["cf"][keep], label=lab)
        ax[1].plot(r["Re_theta"][keep], r["cf"][keep], label=lab)
        ax[2].plot(r["Re_x"][keep], r["Re_theta"][keep], label=lab)
        ax[3].plot(r["Re_theta"][keep], r["H"][keep], label=lab)
    ax[0].set(xlabel="x", ylabel="c_f", xlim=(0, 550), ylim=(0, 0.008))
    ax[2].set(xlim=(0, d["Re_x"].max()), ylim=(0, 1000))
    ax[1].set(xlabel="Re_theta", ylabel="c_f", ylim=(0, 0.008))
    ax[2].set(xlabel="Re_x", ylabel="Re_theta")
    ax[3].set(xlabel="Re_theta", ylabel="H", ylim=(1.3, 2.7))
    for a_ in ax:
        a_.grid(alpha=0.3)
        a_.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(a.plot, dpi=110)
    print("wrote", a.plot)
    if a.profiles:
        profiles(d, runs, labels, a.profiles)
    if a.maps:
        maps(runs[0][1], runs[0][2], a.maps)


def maps(f, r, out, xmax=550.0, ymax=25.0):
    """(x, y) maps of one RANS field: intermittency (log10 omega when the
    transition model is off), k and nu_t/nu, delta_99."""
    import matplotlib.pyplot as plt
    keep = f["x"] <= xmax
    jy = f["y"] <= ymax
    X, Y = f["x"][keep], f["y"][jy]
    first = (("gamma", f["gamma"], "intermittency $\\gamma$", (0, 1)) if "gamma" in f else
             ("omega", np.clip(np.log10(f["omega"]), -2.5, 0.5),
              "$\\log_{10}(\\omega\\, \\delta^*_0 / U_\\infty)$ (clipped to -2.5 .. 0.5)", (-2.5, 0.5)))
    items = [first,
             ("k", f["k"] * 1e3, "$k / U_\\infty^2 \\times 10^3$", None),
             ("nut", f["nut"] * f["re"], "$\\nu_t / \\nu$", None)]
    fig, ax = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    for a_, (_, q, ttl, lim) in zip(ax, items):
        z = q[keep][:, jy].T
        cs = a_.contourf(X, Y, z, 40, cmap="viridis",
                         vmin=lim[0] if lim else None, vmax=lim[1] if lim else None)
        a_.plot(X, r["delta99"][keep], "w--", lw=1)
        fig.colorbar(cs, ax=a_, pad=0.01)
        a_.set(ylabel="y", title=ttl)
    ax[-1].set_xlabel("x  (white dashed: delta_99)")
    fig.tight_layout()
    fig.savefig(out, dpi=100)
    print("wrote", out)


def rans_stresses(f):
    """Boussinesq shear stress -<u'v'> = nu_t (dU/dy + dV/dx) of a RANS field
    (nx, ny); the wall at y = 0 closes the first gradient with U = 0."""
    x, y = f["x"], f["y"]
    dUdy = np.gradient(np.c_[np.zeros(len(x)), f["U"]], np.r_[0.0, y], axis=1, edge_order=2)[:, 1:]
    dVdx = np.gradient(f["V"], x, axis=0, edge_order=2)
    return f["nut"] * (dUdy + dVdx)


def profiles(d, runs, labels, out, stations=(400, 500, 600, 700)):
    """Wall-unit profiles at matched Re_theta: U+, -u'v'+, k+, nu_t/nu."""
    import matplotlib.pyplot as plt
    nu = 1.0 / 450.0
    fig, ax = plt.subplots(len(stations), 4, figsize=(15, 3.2 * len(stations)))
    for row, rt in enumerate(stations):
        # DNS station (x > 100: past the trip / overshoot)
        sel = np.where(d["x"] > 100)[0]
        j = sel[np.argmin(abs(d["Re_theta"][sel] - rt))]
        ut = d["utau"][j]
        yp = d["y"] * ut / nu
        a = ax[row]
        a[0].semilogx(yp, d["U"][:, j] / ut, "k", lw=2, label="DNS")
        a[1].semilogx(yp, -d["uv"][:, j] / ut**2, "k", lw=2)
        a[2].semilogx(yp, d["k"][:, j] / ut**2, "k", lw=2)
        a[3].plot(d["y"] / d["delta99"][j], d["nut_over_nu"][:, j], "k", lw=2)
        for lab, (_, f, r) in zip(labels, runs):
            sel = np.where((r["x"] > 100) & (r["x"] <= 550))[0]   # past the overshoot, before the outlet zone
            i = sel[np.argmin(abs(r["Re_theta"][sel] - rt))]
            ut = r["utau"][i]
            yp = f["y"] * ut / nu
            a[0].semilogx(yp, f["U"][i] / ut, label=f"{lab} (x={r['x'][i]:.0f})")
            a[1].semilogx(yp, rans_stresses(f)[i] / ut**2)
            a[2].semilogx(yp, f["k"][i] / ut**2)
            a[3].plot(f["y"] / r["delta99"][i], f["nut"][i] / nu)
        a[0].set(ylabel=f"Re_theta {rt}\nU+", xlim=(0.2, 2000), ylim=(0, 25))
        a[1].set(ylabel="-u'v'+", xlim=(0.2, 2000))
        a[2].set(ylabel="k+", xlim=(0.2, 2000))
        a[3].set(ylabel="nu_t/nu", xlim=(0, 1.5))
        a[0].legend(fontsize=7)
        for b in a:
            b.grid(alpha=0.3)
    for b, xl in zip(ax[-1], ("y+", "y+", "y+", "y/delta99")):
        b.set_xlabel(xl)
    fig.tight_layout()
    fig.savefig(out, dpi=100)
    print("wrote", out)


if __name__ == "__main__":
    main()
