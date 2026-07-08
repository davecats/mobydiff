#!/usr/bin/env python3
"""RANS T2 channel gate checker (docs/next_session_iddes.md, phase T2).

Reads a solver field file (block-table layout), builds x/z-averaged
profiles of u, k, omega, nut on the y line and checks:

  --mode laminar   gate (a): u(y) vs the exact parabola u = (re/2) y (2-y)
                   (forcing_x = 1, ly = 2); k must have decayed to
                   negligible levels (no sustained turbulence).
  --mode loglaw    gates (b)/(c): U+(y+) vs the log law (kappa 0.41, B 5.0)
                   over y+ in [30, 0.3 Re_tau]. u_tau is measured from the
                   first-cell wall gradient. --wall-lo/--wall-hi shift the
                   walls for the immersed-boundary channel (gate c).
  --mode band      gate (d): multi-level (uniform-y) profiles; k/omega/nut
                   must cross the 2:1 interfaces with no band, and the
                   profiles must match the single-level twin (--reference).

Uniform-y reassembly in band mode uses the block table (origin in
level-l cells, midpoint subdivision), so it needs no per-block coordinate
output.
"""
import argparse
import sys

import h5py
import numpy as np


def load_raw(path):
    with h5py.File(path, "r") as f:
        blocks = f["blocks"][...]            # (n, 4): origin x,y,z + level
        data = {"u": f["un"][...]}
        for name, ds in (("k", "k"), ("om", "omega"), ("nut", "nut")):
            data[name] = f[ds][...] if ds in f else None
        y_nodes = f["y"][...]
        re = float(f.attrs["re"])
        ly = float(f.attrs["ly"])
        ny = int(f.attrs["ny"])
    return blocks, data, y_nodes, re, ly, ny


def profiles_single_level(blocks, data, y_nodes):
    if (blocks[:, 3] != 0).any():
        raise SystemExit("multi-level field: use --mode band")
    nb = data["u"].shape[-1]
    ny = int(blocks[:, 1].max()) + nb
    yc = 0.5*(y_nodes[1:] + y_nodes[:-1])
    cnt = np.zeros(ny)
    prof = {name: (np.zeros(ny) if data[name] is not None else None) for name in data}
    for b in range(blocks.shape[0]):
        oy = blocks[b, 1]
        cnt[oy:oy+nb] += nb*nb
        for name, arr in data.items():
            if arr is not None:
                prof[name][oy:oy+nb] += arr[b].sum(axis=(0, 2))
    for name in data:
        if data[name] is not None:
            prof[name] = prof[name]/cnt
    return yc, prof


def profiles_per_level(blocks, data, ly, ny_glob):
    """Uniform-y multi-level case: per-level x/z-averaged profiles keyed by
    level, sized on the GLOBAL level lattice (ny_glob * 2^lev) -- sizing
    from the highest occupied origin mislocates every row."""
    nb = data["u"].shape[-1]
    levels = sorted(set(blocks[:, 3].tolist()))
    out = {}
    for lev in levels:
        sel = np.where(blocks[:, 3] == lev)[0]
        ny_l = ny_glob*(2**lev)
        sums = {n: np.zeros(ny_l) for n in data if data[n] is not None}
        cnt = np.zeros(ny_l)
        for b in sel:
            o = blocks[b, 1]
            cnt[o:o+nb] += nb*nb
            for n in sums:
                sums[n][o:o+nb] += data[n][b].sum(axis=(0, 2))
        rows = cnt > 0
        prof = {n: np.where(rows, sums[n]/np.maximum(cnt, 1), np.nan) for n in sums}
        # y centre of level-l cell j; the level lattice spans the domain.
        # blocks report origins in level-l cells of the FINEST lattice
        # numbering? No: origin is in level-l cells (leaf table convention).
        out[lev] = (rows, prof, cnt)
    return out


def check_laminar(yc, prof, re, tol):
    exact = 0.5*re*yc*(2.0 - yc)
    err = np.abs(prof["u"] - exact).max()/exact.max()
    kmax = np.nanmax(prof["k"])
    print(f"laminar: max |u - parabola| / u_max = {err:.3e}  (tol {tol})")
    print(f"k max = {kmax:.3e} (seeded ~ 4e-3; must decay by orders of magnitude)")
    return err <= tol and kmax < 1.0e-6


def check_loglaw(yc, prof, re, tol, wall_lo, wall_hi, uplus_center=None):
    nu = 1.0/re
    half_h = 0.5*(wall_hi - wall_lo)
    y_rel = yc - wall_lo
    fluid = (y_rel > 0.0) & (yc < wall_hi)
    yr = y_rel[fluid]
    u = prof["u"][fluid]
    i_lo = np.argmin(yr)
    i_hi = np.argmax(yr)
    dudy_lo = u[i_lo]/yr[i_lo]
    dudy_hi = u[i_hi]/(wall_hi - wall_lo - yr[i_hi])
    u_tau = np.sqrt(nu*0.5*(abs(dudy_lo) + abs(dudy_hi)))
    re_tau = u_tau*half_h*re
    print(f"u_tau (wall gradient) = {u_tau:.4f}, Re_tau = {re_tau:.1f}")

    lower = yr <= half_h
    yp = yr[lower]*u_tau*re
    up = u[lower]/u_tau
    order = np.argsort(yp)
    yp, up = yp[order], up[order]
    log_range = (yp >= 30.0) & (yp <= 0.3*re_tau)
    if log_range.sum() < 3:
        print("log region too short -- FAIL")
        return False
    ref = np.log(yp[log_range])/0.41 + 5.0
    rel = np.abs(up[log_range] - ref)/ref
    print("log region y+ in [%.0f, %.0f], %d points; max |U+ - loglaw|/loglaw = %.3f (tol %g)"
          % (yp[log_range].min(), yp[log_range].max(), log_range.sum(), rel.max(), tol))
    print(f"U+ centreline = {up[-1]:.2f}, bulk U+ = {np.trapezoid(up, yp/(u_tau*re)):.2f}")
    for name in ("k", "om", "nut"):
        if prof[name] is not None:
            p = prof[name][fluid]
            print(f"{name}: min {np.nanmin(p):.3e} max {np.nanmax(p):.3e}")
    ok = rel.max() <= tol and abs(u_tau - 1.0) < 0.06
    if uplus_center is not None:
        # DNS centreline anchor: the pure kappa/B log line deviates from
        # real (and SST) profiles by several % in the overlap region, so
        # the centreline U+ vs DNS is the sharper RANS-quality measure.
        dev_c = abs(up[-1] - uplus_center)/uplus_center
        print(f"U+ centreline vs DNS {uplus_center}: dev {dev_c:.3f} (tol 0.02)")
        ok = ok and dev_c <= 0.02
    return ok


def check_band(path, ref_path, ly, tol_band, tol_ref):
    blocks, data, y_nodes, re, _, ny_glob = load_raw(path)
    nb = data["u"].shape[-1]
    per = profiles_per_level(blocks, data, ly, ny_glob)
    if sorted(per) != [0, 1]:
        raise SystemExit("band mode expects levels {0, 1}")

    rows0, prof0, _ = per[0]
    rows1, prof1, _ = per[1]
    ny0 = rows0.size
    ny1 = rows1.size
    h0 = ly/ny0
    h1 = 0.5*h0
    yc0 = (np.arange(ny0) + 0.5)*h0
    yc1 = (np.arange(ny1) + 0.5)*h1

    ok = True
    # Interface band metric at the LOWER interface (the top one mirrors
    # it). The coarse-pair average of the last two LOWER-band fine rows
    # sits one coarse row below the first coarse row, so for a smooth
    # profile the jump equals ~one local row-to-row increment (ratio ~1);
    # an interface band shows up as ratio >> 1.
    coarse_rows = np.where(rows0)[0]
    jlo_c = int(coarse_rows.min())            # first coarse row above the band
    for name in ("k", "om", "nut", "u"):
        pf = prof1[name]
        pc = prof0[name]
        # last coarse-pair of the LOWER fine band (fine rows 2 jlo_c - 2/-1)
        fine_last = 0.5*(pf[2*jlo_c - 2] + pf[2*jlo_c - 1])
        coarse_first = pc[jlo_c]
        coarse_second = pc[jlo_c + 1]
        local = abs(coarse_second - coarse_first)
        jump = abs(coarse_first - fine_last)
        ratio = jump/max(local, 0.02*abs(coarse_first) + 1e-30)
        print(f"{name}: interface jump {jump:.3e} vs local variation {local:.3e} "
              f"(ratio {ratio:.2f}, tol {tol_band})")
        ok = ok and ratio <= tol_band

    if ref_path:
        rblocks, rdata, ry_nodes, rre, rly, _rny = load_raw(ref_path)
        ryc, rprof = profiles_single_level(rblocks, rdata, ry_nodes)
        # compare coarse-region profiles (y in [0.35, 1.0]) refined vs twin
        selr = (ryc >= 0.35) & (ryc <= 1.0)
        for name in ("u", "k", "om", "nut"):
            pr = np.interp(ryc[selr], yc0[rows0], prof0[name][rows0])
            rel = np.abs(pr - rprof[name][selr]).max()/(np.abs(rprof[name][selr]).max() + 1e-30)
            print(f"{name}: refined-vs-single-level core profile max rel dev = {rel:.3f} (tol {tol_ref})")
            ok = ok and rel <= tol_ref
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("field")
    ap.add_argument("--mode", choices=("laminar", "loglaw", "band"), required=True)
    ap.add_argument("--tolerance", type=float, default=None)
    ap.add_argument("--wall-lo", type=float, default=0.0)
    ap.add_argument("--wall-hi", type=float, default=None)
    ap.add_argument("--reference", default=None,
                    help="resolved single-level reference field for band mode")
    ap.add_argument("--uplus-center", type=float, default=None,
                    help="DNS centreline U+ anchor (2%% tolerance)")
    args = ap.parse_args()

    if args.mode == "band":
        blocks, data, y_nodes, re, ly, ny_glob = load_raw(args.field)
        ok = check_band(args.field, args.reference, ly,
                        args.tolerance if args.tolerance is not None else 1.5, 0.08)
    else:
        blocks, data, y_nodes, re, ly, ny_glob = load_raw(args.field)
        yc, prof = profiles_single_level(blocks, data, y_nodes)
        if args.mode == "laminar":
            ok = check_laminar(yc, prof, re,
                               args.tolerance if args.tolerance is not None else 2.0e-2)
        else:
            wall_hi = args.wall_hi if args.wall_hi is not None else ly
            ok = check_loglaw(yc, prof, re,
                              args.tolerance if args.tolerance is not None else 0.06,
                              args.wall_lo, wall_hi, args.uplus_center)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
