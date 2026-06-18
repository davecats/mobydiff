#!/usr/bin/env python3
"""Post-processing for the 2:1 interface channel validation.

Overlays a band-refined channel run against the uniform reference:
mean velocity, rms, Reynolds stress (from the per-level channel stats
files), streamwise/spanwise u-spectra at y+ ~ 15 and at the interface
height (from field snapshots), and divergence residuals in the
interface band vs the interior. Reports deviations; tunes nothing.

Usage:
  python3 tools/channel_interface_validation.py \
      --reference validation/channel_interface/runs/reference/stats \
      --refined  validation/channel_interface/runs/refined_y110/stats \
      --out plots_y110 [--snapshots 10]
"""

from __future__ import annotations

import argparse
import glob
import os

import h5py
import numpy as np

STAT_U, STAT_V, STAT_W = 0, 1, 2
STAT_UU, STAT_VV, STAT_WW, STAT_UV = 3, 4, 5, 6


# ---------------------------------------------------------------------------
# Profiles

def read_stats(path):
    with h5py.File(path, "r") as h5:
        return {"coord": h5["coord"][...], "count": h5["count"][...],
                "profile": h5["profile"][...], "re": float(h5.attrs["re"])}


def composite_profiles(statsdir, base="channel_stats"):
    """Merge per-level stats files into one profile table sorted by y,
    keeping only rows that actually accumulated samples."""
    parts = []
    for path in sorted(glob.glob(os.path.join(statsdir, base + "*.h5"))):
        s = read_stats(path)
        mask = s["count"] > 0
        parts.append((s["coord"][mask], s["profile"][mask], s["re"]))
    coord = np.concatenate([p[0] for p in parts])
    prof = np.concatenate([p[1] for p in parts])
    order = np.argsort(coord)
    return coord[order], prof[order], parts[0][2]


def utau_of(coord, prof, re):
    """u_tau from the wall gradient of the mean profile (both walls)."""
    du0 = prof[0, STAT_U]/coord[0]
    du1 = prof[-1, STAT_U]/(2.0 - coord[-1])
    return np.sqrt(0.5*(du0 + du1)/re)


def wall_halves(coord, prof):
    """Split into the lower (y in (0, 1]) and upper ([1, 2)) halves, each mapped
    to wall distance, WITHOUT averaging the two together: the two interface
    orientations differ, so they are kept as separate curves. <uv> is negated in
    the upper half so -<u'v'> reads positive against wall distance in both halves.

    Returns [(y_lower, prof_lower), (y_upper, prof_upper)]."""
    halves = []
    for mask, ywall, uv_sign in ((coord <= 1.0, coord, +1.0),
                                 (coord > 1.0, 2.0 - coord, -1.0)):
        y = ywall[mask]
        p = prof[mask].copy()
        p[:, STAT_UV] *= uv_sign
        order = np.argsort(y)
        halves.append((y[order], p[order]))
    return halves


def plot_profiles(ref, refined, label, outdir, iface_y):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    (yr, pr, re_r), (yf, pf, re_f) = ref, refined
    ut_r = utau_of(yr, pr, re_r)
    ut_f = utau_of(yf, pf, re_f)
    print(f"u_tau: reference {ut_r:.5f}, {label} {ut_f:.5f} "
          f"(dev {100*(ut_f/ut_r-1):+.2f}%)")

    # Keep the two channel halves as separate curves so the bottom- and
    # top-wall interfaces can be compared (the interface error is wall-asymmetric).
    ref_h = wall_halves(yr, pr)
    fin_h = wall_halves(yf, pf)
    HALVES = ("lower", "upper")
    REF_LS = {"lower": "-", "upper": "--"}   # reference line style per half
    FIN_MK = {"lower": ".", "upper": "x"}    # refined marker per half

    def yplus(y, ut, re): return y*ut*re

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    ax = axes[0, 0]
    for hi, hn in enumerate(HALVES):
        yrh, prh = ref_h[hi]
        yfh, pfh = fin_h[hi]
        ax.semilogx(yplus(yrh, ut_r, re_r), prh[:, STAT_U]/ut_r, "k"+REF_LS[hn],
                    label=f"reference {hn}")
        ax.semilogx(yplus(yfh, ut_f, re_f), pfh[:, STAT_U]/ut_f, "r"+FIN_MK[hn],
                    ms=3, label=f"{label} {hn}")
    ax.axvline(iface_y*ut_f*re_f, color="b", ls=":", lw=1, label="interface")
    ax.set_xlabel("y+"); ax.set_ylabel("U+"); ax.legend(fontsize=8); ax.set_title("mean velocity")

    ax = axes[0, 1]
    for hi, hn in enumerate(HALVES):
        yrh, prh = ref_h[hi]
        yfh, pfh = fin_h[hi]
        for s, name, c in ((STAT_UU, "u", "C0"), (STAT_VV, "v", "C1"), (STAT_WW, "w", "C2")):
            rms_r = np.sqrt(np.maximum(prh[:, s] - prh[:, s-3]**2, 0))/ut_r
            rms_f = np.sqrt(np.maximum(pfh[:, s] - pfh[:, s-3]**2, 0))/ut_f
            ax.plot(yplus(yrh, ut_r, re_r), rms_r, c+REF_LS[hn],
                    label=(f"{name}rms" if hi == 0 else None))
            ax.plot(yplus(yfh, ut_f, re_f), rms_f, c+FIN_MK[hn], ms=3)
    ax.axvline(iface_y*ut_f*re_f, color="b", ls=":", lw=1)
    ax.set_xlabel("y+"); ax.set_ylabel("rms+"); ax.legend(fontsize=8)
    ax.set_title("rms (line/dot: lower, dash/x: upper)")

    ax = axes[1, 0]
    for hi, hn in enumerate(HALVES):
        yrh, prh = ref_h[hi]
        yfh, pfh = fin_h[hi]
        uv_r = -(prh[:, STAT_UV] - prh[:, STAT_U]*prh[:, STAT_V])/ut_r**2
        uv_f = -(pfh[:, STAT_UV] - pfh[:, STAT_U]*pfh[:, STAT_V])/ut_f**2
        ax.plot(yplus(yrh, ut_r, re_r), uv_r, "k"+REF_LS[hn], label=f"reference {hn}")
        ax.plot(yplus(yfh, ut_f, re_f), uv_f, "r"+FIN_MK[hn], ms=3, label=f"{label} {hn}")
    ax.axvline(iface_y*ut_f*re_f, color="b", ls=":", lw=1)
    ax.set_xlabel("y+"); ax.set_ylabel("-<u'v'>+"); ax.legend(fontsize=8); ax.set_title("Reynolds stress")

    ax = axes[1, 1]
    Umax = max(abs(pr[:, STAT_U]).max(), 1e-30)
    dev_summary = []
    for hi, hn in enumerate(HALVES):
        yrh, prh = ref_h[hi]
        yfh, pfh = fin_h[hi]
        Ui = np.interp(yfh, yrh, prh[:, STAT_U])      # reference of the SAME half
        dev = 100*(pfh[:, STAT_U] - Ui)/Umax
        ax.plot(yplus(yfh, ut_f, re_f), dev, "r"+REF_LS[hn], label=hn)
        k = np.abs(dev).argmax()
        dev_summary.append((hn, np.abs(dev[k]), yplus(yfh, ut_f, re_f)[k]))
    ax.axhline(0.0, color="k", lw=0.5)
    ax.axvline(iface_y*ut_f*re_f, color="b", ls=":", lw=1)
    ax.set_xlabel("y+"); ax.set_ylabel("dU / max(U) [%]")
    ax.legend(fontsize=8); ax.set_title("mean velocity deviation (per half)")

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "profiles.png"), dpi=150)
    print(f"wrote {outdir}/profiles.png")
    for hn, dmax, yloc in dev_summary:
        print(f"mean-U deviation ({hn} wall): max {dmax:.3f}% of centreline, "
              f"at y+ = {yloc:.1f}")
    return ut_r, ut_f


# ---------------------------------------------------------------------------
# Snapshots: planes, spectra, divergence

def snapshot_files(rundir, prefix="channel_field", last=10):
    files = sorted(glob.glob(os.path.join(rundir, prefix + "_*.h5")),
                   key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
    return files[-last:]


class Snapshot:
    """Field snapshot: legacy global 3D or block-table layout."""

    def __init__(self, path):
        self.h5 = h5py.File(path, "r")
        self.block = self.h5["un"].ndim == 4
        if self.block:
            self.blocks = self.h5["blocks"][...]
            self.nb = int(self.h5.attrs["block_nb_x"])
            self.index = {}
            for bid, (ox, oy, oz, lev) in enumerate(self.blocks):
                self.index[(int(lev), ox//self.nb, oy//self.nb, oz//self.nb)] = bid
        self.ny = int(self.h5.attrs["ny"])
        self.y = self.h5["y"][...]

    def uplane(self, var, level, jrow):
        """Global (z, x) plane of `var` at level-l y-cell row jrow (1-based),
        levels as stored (reference: level irrelevant)."""
        if not self.block:
            return self.h5[var][:, jrow-1, :]
        nb = self.nb
        attrs = self.h5.attrs
        nxl = int(attrs["nx"])*2**level
        nzl = int(attrs["nz"])*2**level
        plane = np.full((nzl, nxl), np.nan)
        cy, j = (jrow-1)//nb, (jrow-1) % nb
        data = self.h5[var]
        for (lev, bx, by, bz), bid in self.index.items():
            if lev != level or by != cy:
                continue
            row = data[bid][:, j, :]
            plane[bz*nb:(bz+1)*nb, bx*nb:(bx+1)*nb] = row
        return plane

    def close(self):
        self.h5.close()


def spectra_at(files, level, jrow, lx, lz):
    """Snapshot+homogeneous-direction averaged 1D u-spectra at a y row."""
    ex_acc = ez_acc = None
    n = 0
    for path in files:
        snap = Snapshot(path)
        plane = snap.uplane("un", level, jrow)
        snap.close()
        if np.isnan(plane).any():
            continue
        nz, nx = plane.shape
        up = plane - plane.mean()
        ex = np.mean(np.abs(np.fft.rfft(up, axis=1))**2, axis=0)/nx**2
        ez = np.mean(np.abs(np.fft.rfft(up, axis=0))**2, axis=1)/nz**2
        ex_acc = ex if ex_acc is None else ex_acc + ex
        ez_acc = ez if ez_acc is None else ez_acc + ez
        n += 1
    if n == 0:
        return None
    nz, nx = plane.shape
    kx = 2*np.pi/lx*np.arange(ex_acc.size)
    kz = 2*np.pi/lz*np.arange(ez_acc.size)
    return (kx, ex_acc/n), (kz, ez_acc/n)


def row_of_y(line, y):
    """1-based cell row whose centre is nearest to y on node line `line`."""
    cent = 0.5*(line[:-1] + line[1:])
    return int(np.argmin(np.abs(cent - y))) + 1


def plot_spectra(ref_files, fine_files, y_targets, outdir, lx, lz,
                 base_y, fine_y, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(y_targets), 2, figsize=(11, 4*len(y_targets)),
                             squeeze=False)
    for row, (name, y, level) in enumerate(y_targets):
        line = fine_y if level == 1 else base_y
        jrow = row_of_y(line, y)
        ref_jrow = row_of_y(fine_y, y)
        sp_ref = spectra_at(ref_files, 0, ref_jrow, lx, lz)
        sp_fin = spectra_at(fine_files, level, jrow, lx, lz)
        for col, comp in enumerate(("kx", "kz")):
            ax = axes[row][col]
            if sp_ref:
                k, e = sp_ref[col]
                ax.loglog(k[1:], e[1:], "k-", label="reference")
            if sp_fin:
                k, e = sp_fin[col]
                ax.loglog(k[1:], e[1:], "r--", label=f"{label} (lvl {level})")
            ax.set_xlabel(comp); ax.set_ylabel("E_uu")
            ax.set_title(f"{name}: y = {y:.4f}")
            ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "spectra.png"), dpi=150)
    print(f"wrote {outdir}/spectra.png")


# ---------------------------------------------------------------------------
# Divergence residuals (block-format snapshots)

def divergence_zones(files, base_y, fine_y, iface_lo, iface_hi, lx, lz):
    """rms/max divergence in the interface band vs interiors, averaged
    over snapshots. Interface band: cells within 2 fine cells of either
    interface height."""
    stats = {"interface": [], "fine interior": [], "coarse interior": []}
    for path in files:
        snap = Snapshot(path)
        if not snap.block:
            snap.close()
            continue
        h5 = snap.h5
        nb = snap.nb
        u, v, w = (h5[k] for k in ("un", "vn", "wn"))
        lines = {0: base_y, 1: fine_y}
        nxl = {l: int(h5.attrs["nx"])*2**l for l in (0, 1)}
        dxl = {l: lx/nxl[l] for l in (0, 1)}
        dzl = {l: lz/(int(h5.attrs["nz"])*2**l) for l in (0, 1)}
        # +-2 local fine cells around each interface height
        j_if = np.searchsorted(fine_y, iface_lo)
        band = 2*(fine_y[j_if] - fine_y[j_if-1])

        def face(varr, key, comp, sub):
            """East/north/top closure face plane for a block (z,y,x order)."""
            lev, bx, by, bz = key
            if key in snap.index:
                arr = varr[snap.index[key]]
                return arr[:, :, 0] if comp == 0 else (arr[:, 0, :] if comp == 1 else arr[0, :, :])
            # finer children across the face: 2x2 restriction (the solver's
            # conservative owner value)
            fkey = lambda sx, sy, sz: (lev+1, 2*bx+sx, 2*by+sy, 2*bz+sz)
            if lev + 1 <= 1 and fkey(*sub[0]) in snap.index:
                quad = np.zeros((2*nb, 2*nb))
                for sx, sy, sz in sub:
                    arr = varr[snap.index[fkey(sx, sy, sz)]]
                    p = arr[:, :, 0] if comp == 0 else (arr[:, 0, :] if comp == 1 else arr[0, :, :])
                    # plane axes per component: u (z,y), v (z,x), w (y,x)
                    a0 = sz if comp != 2 else sy
                    a1 = sy if comp == 0 else sx
                    quad[a0*nb:(a0+1)*nb, a1*nb:(a1+1)*nb] = p
                return 0.25*(quad[0::2, 0::2] + quad[1::2, 0::2]
                             + quad[0::2, 1::2] + quad[1::2, 1::2])
            # coarser parent: inject the covering coarse values
            pkey = (lev-1, bx//2, by//2, bz//2)
            if lev - 1 >= 0 and pkey in snap.index:
                arr = varr[snap.index[pkey]]
                p = arr[:, :, 0] if comp == 0 else (arr[:, 0, :] if comp == 1 else arr[0, :, :])
                # take the covering half and repeat
                ta = (bz % 2)*nb//2 if comp != 2 else (by % 2)*nb//2
                tb = (by % 2)*nb//2 if comp == 0 else (bx % 2)*nb//2
                sub_p = p[ta:ta+nb//2, tb:tb+nb//2]
                return np.repeat(np.repeat(sub_p, 2, axis=0), 2, axis=1)
            return None  # wall (v at the top wall: zero flux)

        for (lev, bx, by, bz), bid in snap.index.items():
            uu, vv, ww = u[bid], v[bid], w[bid]
            ue = face(u, (lev, bx+1 if bx+1 < nxl[lev]//nb else 0, by, bz), 0,
                      [(0, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1)])
            nzb = int(h5.attrs["nz"])*2**lev//nb
            wt = face(w, (lev, bx, by, (bz+1) % nzb), 2,
                      [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)])
            vkey = (lev, bx, by+1, bz)
            nyb = int(h5.attrs["ny"])*2**lev//nb
            vn_ = face(v, vkey, 1, [(0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 1)]) \
                if by + 1 < nyb else np.zeros((nb, nb))
            if ue is None or wt is None or vn_ is None:
                continue
            ufull = np.concatenate((uu, ue[:, :, None]), axis=2)
            vfull = np.concatenate((vv, vn_[:, None, :]), axis=1)
            wfull = np.concatenate((ww, wt[None, :, :]), axis=0)
            yline = lines[lev]
            oy = by*nb
            dy = np.diff(yline)[oy:oy+nb][None, :, None]
            div = (np.diff(ufull, axis=2)/dxl[lev]
                   + np.diff(vfull, axis=1)/dy
                   + np.diff(wfull, axis=0)/dzl[lev])
            ycent = 0.5*(yline[oy:oy+nb] + yline[oy+1:oy+nb+1])
            in_band = (np.abs(ycent - iface_lo) < band) | (np.abs(ycent - iface_hi) < band)
            for j in range(nb):
                zone = "interface" if in_band[j] else \
                    ("fine interior" if lev == 1 else "coarse interior")
                stats[zone].append(div[:, j, :].ravel())
        snap.close()
    out = {}
    for zone, chunks in stats.items():
        if chunks:
            allv = np.concatenate(chunks)
            out[zone] = (np.sqrt(np.mean(allv**2)), np.abs(allv).max())
    return out


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, help="reference stats-leg dir")
    parser.add_argument("--refined", required=True, help="refined stats-leg dir")
    parser.add_argument("--out", default="plots")
    parser.add_argument("--snapshots", type=int, default=10)
    parser.add_argument("--label", default="refined")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ref = composite_profiles(args.reference)
    refined = composite_profiles(args.refined)

    # grid lines and interface height from a refined snapshot
    fin_files = snapshot_files(args.refined, last=args.snapshots)
    ref_files = snapshot_files(args.reference, last=args.snapshots)
    with h5py.File(fin_files[-1], "r") as h5:
        base_y = h5["y"][...]
        lx, lz = float(h5.attrs["lx"]), float(h5.attrs["lz"])
        bl = h5["blocks"][...]
        nb = int(h5.attrs["block_nb_x"])
        fine_rows = np.unique(bl[bl[:, 3] == 1][:, 1])
        ny_base = int(h5.attrs["ny"])
        bottom = fine_rows[fine_rows < ny_base]  # level-1 origins of the bottom band
        band_cells = (bottom.max() + nb)//2      # band height in base cells
    fine_y = np.empty(2*(base_y.size-1)+1)
    fine_y[0::2] = base_y
    fine_y[1::2] = 0.5*(base_y[:-1] + base_y[1:])
    iface_lo = base_y[band_cells]
    iface_hi = 2.0 - iface_lo
    print(f"interface heights: y = {iface_lo:.5f} / {iface_hi:.5f} "
          f"(y+ = {iface_lo*180:.1f})")

    ut_r, ut_f = plot_profiles(ref, refined, args.label, args.out, iface_lo)

    y15 = 15.0/(ut_f*180.0)
    dyf = np.diff(fine_y)
    targets = [
        ("y+ ~ 15 (fine)", y15, 1),
        ("below interface (fine)", iface_lo - 0.5*dyf[2*band_cells-1], 1),
        ("above interface (coarse)", iface_lo + 0.5*(base_y[band_cells+1]-base_y[band_cells]), 0),
    ]
    plot_spectra(ref_files, fin_files, targets, args.out, lx, lz,
                 base_y, fine_y, args.label)

    print("divergence residuals (refined, snapshot-averaged):")
    zones = divergence_zones(fin_files, base_y, fine_y, iface_lo, iface_hi, lx, lz)
    for zone, (rms, mx) in zones.items():
        print(f"  {zone:16s} rms {rms:.3e}  max {mx:.3e}")


if __name__ == "__main__":
    main()
