#!/usr/bin/env python3
"""Checkers for the passive-scalar S2 (turbulent closure) gates.

Subcommands
  channel  developed turbulent channel with isothermal antisymmetric walls:
           time- and (x,z)-averaged theta+(y+) against Kader's correlation,
           the resolved/SGS/molecular wall-normal heat-flux split, and the
           CONSTANT-total-flux consistency check that the setup guarantees in
           statistical steady state.  Single level (the uniform LES control
           and any unrefined case).
  band     the interface-band metric of a 2:1-refined case against the
           matched unrefined control: the theta'_rms ratio on the coarse rows
           adjacent to the interface vs the core rows (the method of
           tools/patch_interface_stats.py, and the scalar analogue of the
           u'/v' band lesson in ../channel_interface/README.md).
  rans     steady 1D RANS profile: theta+(y+), its log-layer slope (which
           must be Pr_t/kappa) and intercept, against Kader.

The scalar setup these all assume (turbles.ini / turbslab.ini / turbsst.ini):
theta = +1 at y = 0, -1 at y = ly, no source, so in steady state the total
wall-normal flux J = <v'theta'> - (D + nut/Pr_t) d<theta>/dy is constant and
theta_tau = J/u_tau = J (u_tau = 1 at Re_tau 180 with forcing_x = 1).
"""

from __future__ import annotations

import argparse
import glob
import sys

import h5py
import numpy as np

from scalar_tools import BlockGeometry


# ------------------------------------------------------- Pr_t correlation ---
def prt_kays(pet, prtinf):
    """Kays-Crawford Pr_t -- the Python mirror of scalar.f90's prt_kays
    (same branch structure, so the checker models what the kernel did)."""
    pet = np.maximum(np.asarray(pet, dtype=float), 0.0)
    b = np.sqrt(prtinf)
    a = 0.3 * pet
    ab = np.maximum(a * b, 1e-300)
    x = 1.0 / ab
    small = x < 0.5
    inv = np.empty_like(x)
    # series branch, 1/Prt_inf [1 - x/3! + x^2/4! - ...]
    xs = x[small]
    tot = np.ones_like(xs)
    u = -xs / 6.0
    for m in range(1, 17):
        tot = tot + u
        u = -u * xs / (m + 3)
    inv[small] = tot / prtinf
    xd, ad = x[~small], a[~small]
    inv[~small] = 0.5 / prtinf + ad / b - ad * ad * (1.0 - np.exp(-xd))
    out = 1.0 / inv
    return np.where(a * b < 1e-300, 2.0 * prtinf, out)


def face_prt(nutf, pr, prt, model):
    """Pr_t at a face, the kernel's eddy_diffusivity() choice."""
    if model == "kays":
        return prt_kays(np.maximum(nutf, 0.0) * PRT_RE[0] * pr, prt)
    return prt


PRT_RE = [180.0]   # set by main(); the correlation needs Pe_t = nut Re Pr


# ----------------------------------------------------------------- Kader ---
def kader(yplus: np.ndarray, pr: float) -> np.ndarray:
    """Kader (1981) theta+(y+) for a constant-flux thermal layer."""
    beta = (3.85 * pr ** (1.0 / 3.0) - 1.3) ** 2 + 2.12 * np.log(pr)
    gam = 0.01 * (pr * yplus) ** 4 / (1.0 + 5.0 * pr ** 3 * yplus)
    with np.errstate(divide="ignore", over="ignore"):
        return pr * yplus * np.exp(-gam) + (2.12 * np.log1p(yplus) + beta) * np.exp(-1.0 / gam)


# --------------------------------------------------------------- geometry ---
def level_bands(tab, nb, level):
    """Leaf ids of one level, grouped into CONTIGUOUS y-bands.

    A level need not be a single box: the wall-band-refined channel puts its
    level-1 leaves in two disjoint slabs, one per wall.  Each band IS a filled
    box, which is what the assembly below needs.
    """
    ids = np.flatnonzero(tab[:, 3] == level)
    if ids.size == 0:
        return []
    oy = np.unique(tab[ids, 1])
    groups, cur = [], [oy[0]]
    for prev, nxt in zip(oy[:-1], oy[1:]):
        if nxt - prev == nb[1]:
            cur.append(nxt)
        else:
            groups.append(cur)
            cur = [nxt]
    groups.append(cur)
    return [ids[np.isin(tab[ids, 1], g)] for g in groups]



class LevelGrid:
    """Global cell arrays of ONE refinement level of a block-table file."""

    def __init__(self, f: h5py.File, level: int, band: int = 0):
        self.geo = BlockGeometry(f)
        self.level = level
        self.nb = self.geo.nb
        tab = self.geo.blocks
        bands = level_bands(tab, self.nb, level)
        if not bands:
            raise SystemExit(f"no leaves at level {level}")
        if band >= len(bands):
            raise SystemExit(f"level {level} has {len(bands)} y-band(s), asked for {band}")
        self.ids = bands[band]
        self.nbands = len(bands)
        origins = tab[self.ids, :3]
        self.lo = origins.min(axis=0)
        span = (origins.max(axis=0) - self.lo) // np.array(self.nb) + 1
        self.shape = (span[2] * self.nb[2], span[1] * self.nb[1], span[0] * self.nb[0])
        self.slots = [(int(b), tuple((origins[n] - self.lo) // np.array(self.nb)))
                      for n, b in enumerate(self.ids)]
        # y cell centres / node line of this level, restricted to the covered rows
        line = self.geo.lines[1][level]
        j0 = int(self.lo[1])
        self.ynode = line[j0:j0 + self.shape[1] + 1]
        self.y = 0.5 * (self.ynode[:-1] + self.ynode[1:])
        self.ly = float(f.attrs["ly"])
        # does this level touch the domain walls?
        self.at_lo_wall = abs(self.ynode[0]) < 1e-12
        self.at_hi_wall = abs(self.ynode[-1] - self.ly) < 1e-12

    def assemble(self, dset) -> np.ndarray:
        out = np.full(self.shape, np.nan)
        for bid, (bx, by, bz) in self.slots:
            block = dset[bid]
            out[bz * self.nb[2]:(bz + 1) * self.nb[2],
                by * self.nb[1]:(by + 1) * self.nb[1],
                bx * self.nb[0]:(bx + 1) * self.nb[0]] = block
        if np.isnan(out).any():
            raise SystemExit(f"level {self.level} is not a filled box "
                             f"({int(np.isnan(out).sum())} empty cells)")
        return out


# ------------------------------------------------------------ statistics ---
def row_stats(paths, scalar):
    """(y, <theta>, theta'_rms, level) per stored cell ROW, x/z-averaged and
    averaged over snapshots.

    Works on ANY leaf layout -- it accumulates per row instead of assembling a
    global box, which the 2:1 wall-band case needs: its level-1 leaves cover
    TWO disjoint slabs, so no single level is a filled box.
    """
    acc = {}
    for path in paths:
        with h5py.File(path, "r") as f:
            geo = BlockGeometry(f)
            data = f[scalar][...]
            lev = geo.blocks[:, 3]
            for bid in range(geo.n_blocks):
                (_, _), (yc, _), (_, _) = geo.block_axes(bid)
                block = data[bid]
                n = block.shape[0] * block.shape[2]
                for jj, yv in enumerate(yc):
                    key = round(float(yv), 12)
                    plane = block[:, jj, :]
                    s1, s2, cnt, _ = acc.get(key, (0.0, 0.0, 0, 0))
                    acc[key] = (s1 + float(plane.sum()),
                                s2 + float((plane ** 2).sum()), cnt + n, int(lev[bid]))
    ys = np.array(sorted(acc))
    mean = np.array([acc[y][0] / acc[y][2] for y in ys])
    var = np.array([acc[y][1] / acc[y][2] for y in ys]) - mean ** 2
    lev = np.array([acc[y][3] for y in ys])
    return ys, mean, np.sqrt(np.maximum(var, 0.0)), lev



def snapshot_stats(path, level, scalar, pr, prt, re, walls, band=0, model="constant"):
    """Per-snapshot (x,z)-averaged profiles on one level.

    Returns a dict of y-profiles: theta, theta_rms, and (when the level's
    v/nut stencils are complete) the face-based flux split.
    """
    with h5py.File(path, "r") as f:
        g = LevelGrid(f, level, band)
        th = g.assemble(f[scalar])
        vn = g.assemble(f["vn"])
        nut = g.assemble(f["nut"]) if "nut" in f else np.zeros_like(th)
        t = float(f.attrs["t_current"])

    ny = th.shape[1]
    thm = th.mean(axis=(0, 2))
    thp = th - thm[None, :, None]
    rms = np.sqrt((thp ** 2).mean(axis=(0, 2)))

    # Face quantities at the ny+1 y-faces of this level's box.  vn[:, j, :] is
    # the LOW face of cell j; the high face of the last cell is the domain
    # wall (v = 0) when the level reaches it, and the 2:1 interface otherwise
    # (owned by the block below -- not stored here, hence the nan).
    d = 1.0 / (re * pr)
    ycent = g.y
    yface = g.ynode
    flux = np.full(ny + 1, np.nan)
    turbf = np.full(ny + 1, np.nan)
    nutface = np.zeros(ny + 1)
    for j in range(ny + 1):
        if j == 0:
            if not g.at_lo_wall:
                continue
            # wall face: v = 0 and the mirrored ghost gives d(theta)/dy =
            # (theta_0 - theta_w)/y_0
            grad = (thm[0] - walls[0]) / (ycent[0] - yface[0])
            nutf = 0.0
            turbf[j] = 0.0
        elif j == ny:
            if not g.at_hi_wall:
                continue
            grad = (walls[1] - thm[-1]) / (yface[-1] - ycent[-1])
            nutf = 0.0
            turbf[j] = 0.0
        else:
            grad = (thm[j] - thm[j - 1]) / (ycent[j] - ycent[j - 1])
            nutf = 0.5 * (nut[:, j - 1, :] + nut[:, j, :])
            thface = 0.5 * (th[:, j - 1, :] + th[:, j, :])
            turbf[j] = float((vn[:, j, :] * thface).mean())
            nutf = float(nutf.mean())
        flux[j] = turbf[j] - (d + nutf / face_prt(np.array(nutf), pr, prt, model)) * grad
        nutface[j] = nutf
    return dict(t=t, y=ycent, yface=yface, theta=thm, rms=rms,
                flux=flux, turb=turbf, nutf=nutface)


def average_stats(paths, level, scalar, pr, prt, re, walls, band=0, model="constant"):
    acc = None
    for p in paths:
        s = snapshot_stats(p, level, scalar, pr, prt, re, walls, band, model)
        if acc is None:
            acc = {k: (np.array(v, dtype=float) if k != "t" else [v]) for k, v in s.items()}
            acc["y"] = s["y"]
            acc["yface"] = s["yface"]
            acc["n"] = 1
        else:
            for k in ("theta", "rms", "flux", "turb", "nutf"):
                acc[k] = acc[k] + s[k]
            acc["t"].append(s["t"])
            acc["n"] += 1
    for k in ("theta", "rms", "flux", "turb", "nutf"):
        acc[k] = acc[k] / acc["n"]
    return acc


# ------------------------------------------------------------ subcommands ---
def cmd_channel(a) -> int:
    """Developed turbulent channel: theta+(y+), the flux split, convergence.

    WINDOW: Kader's correlation describes the CONSTANT-FLUX wall layer.  Our
    scalar carries a constant flux everywhere, but the momentum field driving
    it does not -- the stress falls linearly to zero at the centreline -- so
    beyond y/h ~ 0.2 (y+ ~ 35 at Re_tau 180) the two profiles diverge for a
    reason that has nothing to do with the closure: the thermal profile keeps
    a log-like slope where U+ is already flattening into the wake.  The Kader
    comparison is therefore made over the wall layer, and the outer profile
    is reported, not gated.
    """
    paths = sorted_snapshots(a.snapshots, a.skip)
    st = average_stats(paths, a.level, a.scalar, a.pr, a.prt, a.re, a.walls)
    with h5py.File(paths[0], "r") as f:
        ly = float(f.attrs["ly"])

    # J(y) is the SAME sign at both walls (one constant flux crossing the
    # channel from the hot wall to the cold one), so the two wall values are
    # averaged, not differenced.
    jw = 0.5 * (st["flux"][0] + st["flux"][-1])
    theta_tau = abs(jw)
    y = st["y"]
    thp = np.where(y <= 0.5 * ly, (a.walls[0] - st["theta"]) / theta_tau,
                   (st["theta"] - a.walls[1]) / theta_tau)
    ywall = np.where(y <= 0.5 * ly, y, ly - y) * a.re
    order = np.argsort(ywall)
    yw, thp = ywall[order], thp[order]

    print(f"snapshots {len(paths)}  t = {min(st['t']):.3f} .. {max(st['t']):.3f}")
    print(f"theta_tau = {theta_tau:.6f}  ->  theta+ at the centreline "
          f"{1.0/theta_tau*abs(a.walls[0] - 0.5*(a.walls[0]+a.walls[1])):.2f}")

    # (1) constant-flux consistency: the setup makes J(y) constant in steady
    # state, so this measures statistical convergence AND the transport
    # operator's conservation in one number.
    jj = st["flux"][np.isfinite(st["flux"])]
    dev = float(np.max(np.abs(jj - jw))) / theta_tau
    print(f"total flux J(y): {jj.min():.6f} .. {jj.max():.6f}, "
          f"max|J - J_wall|/theta_tau = {dev:.4f}")

    # (2) stationarity: the same quantity from the two halves of the window
    if len(paths) >= 4:
        h = len(paths) // 2
        for tag, sub in (("first half", paths[:h]), ("second half", paths[h:])):
            s2 = average_stats(sub, a.level, a.scalar, a.pr, a.prt, a.re, a.walls)
            j2 = 0.5 * (s2["flux"][0] + s2["flux"][-1])
            print(f"   {tag}: theta_tau = {abs(j2):.6f}")

    # (3) theta+ against Kader over the wall layer
    ref = kader(yw, a.pr)
    win = (yw >= a.ylo) & (yw <= a.yhi)
    rel = np.abs(thp[win] - ref[win]) / ref[win]
    print(f"theta+ vs Kader over y+ in [{a.ylo:g}, {a.yhi:g}]: "
          f"max rel dev = {rel.max():.4f} (at y+ {yw[win][np.argmax(rel)]:.1f}), "
          f"mean = {rel.mean():.4f}")
    sub = yw <= 3.0
    r = thp[sub] / (a.pr * yw[sub])
    print(f"viscous sublayer theta+/(Pr y+) at y+ <= 3: {r.min():.4f} .. {r.max():.4f}")

    # (4) the Reynolds-analogy limit: theta+/U+ -> Pr at the wall (the
    # molecular half of the face diffusivity) and rises through Pr_t.
    _, um, _, _ = row_stats(paths, "un")
    up = um[order]
    print(f"theta+/U+ : first cell {thp[0]/up[0]:.4f} (Pr = {a.pr}), "
          f"y+ 20-40 {np.mean(thp[(yw>20)&(yw<40)]/up[(yw>20)&(yw<40)]):.4f}")

    # (5) the turbulent heat flux: fraction of J carried by the resolved
    # <v'theta'> (the rest is molecular + SGS).
    frac = st["turb"] / jw
    yf = np.where(st["yface"] <= 0.5 * ly, st["yface"], ly - st["yface"]) * a.re
    print("resolved <v'theta'>/J at y+ = "
          + ", ".join(f"{t:.0f}:{np.interp(t, yf[:yf.size//2], frac[:yf.size//2]):.2f}"
                      for t in (5, 15, 30, 60, 120)))
    kres = int(np.nanargmax(np.abs(st["turb"])))
    print(f"  peak {st['turb'][kres]:+.6f} = {abs(st['turb'][kres])/theta_tau:.3f} "
          f"theta_tau u_tau at y+ {yf[kres]:.0f}")

    if a.dump:
        np.savetxt(a.dump, np.column_stack([yw, thp, ref, up]),
                   header="y+ theta+ kader U+")
        print(f"profile written to {a.dump}")

    ok = dev <= a.flux_tolerance and rel.max() <= a.tolerance
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_band(a) -> int:
    """Interface-band metric, the tools/patch_interface_stats.py method.

    A 2:1 interface BAND is a LOCALIZED excess of the fluctuation rms in the
    cells adjacent to the interface, measured against a matched unrefined
    control (which controls for the strong y-dependence of the statistics).
    The comparison is made on the COARSE rows only, where the refined run and
    the control share the same cells exactly: the fine rows carry more
    resolved scales by construction, so a fine/coarse rms difference there is
    physics, not a band.
    """
    ref_paths = sorted_snapshots(a.reference, a.skip)
    paths = sorted_snapshots(a.snapshots, a.skip)

    ys, th, rms, lev = row_stats(paths, a.scalar)
    yc, thc, rmsc, _ = row_stats(ref_paths, a.scalar)
    levels = sorted(set(lev.tolist()))
    print(f"{len(paths)} refined / {len(ref_paths)} control snapshots, levels {levels}")

    changes = np.flatnonzero(lev[1:] != lev[:-1])
    print("2:1 interfaces at y+ = "
          + ", ".join("%.1f" % (0.5 * (ys[i] + ys[i + 1]) * a.re) for i in changes))

    # coarse rows, matched to the control cell by cell
    coarse = np.flatnonzero(lev == min(levels))
    match = np.array([int(np.argmin(np.abs(yc - ys[i]))) for i in coarse])
    dy = np.max(np.abs(yc[match] - ys[coarse]))
    if dy > 1e-12:
        raise SystemExit(f"coarse rows do not match the control (max dy = {dy:.3e})")
    ratio = rms[coarse] / rmsc[match]
    dmean = (th[coarse] - thc[match])

    # distance (in coarse rows) from the nearest interface
    edges = [0, coarse.size - 1]
    dist = np.minimum.reduce([np.abs(np.arange(coarse.size) - e) for e in edges])

    near = dist <= 1
    far = dist >= 4
    r_near, r_far = float(ratio[near].mean()), float(ratio[far].mean())
    m_near = float(np.max(np.abs(dmean[near])))
    m_far = float(np.max(np.abs(dmean[far])))
    print(f"theta'_rms band ratio (refined/control): adjacent to the interface "
          f"{r_near:.4f}, core {r_far:.4f}, excess {r_near/r_far - 1.0:+.4f}")
    print(f"  per-row ratios from the interface outward: "
          + " ".join(f"{ratio[dist == d].mean():.3f}" for d in range(6)))
    print(f"<theta> footprint |refined - control|: adjacent {m_near:.5f}, "
          f"core {m_far:.5f}  (wall difference 1.0)")

    if a.dump:
        np.savetxt(a.dump, np.column_stack([ys * a.re, th, rms, lev]),
                   header="y+ theta rms level (refined run)")
        np.savetxt(a.dump.replace(".dat", "_control.dat"),
                   np.column_stack([yc * a.re, thc, rmsc]), header="y+ theta rms (control)")
        print(f"profiles written to {a.dump}")

    ok = abs(r_near / r_far - 1.0) <= a.tolerance
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_rans(a) -> int:
    """Steady 1D RANS: the scalar profile against the EXACT prediction of the
    model's own eddy viscosity.

    RANS collapses to a 1D fixed point with no resolved fluctuations, so the
    scalar equation reduces to one ODE whose coefficient is the nut the
    solver wrote into the snapshot:

        d(theta)/dy = -J/(1/(Re Pr) + nut(y)/Pr_t(y)),   J from the walls

    Integrating that with the file's own nut and comparing to the file's own
    theta tests EXACTLY what increment S2 added -- the face diffusivity and
    the Pr_t model -- with the turbulence model divided out.  Discretisation
    (the trapezoid here vs the solver's face-flux operator) is the only
    expected difference.

    The log-layer slope is reported for physical context, with the reference
    Pr_t/(kappa (1 - y/h)): the scalar carries a CONSTANT flux while the
    channel's momentum flux falls linearly, so the classic Pr_t/kappa applies
    only as y/h -> 0.
    """
    st = average_stats([a.snapshots[0]], 0, a.scalar, a.pr, a.prt, a.re,
                       a.walls, model=a.prt_model)
    with h5py.File(a.snapshots[0], "r") as f:
        ly = float(f.attrs["ly"])
    jw = 0.5 * (st["flux"][0] + st["flux"][-1])
    theta_tau = abs(jw)
    y, yf = st["y"], st["yface"]

    # the ODE prediction on the face grid (nut at the walls is 0)
    d = 1.0 / (a.re * a.pr)
    nutf = st["nutf"].copy()
    nutf[0] = 0.0
    nutf[-1] = 0.0
    dtot = d + nutf / face_prt(nutf, a.pr, a.prt, a.prt_model)
    integ = np.concatenate([[0.0], np.cumsum(np.diff(yf) / (0.5 * (dtot[1:] + dtot[:-1])))])
    # J that satisfies both wall values, then the profile at the cell centres
    jpred = (a.walls[0] - a.walls[1]) / integ[-1]
    thpred = np.interp(y, yf, a.walls[0] - jpred * integ)

    print(f"{a.snapshots[0]}: scalar '{a.scalar}', Pr = {a.pr}, "
          f"Pr_t = {a.prt} ({a.prt_model})")
    print(f"theta_tau = {theta_tau:.6f}  (predicted {jpred:.6f}, "
          f"dev {abs(jpred-theta_tau)/theta_tau*100:.3f}%)")
    jj = st["flux"][np.isfinite(st["flux"])]
    print(f"total flux J(y) constancy: max|J - J_wall|/theta_tau = "
          f"{float(np.max(np.abs(jj - jw)))/theta_tau:.3e}")
    err = float(np.max(np.abs(st["theta"] - thpred)))
    print(f"theta vs the nut-integral prediction: max|dev| = {err:.3e} "
          f"({err/abs(a.walls[0]-a.walls[1])*100:.4f}% of the wall difference)")
    prtf = np.broadcast_to(face_prt(nutf, a.pr, a.prt, a.prt_model), nutf.shape)
    print(f"Pr_t over the faces: {prtf.min():.4f} .. {prtf.max():.4f} "
          f"(Pe_t = nut Re Pr up to {nutf.max()*a.re*a.pr:.3g})")

    ywall = np.where(y <= 0.5 * ly, y, ly - y) * a.re
    thp = np.where(y <= 0.5 * ly, (a.walls[0] - st["theta"]) / theta_tau,
                   (st["theta"] - a.walls[1]) / theta_tau)
    o = np.argsort(ywall)
    ywall, thp = ywall[o], thp[o]
    sub = ywall <= 2.0
    print(f"viscous sublayer theta+/(Pr y+) at y+ <= 2: "
          f"{(thp[sub]/(a.pr*ywall[sub])).min():.4f} .. "
          f"{(thp[sub]/(a.pr*ywall[sub])).max():.4f}")
    log = (ywall >= a.ylo) & (ywall <= a.yhi)
    slope, icept = np.polyfit(np.log(ywall[log]), thp[log], 1)
    ymid = np.exp(np.mean(np.log(ywall[log])))
    want = a.prt / 0.41 / (1.0 - ymid / (0.5 * ly * a.re))
    print(f"log fit y+ in [{a.ylo:g},{a.yhi:g}]: theta+ = {slope:.4f} ln y+ + {icept:.4f}")
    print(f"   constant-flux reference Pr_t/(kappa(1 - y/h)) = {want:.4f} "
          f"(dev {abs(slope-want)/want*100:.2f}%); the y/h -> 0 limit is "
          f"{a.prt/0.41:.4f}, Kader 2.12")
    if a.dump:
        np.savetxt(a.dump, np.column_stack([ywall, thp, kader(ywall, a.pr)]),
                   header="y+ theta+ kader")
    ok = err <= a.tolerance * abs(a.walls[0] - a.walls[1])
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def sorted_snapshots(patterns, skip):
    out = []
    for p in patterns:
        out.extend(glob.glob(p) if any(c in p for c in "*?[") else [p])
    if not out:
        raise SystemExit("no snapshots matched")
    out = sorted(set(out), key=lambda p: h5_time(p))
    return out[skip:] if skip else out


def h5_time(path):
    with h5py.File(path, "r") as f:
        return float(f.attrs["t_current"])


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("snapshots", nargs="+")
        p.add_argument("--scalar", default="theta")
        p.add_argument("--pr", type=float, default=0.71)
        p.add_argument("--prt", type=float, default=0.85)
        p.add_argument("--re", type=float, default=180.0)
        p.add_argument("--walls", type=float, nargs=2, default=[1.0, -1.0])
        p.add_argument("--skip", type=int, default=0,
                       help="drop the first N snapshots (thermal transient)")
        p.add_argument("--level", type=int, default=0)
        p.add_argument("--prt-model", default="constant", choices=["constant", "kays"])
        p.add_argument("--dump", default=None)

    p = sub.add_parser("channel")
    common(p)
    p.add_argument("--ylo", type=float, default=1.0)
    p.add_argument("--yhi", type=float, default=35.0)
    p.add_argument("--tolerance", type=float, default=0.25)
    p.add_argument("--flux-tolerance", type=float, default=0.10)
    p.set_defaults(func=cmd_channel)

    p = sub.add_parser("band")
    common(p)
    p.add_argument("--reference", nargs="+", required=True)
    p.add_argument("--tolerance", type=float, default=0.10)
    p.set_defaults(func=cmd_band)

    p = sub.add_parser("rans")
    common(p)
    p.add_argument("--ylo", type=float, default=30.0)
    p.add_argument("--yhi", type=float, default=60.0)
    p.add_argument("--tolerance", type=float, default=0.01,
                   help="max |theta - prediction| as a fraction of the wall difference")
    p.set_defaults(func=cmd_rans)

    a = ap.parse_args()
    PRT_RE[0] = a.re
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
