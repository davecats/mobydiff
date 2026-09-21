#!/usr/bin/env python3
"""Radial statistics of the conjugate pipe, out of block-table snapshots.

The solver's own statistics module ([scalar] stats_layout = plane) gives the
z-averaged (x,y) cross-section of every SCALAR quantity, with the wall fluxes
built from the transport kernel's own face diffusivity -- that is the exact
discrete flux and it is what scalar_stats.py reads. It does NOT give the
velocity statistics (the stats module is scalar-only) and it does not give the
axial or azimuthal turbulent heat flux. This script does all of it from
snapshots instead, so one accumulator serves both sides of the comparison.

WHAT IT DOES. Every cell is binned by the radius of its CENTRE. Because the
pipe is homogeneous in z and (at bccode = 0) in azimuth, a bin averages over
both directions at once: at the coarse grid one bin at the wall holds ~6e4
samples per snapshot. Velocities are interpolated to the cell centre from the
staggered faces (q(i) is the LOW face of cell i, so u_c = (q(i)+q(i+1))/2 --
scalar_stats.f90:388) and rotated into the pipe frame,
    u_r =  u cos(phi) + v sin(phi)      u_phi = -u sin(phi) + v cos(phi)
with u_z = w the axial component.

FLUID AND SOLID ARE ACCUMULATED SEPARATELY in every bin, because the immersed
wall is a STAIRCASE: a bin that straddles r = R holds cells of both materials,
and mixing them would smear the fluid profile into the solid one. The material
comes from the case file's own marker (|coef_p| > 1e20 is exactly "this centre
is in the solid" -- the test the solver uses), falling back to r > R.

THE BIN MEAN IS THE AVERAGING OPERATOR, so fluctuations are about it: the rms
reported is sqrt(<q^2> - <q>^2) over the bin. A bin has a finite width, so a
steep mean profile inside it inflates the rms slightly; the default width is
half a cell, which keeps that below the difference between the two grids.

Accumulators are saved to .npz and are ADDITIVE: run it on each batch of
snapshots and combine with --add, so a statistics leg does not have to keep
every snapshot on disk.

    ./pipe_stats.py pipeD_*.h5 --case pipe_coarse.h5 --out pipeD.npz
    ./pipe_stats.py --add pipeD_a.npz pipeD_b.npz --out pipeD.npz
    ./pipe_stats.py --report pipeD.npz
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np

SOLID_THRESHOLD = 1e20
VEL = ("un", "vn", "wn")


# --------------------------------------------------------------------------
# block-table io (the layout tools/compare_fields.py defines)
# --------------------------------------------------------------------------

def block_geometry(h5):
    """Block table, block size and the global lattice shape. A snapshot
    carries block_nb_x/y/z, a case file the single cubic block_nb."""
    blocks = h5["blocks"][...]
    if "block_nb_x" in h5.attrs:
        nb = (int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]),
              int(h5.attrs["block_nb_z"]))
    else:
        nb = (int(h5.attrs["block_nb"]),) * 3
    if int(blocks[:, 3].max()) != 0:
        raise SystemExit("pipe_stats: single-level files only")
    shape = (int(h5.attrs["nz"]), int(h5.attrs["ny"]), int(h5.attrs["nx"]))
    return blocks, nb, shape


def load_field(h5, name):
    """Reassemble one interior field onto the global (z, y, x) lattice."""
    blocks, nb, shape = block_geometry(h5)
    arr = np.zeros(shape, dtype=np.float64)
    data = h5[name]
    for bid, (ox, oy, oz, _lev) in enumerate(blocks):
        arr[oz:oz + nb[2], oy:oy + nb[1], ox:ox + nb[0]] = data[bid]
    return arr


def load_solid_mask(case):
    """The solver's own material marker, from the case file's INTERIOR cells.

    Case-file tiles are ghost-inclusive and (i, j, k) ordered -- the opposite
    of the field datasets, a trap check_oblique.py paid for first and
    check_annulus.py documents. The interior slice [1:-1] in each direction is
    what lines up with a snapshot.
    """
    with h5py.File(case, "r") as h5:
        blocks, nb, shape = block_geometry(h5)
        coef = h5["coef_p_blocks"]
        mask = np.zeros(shape, dtype=bool)
        for bid, (ox, oy, oz, _lev) in enumerate(blocks):
            tile = coef[bid][1:-1, 1:-1, 1:-1]          # (i, j, k)
            mask[oz:oz + nb[2], oy:oy + nb[1], ox:ox + nb[0]] = \
                (np.abs(tile) > SOLID_THRESHOLD).transpose(2, 1, 0)
    return mask


def cell_centres(h5):
    """Cell centres from the file's node lines. A snapshot names them x/y/z,
    a case file x_nodes/y_nodes/z_nodes -- the same lines either way."""
    node = {d: h5[d if d in h5 else f"{d}_nodes"][...] for d in "xyz"}
    return {d: 0.5 * (node[d][:-1] + node[d][1:]) for d in "xyz"}


def scalar_names(h5):
    """Field datasets that are not velocity/pressure/turbulence quantities."""
    known = {"un", "vn", "wn", "pn", "nut", "k", "omega", "gamma", "rethetat",
             "fd", "vfrac"}
    return [k for k, v in h5.items()
            if isinstance(v, h5py.Dataset) and v.ndim == 4 and k not in known
            and k != "blocks"]


def azimuthal_mean(field, r, dx, min_cells=64, w_min=None):
    """The azimuthal mean of a cross-section field, as a function of radius,
    evaluated at each cell's own radius.

    THE BIN WIDTH MUST ADAPT. A ring of radius r and width w holds
    2*pi*r*w/dx^2 cells, so a FIXED width narrower than a cell samples only an
    ALIASED subset of azimuths -- and because the cells at a given radius have
    four-fold symmetry, that aliasing lands on m = 4, 8, 12, exactly where a
    Cartesian grid artifact would. It is harmless where the mean profile is
    flat and ruinous where it is steep: near the axis the bins hold no cells
    at all, and in a near-insulating shell (kappa_s/kappa_f = 0.01, whose mean
    drops by thousands of units across it) it manufactures a striation that
    looks like a numerical instability and is not.

    So each ring is widened until it holds at least `min_cells` cells:
        w(r) = max(w_min, min_cells * dx^2 / (2 pi r)).
    """
    w_min = w_min or 0.25 * dx
    edges = [0.0]
    while edges[-1] < r.max():
        rc = edges[-1]
        w = max(w_min, min_cells * dx * dx / (2.0 * np.pi * max(rc, 0.5 * dx)))
        edges.append(rc + w)
    edges = np.asarray(edges)
    idx = np.clip(np.digitize(r.ravel(), edges) - 1, 0, len(edges) - 2)
    v = np.nan_to_num(field.ravel())
    good = np.isfinite(field.ravel())
    n = np.bincount(idx[good], minlength=len(edges) - 1)
    sm = np.bincount(idx[good], weights=v[good], minlength=len(edges) - 1)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    m = n > 0
    return np.interp(r, ctr[m], (sm[m] / n[m]))


# --------------------------------------------------------------------------
# accumulation
# --------------------------------------------------------------------------

# Per-bin sums, in the order they are stored. Velocity moments first, then
# five moments per scalar.
VEL_KEYS = ("uz", "ur", "up", "uzz", "urr", "upp", "urz")
SC_KEYS = ("s", "ss", "sur", "suz", "sup")


def accumulate(files, case, centre, radius, dr, nsector, verbose=True):
    cx, cy = centre
    solid_ref = load_solid_mask(case) if case else None

    acc = None
    for path in files:
        with h5py.File(path, "r") as h5:
            cent = cell_centres(h5)
            names = scalar_names(h5)
            shape = (len(cent["z"]), len(cent["y"]), len(cent["x"]))

            if acc is None:
                # Geometry is the same for every snapshot of a campaign, so
                # build the bin map once.
                zc, yc, xc = np.meshgrid(cent["z"], cent["y"], cent["x"],
                                         indexing="ij")
                del zc
                px, py = xc - cx, yc - cy
                del xc, yc
                r = np.hypot(px, py)
                cphi, sphi = px / np.maximum(r, 1e-300), py / np.maximum(r, 1e-300)
                del px, py
                nbin = int(np.ceil(r.max() / dr))
                ibin = np.minimum((r / dr).astype(np.int32), nbin - 1)
                solid = solid_ref if solid_ref is not None else (r > radius)
                if solid.shape != shape:
                    raise SystemExit("pipe_stats: case file and snapshot disagree"
                                     f" on the grid: {solid.shape} vs {shape}")
                # Azimuthal sectors, for the uniformity check (Section 6.5):
                # a second, coarser binning of the SAME cells.
                phi = np.arctan2(sphi, cphi)
                isec = np.minimum(((phi + np.pi) / (2 * np.pi) * nsector)
                                  .astype(np.int32), nsector - 1)
                del phi
                idx_fluid = (ibin * 2 + 0).ravel()
                idx_solid = (ibin * 2 + 1).ravel()
                idx = np.where(solid.ravel(), idx_solid, idx_fluid)
                del idx_fluid, idx_solid
                idx_sec = (ibin * nsector + isec).ravel()
                acc = {
                    "nbin": nbin, "dr": dr, "nsector": nsector,
                    "radius": radius, "centre": np.array([cx, cy]),
                    "r_edges": np.arange(nbin + 1) * dr,
                    "names": np.array(names, dtype=object),
                    "nfile": 0,
                    "count": np.zeros(2 * nbin),
                    "count_sec": np.zeros(nbin * nsector),
                }
                for key in VEL_KEYS:
                    acc[key] = np.zeros(2 * nbin)
                for nm in names:
                    for key in SC_KEYS:
                        acc[f"{nm}_{key}"] = np.zeros(2 * nbin)
                    # sector sums of the scalar and its square, fluid+solid
                    # together (a sector at fixed r is one material almost
                    # everywhere, and the staircase is exactly what this
                    # check is for).
                    acc[f"{nm}_sec"] = np.zeros(nbin * nsector)
                    acc[f"{nm}_secss"] = np.zeros(nbin * nsector)

            # --- velocities at the cell centre -------------------------
            u = load_field(h5, "un")
            uc = 0.5 * (u + np.roll(u, -1, axis=2))
            uc[:, :, -1] = u[:, :, -1]        # x is not periodic; jacket only
            del u
            v = load_field(h5, "vn")
            vc = 0.5 * (v + np.roll(v, -1, axis=1))
            vc[:, -1, :] = v[:, -1, :]
            del v
            w = load_field(h5, "wn")
            wc = 0.5 * (w + np.roll(w, -1, axis=0))   # z IS periodic
            del w

            ur = (uc * cphi + vc * sphi).ravel()
            up = (-uc * sphi + vc * cphi).ravel()
            uz = wc.ravel()
            del uc, vc, wc

            nb2 = 2 * acc["nbin"]
            for key, val in (("uz", uz), ("ur", ur), ("up", up),
                             ("uzz", uz * uz), ("urr", ur * ur),
                             ("upp", up * up), ("urz", ur * uz)):
                acc[key] += np.bincount(idx, weights=val, minlength=nb2)

            # --- scalars ------------------------------------------------
            for nm in acc["names"]:
                s = load_field(h5, str(nm)).ravel()
                acc[f"{nm}_s"] += np.bincount(idx, weights=s, minlength=nb2)
                acc[f"{nm}_ss"] += np.bincount(idx, weights=s * s, minlength=nb2)
                acc[f"{nm}_sur"] += np.bincount(idx, weights=s * ur, minlength=nb2)
                acc[f"{nm}_suz"] += np.bincount(idx, weights=s * uz, minlength=nb2)
                acc[f"{nm}_sup"] += np.bincount(idx, weights=s * up, minlength=nb2)
                nsec_tot = acc["nbin"] * acc["nsector"]
                acc[f"{nm}_sec"] += np.bincount(idx_sec, weights=s,
                                                minlength=nsec_tot)
                acc[f"{nm}_secss"] += np.bincount(idx_sec, weights=s * s,
                                                  minlength=nsec_tot)
                del s
            del ur, up, uz

            acc["count"] += np.bincount(idx, minlength=2 * acc["nbin"])
            acc["count_sec"] += np.bincount(
                idx_sec, minlength=acc["nbin"] * acc["nsector"])
            acc["nfile"] += 1
            if verbose:
                print(f"   {path}: binned", flush=True)
    return acc


def bin_plane_stats(path, case, centre, radius, dr, names):
    """Radially bin the solver's own PLANE statistics ([scalar] stats_layout =
    plane), which are z-averaged (x,y) rows sampled every
    `stats_sample_interval` steps -- far better averaged than snapshots, and
    they cover the SOLID as well, which is where the conjugate signature
    lives.

    Row ordering is the solver's: `row = (gx - 1)*ny + gy`, i.e. y fastest
    (scalar_stats.f90). Columns are 7 per scalar, in the order
    <s> <s^2> <u_c s> <v s>|lo <J>|lo <v s>|hi <J>|hi.

    THE RADIAL CONVECTIVE FLUX built here is
        <u_r s> ~ <u_x s> cos(phi) + <v s> sin(phi),
    with <v s> the mean of the cell's two y faces. It is the TURBULENT flux
    only to the extent that <u_x> and <v> vanish, which they do in a pipe up
    to the statistical error -- the stats file stores no mean velocity to
    subtract. Use the snapshot path when that matters, and for the AXIAL flux
    <u_z s>, which the plane layout does not carry at all.
    """
    NSTAT = 7
    with h5py.File(path, "r") as h5:
        prof = h5["profile"][...]
        x, y = h5["xcoord"][...], h5["ycoord"][...]
        nx, ny = int(h5.attrs["nx"]), int(h5.attrs["ny"])
        nscalar = int(h5.attrs["nstat"]) // NSTAT
    if len(names) != nscalar:
        raise SystemExit(f"pipe_stats: file has {nscalar} scalars, "
                         f"{len(names)} names given (--names)")

    cx, cy = centre
    X, Y = np.meshgrid(x, y, indexing="ij")              # (nx, ny), row order
    r = np.hypot(X - cx, Y - cy)
    c, sn = (X - cx) / np.maximum(r, 1e-300), (Y - cy) / np.maximum(r, 1e-300)
    if case:
        solid = load_solid_mask(case)[0].T               # (z,y,x) -> (x,y)
    else:
        solid = (r > radius)

    nbin = int(np.ceil(r.max() / dr))
    ibin = np.minimum((r / dr).astype(np.int32), nbin - 1)
    idx = (ibin * 2 + solid.astype(np.int32)).ravel()
    nb2 = 2 * nbin

    acc = {"nbin": nbin, "dr": dr, "nsector": 1, "radius": radius,
           "centre": np.array([cx, cy]), "r_edges": np.arange(nbin + 1) * dr,
           "names": np.array(names, dtype=object), "nfile": 1,
           "count": np.bincount(idx, minlength=nb2).astype(float),
           "count_sec": np.zeros(nbin)}
    for key in VEL_KEYS:
        acc[key] = np.zeros(nb2)
    for k, nm in enumerate(names):
        base = NSTAT * k

        def col(j, base=base):
            return prof[:, base + j].reshape(nx, ny).ravel()

        # THE PER-CELL TEMPORAL VARIANCE, kept separately. Forming the rms
        # from bin-averaged moments instead adds the azimuthal variance of the
        # steady (time-mean) field to the temporal fluctuation -- and in the
        # solid that steady part is the STAIRCASE IMPRINT, which does not
        # decay with depth the way a fluctuation does. Measured inflation of
        # the solid-side rms: 22-36 % (c0), 11-17 % (c3), 72-113 % (isoflux).
        # The temporal definition is also the one the reference data uses.
        var_cell = np.maximum(col(1) - col(0) ** 2, 0.0)
        vs = 0.5 * (col(3) + col(5))
        acc[f"{nm}_s"] = np.bincount(idx, weights=col(0), minlength=nb2)
        acc[f"{nm}_ss"] = np.bincount(idx, weights=col(1), minlength=nb2)
        acc[f"{nm}_var"] = np.bincount(idx, weights=var_cell, minlength=nb2)
        acc[f"{nm}_sur"] = np.bincount(
            idx, weights=col(2) * c.ravel() + vs * sn.ravel(), minlength=nb2)
        acc[f"{nm}_suz"] = np.zeros(nb2)
        acc[f"{nm}_sup"] = np.zeros(nb2)
        acc[f"{nm}_sec"] = np.zeros(nbin)
        acc[f"{nm}_secss"] = np.zeros(nbin)
    return acc


def combine(paths):
    out = None
    for p in paths:
        d = dict(np.load(p, allow_pickle=True))
        if out is None:
            out = d
            continue
        for k, v in d.items():
            if k in ("nbin", "dr", "nsector", "radius", "centre", "r_edges", "names"):
                continue
            out[k] = out[k] + v
    return out


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def profiles(acc):
    """Bin means and rms, fluid and solid rows kept apart."""
    nbin = int(acc["nbin"])
    r = 0.5 * (acc["r_edges"][:-1] + acc["r_edges"][1:])
    cnt = np.asarray(acc["count"]).reshape(nbin, 2)
    out = {"r": r, "count": cnt}
    safe = np.where(cnt > 0, cnt, 1.0)

    for key in VEL_KEYS:
        out[key] = np.asarray(acc[key]).reshape(nbin, 2) / safe
    out["ur_rms"] = np.sqrt(np.maximum(out["urr"] - out["ur"] ** 2, 0.0))
    out["up_rms"] = np.sqrt(np.maximum(out["upp"] - out["up"] ** 2, 0.0))
    out["uz_rms"] = np.sqrt(np.maximum(out["uzz"] - out["uz"] ** 2, 0.0))
    out["urz_cov"] = out["urz"] - out["ur"] * out["uz"]

    names = [str(n) for n in acc["names"]]
    out["names"] = names
    for nm in names:
        s = np.asarray(acc[f"{nm}_s"]).reshape(nbin, 2) / safe
        ss = np.asarray(acc[f"{nm}_ss"]).reshape(nbin, 2) / safe
        out[f"{nm}_mean"] = s
        if f"{nm}_var" in acc:
            # the temporal rms, averaged over the bin (see bin_plane_stats)
            out[f"{nm}_rms"] = np.sqrt(
                np.asarray(acc[f"{nm}_var"]).reshape(nbin, 2) / safe)
            out[f"{nm}_rms_binmoment"] = np.sqrt(np.maximum(ss - s ** 2, 0.0))
        else:
            out[f"{nm}_rms"] = np.sqrt(np.maximum(ss - s ** 2, 0.0))
        for comp in ("ur", "uz", "up"):
            f = np.asarray(acc[f"{nm}_s{comp}"]).reshape(nbin, 2) / safe
            out[f"{nm}_flux_{comp}"] = f - s * out[comp]
    return out


def report(acc, radius):
    p = profiles(acc)
    nbin = int(acc["nbin"])
    r, cnt = p["r"], p["count"]
    print(f"# snapshots = {int(acc['nfile'])}, bins = {nbin}, "
          f"dr = {float(acc['dr']):.5g}")
    print("#     r       cells      <uz>      uz'      ur'      up'"
          "      <ur uz>")
    for i in range(nbin):
        if cnt[i, 0] < 1:
            continue
        if r[i] > radius:
            break
        print(f"{r[i]:8.5f} {int(cnt[i,0]):9d} "
              f"{p['uz'][i,0]:9.5f} {p['uz_rms'][i,0]:8.5f} "
              f"{p['ur_rms'][i,0]:8.5f} {p['up_rms'][i,0]:8.5f} "
              f"{p['urz_cov'][i,0]:10.3e}")
    for nm in p["names"]:
        print(f"\n# scalar {nm}:   r    material   <theta>      theta'   "
              "  <ur theta>'   <uz theta>'")
        for i in range(nbin):
            for m, lab in ((0, "fluid"), (1, "solid")):
                if cnt[i, m] < 1:
                    continue
                print(f"{r[i]:8.5f}  {lab:6} {p[f'{nm}_mean'][i,m]:12.5e} "
                      f"{p[f'{nm}_rms'][i,m]:12.5e} "
                      f"{p[f'{nm}_flux_ur'][i,m]:12.4e} "
                      f"{p[f'{nm}_flux_uz'][i,m]:12.4e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", help="snapshots to bin")
    ap.add_argument("--case", help="case file, for the solver's own material marker")
    ap.add_argument("--centre", type=float, nargs=2, default=(0.65, 0.65))
    ap.add_argument("--radius", type=float, default=0.5)
    ap.add_argument("--dr", type=float, default=None,
                    help="radial bin width (default: half a cell)")
    ap.add_argument("--sectors", type=int, default=16)
    ap.add_argument("--out")
    ap.add_argument("--add", nargs="*", help="combine existing .npz accumulators")
    ap.add_argument("--report", help="print the profiles of an .npz")
    ap.add_argument("--stats", help="bin the solver's plane statistics file"
                                    " instead of snapshots")
    ap.add_argument("--names", default="c0,c1,c2,c3,c4,mbc,isof",
                    help="scalar names, in configuration order (--stats)")
    a = ap.parse_args()

    if a.report:
        acc = dict(np.load(a.report, allow_pickle=True))
        report(acc, a.radius)
        return 0

    if a.stats:
        if a.dr is None:
            with h5py.File(a.stats, "r") as h5:
                xc = h5["xcoord"][...]
            a.dr = 0.5 * float(xc[1] - xc[0])
        acc = bin_plane_stats(a.stats, a.case, a.centre, a.radius, a.dr,
                              a.names.split(","))
    elif a.add:
        acc = combine(a.add)
    else:
        if not a.files:
            ap.error("nothing to do: give snapshots, --add or --report")
        if a.dr is None:
            with h5py.File(a.files[0], "r") as h5:
                cent = cell_centres(h5)
                a.dr = 0.5 * float(cent["x"][1] - cent["x"][0])
        acc = accumulate(a.files, a.case, a.centre, a.radius, a.dr, a.sectors)

    if a.out:
        np.savez(a.out, **acc)
        print(f"   wrote {a.out}")
    else:
        report(acc, a.radius)
    return 0


if __name__ == "__main__":
    sys.exit(main())
