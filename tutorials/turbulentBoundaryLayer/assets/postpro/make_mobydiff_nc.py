#!/usr/bin/env python3
"""Export a mobydiff boundary-layer statistics file to the res-study NetCDF format
used by the CaNS / AMPHIBIOUS comparison data
(~/.../tbl-dns/res_study/data_processed/<code>/xyz_<nx>_<ny>_<nz>/data.nc).

    make_mobydiff_nc.py <bl_stats.h5> <field_with_nodes.h5> <out.nc> [--ini production_stats.ini]

The stats file (fdm_h5_write_bl_stats) stores span+time-averaged profiles on the
(x, y) plane, flattened y-fastest, with the 10 accumulators
U V W UU VV WW UV UW VW P (all velocities interpolated to the cell centre, so the
means are collocated -- like the reference codes). The `profile` dataset is already
sum/count (the mean of each accumulator); Reynolds stresses are the mean-subtracted
covariances, e.g. uu_stress = <uu> - <u><u>.

The grid NODE lines (x,y,z of length n+1, starting at 0) come from any field
snapshot of the same run; cell centres are 0.5*(node[i]+node[i+1]) (== the stats
xcoord/ycoord) and the staggered "upper face" coordinates xu/yv/zw are node[1:].
Coordinates are already in delta*_in units (the case is nondimensional, delta*_0=1).

Output layout matches the CaNS files: 2D variables dimensioned (x, y), plus the
1D coordinate variables and the snapshot bookkeeping, and the same global
attributes (n, l, time, step, data_dim, git_*, nu, re_dstar, raw_input_nml).
"""
import argparse
import subprocess

import h5py
import numpy as np
import netCDF4 as nc

# accumulator order in the stats `profile` (0-indexed), matching
# boundarylayer_stats.f90 STAT_* (1-indexed there).
U, V, W, UU, VV, WW, UV, UW, VW, P = range(10)


def git_info(repo):
    def g(*a):
        try:
            return subprocess.check_output(["git", "-C", repo, *a], text=True).strip()
        except Exception:
            return ""
    return g("rev-parse", "HEAD"), g("log", "-1", "--format=%cI")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stats")
    ap.add_argument("field", help="a field snapshot of the same run (for the node lines)")
    ap.add_argument("out")
    ap.add_argument("--ini", default=None, help="run ini, stored verbatim as raw_input_nml")
    ap.add_argument("--repo", default=".", help="repo path for git metadata")
    a = ap.parse_args()

    with h5py.File(a.stats, "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        re = float(f.attrs["re"]); step = int(f.attrs["step"]); t = float(f.attrs["t_current"])
        prof = f["profile"][...].reshape(nx, ny, -1)
        xc = f["xcoord"][...]; yc = f["ycoord"][...]
    with h5py.File(a.field, "r") as f:
        xn, yn, zn = f["x"][...], f["y"][...], f["z"][...]      # node lines, length n+1
        lx, ly, lz = (float(f.attrs[k]) for k in ("lx", "ly", "lz"))
        nz = int(f.attrs["nz"])
    assert len(xn) == nx + 1 and len(yn) == ny + 1 and len(zn) == nz + 1
    zc = 0.5 * (zn[:-1] + zn[1:])                              # span cell centres

    nu = 1.0 / re
    Um, Vm, Wm, Pm = prof[:, :, U], prof[:, :, V], prof[:, :, W], prof[:, :, P]
    # Reynolds stresses = mean-subtracted covariances.
    stress = {
        "uu_stress": prof[:, :, UU] - Um * Um,
        "vv_stress": prof[:, :, VV] - Vm * Vm,
        "ww_stress": prof[:, :, WW] - Wm * Wm,
        "uv_stress": prof[:, :, UV] - Um * Vm,
        "uw_stress": prof[:, :, UW] - Um * Wm,
        "vw_stress": prof[:, :, VW] - Vm * Wm,
    }
    means = {"u_mean": Um, "v_mean": Vm, "w_mean": Wm, "p_mean": Pm}

    commit, date = git_info(a.repo)
    raw_ini = ""
    if a.ini:
        with open(a.ini) as fh:
            raw_ini = fh.read()

    ds = nc.Dataset(a.out, "w", format="NETCDF4")
    ds.n = np.array([nx, ny, nz], dtype=np.int32)
    ds.l = np.array([lx, ly, lz], dtype=np.float64)
    ds.time = t
    ds.step = step
    ds.data_dim = "2d"
    ds.git_repo = "mobydiff"
    ds.git_commit = commit
    ds.Host = "HoreKa"
    ds.Date = date
    ds.nu = nu
    ds.re_dstar = re
    ds.statistics_grid = "cell centres (velocities interpolated); means collocated"
    ds.statistics_xy_interpolation = "none"
    ds.raw_input_nml = raw_ini

    ds.createDimension("x", nx); ds.createDimension("y", ny); ds.createDimension("z", nz)
    ds.createDimension("snapshot_indices", 1)

    def coord(name, data, long_name):
        v = ds.createVariable(name, "f8", (name[0],), fill_value=np.nan)
        v[:] = data
        if long_name:
            v.long_name = long_name

    coord("x", xc, r"$x / \delta^*_{in}$")
    coord("y", yc, r"$y / \delta^*_{in}$")
    coord("z", zc, r"$z / \delta^*_{in}$")
    # staggered "upper face" coordinates (node[1:]); named like the reference files.
    for name, data in (("xu", xn[1:]), ("yv", yn[1:]), ("zw", zn[1:])):
        v = ds.createVariable(name, "f8", (name[0],), fill_value=np.nan); v[:] = data

    for name, fld in {**means, **stress}.items():
        v = ds.createVariable(name, "f8", ("x", "y"), fill_value=np.nan)
        v[:] = fld
        v.coordinates = "x y"

    si = ds.createVariable("snapshot_indices", "i8", ("snapshot_indices",)); si[:] = [0]
    ss = ds.createVariable("snapshot_steps", "i8", ("snapshot_indices",)); ss[:] = [step]
    st = ds.createVariable("snapshot_times", "f8", ("snapshot_indices",), fill_value=np.nan)
    st[:] = [t]
    ds.close()
    print(f"wrote {a.out}: n={[nx, ny, nz]} l={[lx, ly, lz]} step={step} t={t:.0f}")


if __name__ == "__main__":
    main()
