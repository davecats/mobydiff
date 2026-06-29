#!/usr/bin/env python3
"""Generate initial conditions for the 2:1 interface channel validation.

Interpolates the fields of an existing channel restart (e.g.
tutorials/channel_kmm180/channel_kmm180_restart.h5, 288x136x288) onto the
validation grids and writes ready-to-restart files:

- reference: uniform 256x128x256 run, legacy global-3D layout. Its grid
  lines are the midpoint subdivision of the 128x64x128 base lines
  ([grid.*] subdivided = true in the matching input), so the refined
  case's fine level shares them bitwise.
- refined: 128x64x128 base with both wall bands refined to level 1,
  block-table layout. The leaf table reproduces the solver's enumeration
  for [blocks] nb=8 with two refine boxes covering `band` base cells at
  each wall.

The interpolated fields are not divergence-free on the new grids; the
first projection absorbs that, which is why the validation discards a
transient before sampling statistics.
"""

from __future__ import annotations

import argparse

import h5py
import numpy as np

NB = 8


# ---------------------------------------------------------------------------
# Grid lines (ports of src/modules/init.f90; keep in sync)

def natural_wall_coordinate(j, blend_index, dy_wall_plus, alpha=1.25, c_eta=0.8):
    j = np.asarray(j, dtype=np.float64)
    jb = blend_index if blend_index > 0.0 else 40.0
    dy_wall = dy_wall_plus if dy_wall_plus > 0.0 else 0.05
    blend = (j/jb)**2
    outer = (0.75*alpha*c_eta*j)**(4.0/3.0)
    return np.where(j <= 0.0, 0.0, (dy_wall*j + outer*blend)/(1.0 + blend))


def natural_line(n, length, blend_index, dy_wall_plus):
    i = np.arange(n + 1, dtype=np.float64)
    s = i/n
    j = np.minimum(s, 1.0 - s)*n
    yp = natural_wall_coordinate(j, blend_index, dy_wall_plus)
    ypm = natural_wall_coordinate(0.5*n, blend_index, dy_wall_plus)
    half = 0.5*length
    node = np.where(s <= 0.5, half*yp/ypm, length - half*yp/ypm)
    node[0] = 0.0
    node[-1] = length
    return node


def uniform_line(n, length):
    i = np.arange(n + 1, dtype=np.float64)
    node = length*(i/n)
    node[0] = 0.0
    node[-1] = length
    return node


def subdivide(line):
    """Midpoint subdivision, bitwise as blocks.f90 subdivide_node_line."""
    fine = np.empty(2*(line.size - 1) + 1, dtype=np.float64)
    fine[0::2] = line
    fine[1::2] = 0.5*(line[:-1] + line[1:])
    return fine


# ---------------------------------------------------------------------------
# Tensor-product linear interpolation of the staggered fields

def axis_interp(arr, axis, src, dst, periodic, length):
    """Linear interpolation along one axis of a 3D array."""
    src = np.asarray(src)
    if periodic:
        srcx = np.concatenate((src - length, src, src + length))
        arr = np.concatenate((arr,)*3, axis=axis)
    else:
        srcx = src
    idx = np.clip(np.searchsorted(srcx, dst) - 1, 0, srcx.size - 2)
    w = (dst - srcx[idx])/(srcx[idx+1] - srcx[idx])
    w = np.clip(w, 0.0, 1.0)
    a0 = np.take(arr, idx, axis=axis)
    a1 = np.take(arr, idx + 1, axis=axis)
    shape = [1, 1, 1]
    shape[axis] = w.size
    w = w.reshape(shape)
    return a0*(1.0 - w) + a1*w


def interp_field(field, src_pos, dst_pos, periodic, lengths):
    out = field
    for axis in range(3):  # array axes are (z, y, x) -> dims (3, 2, 1)
        d = 2 - axis
        out = axis_interp(out, axis, src_pos[d], dst_pos[d], periodic[d], lengths[d])
    return out


def staggered_positions(nodes, var):
    """Positions of the stored values per dimension. nodes = (x, y, z)."""
    cent = [0.5*(n[:-1] + n[1:]) for n in nodes]
    face = [n[:-1] for n in nodes]
    pos = list(cent)
    if var < 3:
        pos[var] = face[var]
    return pos  # (x positions, y positions, z positions)


# ---------------------------------------------------------------------------
# Leaf table for the band-refined case (mirror of blocks.f90 for this
# specific layout: whole base block rows refined at both walls, which is
# already 2:1 smooth)

def morton_key(cx, cy, cz):
    key = 0
    for bit in range(21):
        key |= ((cx >> bit) & 1) << (3*bit)
        key |= ((cy >> bit) & 1) << (3*bit + 1)
        key |= ((cz >> bit) & 1) << (3*bit + 2)
    return key


def band_leaf_table(gnbt, band_rows):
    """(level, origin-in-level-cells) leaves: y block rows < band_rows or
    >= gnbt[1]-band_rows refined to level 1; Morton ids on the finest
    lattice."""
    leaves = []
    for cz in range(gnbt[2]):
        for cy in range(gnbt[1]):
            for cx in range(gnbt[0]):
                refined = cy < band_rows or cy >= gnbt[1] - band_rows
                if refined:
                    for sz in (0, 1):
                        for sy in (0, 1):
                            for sx in (0, 1):
                                leaves.append((1, 2*cx + sx, 2*cy + sy, 2*cz + sz))
                else:
                    leaves.append((0, cx, cy, cz))
    lmax = 1
    def key(leaf):
        lev, cx, cy, cz = leaf
        f = 2**(lmax - lev)
        return morton_key(cx*f, cy*f, cz*f)
    leaves.sort(key=key)
    return leaves


# ---------------------------------------------------------------------------

def load_source(path):
    with h5py.File(path, "r") as h5:
        src = {
            "fields": {k: h5[k][...] for k in ("un", "vn", "wn", "pn")},
            "nodes": (h5["x"][...], h5["y"][...], h5["z"][...]),
            "attrs": dict(h5.attrs),
        }
    return src


def common_attrs(src_attrs, nx, ny, nz, dyw_plus):
    a = dict(src_attrs)
    a.update({
        "nx": np.int32(nx), "ny": np.int32(ny), "nz": np.int32(nz),
        "grid_distribution": np.array([1, 4, 1], dtype=np.int32),
        "grid_stretch": np.array([0.0, 16.0, 0.0]),
        "grid_natural_dyw_plus": np.array([0.05, dyw_plus, 0.05]),
        "step": np.int32(0), "t_current": 0.0, "nsteps": np.int32(0),
        "dt": 1.0e-5,
        "nranks": np.int32(1),
    })
    return a


def write_attrs(h5, attrs):
    for k, v in attrs.items():
        h5.attrs[k] = v


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", default="tutorials/channel_kmm180/channel_kmm180_restart.h5")
    parser.add_argument("--mode", choices=("reference", "refined", "base"), required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--band-cells", type=int, default=24,
                        help="refined band height in base y-cells per wall (refined mode)")
    parser.add_argument("--dyw-plus", type=float, default=0.5)
    args = parser.parse_args()

    src = load_source(args.source)
    lx = float(src["attrs"]["lx"]); ly = float(src["attrs"]["ly"]); lz = float(src["attrs"]["lz"])
    lengths = (lx, ly, lz)
    periodic = (True, False, True)

    base = (128, 64, 128)
    base_nodes = (uniform_line(base[0], lx),
                  natural_line(base[1], ly, 16.0, args.dyw_plus),
                  uniform_line(base[2], lz))
    fine_nodes = tuple(subdivide(n) for n in base_nodes)

    # Interpolate each staggered field onto base and fine lattices.
    src_nodes = src["nodes"]
    names = ["un", "vn", "wn", "pn"]
    fine = {}
    coarse = {}
    for var, name in enumerate(names):
        spos = staggered_positions(src_nodes, var)
        f = src["fields"][name]
        fine[name] = interp_field(f, spos, staggered_positions(fine_nodes, var),
                                  periodic, lengths)
        if args.mode in ("refined", "base"):
            coarse[name] = interp_field(f, spos, staggered_positions(base_nodes, var),
                                        periodic, lengths)

    if args.mode == "base":
        # Uniform 128x64x128 at the base resolution (no refinement), legacy
        # global-3D layout. Its grid lines ARE the refined case's level-0
        # (coarse core) lines bitwise, so refined-core vs base isolates the
        # interface effect from the coarse-resolution deficit.
        attrs = common_attrs(src["attrs"], base[0], base[1], base[2], args.dyw_plus)
        with h5py.File(args.out, "w") as h5:
            write_attrs(h5, attrs)
            for name in names:
                h5.create_dataset(name, data=coarse[name])
            h5.create_dataset("x", data=base_nodes[0])
            h5.create_dataset("y", data=base_nodes[1])
            h5.create_dataset("z", data=base_nodes[2])
            h5.create_dataset("rank_local_range", data=np.array(
                [[1, base[0], 1, base[1], 1, base[2]]], dtype=np.int32))
        iface = None
    elif args.mode == "reference":
        attrs = common_attrs(src["attrs"], 2*base[0], 2*base[1], 2*base[2], args.dyw_plus)
        with h5py.File(args.out, "w") as h5:
            write_attrs(h5, attrs)
            for name in names:
                h5.create_dataset(name, data=fine[name])
            h5.create_dataset("x", data=fine_nodes[0])
            h5.create_dataset("y", data=fine_nodes[1])
            h5.create_dataset("z", data=fine_nodes[2])
            h5.create_dataset("rank_local_range", data=np.array(
                [[1, 2*base[0], 1, 2*base[1], 1, 2*base[2]]], dtype=np.int32))
        iface = None
    else:
        gnbt = tuple(n//NB for n in base)
        band_rows = args.band_cells//NB
        assert band_rows*NB == args.band_cells, "band must be whole block rows"
        leaves = band_leaf_table(gnbt, band_rows)
        nleaf = len(leaves)
        blocks = np.zeros((nleaf, 4), dtype=np.int32)
        rows = {name: np.zeros((nleaf, NB, NB, NB)) for name in names}
        grids = {0: coarse, 1: fine}
        for bid, (lev, cx, cy, cz) in enumerate(leaves):
            blocks[bid] = (cx*NB, cy*NB, cz*NB, lev)
            g = grids[lev]
            ox, oy, oz = cx*NB, cy*NB, cz*NB
            for name in names:
                rows[name][bid] = g[name][oz:oz+NB, oy:oy+NB, ox:ox+NB]
        attrs = common_attrs(src["attrs"], base[0], base[1], base[2], args.dyw_plus)
        attrs.update({"block_nb_x": np.int32(NB), "block_nb_y": np.int32(NB),
                      "block_nb_z": np.int32(NB), "n_blocks": np.int32(nleaf)})
        with h5py.File(args.out, "w") as h5:
            write_attrs(h5, attrs)
            h5.create_dataset("blocks", data=blocks)
            for name in names:
                h5.create_dataset(name, data=rows[name])
            h5.create_dataset("x", data=base_nodes[0])
            h5.create_dataset("y", data=base_nodes[1])
            h5.create_dataset("z", data=base_nodes[2])
        iface = 0.5*(base_nodes[1][args.band_cells] + base_nodes[1][args.band_cells + 1])
        nref = int((blocks[:, 3] == 1).sum())
        print(f"refined leaf table: {nleaf} leaves, {nref} fine")

    y = base_nodes[1]
    print(f"wrote {args.out}")
    print(f"base y-line: dy_wall+ = {y[1]*180:.3f} (fine {y[1]/2*180:.3f}), "
          f"interface candidates: node16 y+={y[16]*180:.1f}, node24 y+={y[24]*180:.1f}")
    if iface is not None:
        print(f"refine box edge for band={args.band_cells}: y1 = {iface:.16e} "
              f"(y+ = {y[args.band_cells]*180:.1f})")


if __name__ == "__main__":
    main()
