#!/usr/bin/env python3
"""External geometry-to-IBM coefficient utility for mobyDiff.

RETIRED FOR PRODUCTION (prepare/solve split P1b): the geometry subcommands
(stl-ibm-coeff, block-active, block-table) are superseded by the Fortran
moby_prepare executable ([ibm] stl_file; see docs/prepare_solve_strategy.md
and tools/README_mobygeom.md). Kept working as the independent
cross-implementation reference for the validation/prepare/ gates; the STL
generators and checkers remain in normal use.

The Python front-end owns STL mesh classification, HDF5 output, and analytic tests.
Grid coordinates are imported from a mobygrid HDF5 file written by the Fortran solver setup.
"""
from __future__ import annotations

import argparse
import csv
import concurrent.futures
import json
import math
import multiprocessing
import os
from pathlib import Path
from xml.sax.saxutils import escape

import h5py
import numpy as np

VAR_U = 1
VAR_V = 2
VAR_W = 3
SOLID = 1.0e30
AXIS_NAMES = {1: "x", 2: "y", 3: "z"}
GRID_DISTRIBUTION_IDS = {1: "uniform", 2: "cosine", 3: "tanh", 4: "natural"}

ROOT = Path(__file__).resolve().parents[1]


def coefficient_xdmf_path(h5_path: Path) -> Path:
    """Return the sidecar XDMF path associated with one coefficient HDF5 file."""
    return Path(h5_path).with_suffix(".xdmf")


def xdmf_hdf_reference(xdmf_path: Path, h5_path: Path, dataset: str) -> str:
    """Build a relative HDF5 reference suitable for an XDMF DataItem."""
    relative = Path(os.path.relpath(h5_path, xdmf_path.parent)).as_posix()
    return f"{escape(relative)}:{escape(dataset)}"


def write_coefficient_xdmf(xdmf_path: Path, h5_path: Path, args: argparse.Namespace) -> None:
    """Write a ParaView-readable XDMF wrapper for staggered IBM coefficients."""
    finalize_grid_args(args)
    xdmf_path = Path(xdmf_path).resolve()
    h5_path = Path(h5_path).resolve()
    dims = (int(args.nz) + 2, int(args.ny) + 2, int(args.nx) + 2)
    grids = []
    for var_name in ("u", "v", "w"):
        grid_name = f"coef_{var_name}"
        group = f"/xdmf/{var_name}"
        grids.append(f"""      <Grid Name=\"{grid_name}\" GridType=\"Uniform\">
        <Topology TopologyType=\"3DRectMesh\" Dimensions=\"{dims[0]} {dims[1]} {dims[2]}\"/>
        <Geometry GeometryType=\"VXVYVZ\">
          <DataItem Name=\"X\" Dimensions=\"{dims[2]}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{xdmf_hdf_reference(xdmf_path, h5_path, group + "/x")}</DataItem>
          <DataItem Name=\"Y\" Dimensions=\"{dims[1]}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{xdmf_hdf_reference(xdmf_path, h5_path, group + "/y")}</DataItem>
          <DataItem Name=\"Z\" Dimensions=\"{dims[0]}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{xdmf_hdf_reference(xdmf_path, h5_path, group + "/z")}</DataItem>
        </Geometry>
        <Attribute Name=\"{grid_name}\" AttributeType=\"Scalar\" Center=\"Node\">
          <DataItem Dimensions=\"{dims[0]} {dims[1]} {dims[2]}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{xdmf_hdf_reference(xdmf_path, h5_path, group + "/coef")}</DataItem>
        </Attribute>
      </Grid>
""")
    xdmf = f"""<?xml version=\"1.0\" ?>
<Xdmf Version=\"2.0\">
  <Domain>
    <Grid Name=\"mobyDiff IBM coefficients\" GridType=\"Collection\" CollectionType=\"Spatial\">
{''.join(grids)}    </Grid>
  </Domain>
</Xdmf>
"""
    xdmf_path.write_text(xdmf)


def write_grid_metadata(h5: h5py.File, args: argparse.Namespace) -> None:
    """Record grid generation options used by the coefficient preprocessor."""
    finalize_grid_args(args)
    h5.attrs["grid_nonuniform"] = int(grid_is_nonuniform(args))
    h5.attrs["grid_file"] = "" if getattr(args, "grid_file", None) is None else str(Path(args.grid_file).resolve())
    for direction in (1, 2, 3):
        axis = axis_name(direction)
        h5.attrs[f"grid_{axis}_distribution"] = getattr(args, f"{axis}_distribution", "uniform")
        h5.attrs[f"grid_{axis}_stretch"] = float(getattr(args, f"{axis}_stretch", 0.0))
        h5.attrs[f"grid_{axis}_natural_one_sided"] = int(bool(getattr(args, f"{axis}_natural_one_sided", False)))
        h5.attrs[f"grid_{axis}_nodes_file"] = "" if getattr(args, f"{axis}_nodes", None) is None else str(Path(getattr(args, f"{axis}_nodes")).resolve())
        h5.attrs[f"periodic_{axis}"] = int(grid_periodic(direction, args))


def write_hdf5(path: Path, coef: np.ndarray, args: argparse.Namespace, geometry: Path,
               extra_attrs: dict[str, object] | None = None) -> Path:
    """Write coefficients in mobyDiff format and add a ParaView XDMF sidecar."""
    finalize_grid_args(args)
    path = Path(path).resolve()
    with h5py.File(path, "w") as h5:
        h5.create_dataset("coef", data=coef)
        h5.create_dataset("coef_u", data=coef[..., 0])
        h5.create_dataset("coef_v", data=coef[..., 1])
        h5.create_dataset("coef_w", data=coef[..., 2])
        xdmf_group = h5.create_group("xdmf")
        indices = {
            1: np.arange(0, int(args.nx) + 2, dtype=np.float64),
            2: np.arange(0, int(args.ny) + 2, dtype=np.float64),
            3: np.arange(0, int(args.nz) + 2, dtype=np.float64),
        }
        extents = {1: (int(args.nx), float(args.lx)),
                   2: (int(args.ny), float(args.ly)),
                   3: (int(args.nz), float(args.lz))}
        for var_name, var in (("u", VAR_U), ("v", VAR_V), ("w", VAR_W)):
            group = xdmf_group.create_group(var_name)
            for direction, coord_name in ((1, "x"), (2, "y"), (3, "z")):
                group.create_dataset(coord_name, data=stl_axis_coords(indices[direction], direction, var, args))
            group.create_dataset("coef", data=np.transpose(coef[..., var - 1], (2, 1, 0)))
        h5.attrs["nx"] = int(args.nx)
        h5.attrs["ny"] = int(args.ny)
        h5.attrs["nz"] = int(args.nz)
        h5.attrs["lx"] = float(args.lx)
        h5.attrs["ly"] = float(args.ly)
        h5.attrs["lz"] = float(args.lz)
        h5.attrs["re"] = float(args.re)
        h5.attrs["source_geometry"] = str(geometry)
        h5.attrs["convention"] = "mobyDiff staggered u/v/w coefficients, shape (nx+2,ny+2,nz+2,3)"
        h5.attrs["xdmf_file"] = coefficient_xdmf_path(path).name
        write_grid_metadata(h5, args)
        if extra_attrs:
            for key, value in extra_attrs.items():
                h5.attrs[key] = value
    xdmf_path = coefficient_xdmf_path(path)
    write_coefficient_xdmf(xdmf_path, path, args)
    return xdmf_path


def is_face_staggered(direction: int, var: int) -> bool:
    return direction == var


def sphere_inside(x: np.ndarray, center: np.ndarray, radius: float) -> bool:
    return float(np.dot(x - center, x - center)) < radius * radius


def sphere_segment_distance(xa: np.ndarray, xb: np.ndarray, center: np.ndarray, radius: float) -> float:
    d = xb - xa
    q = xa - center
    a = float(np.dot(d, d))
    b = 2.0 * float(np.dot(q, d))
    c = float(np.dot(q, q) - radius * radius)
    disc = b * b - 4.0 * a * c
    if disc < 0.0 and disc > -1.0e-14:
        disc = 0.0
    if disc < 0.0:
        raise RuntimeError("segment does not intersect analytic sphere")
    sqrt_disc = math.sqrt(disc)
    roots = [(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)]
    roots = [t for t in roots if -1.0e-12 <= t <= 1.0 + 1.0e-12]
    if not roots:
        raise RuntimeError("sphere intersection falls outside neighbor segment")
    t = min(max(0.0, min(1.0, t)) for t in roots)
    return t * math.sqrt(a)


def analytic_sphere_coeff(args: argparse.Namespace, center: np.ndarray, radius: float) -> np.ndarray:
    """Reference implementation for analytic sphere coefficients used to validate STL paths."""
    finalize_grid_args(args)
    nx, ny, nz = int(args.nx), int(args.ny), int(args.nz)
    coef = np.zeros((nx + 2, ny + 2, nz + 2, 3), dtype=np.float64)
    re_inv = 1.0 / float(args.re)
    solid_coef = SOLID * re_inv
    offsets = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]
    for var in (VAR_U, VAR_V, VAR_W):
        out_var = var - 1
        x = stl_axis_coords(np.arange(-1, nx + 3, dtype=np.int64), 1, var, args)
        y = stl_axis_coords(np.arange(-1, ny + 3, dtype=np.int64), 2, var, args)
        z = stl_axis_coords(np.arange(-1, nz + 3, dtype=np.int64), 3, var, args)
        for k in range(nz + 2):
            for j in range(ny + 2):
                for i in range(nx + 2):
                    xa = np.array([x[i + 1], y[j + 1], z[k + 1]], dtype=np.float64)
                    if sphere_inside(xa, center, radius):
                        coef[i, j, k, out_var] = solid_coef
                        continue
                    value = 0.0
                    for di, dj, dk in offsets:
                        xb = np.array([x[i + di + 1], y[j + dj + 1], z[k + dk + 1]], dtype=np.float64)
                        if not sphere_inside(xb, center, radius):
                            continue
                        d0 = float(np.linalg.norm(xb - xa))
                        d = sphere_segment_distance(xa, xb, center, radius)
                        value += ((d0 - d) / d) / (d0 * d0)
                    coef[i, j, k, out_var] = value * re_inv
    return coef


def require_stl_tools():
    """Import optional STL dependencies lazily so non-STL commands stay lightweight."""
    try:
        import igl
        import trimesh
        from trimesh.ray.ray_triangle import RayMeshIntersector
    except ImportError as exc:
        raise SystemExit(
            "STL support requires the geometry venv, e.g.\n"
            "  /home/davide/ibmc/bin/python tools/mobygeom.py ...\n"
            "with numpy, h5py, scipy, trimesh, rtree, and libigl installed."
        ) from exc
    return igl, trimesh, RayMeshIntersector



def load_points(path: Path, *, label: str = "point cloud",
                dataset_names: tuple[str, ...] = ("points",)) -> np.ndarray:
    """Load a generic 3D point cloud from text, CSV, NumPy, or HDF5."""
    path = Path(path).resolve()
    if path.suffix.lower() == ".npy":
        points = np.load(path)
    elif path.suffix.lower() in (".h5", ".hdf5"):
        with h5py.File(path, "r") as h5:
            for name in dataset_names:
                if name in h5:
                    points = h5[name][...]
                    break
            else:
                quoted = ", ".join(f"'{name}'" for name in dataset_names)
                raise ValueError(f"{path} must contain one of: {quoted}")
    else:
        try:
            points = np.loadtxt(path, delimiter=",")
        except ValueError:
            points = np.loadtxt(path)
    points = np.asarray(points, dtype=np.float64)
    if points.ndim == 1:
        points = points.reshape(1, -1)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"{label} {path} must have shape (n,3)")
    if points.shape[0] == 0:
        raise ValueError(f"{label} {path} is empty")
    return np.ascontiguousarray(points, dtype=np.float64)


def load_fluid_points(path: Path) -> np.ndarray:
    return load_points(path, label="fluid point cloud", dataset_names=("fluid_points", "points"))


def sphere_validation_fluid_points(args: argparse.Namespace, center: np.ndarray, radius: float) -> np.ndarray:
    coords = []
    xs = [0.0, 0.5 * args.lx, args.lx]
    ys = [0.0, 0.5 * args.ly, args.ly]
    zs = [0.0, 0.5 * args.lz, args.lz]
    for x in xs:
        for y in ys:
            for z in zs:
                on_boundary = x in (0.0, args.lx) or y in (0.0, args.ly) or z in (0.0, args.lz)
                if on_boundary:
                    p = np.array([x, y, z], dtype=np.float64)
                    if np.linalg.norm(p - center) > 1.25 * radius:
                        coords.append(p)
    points = np.unique(np.asarray(coords, dtype=np.float64), axis=0)
    if points.shape[0] == 0:
        raise ValueError("could not create validation fluid point cloud outside the analytic sphere")
    return np.ascontiguousarray(points, dtype=np.float64)


def prepare_fluid_points(args: argparse.Namespace, *, center: np.ndarray | None = None,
                         radius: float | None = None, write_path: Path | None = None) -> np.ndarray | None:
    """Load or build known-fluid probes used to reinforce STL inside/outside classification."""
    if hasattr(args, "_fluid_points"):
        return args._fluid_points
    points = None
    source = getattr(args, "fluid_points", None)
    if source:
        points = load_fluid_points(Path(source))
    elif getattr(args, "auto_fluid_points", False) and center is not None and radius is not None:
        points = sphere_validation_fluid_points(args, center, radius)
        if write_path is not None:
            np.savetxt(write_path, points, fmt="%.17e")
            args.fluid_points = str(write_path)
    if points is not None and getattr(args, "fluid_ray_max_seeds", 0) > 0 and points.shape[0] > args.fluid_ray_max_seeds:
        ids = np.linspace(0, points.shape[0] - 1, args.fluid_ray_max_seeds, dtype=np.int64)
        points = points[ids]
    args._fluid_points = points
    return points

def count_ray_hits_to_fluid(mesh, origins: np.ndarray, directions: np.ndarray,
                            lengths: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    """Count unique triangle crossings along rays from query points to known-fluid seeds."""
    counts = np.zeros(origins.shape[0], dtype=np.int32)
    if origins.shape[0] == 0:
        return counts
    intersector = getattr(args, "_ray_intersector", None)
    if intersector is None:
        _, _, RayMeshIntersector = require_stl_tools()
        intersector = RayMeshIntersector(mesh)
        args._ray_intersector = intersector
    locations, index_ray, _ = intersector.intersects_location(origins, directions, multiple_hits=True)
    if len(locations) == 0:
        return counts
    projected = np.einsum("ij,ij->i", locations - origins[index_ray], directions[index_ray])
    base_tol = max(1.0e-14, float(args.fluid_ray_hit_tol))
    tol = np.maximum(base_tol, base_tol * lengths[index_ray])
    valid = (projected > tol) & (projected < lengths[index_ray] - tol)
    if not np.any(valid):
        return counts
    ray = index_ray[valid]
    dist = projected[valid]
    order = np.lexsort((dist, ray))
    ray = ray[order]
    dist = dist[order]
    current_ray = -1
    last_dist = -np.inf
    for ray_id, hit_dist in zip(ray, dist):
        unique_tol = max(base_tol, base_tol * lengths[ray_id])
        if ray_id != current_ray:
            counts[ray_id] += 1
            current_ray = int(ray_id)
            last_dist = float(hit_dist)
        elif abs(float(hit_dist) - last_dist) > unique_tol:
            counts[ray_id] += 1
            last_dist = float(hit_dist)
    return counts


def classify_stl_points_by_fluid_rays(mesh, points: np.ndarray,
                                       fluid_points: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, int, np.ndarray]:
    """Classify points by majority ray parity to known-fluid probes."""
    points = np.ascontiguousarray(points, dtype=np.float64)
    fluid_points = np.ascontiguousarray(fluid_points, dtype=np.float64)
    inside = np.empty(points.shape[0], dtype=bool)
    vote_fraction_all = np.empty(points.shape[0], dtype=np.float64)
    ambiguous = 0
    chunk_size = max(1, int(args.chunk_size))
    majority = float(args.fluid_ray_majority)
    margin = float(args.fluid_ray_margin)
    for start in range(0, points.shape[0], chunk_size):
        stop = min(points.shape[0], start + chunk_size)
        query = points[start:stop]
        inside_votes = np.zeros(query.shape[0], dtype=np.int32)
        valid_votes = np.zeros(query.shape[0], dtype=np.int32)
        for seed in fluid_points:
            delta = seed.reshape(1, 3) - query
            lengths = np.linalg.norm(delta, axis=1)
            active = lengths > 1.0e-14
            valid_votes += 1
            if not np.any(active):
                continue
            dirs = np.zeros_like(delta)
            dirs[active] = delta[active] / lengths[active, None]
            counts = count_ray_hits_to_fluid(mesh, query[active], dirs[active], lengths[active], args)
            odd = (counts % 2 == 1)
            if getattr(args, "inside_is_fluid", False):
                inside_votes[active] += (~odd).astype(np.int32)
            else:
                inside_votes[active] += odd.astype(np.int32)
        vote_fraction = inside_votes.astype(np.float64) / np.maximum(1, valid_votes)
        vote_fraction_all[start:stop] = vote_fraction
        inside[start:stop] = vote_fraction > majority
        ambiguous += int(np.count_nonzero(np.abs(vote_fraction - majority) <= margin))
    return inside, ambiguous, vote_fraction_all


def axis_name(direction: int) -> str:
    return AXIS_NAMES[int(direction)]


def grid_size_and_length(direction: int, args: argparse.Namespace) -> tuple[int, float]:
    axis = axis_name(direction)
    if getattr(args, f"n{axis}", None) is None or getattr(args, f"l{axis}", None) is None:
        finalize_grid_args(args)
    return int(getattr(args, f"n{axis}")), float(getattr(args, f"l{axis}"))


def grid_periodic(direction: int, args: argparse.Namespace) -> bool:
    return bool(getattr(args, f"periodic_{axis_name(direction)}", False))


def validate_axis_nodes(nodes: np.ndarray, *, source: Path, axis: str,
                        expected_size: int, length: float) -> np.ndarray:
    """Validate one solver face/node coordinate line."""
    nodes = np.asarray(nodes, dtype=np.float64).reshape(-1)
    if nodes.size != expected_size:
        raise ValueError(f"{source} {axis}-nodes has {nodes.size} nodes, expected {expected_size}")
    if np.any(np.diff(nodes) <= 0.0):
        raise ValueError(f"{source} {axis}-nodes must be strictly increasing")
    tol = max(1.0e-12, 1.0e-10 * abs(length))
    if abs(nodes[0]) > tol or abs(nodes[-1] - length) > tol:
        raise ValueError(f"{source} {axis}-node endpoints must be 0 and {length}")
    return nodes


def h5_scalar_attr(attrs, name: str):
    """Return a scalar HDF5 attribute as a Python value, or None if missing/non-scalar."""
    if name not in attrs:
        return None
    value = np.asarray(attrs[name])
    if value.shape != ():
        return None
    return value.item()


def assign_or_validate_grid_arg(args: argparse.Namespace, name: str, value, source: Path) -> None:
    """Fill a grid argument from a grid file, or reject an explicit mismatch."""
    current = getattr(args, name, None)
    if current is None:
        setattr(args, name, value)
        return
    if isinstance(value, bool):
        if bool(current) != value:
            raise ValueError(f"{source} has {name}={value}, but command line has {current}")
        return
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        if int(current) != int(value):
            raise ValueError(f"{source} has {name}={int(value)}, but command line has {current}")
        return
    if isinstance(value, (float, np.floating)):
        current_f = float(current)
        value_f = float(value)
        tol = max(1.0e-12, 1.0e-10 * abs(value_f))
        if abs(current_f - value_f) > tol:
            raise ValueError(f"{source} has {name}={value_f}, but command line has {current_f}")
        return
    if str(current) != str(value):
        raise ValueError(f"{source} has {name}={value}, but command line has {current}")


def load_grid_file(args: argparse.Namespace) -> None:
    """Load solver-exported grid nodes and metadata into the argparse namespace."""
    grid_file = getattr(args, "grid_file", None)
    if not grid_file:
        raise ValueError("provide --grid-file from the mobygrid executable")
    path = Path(grid_file).resolve()
    if getattr(args, "_grid_file_loaded", None) == str(path):
        return
    for axis in ("x", "y", "z"):
        if getattr(args, f"{axis}_nodes", None) is not None:
            raise ValueError(f"use either --grid-file or --{axis}-nodes, not both")

    axis_nodes: dict[int, np.ndarray] = {}
    with h5py.File(path, "r") as h5:
        if int(np.asarray(h5.attrs.get("mobygrid_format", 0)).item()) != 1:
            raise ValueError(f"{path} is not a mobygrid HDF5 file")
        for direction, axis in AXIS_NAMES.items():
            dataset = f"{axis}_nodes"
            if dataset not in h5:
                raise ValueError(f"{path} does not contain /{dataset}")
            nodes = np.asarray(h5[dataset][...], dtype=np.float64).reshape(-1)
            if nodes.size < 2:
                raise ValueError(f"{path} {axis}-nodes must contain at least two values")
            n = nodes.size - 1
            length_attr = h5_scalar_attr(h5.attrs, f"l{axis}")
            length = float(length_attr) if length_attr is not None else float(nodes[-1])
            assign_or_validate_grid_arg(args, f"n{axis}", int(n), path)
            assign_or_validate_grid_arg(args, f"l{axis}", length, path)
            axis_nodes[direction] = validate_axis_nodes(
                nodes, source=path, axis=axis, expected_size=n + 1, length=length
            )

        periodic = np.asarray(h5.attrs["periodic"], dtype=np.int64).reshape(-1) if "periodic" in h5.attrs else None
        if periodic is not None and periodic.size == 3:
            for direction, axis in AXIS_NAMES.items():
                assign_or_validate_grid_arg(args, f"periodic_{axis}", bool(periodic[direction - 1]), path)

        distribution = (
            np.asarray(h5.attrs["grid_distribution"], dtype=np.int64).reshape(-1)
            if "grid_distribution" in h5.attrs else None
        )
        if distribution is not None and distribution.size == 3:
            for direction, axis in AXIS_NAMES.items():
                value = GRID_DISTRIBUTION_IDS.get(int(distribution[direction - 1]), "uniform")
                assign_or_validate_grid_arg(args, f"{axis}_distribution", value, path)

        stretch = np.asarray(h5.attrs["grid_stretch"], dtype=np.float64).reshape(-1) if "grid_stretch" in h5.attrs else None
        if stretch is not None and stretch.size == 3:
            for direction, axis in AXIS_NAMES.items():
                assign_or_validate_grid_arg(args, f"{axis}_stretch", float(stretch[direction - 1]), path)

        natural_one_sided = (
            np.asarray(h5.attrs["grid_natural_one_sided"], dtype=np.int64).reshape(-1)
            if "grid_natural_one_sided" in h5.attrs else None
        )
        if natural_one_sided is not None and natural_one_sided.size == 3:
            for direction, axis in AXIS_NAMES.items():
                assign_or_validate_grid_arg(
                    args, f"{axis}_natural_one_sided", bool(natural_one_sided[direction - 1]), path
                )

    cache = getattr(args, "_grid_axis_nodes", None)
    if cache is None:
        cache = {}
        args._grid_axis_nodes = cache
    cache.update(axis_nodes)
    args.grid_file = str(path)
    args._grid_file_loaded = str(path)


def finalize_grid_args(args: argparse.Namespace) -> None:
    """Resolve grid arguments, optionally from --grid-file, before coordinates are used."""
    if getattr(args, "_grid_finalized", False):
        return
    load_grid_file(args)
    for axis in ("x", "y", "z"):
        n_name = f"n{axis}"
        l_name = f"l{axis}"
        if getattr(args, n_name, None) is None:
            raise ValueError(f"grid file did not provide --{n_name}")
        if getattr(args, l_name, None) is None:
            setattr(args, l_name, 1.0)
        if getattr(args, f"{axis}_distribution", None) is None:
            setattr(args, f"{axis}_distribution", "uniform")
        if getattr(args, f"{axis}_stretch", None) is None:
            setattr(args, f"{axis}_stretch", 0.0)
        if getattr(args, f"{axis}_natural_one_sided", None) is None:
            setattr(args, f"{axis}_natural_one_sided", False)
        if getattr(args, f"periodic_{axis}", None) is None:
            setattr(args, f"periodic_{axis}", False)
    args.nx = int(args.nx)
    args.ny = int(args.ny)
    args.nz = int(args.nz)
    args.lx = float(args.lx)
    args.ly = float(args.ly)
    args.lz = float(args.lz)
    args.re = float(getattr(args, "re", 100.0))
    args._grid_finalized = True


def grid_axis_nodes(direction: int, args: argparse.Namespace) -> np.ndarray:
    """Return solver-compatible face/node coordinates for one global direction."""
    finalize_grid_args(args)
    cache = getattr(args, "_grid_axis_nodes", None)
    if cache is not None and direction in cache:
        return cache[direction]
    raise ValueError(f"grid file did not load {axis_name(direction)}-nodes")


def face_at(nodes: np.ndarray, n: int, length: float, indices: np.ndarray, periodic: bool) -> np.ndarray:
    """Vectorized version of the solver's face_at helper, including one/two-cell halos."""
    idx = np.asarray(indices, dtype=np.int64)
    out = np.empty(idx.shape, dtype=np.float64)

    mid = (idx >= 0) & (idx <= n)
    out[mid] = nodes[idx[mid]]

    low = idx < 0
    if np.any(low):
        if periodic:
            src = idx[low] + n
            if np.any((src < 0) | (src > n)):
                raise ValueError("periodic lower halo index exceeds available node line")
            out[low] = nodes[src] - length
        else:
            src = -idx[low]
            if np.any(src > n):
                raise ValueError("lower halo index exceeds available node line")
            out[low] = 2.0 * nodes[0] - nodes[src]

    high = idx > n
    if np.any(high):
        if periodic:
            src = idx[high] - n
            if np.any((src < 0) | (src > n)):
                raise ValueError("periodic upper halo index exceeds available node line")
            out[high] = nodes[src] + length
        else:
            src = 2 * n - idx[high]
            if np.any(src < 0):
                raise ValueError("upper halo index exceeds available node line")
            out[high] = 2.0 * nodes[n] - nodes[src]
    return out


def cell_center_at(nodes: np.ndarray, n: int, length: float, indices: np.ndarray, periodic: bool) -> np.ndarray:
    return 0.5 * (
        face_at(nodes, n, length, np.asarray(indices, dtype=np.int64) - 1, periodic) +
        face_at(nodes, n, length, np.asarray(indices, dtype=np.int64), periodic)
    )


def stl_axis_coords(indices: np.ndarray, direction: int, var: int, args: argparse.Namespace) -> np.ndarray:
    """Return staggered-grid coordinates from mobygrid node lines."""
    indices_i = np.asarray(indices, dtype=np.int64)
    n, axis_length = grid_size_and_length(direction, args)
    nodes = grid_axis_nodes(direction, args)
    periodic = grid_periodic(direction, args)
    if is_face_staggered(direction, var):
        return face_at(nodes, n, axis_length, indices_i - 1, periodic)
    return cell_center_at(nodes, n, axis_length, indices_i, periodic)


def max_axis_spacing(direction: int, args: argparse.Namespace) -> float:
    nodes = grid_axis_nodes(direction, args)
    return float(np.max(np.diff(nodes)))


def bbox_padding(args: argparse.Namespace) -> np.ndarray:
    padding_cells = max(0.0, float(getattr(args, "bbox_padding_cells", 2.0)))
    return padding_cells * np.array([
        max_axis_spacing(1, args),
        max_axis_spacing(2, args),
        max_axis_spacing(3, args),
    ], dtype=np.float64)


def grid_is_nonuniform(args: argparse.Namespace) -> bool:
    cache = getattr(args, "_grid_axis_nodes", None)
    for direction in (1, 2, 3):
        axis = axis_name(direction)
        if cache is not None and direction in cache:
            dx = np.diff(np.asarray(cache[direction], dtype=np.float64))
            if dx.size > 0:
                tol = max(1.0e-12, 1.0e-10 * float(np.max(np.abs(dx))))
                if np.max(np.abs(dx - dx[0])) > tol:
                    return True
        if getattr(args, f"{axis}_nodes", None):
            return True
        if (getattr(args, f"{axis}_distribution", None) or "uniform") != "uniform":
            return True
    return False


def stl_points_for_indices(ii: np.ndarray, jj: np.ndarray, kk: np.ndarray, var: int,
                           args: argparse.Namespace) -> np.ndarray:
    """Vectorized staggered-grid coordinate builder for selected index triplets."""
    return np.column_stack((
        stl_axis_coords(ii, 1, var, args),
        stl_axis_coords(jj, 2, var, args),
        stl_axis_coords(kk, 3, var, args),
    ))


def stl_extended_grid_points(var: int, args: argparse.Namespace) -> tuple[np.ndarray, tuple[int, int, int]]:
    """Build a one-cell-expanded staggered grid used to classify neighbor crossings."""
    _, _, _, x, y, z = stl_extended_grid_axes(var, args)
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    return points, xx.shape


def stl_extended_grid_axes(var: int, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                                                        np.ndarray, np.ndarray, np.ndarray]:
    """Return extended staggered indices and coordinates used by the STL coefficient path."""
    ix = np.arange(-1, args.nx + 3, dtype=np.int64)
    iy = np.arange(-1, args.ny + 3, dtype=np.int64)
    iz = np.arange(-1, args.nz + 3, dtype=np.int64)
    x = stl_axis_coords(ix, 1, var, args)
    y = stl_axis_coords(iy, 2, var, args)
    z = stl_axis_coords(iz, 3, var, args)
    return ix, iy, iz, x, y, z


def stl_axis_window(coords: np.ndarray, lower: float, upper: float) -> slice:
    """Return the smallest contiguous index slice covering coordinates inside [lower, upper]."""
    active = np.flatnonzero((coords >= lower) & (coords <= upper))
    if active.size == 0:
        return slice(0, 0)
    return slice(int(active[0]), int(active[-1]) + 1)


def stl_bbox_classification_window(mesh, var: int, args: argparse.Namespace) -> tuple[slice, slice, slice,
                                                                                      tuple[int, int, int],
                                                                                      tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Build a conservative extended-grid window around the transformed STL bounding box."""
    ix, iy, iz, x, y, z = stl_extended_grid_axes(var, args)
    lower, upper = np.asarray(mesh.bounds, dtype=np.float64)
    padding = bbox_padding(args)

    sx = stl_axis_window(x, lower[0] - padding[0], upper[0] + padding[0])
    sy = stl_axis_window(y, lower[1] - padding[1], upper[1] + padding[1])
    sz = stl_axis_window(z, lower[2] - padding[2], upper[2] + padding[2])
    shape = (x.size, y.size, z.size)
    return sx, sy, sz, shape, (ix, iy, iz)


def stl_window_points(ix: np.ndarray, iy: np.ndarray, iz: np.ndarray,
                      sx: slice, sy: slice, sz: slice, var: int,
                      args: argparse.Namespace) -> tuple[np.ndarray, tuple[int, int, int]]:
    """Build staggered physical points for one rectangular extended-grid subwindow."""
    wx = ix[sx]
    wy = iy[sy]
    wz = iz[sz]
    if wx.size == 0 or wy.size == 0 or wz.size == 0:
        return np.empty((0, 3), dtype=np.float64), (wx.size, wy.size, wz.size)
    xx, yy, zz = np.meshgrid(
        stl_axis_coords(wx, 1, var, args),
        stl_axis_coords(wy, 2, var, args),
        stl_axis_coords(wz, 3, var, args),
        indexing="ij",
    )
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    return points, xx.shape


def classify_stl_extended_grid(mesh, vertices: np.ndarray, faces: np.ndarray, var: int,
                               args: argparse.Namespace) -> tuple[np.ndarray, int, int, int]:
    """Classify the extended staggered grid, optionally only near the STL bounding box."""
    if getattr(args, "inside_is_fluid", False) or getattr(args, "no_bbox_cull", False):
        points, shape = stl_extended_grid_points(var, args)
        inside_flat, ambiguous = classify_stl_points(mesh, vertices, faces, points, args)
        return inside_flat.reshape(shape), ambiguous, int(points.shape[0]), int(points.shape[0])

    sx, sy, sz, shape, indices = stl_bbox_classification_window(mesh, var, args)
    points, window_shape = stl_window_points(indices[0], indices[1], indices[2], sx, sy, sz, var, args)
    inside_ext = np.zeros(shape, dtype=bool)
    if points.shape[0] == 0:
        return inside_ext, 0, 0, int(np.prod(shape))

    inside_flat, ambiguous = classify_stl_points(mesh, vertices, faces, points, args)
    inside_ext[sx, sy, sz] = inside_flat.reshape(window_shape)
    return inside_ext, ambiguous, int(points.shape[0]), int(np.prod(shape))


def stl_tile_extended_indices(i0: int, i1: int, j0: int, j1: int,
                              k0: int, k1: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return output-tile indices plus the one neighbor layer needed for IBM coefficients."""
    return (
        np.arange(i0 - 1, i1 + 1, dtype=np.int64),
        np.arange(j0 - 1, j1 + 1, dtype=np.int64),
        np.arange(k0 - 1, k1 + 1, dtype=np.int64),
    )


def classify_stl_tile_extended_grid(mesh, vertices: np.ndarray, faces: np.ndarray, var: int,
                                    tile: tuple[int, int, int, int, int, int],
                                    args: argparse.Namespace) -> tuple[np.ndarray, int, int, int]:
    """Classify the extended grid for one output tile, with optional bbox culling."""
    i0, i1, j0, j1, k0, k1 = tile
    ix, iy, iz = stl_tile_extended_indices(i0, i1, j0, j1, k0, k1)
    x = stl_axis_coords(ix, 1, var, args)
    y = stl_axis_coords(iy, 2, var, args)
    z = stl_axis_coords(iz, 3, var, args)
    shape = (x.size, y.size, z.size)
    total_points = int(np.prod(shape))
    inside_ext = np.zeros(shape, dtype=bool)

    if getattr(args, "inside_is_fluid", False) or getattr(args, "no_bbox_cull", False):
        sx, sy, sz = slice(None), slice(None), slice(None)
    else:
        lower, upper = np.asarray(mesh.bounds, dtype=np.float64)
        padding = bbox_padding(args)
        sx = stl_axis_window(x, lower[0] - padding[0], upper[0] + padding[0])
        sy = stl_axis_window(y, lower[1] - padding[1], upper[1] + padding[1])
        sz = stl_axis_window(z, lower[2] - padding[2], upper[2] + padding[2])

    points, window_shape = stl_window_points(ix, iy, iz, sx, sy, sz, var, args)
    if points.shape[0] == 0:
        return inside_ext, 0, 0, total_points

    inside_flat, ambiguous = classify_stl_points(mesh, vertices, faces, points, args)
    inside_ext[sx, sy, sz] = inside_flat.reshape(window_shape)
    return inside_ext, ambiguous, int(points.shape[0]), total_points


def load_stl_mesh(path: Path, *, repair: bool = True):
    """Load one STL mesh and apply conservative repairs that do not assume a specific geometry."""
    _, trimesh, _ = require_stl_tools()
    mesh = trimesh.load_mesh(path, force="mesh", process=True)
    if hasattr(mesh, "geometry") and not hasattr(mesh, "faces"):
        parts = [g for g in mesh.geometry.values() if hasattr(g, "faces") and len(g.faces) > 0]
        if not parts:
            raise ValueError(f"no triangle meshes found in {path}")
        mesh = trimesh.util.concatenate(parts)
    if repair:
        for method in ("merge_vertices", "remove_duplicate_faces", "remove_degenerate_faces", "remove_unreferenced_vertices"):
            if hasattr(mesh, method):
                try:
                    getattr(mesh, method)()
                except Exception:
                    pass
        try:
            trimesh.repair.fix_normals(mesh, multibody=True)
        except Exception:
            pass
    if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
        raise ValueError(f"empty STL mesh in {path}")
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    return mesh, vertices, faces


def load_stl_meshes(paths: list[Path], *, repair: bool = True):
    """Load and concatenate multiple STL shells, for example a body plus optional lids."""
    _, trimesh, _ = require_stl_tools()
    if len(paths) == 1:
        return load_stl_mesh(paths[0], repair=repair)
    meshes = [load_stl_mesh(path, repair=repair)[0] for path in paths]
    mesh = trimesh.util.concatenate(meshes)
    if repair:
        try:
            mesh.merge_vertices()
            mesh.remove_unreferenced_vertices()
            trimesh.repair.fix_normals(mesh, multibody=True)
        except Exception:
            pass
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    return mesh, vertices, faces


def stl_transform_values(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    """Return per-axis scale and translation vectors for STL geometry."""
    scale_values = getattr(args, "scale", [1.0])
    if isinstance(scale_values, (float, int)):
        scale_values = [float(scale_values)]
    if len(scale_values) == 1:
        scale = np.full(3, float(scale_values[0]), dtype=np.float64)
    elif len(scale_values) == 3:
        scale = np.asarray(scale_values, dtype=np.float64)
    else:
        raise SystemExit("--scale expects either one value or three values")

    translate = np.asarray(getattr(args, "translate", [0.0, 0.0, 0.0]), dtype=np.float64)
    if translate.shape != (3,):
        raise SystemExit("--translate expects three values")
    return scale, translate


def apply_stl_transform(mesh, args: argparse.Namespace):
    """Apply optional scale/translation to an STL mesh and return updated arrays."""
    scale, translate = stl_transform_values(args)
    if np.allclose(scale, 1.0) and np.allclose(translate, 0.0):
        vertices = np.asarray(mesh.vertices, dtype=np.float64)
        faces = np.asarray(mesh.faces, dtype=np.int64)
        return mesh, vertices, faces

    mesh = mesh.copy()
    mesh.vertices = np.asarray(mesh.vertices, dtype=np.float64)*scale + translate
    try:
        mesh.remove_unreferenced_vertices()
    except Exception:
        pass
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    return mesh, vertices, faces


def classify_stl_points(mesh, vertices: np.ndarray, faces: np.ndarray, points: np.ndarray,
                        args: argparse.Namespace) -> tuple[np.ndarray, int]:
    """Classify STL query points with winding numbers and optional known-fluid reinforcement."""
    igl, _, _ = require_stl_tools()
    chunk_size = max(1, int(args.chunk_size))
    inside = np.empty(points.shape[0], dtype=bool)
    winding_margin = float(args.winding_margin)
    winding_ambiguous = 0
    magnitudes = np.empty(points.shape[0], dtype=np.float64)

    for start in range(0, points.shape[0], chunk_size):
        stop = min(points.shape[0], start + chunk_size)
        query = np.ascontiguousarray(points[start:stop], dtype=np.float64)
        if args.classification == "winding":
            winding = igl.winding_number(vertices, faces, query)
        else:
            winding = igl.fast_winding_number(vertices, faces, query)
        magnitude = np.abs(winding)
        magnitudes[start:stop] = magnitude
        inside[start:stop] = magnitude >= args.winding_threshold
        winding_ambiguous += int(np.count_nonzero(np.abs(magnitude - args.winding_threshold) <= winding_margin))

    fluid_points = prepare_fluid_points(args)
    if fluid_points is None:
        return inside, winding_ambiguous

    winding_near = np.abs(magnitudes - args.winding_threshold) <= winding_margin
    if getattr(args, "fluid_ray_scope", "ambiguous") == "all":
        ray_mask = np.ones(points.shape[0], dtype=bool)
    else:
        ray_mask = winding_near

    if not np.any(ray_mask):
        return inside, winding_ambiguous

    ray_inside, ray_ambiguous, vote_fraction = classify_stl_points_by_fluid_rays(
        mesh, points[ray_mask], fluid_points, args
    )
    ray_confident = np.abs(vote_fraction - args.fluid_ray_majority) > args.fluid_ray_margin
    inside_ray_subset = inside[ray_mask]
    disagreements_subset = ray_confident & (ray_inside != inside_ray_subset)
    if args.fluid_ray_policy == "check":
        override_subset = np.zeros_like(ray_inside, dtype=bool)
    elif args.fluid_ray_policy == "override":
        override_subset = disagreements_subset
    else:
        override_subset = disagreements_subset & winding_near[ray_mask]

    ray_indices = np.flatnonzero(ray_mask)
    inside[ray_indices[override_subset]] = ray_inside[override_subset]

    args._fluid_ray_checks = getattr(args, "_fluid_ray_checks", 0) + int(ray_indices.size)
    args._fluid_ray_ambiguous = getattr(args, "_fluid_ray_ambiguous", 0) + int(ray_ambiguous)
    args._fluid_ray_disagreements = getattr(args, "_fluid_ray_disagreements", 0) + int(np.count_nonzero(disagreements_subset))
    args._fluid_ray_overrides = getattr(args, "_fluid_ray_overrides", 0) + int(np.count_nonzero(override_subset))
    return inside, winding_ambiguous + ray_ambiguous

def stl_inside_point(mesh, vertices: np.ndarray, faces: np.ndarray, point: np.ndarray, args: argparse.Namespace) -> bool:
    query = np.asarray(point, dtype=np.float64).reshape(1, 3)
    inside, _ = classify_stl_points(mesh, vertices, faces, query, args)
    return bool(inside[0])

def winding_root_distance(mesh, vertices: np.ndarray, faces: np.ndarray, xa: np.ndarray, xb: np.ndarray,
                          length: float, args: argparse.Namespace) -> float:
    """Fallback distance estimate by bisecting the winding-number solid/fluid transition."""
    lo = 0.0
    hi = 1.0
    inside_hi = stl_inside_point(mesh, vertices, faces, xb, args)
    if stl_inside_point(mesh, vertices, faces, xa, args) == inside_hi:
        return 0.5 * length
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        xm = xa + mid * (xb - xa)
        if stl_inside_point(mesh, vertices, faces, xm, args) == inside_hi:
            hi = mid
        else:
            lo = mid
    return hi * length


def stl_segment_distances(mesh, vertices: np.ndarray, faces: np.ndarray, xa: np.ndarray, xb: np.ndarray,
                          args: argparse.Namespace) -> tuple[np.ndarray, int]:
    """Find distances from fluid points to crossed STL surfaces, with winding fallback."""
    if xa.shape[0] == 0:
        return np.empty(0, dtype=np.float64), 0
    delta = xb - xa
    length = np.linalg.norm(delta, axis=1)
    direction = delta / length[:, None]
    distance = np.full(xa.shape[0], np.nan, dtype=np.float64)

    if args.distance_mode != "winding-root":
        intersector = getattr(args, "_ray_intersector", None)
        if intersector is None:
            _, _, RayMeshIntersector = require_stl_tools()
            intersector = RayMeshIntersector(mesh)
            args._ray_intersector = intersector
        locations, index_ray, _ = intersector.intersects_location(xa, direction, multiple_hits=True)
        if len(locations) > 0:
            projected = np.einsum("ij,ij->i", locations - xa[index_ray], direction[index_ray])
            tol = np.maximum(1.0e-12, 1.0e-10 * length[index_ray])
            valid = (projected > tol) & (projected <= length[index_ray] + tol)
            for ray, dist in zip(index_ray[valid], projected[valid]):
                if np.isnan(distance[ray]) or dist < distance[ray]:
                    distance[ray] = dist

    missing = np.flatnonzero(~np.isfinite(distance))
    for idx in missing:
        distance[idx] = winding_root_distance(mesh, vertices, faces, xa[idx], xb[idx], length[idx], args)
    return np.maximum(distance, 1.0e-14), int(missing.size)


def stl_coeff_from_mesh(mesh, vertices: np.ndarray, faces: np.ndarray,
                        args: argparse.Namespace) -> tuple[np.ndarray, dict[str, int]]:
    """Assemble staggered IBM coefficients by classifying points and checking solid neighbors."""
    finalize_grid_args(args)
    coef = np.zeros((args.nx + 2, args.ny + 2, args.nz + 2, 3), dtype=np.float64)
    re_inv = 1.0 / args.re
    solid_coef = SOLID * re_inv
    offsets = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]
    stats = {
        "ambiguous_points": 0, "crossing_segments": 0, "fallback_segments": 0,
        "fluid_ray_checks": 0, "fluid_ray_ambiguous": 0,
        "fluid_ray_disagreements": 0, "fluid_ray_overrides": 0,
        "classified_points": 0, "total_points": 0, "bbox_culled_points": 0,
    }
    for name in ("_fluid_ray_checks", "_fluid_ray_ambiguous", "_fluid_ray_disagreements", "_fluid_ray_overrides"):
        setattr(args, name, 0)
    args._ray_intersector = None

    for var in (VAR_U, VAR_V, VAR_W):
        inside_ext, ambiguous, classified_points, total_points = classify_stl_extended_grid(
            mesh, vertices, faces, var, args
        )
        solid_ext = ~inside_ext if getattr(args, "inside_is_fluid", False) else inside_ext
        stats["ambiguous_points"] += ambiguous
        stats["classified_points"] += classified_points
        stats["total_points"] += total_points
        stats["bbox_culled_points"] += total_points - classified_points

        solid = solid_ext[1:args.nx + 3, 1:args.ny + 3, 1:args.nz + 3]
        value = np.zeros_like(solid, dtype=np.float64)
        for di, dj, dk in offsets:
            neighbor = solid_ext[
                1 + di:1 + di + args.nx + 2,
                1 + dj:1 + dj + args.ny + 2,
                1 + dk:1 + dk + args.nz + 2,
            ]
            crossing = (~solid) & neighbor
            if not np.any(crossing):
                continue
            ii, jj, kk = np.nonzero(crossing)
            xa = stl_points_for_indices(ii, jj, kk, var, args)
            xb = stl_points_for_indices(ii + di, jj + dj, kk + dk, var, args)
            d0 = np.linalg.norm(xb - xa, axis=1)
            d, fallback_count = stl_segment_distances(mesh, vertices, faces, xa, xb, args)
            value[crossing] += ((d0 - d) / d) / (d0 * d0)
            stats["crossing_segments"] += int(ii.size)
            stats["fallback_segments"] += fallback_count

        coef[..., var - 1] = np.where(solid, solid_coef, value * re_inv)

    stats["fluid_ray_checks"] = int(getattr(args, "_fluid_ray_checks", 0))
    stats["fluid_ray_ambiguous"] = int(getattr(args, "_fluid_ray_ambiguous", 0))
    stats["fluid_ray_disagreements"] = int(getattr(args, "_fluid_ray_disagreements", 0))
    stats["fluid_ray_overrides"] = int(getattr(args, "_fluid_ray_overrides", 0))
    return coef, stats


_TILE_CONTEXT = None


def stl_coeff_tile_from_mesh(mesh, vertices: np.ndarray, faces: np.ndarray,
                             args: argparse.Namespace,
                             tile: tuple[int, int, int, int, int, int]) -> tuple[np.ndarray, dict[str, int]]:
    """Compute one output coefficient tile using a one-point local classification overlap."""
    finalize_grid_args(args)
    i0, i1, j0, j1, k0, k1 = tile
    tile_shape = (i1 - i0, j1 - j0, k1 - k0)
    coef = np.zeros((*tile_shape, 3), dtype=np.float64)
    re_inv = 1.0 / args.re
    solid_coef = SOLID * re_inv
    offsets = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]
    stats = {
        "ambiguous_points": 0, "crossing_segments": 0, "fallback_segments": 0,
        "fluid_ray_checks": 0, "fluid_ray_ambiguous": 0,
        "fluid_ray_disagreements": 0, "fluid_ray_overrides": 0,
        "classified_points": 0, "total_points": 0, "bbox_culled_points": 0,
    }

    for name in ("_fluid_ray_checks", "_fluid_ray_ambiguous", "_fluid_ray_disagreements", "_fluid_ray_overrides"):
        setattr(args, name, 0)
    args._ray_intersector = None

    for var in (VAR_U, VAR_V, VAR_W):
        inside_ext, ambiguous, classified_points, total_points = classify_stl_tile_extended_grid(
            mesh, vertices, faces, var, tile, args
        )
        solid_ext = ~inside_ext if getattr(args, "inside_is_fluid", False) else inside_ext
        stats["ambiguous_points"] += ambiguous
        stats["classified_points"] += classified_points
        stats["total_points"] += total_points
        stats["bbox_culled_points"] += total_points - classified_points

        solid = solid_ext[1:-1, 1:-1, 1:-1]
        value = np.zeros_like(solid, dtype=np.float64)
        for di, dj, dk in offsets:
            neighbor = solid_ext[
                1 + di:1 + di + tile_shape[0],
                1 + dj:1 + dj + tile_shape[1],
                1 + dk:1 + dk + tile_shape[2],
            ]
            crossing = (~solid) & neighbor
            if not np.any(crossing):
                continue
            ii, jj, kk = np.nonzero(crossing)
            xa = stl_points_for_indices(ii + i0, jj + j0, kk + k0, var, args)
            xb = stl_points_for_indices(ii + i0 + di, jj + j0 + dj, kk + k0 + dk, var, args)
            d0 = np.linalg.norm(xb - xa, axis=1)
            d, fallback_count = stl_segment_distances(mesh, vertices, faces, xa, xb, args)
            value[crossing] += ((d0 - d) / d) / (d0 * d0)
            stats["crossing_segments"] += int(ii.size)
            stats["fallback_segments"] += fallback_count

        coef[..., var - 1] = np.where(solid, solid_coef, value * re_inv)

    stats["fluid_ray_checks"] = int(getattr(args, "_fluid_ray_checks", 0))
    stats["fluid_ray_ambiguous"] = int(getattr(args, "_fluid_ray_ambiguous", 0))
    stats["fluid_ray_disagreements"] = int(getattr(args, "_fluid_ray_disagreements", 0))
    stats["fluid_ray_overrides"] = int(getattr(args, "_fluid_ray_overrides", 0))
    return coef, stats


def init_tile_worker(mesh, vertices: np.ndarray, faces: np.ndarray, args: argparse.Namespace) -> None:
    """Store read-only STL coefficient context once per worker process."""
    global _TILE_CONTEXT
    _TILE_CONTEXT = (mesh, vertices, faces, args)


def compute_stl_coeff_tile_worker(tile: tuple[int, int, int, int, int, int]):
    """Multiprocessing worker entry point for one coefficient tile."""
    mesh, vertices, faces, args = _TILE_CONTEXT
    coef, stats = stl_coeff_tile_from_mesh(mesh, vertices, faces, args, tile)
    return tile, coef, stats


def stl_output_tiles(args: argparse.Namespace, tile_size: tuple[int, int, int]) -> list[tuple[int, int, int, int, int, int]]:
    """Enumerate coefficient output tiles over the full solver coefficient domain."""
    nx = int(args.nx) + 2
    ny = int(args.ny) + 2
    nz = int(args.nz) + 2
    tx, ty, tz = tile_size
    tiles = []
    for i0 in range(0, nx, tx):
        i1 = min(nx, i0 + tx)
        for j0 in range(0, ny, ty):
            j1 = min(ny, j0 + ty)
            for k0 in range(0, nz, tz):
                k1 = min(nz, k0 + tz)
                tiles.append((i0, i1, j0, j1, k0, k1))
    return tiles


def stl_tile_size(args: argparse.Namespace) -> tuple[int, int, int]:
    """Return tile dimensions for chunked STL coefficient generation."""
    values = getattr(args, "tile_size", None)
    if values is None:
        return (
            min(int(args.nx) + 2, 64),
            min(int(args.ny) + 2, 64),
            min(int(args.nz) + 2, 64),
        )
    if len(values) != 3:
        raise SystemExit("--tile-size expects three integers")
    tile = tuple(max(1, int(v)) for v in values)
    return tile


def hdf5_compression_kwargs(args: argparse.Namespace) -> dict[str, object]:
    """Translate CLI compression options to h5py.create_dataset keyword arguments."""
    compression = getattr(args, "h5_compression", "none")
    if compression == "none":
        return {}
    if compression == "gzip":
        return {"compression": "gzip", "compression_opts": int(getattr(args, "gzip_level", 1))}
    return {"compression": "lzf"}


def create_component_virtual_datasets(h5: h5py.File, path: Path, shape: tuple[int, int, int]) -> None:
    """Expose coef_u/v/w as no-copy virtual datasets backed by the 4D coef dataset."""
    labels = ("u", "v", "w")
    source = h5py.VirtualSource(path.name, "coef", shape=(*shape, 3), dtype=np.float64)
    for var, label in enumerate(labels):
        layout = h5py.VirtualLayout(shape=shape, dtype=np.float64)
        layout[:, :, :] = source[:, :, :, var]
        h5.create_virtual_dataset(f"coef_{label}", layout)


def write_chunked_coefficient_xdmf(xdmf_path: Path, h5_path: Path, args: argparse.Namespace) -> None:
    """Write a lightweight XDMF sidecar for chunked coefficient files."""
    xdmf_path = Path(xdmf_path).resolve()
    h5_path = Path(h5_path).resolve()
    nx, ny, nz = int(args.nx) + 2, int(args.ny) + 2, int(args.nz) + 2
    grids = []
    for var_name in ("u", "v", "w"):
        grids.append(f"""      <Grid Name=\"coef_{var_name}\" GridType=\"Uniform\">
        <Topology TopologyType=\"3DRectMesh\" Dimensions=\"{nx} {ny} {nz}\"/>
        <Geometry GeometryType=\"VXVYVZ\">
          <DataItem Name=\"X\" Dimensions=\"{nx}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{xdmf_hdf_reference(xdmf_path, h5_path, f"/xdmf/{var_name}/x")}</DataItem>
          <DataItem Name=\"Y\" Dimensions=\"{ny}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{xdmf_hdf_reference(xdmf_path, h5_path, f"/xdmf/{var_name}/y")}</DataItem>
          <DataItem Name=\"Z\" Dimensions=\"{nz}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{xdmf_hdf_reference(xdmf_path, h5_path, f"/xdmf/{var_name}/z")}</DataItem>
        </Geometry>
        <Attribute Name=\"coef_{var_name}\" AttributeType=\"Scalar\" Center=\"Node\">
          <DataItem Dimensions=\"{nx} {ny} {nz}\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">{xdmf_hdf_reference(xdmf_path, h5_path, f"/coef_{var_name}")}</DataItem>
        </Attribute>
      </Grid>
""")
    xdmf = f"""<?xml version=\"1.0\" ?>
<Xdmf Version=\"2.0\">
  <Domain>
    <Grid Name=\"mobyDiff IBM coefficients\" GridType=\"Collection\" CollectionType=\"Spatial\">
{''.join(grids)}    </Grid>
  </Domain>
</Xdmf>
"""
    xdmf_path.write_text(xdmf)


def write_chunked_stl_hdf5(path: Path, mesh, vertices: np.ndarray, faces: np.ndarray,
                           args: argparse.Namespace, geometry: str) -> tuple[dict[str, int], Path]:
    """Compute STL coefficients in independent tiles and write them to one chunked HDF5 file."""
    finalize_grid_args(args)
    path = Path(path).resolve()
    tile_size = stl_tile_size(args)
    chunks = (
        min(tile_size[0], int(args.nx) + 2),
        min(tile_size[1], int(args.ny) + 2),
        min(tile_size[2], int(args.nz) + 2),
        3,
    )
    shape = (int(args.nx) + 2, int(args.ny) + 2, int(args.nz) + 2)
    stats = {
        "ambiguous_points": 0, "crossing_segments": 0, "fallback_segments": 0,
        "fluid_ray_checks": 0, "fluid_ray_ambiguous": 0,
        "fluid_ray_disagreements": 0, "fluid_ray_overrides": 0,
        "classified_points": 0, "total_points": 0, "bbox_culled_points": 0,
    }
    tiles = stl_output_tiles(args, tile_size)
    jobs = max(1, int(getattr(args, "jobs", 1)))
    kwargs = hdf5_compression_kwargs(args)

    with h5py.File(path, "w") as h5:
        coef_dset = h5.create_dataset("coef", shape=(*shape, 3), dtype=np.float64,
                                      chunks=chunks, fillvalue=0.0, **kwargs)
        xdmf_group = h5.create_group("xdmf")
        indices = {
            1: np.arange(0, int(args.nx) + 2, dtype=np.float64),
            2: np.arange(0, int(args.ny) + 2, dtype=np.float64),
            3: np.arange(0, int(args.nz) + 2, dtype=np.float64),
        }
        extents = {1: (int(args.nx), float(args.lx)),
                   2: (int(args.ny), float(args.ly)),
                   3: (int(args.nz), float(args.lz))}
        for var_name, var in (("u", VAR_U), ("v", VAR_V), ("w", VAR_W)):
            group = xdmf_group.create_group(var_name)
            for direction, coord_name in ((1, "x"), (2, "y"), (3, "z")):
                group.create_dataset(coord_name, data=stl_axis_coords(indices[direction], direction, var, args))

        h5.attrs["nx"] = int(args.nx)
        h5.attrs["ny"] = int(args.ny)
        h5.attrs["nz"] = int(args.nz)
        h5.attrs["lx"] = float(args.lx)
        h5.attrs["ly"] = float(args.ly)
        h5.attrs["lz"] = float(args.lz)
        h5.attrs["re"] = float(args.re)
        h5.attrs["source_geometry"] = geometry
        h5.attrs["convention"] = "mobyDiff staggered u/v/w coefficients, shape (nx+2,ny+2,nz+2,3)"
        h5.attrs["chunked_preprocessing"] = 1
        h5.attrs["preprocess_jobs"] = jobs
        h5.attrs["tile_size"] = np.array(tile_size, dtype=np.int64)
        h5.attrs["h5_compression"] = getattr(args, "h5_compression", "none")
        h5.attrs["xdmf_file"] = coefficient_xdmf_path(path).name
        write_grid_metadata(h5, args)

        if jobs == 1:
            for tile in tiles:
                coef_tile, tile_stats = stl_coeff_tile_from_mesh(mesh, vertices, faces, args, tile)
                i0, i1, j0, j1, k0, k1 = tile
                if np.any(coef_tile):
                    coef_dset[i0:i1, j0:j1, k0:k1, :] = coef_tile
                for key in stats:
                    stats[key] += int(tile_stats.get(key, 0))
        else:
            try:
                context = multiprocessing.get_context("fork")
            except ValueError:
                context = multiprocessing.get_context()
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=jobs,
                mp_context=context,
                initializer=init_tile_worker,
                initargs=(mesh, vertices, faces, args),
            ) as pool:
                futures = [pool.submit(compute_stl_coeff_tile_worker, tile) for tile in tiles]
                for future in concurrent.futures.as_completed(futures):
                    tile, coef_tile, tile_stats = future.result()
                    i0, i1, j0, j1, k0, k1 = tile
                    if np.any(coef_tile):
                        coef_dset[i0:i1, j0:j1, k0:k1, :] = coef_tile
                    for key in stats:
                        stats[key] += int(tile_stats.get(key, 0))

        create_component_virtual_datasets(h5, path, shape)
        for key, value in stl_extra_attrs(mesh, stats, args).items():
            h5.attrs[key] = value

    xdmf_path = coefficient_xdmf_path(path)
    write_chunked_coefficient_xdmf(xdmf_path, path, args)
    return stats, xdmf_path


def stl_extra_attrs(mesh, stats: dict[str, int], args: argparse.Namespace) -> dict[str, object]:
    """Collect STL diagnostic metadata stored alongside the coefficient field."""
    scale, translate = stl_transform_values(args)
    return {
        "geometry_type": "stl",
        "stl_vertices": int(len(mesh.vertices)),
        "stl_faces": int(len(mesh.faces)),
        "stl_is_watertight": int(bool(mesh.is_watertight)),
        "stl_euler_number": int(mesh.euler_number),
        "stl_classification": args.classification,
        "stl_winding_threshold": float(args.winding_threshold),
        "stl_distance_mode": args.distance_mode,
        "stl_inside_is_fluid": int(bool(getattr(args, "inside_is_fluid", False))),
        "stl_fluid_points": 0 if getattr(args, "_fluid_points", None) is None else int(args._fluid_points.shape[0]),
        "stl_fluid_ray_policy": getattr(args, "fluid_ray_policy", "ambiguous"),
        "stl_fluid_ray_scope": getattr(args, "fluid_ray_scope", "ambiguous"),
        "stl_fluid_ray_checks": int(stats.get("fluid_ray_checks", 0)),
        "stl_fluid_ray_ambiguous": int(stats.get("fluid_ray_ambiguous", 0)),
        "stl_fluid_ray_disagreements": int(stats.get("fluid_ray_disagreements", 0)),
        "stl_fluid_ray_overrides": int(stats.get("fluid_ray_overrides", 0)),
        "stl_ambiguous_points": int(stats["ambiguous_points"]),
        "stl_crossing_segments": int(stats["crossing_segments"]),
        "stl_fallback_segments": int(stats["fallback_segments"]),
        "stl_classified_points": int(stats.get("classified_points", 0)),
        "stl_total_points": int(stats.get("total_points", 0)),
        "stl_bbox_culled_points": int(stats.get("bbox_culled_points", 0)),
        "stl_bbox_culling": int(not bool(getattr(args, "inside_is_fluid", False)) and
                                not bool(getattr(args, "no_bbox_cull", False))),
        "stl_bbox_padding_cells": float(getattr(args, "bbox_padding_cells", 2.0)),
        "stl_grid_nonuniform": int(grid_is_nonuniform(args)),
        "stl_scale": scale,
        "stl_translate": translate,
    }


def stl_solid_mask_for_points(mesh, vertices: np.ndarray, faces: np.ndarray,
                              points: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    """Convert inside/outside classification into solid/fluid meaning for probe checks."""
    inside, _ = classify_stl_points(mesh, vertices, faces, np.asarray(points, dtype=np.float64), args)
    return ~inside if getattr(args, "inside_is_fluid", False) else inside


def load_check_points(args: argparse.Namespace, attr_name: str) -> np.ndarray | None:
    """Load generic must-be-fluid or must-be-solid probe points."""
    embedded = getattr(args, f"_{attr_name}", None)
    if embedded is not None:
        points = np.asarray(embedded, dtype=np.float64)
        if points.ndim == 1:
            points = points.reshape(1, -1)
        return points
    path = getattr(args, attr_name, None)
    if path is None:
        return None
    dataset_names = (attr_name, attr_name.removeprefix("check_"), "points")
    return load_points(Path(path), label=attr_name.replace("_", " "), dataset_names=dataset_names)


def generic_stl_checks(mesh, vertices: np.ndarray, faces: np.ndarray,
                       args: argparse.Namespace) -> dict[str, int | bool]:
    """Run geometry-independent mesh health and probe-classification checks."""
    boundary_edges, nonmanifold_edges = mesh_edge_counts(mesh)
    report: dict[str, int | bool] = {
        "watertight": bool(mesh.is_watertight),
        "boundary_edges": int(boundary_edges),
        "nonmanifold_edges": int(nonmanifold_edges),
        "fluid_check_points": 0,
        "fluid_check_failures": 0,
        "solid_check_points": 0,
        "solid_check_failures": 0,
    }

    fluid_points = load_check_points(args, "check_fluid_points")
    if fluid_points is not None:
        solid = stl_solid_mask_for_points(mesh, vertices, faces, fluid_points, args)
        report["fluid_check_points"] = int(fluid_points.shape[0])
        report["fluid_check_failures"] = int(np.count_nonzero(solid))

    solid_points = load_check_points(args, "check_solid_points")
    if solid_points is not None:
        solid = stl_solid_mask_for_points(mesh, vertices, faces, solid_points, args)
        report["solid_check_points"] = int(solid_points.shape[0])
        report["solid_check_failures"] = int(np.count_nonzero(~solid))
    return report


def print_generic_stl_check_report(report: dict[str, int | bool]) -> None:
    print(f"watertight: {report['watertight']}")
    print(f"boundary_edges: {report['boundary_edges']}")
    print(f"nonmanifold_edges: {report['nonmanifold_edges']}")
    if report["fluid_check_points"]:
        print(f"fluid check points: {report['fluid_check_points']}")
        print(f"fluid check failures: {report['fluid_check_failures']}")
    if report["solid_check_points"]:
        print(f"solid check points: {report['solid_check_points']}")
        print(f"solid check failures: {report['solid_check_failures']}")


def enforce_generic_stl_checks(report: dict[str, int | bool], args: argparse.Namespace,
                               *, label: str = "STL geometry") -> None:
    """Turn generic STL check results into fail-fast CLI errors when limits are exceeded."""
    if getattr(args, "expect_watertight", False) and not report["watertight"]:
        raise SystemExit(f"{label} check failed: mesh is not watertight")
    if report["boundary_edges"] > getattr(args, "max_boundary_edges", 10**18):
        raise SystemExit(
            f"{label} check failed: {report['boundary_edges']} boundary edges > {args.max_boundary_edges}"
        )
    if report["nonmanifold_edges"] > getattr(args, "max_nonmanifold_edges", 10**18):
        raise SystemExit(
            f"{label} check failed: {report['nonmanifold_edges']} non-manifold edges > {args.max_nonmanifold_edges}"
        )
    if report["fluid_check_failures"] > getattr(args, "max_fluid_check_failures", 0):
        raise SystemExit(
            f"{label} check failed: {report['fluid_check_failures']} fluid probe failures "
            f"> {args.max_fluid_check_failures}"
        )
    if report["solid_check_failures"] > getattr(args, "max_solid_check_failures", 0):
        raise SystemExit(
            f"{label} check failed: {report['solid_check_failures']} solid probe failures "
            f"> {args.max_solid_check_failures}"
        )


def check_stl_geometry(args: argparse.Namespace) -> None:
    """Standalone CLI entry point for validating arbitrary STL geometry before coefficient generation."""
    geometries = [Path(value).resolve() for value in args.geometry]
    mesh, vertices, faces = load_stl_meshes(geometries, repair=not args.no_repair)
    mesh, vertices, faces = apply_stl_transform(mesh, args)
    prepare_fluid_points(args)
    report = generic_stl_checks(mesh, vertices, faces, args)
    print(f"STL mesh: {';'.join(str(path) for path in geometries)}")
    print(f"vertices/faces: {len(mesh.vertices)} / {len(mesh.faces)}")
    print(f"inside_is_fluid: {bool(getattr(args, 'inside_is_fluid', False))}")
    print_generic_stl_check_report(report)
    enforce_generic_stl_checks(report, args)


def coeff_from_stl(args: argparse.Namespace) -> np.ndarray:
    """Generate IBM coefficients from one or more STL meshes and run optional generic checks."""
    finalize_grid_args(args)
    geometry_values = args.geometry if isinstance(args.geometry, list) else [args.geometry]
    geometries = [Path(value).resolve() for value in geometry_values]
    output = Path(args.output).resolve()
    mesh, vertices, faces = load_stl_meshes(geometries, repair=not args.no_repair)
    mesh, vertices, faces = apply_stl_transform(mesh, args)
    prepare_fluid_points(args)
    source_geometry = ";".join(str(path) for path in geometries)
    production_command = getattr(args, "command", "") == "stl-ibm-coeff"
    use_chunked = (
        (production_command and not bool(getattr(args, "no_tiled_output", False))) or
        int(getattr(args, "jobs", 1)) > 1 or
        getattr(args, "tile_size", None) is not None
    )
    if use_chunked:
        stats, xdmf_path = write_chunked_stl_hdf5(output, mesh, vertices, faces, args, source_geometry)
    else:
        coef, stats = stl_coeff_from_mesh(mesh, vertices, faces, args)
        xdmf_path = write_hdf5(output, coef, args, source_geometry, stl_extra_attrs(mesh, stats, args))
    print(f"STL mesh: {source_geometry}")
    print(f"vertices/faces: {len(mesh.vertices)} / {len(mesh.faces)}")
    print(f"watertight: {mesh.is_watertight}")
    if stats["total_points"] > 0:
        kept = 100.0*float(stats["classified_points"])/float(stats["total_points"])
        print(f"classified grid points: {stats['classified_points']} / {stats['total_points']} ({kept:.2f}%)")
    print(f"ambiguous classification points: {stats['ambiguous_points']}")
    if stats.get("fluid_ray_checks", 0) > 0:
        print(f"fluid-point ray checks: {stats['fluid_ray_checks']}")
        print(f"fluid-point ray disagreements: {stats['fluid_ray_disagreements']}")
        print(f"fluid-point ray overrides: {stats['fluid_ray_overrides']}")
    print(f"crossing segments: {stats['crossing_segments']}")
    print(f"winding-root fallback segments: {stats['fallback_segments']}")
    check_report = generic_stl_checks(mesh, vertices, faces, args)
    if check_report["fluid_check_points"] or check_report["solid_check_points"] or getattr(args, "expect_watertight", False):
        print_generic_stl_check_report(check_report)
    enforce_generic_stl_checks(check_report, args, label="STL coefficient geometry")
    print(f"HDF5 coefficients: {output}")
    print(f"XDMF coefficients: {xdmf_path}")
    return None if use_chunked else coef


def window_solid_counts(solid_ext: np.ndarray, nb: int,
                        gnbt: tuple[int, int, int]) -> np.ndarray:
    """Per-lattice-block count of solid points in the one-halo dilated window.

    solid_ext covers extended indices -1..n+2 (cell index + 1). The dilated
    window of the block at lattice (bx,by,bz) spans cells [o, o+nb+1] with
    o = b*nb, i.e. extended slice [o+1, o+nb+3). Windows of neighbouring
    blocks overlap, so the reduction uses a 3D integral image — built in
    x-lattice chunks, because the int64 cumsum at a deep refinement level's
    full grid does not fit in memory (the airfoil L4 lattice would need
    ~70 GB; integer arithmetic makes the chunked result identical).
    """
    w = nb + 2
    total = np.empty(gnbt, dtype=np.int64)
    ny3, nz3 = solid_ext.shape[1] + 1, solid_ext.shape[2] + 1
    rows = max(1, int((2 << 30) // (ny3*nz3*8*nb)))   # ~2 GB cnt slabs
    lo_yz = [np.arange(g)*nb + 1 for g in gnbt[1:]]
    hi_yz = [l + w for l in lo_yz]
    for b0 in range(0, gnbt[0], rows):
        b1 = min(b0 + rows, gnbt[0])
        # extended x-cells feeding blocks b0..b1-1: [b0*nb+1 - 1, (b1-1)*nb+w]
        x_from = b0*nb            # one layer before the first window start
        x_to = (b1 - 1)*nb + 1 + w
        sub = solid_ext[x_from:x_to]
        cnt = np.zeros((sub.shape[0] + 1, ny3, nz3), dtype=np.int64)
        cnt[1:, 1:, 1:] = sub
        cnt = cnt.cumsum(axis=0).cumsum(axis=1).cumsum(axis=2)
        lo_x = np.arange(b0, b1)*nb + 1 - x_from
        hi_x = lo_x + w
        ix0, iy0, iz0 = np.meshgrid(lo_x, lo_yz[0], lo_yz[1], indexing="ij")
        ix1, iy1, iz1 = np.meshgrid(hi_x, hi_yz[0], hi_yz[1], indexing="ij")
        total[b0:b1] = (cnt[ix1, iy1, iz1] - cnt[ix0, iy1, iz1]
                        - cnt[ix1, iy0, iz1] - cnt[ix1, iy1, iz0]
                        + cnt[ix0, iy0, iz1] + cnt[ix0, iy1, iz0]
                        + cnt[ix1, iy0, iz0] - cnt[ix0, iy0, iz0])
    return total


def window_all_solid(solid_ext: np.ndarray, nb: int, gnbt: tuple[int, int, int]) -> np.ndarray:
    """Per-lattice-block test: every point of the one-halo dilated window solid."""
    return window_solid_counts(solid_ext, nb, gnbt) == (nb + 2)**3


def refine_mask(args: argparse.Namespace) -> np.ndarray:
    """--refine-dims as a per-direction mask (1 = the direction halves per
    level, 0 = fixed), the [blocks] refine_dims counterpart."""
    dims = getattr(args, "refine_dims", "xyz") or "xyz"
    if dims == "xyz":
        return np.array([1, 1, 1], dtype=np.int64)
    if dims == "xz":
        return np.array([1, 0, 1], dtype=np.int64)
    raise SystemExit("--refine-dims must be xyz or xz")


def subdivided_args(args: argparse.Namespace, level: int) -> argparse.Namespace:
    """Copy of the grid args at refinement level `level`: sizes doubled per
    level and node lines midpoint-subdivided in the refined directions
    (--refine-dims; fixed directions keep the global line), matching the
    solver's build_level_lines so coefficients live at identical
    coordinates."""
    finalize_grid_args(args)
    out = argparse.Namespace(**vars(args))
    mask = refine_mask(args)
    f = 2**(level*mask)
    out.nx, out.ny, out.nz = int(args.nx)*int(f[0]), int(args.ny)*int(f[1]), int(args.nz)*int(f[2])
    nodes = {}
    for d in (1, 2, 3):
        line = grid_axis_nodes(d, args)
        for _ in range(level if mask[d - 1] else 0):
            fine = np.empty(2*(line.size - 1) + 1, dtype=np.float64)
            fine[0::2] = line
            fine[1::2] = 0.5*(line[:-1] + line[1:])
            line = fine
        nodes[d] = line
    out._grid_axis_nodes = nodes
    out._ray_intersector = None
    return out


def window_any_solid(solid_ext: np.ndarray, nb: int, gnbt: tuple[int, int, int]) -> np.ndarray:
    return window_solid_counts(solid_ext, nb, gnbt) > 0


def classify_stl_window(mesh, vertices, faces, var: int,
                        args: argparse.Namespace) -> tuple[np.ndarray, tuple[int, int, int]]:
    """Classify ONLY the bbox window of the extended staggered grid and
    return (window bool array, window offset in extended indices). The
    full-level raster of deep refinements does not fit in memory (the
    airfoil L5/L6 rasters are 69/550 GB as one bool array); everything
    outside the window is fluid by the bbox argument, which the windowed
    block-count reduction (window_solid_counts_win) encodes as zero."""
    sx, sy, sz, shape, indices = stl_bbox_classification_window(mesh, var, args)
    points, window_shape = stl_window_points(indices[0], indices[1], indices[2], sx, sy, sz, var, args)
    win = np.zeros(window_shape, dtype=bool)
    if points.shape[0]:
        inside_flat, _ = classify_stl_points(mesh, vertices, faces, points, args)
        win[...] = inside_flat.reshape(window_shape)
    return win, (sx.start or 0, sy.start or 0, sz.start or 0)


def window_solid_counts_win(win: np.ndarray, offs: tuple[int, int, int], nb: int,
                            gnbt: tuple[int, int, int],
                            blk_lo: tuple[int, int, int] = (0, 0, 0),
                            blk_dims: tuple[int, int, int] | None = None) -> np.ndarray:
    """window_solid_counts on a bbox-windowed solid raster: per-lattice-block
    count of solid points in the one-halo dilated window, everything outside
    the classification window counting as fluid (zero count). Same integer
    arithmetic as the full-raster version (clipped inclusion-exclusion on the
    window's integral image, built in x-lattice chunks like
    window_solid_counts), so the resulting masks are identical.
    blk_lo/blk_dims restrict the OUTPUT to a block window (deep-refinement
    lattices do not fit as dense arrays; blocks outside a padded STL bbox
    are all-fluid by construction)."""
    if blk_dims is None:
        blk_dims = gnbt
    w = nb + 2
    total = np.zeros(blk_dims, dtype=np.int64)
    ny1, nz1 = win.shape[1] + 1, win.shape[2] + 1
    by = blk_lo[1] + np.arange(blk_dims[1])
    bz = blk_lo[2] + np.arange(blk_dims[2])
    # y/z cnt-index ranges (cnt index = window row + 1), clipped
    lo_y = np.clip(by*nb + 1 - offs[1], 0, win.shape[1])
    hi_y = np.clip(by*nb + 1 - offs[1] + w, 0, win.shape[1])
    lo_z = np.clip(bz*nb + 1 - offs[2], 0, win.shape[2])
    hi_z = np.clip(bz*nb + 1 - offs[2] + w, 0, win.shape[2])
    rows = max(1, int((2 << 30) // (ny1*nz1*8*nb)))   # ~2 GB cnt slabs
    for b0 in range(0, blk_dims[0], rows):
        b1 = min(b0 + rows, blk_dims[0])
        bx = blk_lo[0] + np.arange(b0, b1)
        # window rows feeding these blocks (unclipped extended range)
        x_from = np.clip(bx[0]*nb + 1 - offs[0] - 1, 0, win.shape[0])
        x_to = np.clip(bx[-1]*nb + 1 + w - offs[0], 0, win.shape[0])
        if x_to <= x_from:
            continue   # chunk entirely outside the window: all fluid
        sub = win[x_from:x_to]
        cnt = np.zeros((sub.shape[0] + 1, ny1, nz1), dtype=np.int64)
        cnt[1:, 1:, 1:] = sub
        cnt = cnt.cumsum(axis=0).cumsum(axis=1).cumsum(axis=2)
        lo_x = np.clip(bx*nb + 1 - offs[0], x_from, x_to) - x_from
        hi_x = np.clip(bx*nb + 1 - offs[0] + w, x_from, x_to) - x_from
        ix0, iy0, iz0 = np.meshgrid(lo_x, lo_y, lo_z, indexing="ij")
        ix1, iy1, iz1 = np.meshgrid(hi_x, hi_y, hi_z, indexing="ij")
        total[b0:b1] = (cnt[ix1, iy1, iz1] - cnt[ix0, iy1, iz1]
                        - cnt[ix1, iy1, iz0] - cnt[ix1, iy0, iz1]
                        + cnt[ix0, iy0, iz1] + cnt[ix0, iy1, iz0]
                        + cnt[ix1, iy0, iz0] - cnt[ix0, iy0, iz0])
    return total


def touch_block_window(mesh, la, nb: int,
                       pad_blocks: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Block window (lo, dims) at this level that can possibly touch or
    bury against the STL: the mesh bbox in level-l block coords, padded.
    Unrefined (span) directions keep the full range (the extrusion
    crosses them)."""
    b = mesh.bounds
    lo = np.zeros(3, dtype=np.int64)
    hi = np.zeros(3, dtype=np.int64)
    for d in range(3):
        nodes = grid_axis_nodes(d + 1, la)
        nbl = (nodes.size - 1)//nb
        j0 = int(np.searchsorted(nodes, b[0, d], side="right") - 1)//nb
        j1 = int(np.searchsorted(nodes, b[1, d], side="left"))//nb
        lo[d] = max(0, j0 - pad_blocks)
        hi[d] = min(nbl, j1 + pad_blocks + 1)
    return lo, hi - lo


def level_masks(mesh, vertices, faces, args: argparse.Namespace, nb: int,
                levels: int):
    """Per-level touch/buried block masks, WINDOWED to the padded STL bbox
    (deep-refinement lattices do not fit as dense arrays; every block
    outside the window is untouched/unburied by construction). Returns
    (touch, buried, win_lo, win_dims): mask arrays have the window shape.
    The file counterpart of the solver's classify_block_geometry."""
    touch, buried, wlo, wdim = [], [], [], []
    inside_is_fluid = bool(getattr(args, "inside_is_fluid", False))
    for l in range(levels):
        la = subdivided_args(args, l)
        gnbt = (la.nx // nb, la.ny // nb, la.nz // nb)
        if inside_is_fluid:
            # Cavity geometries: outside the bbox window is SOLID, so the
            # windowed shortcut does not apply — full raster (only viable
            # at shallow refinement).
            blo, bdim = np.zeros(3, dtype=np.int64), np.asarray(gnbt, dtype=np.int64)
        else:
            blo, bdim = touch_block_window(mesh, la, nb)
        anySolid = np.zeros(tuple(bdim), dtype=bool)
        allSolid = np.ones(tuple(bdim), dtype=bool)
        anyFluid = np.zeros(tuple(bdim), dtype=bool)
        w3 = (nb + 2)**3
        for var in (VAR_U, VAR_V, VAR_W, 0):
            if inside_is_fluid:
                inside_ext, _, _, _ = classify_stl_extended_grid(mesh, vertices, faces, var, la)
                cnt = window_solid_counts(~inside_ext, nb, gnbt)
            else:
                # Solid-body geometries: everything outside the bbox window
                # is fluid — classify and reduce on the window only (the
                # full raster of deep refinements does not fit in memory).
                win, offs = classify_stl_window(mesh, vertices, faces, var, la)
                cnt = window_solid_counts_win(win, offs, nb, gnbt,
                                              blk_lo=tuple(blo), blk_dims=tuple(bdim))
            anySolid |= cnt > 0
            allSolid &= cnt == w3
            anyFluid |= cnt < w3
        touch.append(anySolid & anyFluid)
        buried.append(allSolid)
        wlo.append(blo)
        wdim.append(bdim)
    return touch, buried, wlo, wdim


def morton_key3(c: np.ndarray) -> np.ndarray:
    key = np.zeros(c.shape[0], dtype=np.int64)
    for bit in range(21):
        key |= ((c[:, 0].astype(np.int64) >> bit) & 1) << (3*bit)
        key |= ((c[:, 1].astype(np.int64) >> bit) & 1) << (3*bit + 1)
        key |= ((c[:, 2].astype(np.int64) >> bit) & 1) << (3*bit + 2)
    return key


def morton_key2(cx: np.ndarray, cz: np.ndarray) -> np.ndarray:
    key = np.zeros(cx.shape[0], dtype=np.int64)
    for bit in range(21):
        key |= ((cx.astype(np.int64) >> bit) & 1) << (2*bit)
        key |= ((cz.astype(np.int64) >> bit) & 1) << (2*bit + 1)
    return key


def leaf_keys(crd: np.ndarray, lev: np.ndarray, lmax: int, mask: np.ndarray) -> np.ndarray:
    """Canonical leaf ordering keys (blocks.f90 leaf_key): the 3D Morton
    key of the finest-lattice coords (xyz octree), or — xz quadtree — the
    y tile index in the high bits (42+) above the 2D (x,z) Morton key."""
    cf = crd*(2**((lmax - lev)[:, None]*mask[None, :]))
    if mask.all():
        return morton_key3(cf)
    return (cf[:, 1].astype(np.int64) << 42) | morton_key2(cf[:, 0], cf[:, 2])


def leaf_level_windows(gnbt, levels, mask, boxes, lines, nb, touch_wlo, touch_wdim,
                       pad=2):
    """Per-level OCCUPANCY windows (lo, dims) of the leaf builder, in
    level-l block coords: dense per-level lattices are impossible at deep
    refinement (the B11 finest lattice is 7.2e9 blocks), but level-l cells
    can only exist where level-(l-1) cells split — inside a refine box
    targeting >= l, near the body (touch + 1-block buffer), or within the
    2:1-smoothing spill of the finer window. Computed fine-to-coarse as
    conservative bbox hulls, padded; the builder ASSERTS every write stays
    inside (an undersized window fails loudly, never silently)."""
    levels_lat = [np.array(gnbt)*2**(l*mask) for l in range(levels)]
    lo = [None]*levels
    hi = [None]*levels
    # split-region hull S[l] (level-l block coords), fine to coarse
    S_lo = [None]*levels
    S_hi = [None]*levels
    for l in range(levels - 2, -1, -1):
        nl = levels_lat[l]
        slo = np.full(3, np.iinfo(np.int64).max, dtype=np.int64)
        shi = np.full(3, np.iinfo(np.int64).min, dtype=np.int64)
        have = False
        # refine boxes whose target level exceeds l split level-l cells
        if boxes:
            for b in boxes:
                tgt = int(b[6]) if len(b) >= 7 else levels - 1
                if tgt <= l:
                    continue
                for d in range(3):
                    ln = lines[d][l]
                    j0 = max(0, int(np.searchsorted(ln, b[2*d], side="right") - 1)//nb)
                    j1 = min(int(nl[d]), int(np.searchsorted(ln, b[2*d+1], side="left"))//nb + 1)
                    slo[d] = min(slo[d], j0)
                    shi[d] = max(shi[d], j1)
                have = True
        # touch + one-block buffer splits level-l cells
        if touch_wlo is not None:
            tlo, tdim = touch_wlo[l], touch_wdim[l]
            slo = np.minimum(slo, np.asarray(tlo) - 1)
            shi = np.maximum(shi, np.asarray(tlo) + np.asarray(tdim) + 1)
            have = True
        # 2:1 smoothing: level-l cells within 1 block of finer splits
        if S_lo[l+1] is not None:
            plo = S_lo[l+1]//(1 + mask) - 1
            phi = (S_hi[l+1] + mask)//(1 + mask) + 1
            slo = np.minimum(slo, plo)
            shi = np.maximum(shi, phi)
            have = True
        if have:
            S_lo[l] = np.maximum(0, slo - pad)
            S_hi[l] = np.minimum(nl, shi + pad)
    lo[0] = np.zeros(3, dtype=np.int64)
    hi[0] = levels_lat[0].astype(np.int64)
    for l in range(1, levels):
        if S_lo[l-1] is None:
            lo[l] = np.zeros(3, dtype=np.int64)
            hi[l] = np.zeros(3, dtype=np.int64)
        else:
            lo[l] = np.maximum(0, S_lo[l-1]*(1 + mask))
            hi[l] = np.minimum(levels_lat[l], S_hi[l-1]*(1 + mask))
    return lo, [h - l_ for l_, h in zip(lo, hi)]


def build_leaf_table_py(gnbt, levels, periodic, touch, buried, refine_box, lines, nb,
                        mask=None, touch_wlo=None, touch_wdim=None):
    """Mirror of the solver's build_leaf_table: root tiling, refinement by
    body-touch (+ one-block 26-neighbour buffer) and/or physical box, 2:1
    smoothing, removal of buried leaves at every level, ids along the
    finest-lattice Morton curve (mixed y-major form in xz quadtree mode).
    Must stay rule-for-rule identical to blocks.f90; the solver verifies
    the resulting table at startup. All per-level lattices are WINDOWED
    (leaf_level_windows); touch/buried may be windowed too (touch_wlo/
    touch_wdim; None = full rasters). Returns (lev, crd, win_lo, win_dims)."""
    M_NONE, M_SPLIT, M_LEAF = -1, -2, 0
    if mask is None:
        mask = np.array([1, 1, 1], dtype=np.int64)
    lmax = levels - 1
    boxes = refine_box
    if boxes and not hasattr(boxes[0], "__len__"):
        boxes = [boxes]
    if touch is not None and touch_wlo is None:
        touch_wlo = [np.zeros(3, dtype=np.int64)]*levels
        touch_wdim = [np.array(t.shape, dtype=np.int64) for t in touch]
    wlo, wdim = leaf_level_windows(gnbt, levels, mask, boxes, lines, nb,
                                   touch_wlo if touch is not None else None,
                                   touch_wdim if touch is not None else None)
    occ = [np.full(tuple(wdim[l]), M_NONE, dtype=np.int64) for l in range(levels)]
    occ[0][:, :, :] = M_LEAF

    def get(l, c):
        r = np.asarray(c) - wlo[l]
        if (r < 0).any() or (r >= wdim[l]).any():
            return M_NONE
        return occ[l][tuple(r)]

    def put(l, c, v):
        r = np.asarray(c) - wlo[l]
        assert (r >= 0).all() and (r < wdim[l]).all(), \
            f"leaf builder window undersized at level {l}: {c} outside " \
            f"{wlo[l]}+{wdim[l]}"
        occ[l][tuple(r)] = v

    def touch_at(l, c):
        if touch is None:
            return False
        r = np.asarray(c) - touch_wlo[l]
        if (r < 0).any() or (r >= touch_wdim[l]).any():
            return False
        return bool(touch[l][tuple(r)])

    def wrap(l, c):
        nl = np.array(gnbt)*2**(l*mask)
        c = np.array(c)
        for d in range(3):
            if c[d] < 0 or c[d] >= nl[d]:
                if not periodic[d]:
                    return None
                c[d] = c[d] % nl[d]
        return tuple(c)

    def split(l, c):
        put(l, c, M_SPLIT)
        o = np.array(c)*(1 + mask)
        for sx in range(1 + mask[0]):
            for sy in range(1 + mask[1]):
                for sz in range(1 + mask[2]):
                    put(l+1, (o[0]+sx, o[1]+sy, o[2]+sz), M_LEAF)

    for l in range(lmax):
        lvl_boxes = ([b for b in boxes if len(b) < 7 or b[6] > l]
                     if boxes else [])
        for cx in range(wlo[l][0], wlo[l][0] + wdim[l][0]):
            for cy in range(wlo[l][1], wlo[l][1] + wdim[l][1]):
                for cz in range(wlo[l][2], wlo[l][2] + wdim[l][2]):
                    if get(l, (cx, cy, cz)) != M_LEAF:
                        continue
                    hit = False
                    if lvl_boxes:
                        # An optional 7th box value is its TARGET LEVEL: it
                        # refines only rounds below it (absent = the finest,
                        # blocks.f90 rule).
                        lo = [lines[d][l][[cx, cy, cz][d]*nb] for d in range(3)]
                        hi = [lines[d][l][([cx, cy, cz][d]+1)*nb] for d in range(3)]
                        hit = any(all(box[2*d] < hi[d] and box[2*d+1] > lo[d]
                                      for d in range(3)) for box in lvl_boxes)
                    if not hit and touch is not None:
                        hit = touch_at(l, (cx, cy, cz))
                        if not hit:
                            for ox in (-1, 0, 1):
                                for oy in (-1, 0, 1):
                                    for oz in (-1, 0, 1):
                                        if ox == oy == oz == 0:
                                            continue
                                        cn = wrap(l, (cx+ox, cy+oy, cz+oz))
                                        if cn is not None and touch_at(l, cn):
                                            hit = True
                                            break
                                    if hit: break
                                if hit: break
                    if hit:
                        split(l, (cx, cy, cz))

    changed = True
    while changed:
        changed = False
        for l in range(lmax - 1):
            for cx in range(wlo[l][0], wlo[l][0] + wdim[l][0]):
                for cy in range(wlo[l][1], wlo[l][1] + wdim[l][1]):
                    for cz in range(wlo[l][2], wlo[l][2] + wdim[l][2]):
                        if get(l, (cx, cy, cz)) != M_LEAF:
                            continue
                        must = False
                        for ox in (-1, 0, 1):
                            for oy in (-1, 0, 1):
                                for oz in (-1, 0, 1):
                                    if ox == oy == oz == 0:
                                        continue
                                    cn = wrap(l, (cx+ox, cy+oy, cz+oz))
                                    if cn is None or get(l, cn) != M_SPLIT:
                                        continue
                                    o = np.array(cn)*(1 + mask)
                                    ch = [get(l+1, (o[0]+sx, o[1]+sy, o[2]+sz))
                                          for sx in range(1 + mask[0])
                                          for sy in range(1 + mask[1])
                                          for sz in range(1 + mask[2])]
                                    if M_SPLIT in ch:
                                        must = True
                                        break
                                if must: break
                            if must: break
                        if must:
                            split(l, (cx, cy, cz))
                            changed = True

    levels_out, coords_out = [], []
    for l in range(levels):
        leaf = np.argwhere(occ[l] == M_LEAF) + wlo[l]
        if buried is not None and leaf.size:
            r = leaf - touch_wlo[l]
            inside = ((r >= 0).all(axis=1) & (r < touch_wdim[l]).all(axis=1))
            bur = np.zeros(leaf.shape[0], dtype=bool)
            if inside.any():
                ri = r[inside]
                bur[inside] = buried[l][ri[:, 0], ri[:, 1], ri[:, 2]]
            leaf = leaf[~bur]
        if leaf.size:
            levels_out.append(np.full(leaf.shape[0], l, dtype=np.int64))
            coords_out.append(leaf)
    lev = np.concatenate(levels_out)
    crd = np.concatenate(coords_out)
    keys = leaf_keys(crd, lev, lmax, mask)
    order = np.argsort(keys, kind="stable")
    return lev[order], crd[order], wlo, wdim


def init_block_tile_worker(mesh, vertices: np.ndarray, faces: np.ndarray,
                           level_args: list, nb: int, want_dwall: bool) -> None:
    """Read-only per-worker context for the block-table tile pool."""
    global _BLOCK_TILE_CONTEXT
    _BLOCK_TILE_CONTEXT = (mesh, vertices, faces, level_args, nb, want_dwall)


def block_tile_worker(item):
    """One near-body leaf: exact coefficient (+ dwall) ghost-window tiles."""
    i, level, o = item
    mesh, vertices, faces, level_args, nb, want_dwall = _BLOCK_TILE_CONTEXT
    la = level_args[level]
    tile = (o[0], o[0]+nb+2, o[1], o[1]+nb+2, o[2], o[2]+nb+2)
    coef_tile, _ = stl_coeff_tile_from_mesh(mesh, vertices, faces, la, tile)
    dw = (dwall_tile_from_mesh(vertices, faces, la, np.asarray(o), nb)
          if want_dwall else None)
    return i, coef_tile, dw


def dwall_tile_from_mesh(vertices: np.ndarray, faces: np.ndarray,
                         la: argparse.Namespace, o: np.ndarray, nb: int) -> np.ndarray:
    """Cell-centred unsigned distance to the immersed surface over one leaf's
    ghost-inclusive (nb+2)^3 window, at the leaf's level (la = level args).
    Global 1-based cell indices o..o+nb+1 mirror the solver's local 0:nb+1.
    igl's AABB query is exact; trimesh's proximity.on_surface was measured
    O(1e-6) high near the quantized surface (A3 INCREMENT 0 gate)."""
    igl, _, _ = require_stl_tools()
    idx = np.arange(nb + 2, dtype=np.int64)
    ii, jj, kk = np.meshgrid(idx + int(o[0]), idx + int(o[1]), idx + int(o[2]), indexing="ij")
    points = stl_points_for_indices(ii.ravel(), jj.ravel(), kk.ravel(), 0, la)
    sq_dist, _, _ = igl.point_mesh_squared_distance(points, vertices, faces)
    return np.sqrt(sq_dist).reshape(nb + 2, nb + 2, nb + 2)


def block_table_from_stl(args: argparse.Namespace) -> None:
    """Write the block-table IBM coefficient file for refined runs: blocks
    (id -> origin, level), per-leaf coefficient windows evaluated at each
    leaf's level-l staggered locations, and the per-level touch/buried
    rasters the solver rebuilds its leaf table from."""
    finalize_grid_args(args)
    nb = int(args.block_nb)
    levels = int(args.levels)
    if nb < 4 or nb % 2:
        raise SystemExit("--block-nb must be even and >= 4")
    for n, name in ((args.nx, "nx"), (args.ny, "ny"), (args.nz, "nz")):
        if n % nb:
            raise SystemExit(f"--block-nb {nb} does not divide {name} = {n}")
    gnbt = (args.nx // nb, args.ny // nb, args.nz // nb)
    periodic = [grid_periodic(d, args) for d in (1, 2, 3)]
    # node lines per direction (0-based) and level for the box-refine rule
    lines = {d: [grid_axis_nodes(d + 1, subdivided_args(args, l)) for l in range(levels)]
             for d in range(3)}

    geometry_values = args.geometry if isinstance(args.geometry, list) else [args.geometry]
    geometries = [Path(value).resolve() for value in geometry_values]
    output = Path(args.output).resolve()
    mesh, vertices, faces = load_stl_meshes(geometries, repair=not args.no_repair)
    mesh, vertices, faces = apply_stl_transform(mesh, args)
    prepare_fluid_points(args)
    for name in ("_fluid_ray_checks", "_fluid_ray_ambiguous", "_fluid_ray_disagreements", "_fluid_ray_overrides"):
        setattr(args, name, 0)
    args._ray_intersector = None

    refine_box = getattr(args, "refine_box", None)
    if refine_box:
        for b in refine_box:
            if len(b) not in (6, 7):
                raise SystemExit("--refine-box takes 6 values or 6 + a target level")
            if len(b) == 7 and not (1 <= b[6] <= levels - 1):
                raise SystemExit(f"--refine-box level {b[6]} outside 1..{levels-1}")
    # Body classification always runs (the command requires --geometry);
    # --refine-box ADDS box refinement on top, mirroring the solver's
    # combined builder ([blocks] refine + refine_body). The solver ini must
    # carry the same [blocks] refine box or its builder cross-check errors.
    touch, buried, mwlo, mwdim = level_masks(mesh, vertices, faces, args, nb, levels)
    if getattr(args, "keep_buried", False):
        # Keep leaves buried inside the body (zeroed masks reach the
        # solver, whose own builder then also removes nothing). REQUIRED
        # for penalization-force cases: a removed core's closed faces
        # absorb the body's pressure loading with no coef bookkeeping,
        # so sum(coef u dV) misses most of the (pressure-dominated) lift.
        buried = [np.zeros_like(b) for b in buried]

    lev, crd, lwlo, lwdim = build_leaf_table_py(
        gnbt, levels, periodic, touch, buried, refine_box, lines, nb,
        mask=refine_mask(args), touch_wlo=mwlo, touch_wdim=mwdim)
    n_leaves = lev.shape[0]
    print(f"block table: {n_leaves} leaves, "
          f"{int((lev > 0).sum())} refined, levels {np.bincount(lev, minlength=levels)}")

    level_args = [subdivided_args(args, l) for l in range(levels)]
    with h5py.File(output, "w") as h5:
        h5.attrs["nx"] = int(args.nx)
        h5.attrs["ny"] = int(args.ny)
        h5.attrs["nz"] = int(args.nz)
        h5.attrs["lx"] = float(args.lx)
        h5.attrs["ly"] = float(args.ly)
        h5.attrs["lz"] = float(args.lz)
        h5.attrs["re"] = float(args.re)
        h5.attrs["block_nb"] = np.int32(nb)
        h5.attrs["block_levels"] = np.int32(levels)
        if not refine_mask(args).all():
            # xz quadtree file variant marker (same convention as the
            # solver's field files: absent = xyz octree).
            h5.attrs["refine_dims"] = refine_mask(args).astype(np.int32)
        h5.attrs["convention"] = ("mobyDiff block-table coefficients: coef_blocks row id = "
                                  "(nb+2)^3 ghost window x 3 staggered vars at the leaf's level")
        blocks = np.empty((n_leaves, 4), dtype=np.int32)
        blocks[:, 0:3] = (crd*nb).astype(np.int32)
        blocks[:, 3] = lev.astype(np.int32)
        h5.create_dataset("blocks", data=blocks)
        if touch is not None:
            # WINDOWED per-level mask rasters (x-fastest within the window)
            # + window attrs; the deep-refinement full rasters (7.2e9
            # blocks at B11's finest level) are unrepresentable. Legacy
            # files without the attrs keep the full-raster meaning.
            for l in range(levels):
                h5.create_dataset(f"block_touch_l{l}",
                                  data=touch[l].transpose(2, 1, 0).ravel().astype(np.int32))
                h5.create_dataset(f"block_buried_l{l}",
                                  data=buried[l].transpose(2, 1, 0).ravel().astype(np.int32))
                h5.attrs[f"mask_win_lo_l{l}"] = np.asarray(mwlo[l], dtype=np.int32)
                h5.attrs[f"mask_win_dims_l{l}"] = np.asarray(mwdim[l], dtype=np.int32)
                # the builder's occupancy windows: the solver adopts them
                # verbatim (windowed lidOf) so the two builders share the
                # exact same lattice bounds by construction.
                h5.attrs[f"lev_win_lo_l{l}"] = np.asarray(lwlo[l], dtype=np.int32)
                h5.attrs[f"lev_win_dims_l{l}"] = np.asarray(lwdim[l], dtype=np.int32)
        coef = h5.create_dataset("coef_blocks", shape=(n_leaves, nb+2, nb+2, nb+2, 3),
                                 dtype=np.float64, chunks=(1, nb+2, nb+2, nb+2, 3))
        want_dwall = not getattr(args, "no_dwall", False)
        dwall = None
        if want_dwall:
            # Per-leaf cell-centred distance to the immersed surface for the
            # RANS wall treatment, at each leaf's level like coef_blocks
            # (ghost-inclusive (nb+2)^3 windows; the solver mins in the
            # domain-wall part and applies the half-cell floor itself).
            dwall = h5.create_dataset("dwall_blocks", shape=(n_leaves, nb+2, nb+2, nb+2),
                                      dtype=np.float64, chunks=(1, nb+2, nb+2, nb+2))

        # Far-field leaf shortcut: a leaf whose ghost window (padded by 3
        # of its own cells) misses the STL bbox in x or z has identically
        # ZERO coefficients (the graded coefficient reaches one cell from
        # the surface), so its coef row keeps the dataset fill value with
        # no classification work; its dwall row still uses the EXACT igl
        # distance, but batched over MANY leaves per query (the per-leaf
        # cost was the AABB rebuild + call overhead, not the query). The
        # span dim y always overlaps (the extrusion crosses the whole
        # span). Deep-refinement far-field layouts (tutorials/naca B11:
        # >90% far leaves) drop from hours to minutes.
        far = np.zeros(n_leaves, dtype=bool)
        if n_leaves:
            bounds = mesh.bounds
            for i in range(n_leaves):
                la = level_args[lev[i]]
                o = crd[i]*nb
                ext = np.zeros((3, 2))
                for d in range(3):
                    nodes = grid_axis_nodes(d + 1, la)
                    n_d = nodes.size - 1
                    lo = nodes[min(max(int(o[d]) - 1, 0), n_d)]
                    hi = nodes[min(int(o[d]) + nb + 2, n_d)]
                    h_d = (nodes[-1] - nodes[0])/n_d
                    ext[d] = (lo - 3.0*h_d, hi + 3.0*h_d)
                far[i] = (ext[0, 1] < bounds[0, 0] or ext[0, 0] > bounds[1, 0] or
                          ext[2, 1] < bounds[0, 2] or ext[2, 0] > bounds[1, 2])
        near_ids = np.nonzero(~far)[0]
        print(f"tiles: {near_ids.size} near-body leaves through the STL machinery, "
              f"{int(far.sum())} far-field leaves shortcut")

        if want_dwall and far.any():
            igl, _, _ = require_stl_tools()
            idx = np.arange(nb + 2, dtype=np.int64)
            far_ids = np.nonzero(far)[0]
            w3 = (nb + 2)**3
            chunk = max(1, 400)
            for c0 in range(0, far_ids.size, chunk):
                ids = far_ids[c0:c0 + chunk]
                pts = np.empty((ids.size*w3, 3))
                for j, i in enumerate(ids):
                    la = level_args[lev[i]]
                    o = crd[i]*nb
                    ii, jj, kk = np.meshgrid(idx + int(o[0]), idx + int(o[1]),
                                             idx + int(o[2]), indexing="ij")
                    pts[j*w3:(j+1)*w3] = stl_points_for_indices(
                        ii.ravel(), jj.ravel(), kk.ravel(), 0, la)
                sq, _, _ = igl.point_mesh_squared_distance(pts, vertices, faces)
                d = np.sqrt(sq)
                for j, i in enumerate(ids):
                    dwall[i] = d[j*w3:(j+1)*w3].reshape(nb+2, nb+2, nb+2)

        # Near-body tiles: the exact per-leaf machinery, optionally over a
        # worker pool (--jobs; the legacy chunked path's pattern).
        jobs = max(1, int(getattr(args, "jobs", 1)))
        items = [(int(i), int(lev[i]), (int(crd[i][0]*nb), int(crd[i][1]*nb),
                                        int(crd[i][2]*nb))) for i in near_ids]
        if jobs == 1:
            init_block_tile_worker(mesh, vertices, faces, level_args, nb, want_dwall)
            for item in items:
                i, coef_tile, dw = block_tile_worker(item)
                coef[i] = coef_tile
                if dw is not None:
                    dwall[i] = dw
        else:
            try:
                context = multiprocessing.get_context("fork")
            except ValueError:
                context = multiprocessing.get_context()
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=jobs,
                mp_context=context,
                initializer=init_block_tile_worker,
                initargs=(mesh, vertices, faces, level_args, nb, want_dwall),
            ) as pool:
                futures = [pool.submit(block_tile_worker, item) for item in items]
                for future in concurrent.futures.as_completed(futures):
                    i, coef_tile, dw = future.result()
                    coef[i] = coef_tile
                    if dw is not None:
                        dwall[i] = dw
    print(f"block-table coefficients written to: {output}")


def block_active_from_stl(args: argparse.Namespace) -> None:
    """Write the per-block keep flags used by the solver's solid-block removal.

    A block is removable iff the block dilated by one halo cell is solid at
    cell centres and all three staggered locations. The flags go into the
    coefficient file as the 1D dataset "block_active" (x-fastest lattice
    raster order, 1 = keep) plus a block_nb attribute.
    """
    finalize_grid_args(args)
    nb = int(args.block_nb)
    if nb < 4 or nb % 2:
        raise SystemExit("--block-nb must be even and >= 4")
    for n, name in ((args.nx, "nx"), (args.ny, "ny"), (args.nz, "nz")):
        if n % nb:
            raise SystemExit(f"--block-nb {nb} does not divide {name} = {n}")
    gnbt = (args.nx // nb, args.ny // nb, args.nz // nb)

    geometry_values = args.geometry if isinstance(args.geometry, list) else [args.geometry]
    geometries = [Path(value).resolve() for value in geometry_values]
    output = Path(args.output).resolve()
    mesh, vertices, faces = load_stl_meshes(geometries, repair=not args.no_repair)
    mesh, vertices, faces = apply_stl_transform(mesh, args)
    prepare_fluid_points(args)
    for name in ("_fluid_ray_checks", "_fluid_ray_ambiguous", "_fluid_ray_disagreements", "_fluid_ray_overrides"):
        setattr(args, name, 0)
    args._ray_intersector = None

    removable = np.ones(gnbt, dtype=bool)
    for var in (VAR_U, VAR_V, VAR_W, 0):  # 0 = cell centres (never face-staggered)
        inside_ext, _, _, _ = classify_stl_extended_grid(mesh, vertices, faces, var, args)
        solid_ext = ~inside_ext if getattr(args, "inside_is_fluid", False) else inside_ext
        removable &= window_all_solid(solid_ext, nb, gnbt)

    active = (~removable).astype(np.int32)
    with h5py.File(output, "r+") as h5:
        if "block_active" in h5:
            del h5["block_active"]
        # x-fastest raster order, matching the solver's lattice indexing
        h5.create_dataset("block_active", data=active.transpose(2, 1, 0).ravel())
        h5.attrs["block_nb"] = np.int32(nb)
    print(f"block lattice: {gnbt[0]} x {gnbt[1]} x {gnbt[2]} (nb = {nb})")
    print(f"removable blocks: {int(removable.sum())} of {int(removable.size)}")
    print(f"block_active written to: {output}")


def make_sphere_stl(args: argparse.Namespace) -> None:
    center = np.asarray(args.center, dtype=np.float64)
    mesh = make_stl_sphere_mesh(center, args.radius, args.subdivisions)
    remove_stl_cap(mesh, center, args.radius, args.open_cap_angle)
    Path(args.output).resolve().parent.mkdir(parents=True, exist_ok=True)
    mesh.export(Path(args.output).resolve())


def test_stl_sphere(args: argparse.Namespace) -> None:
    """STL sphere regression against the analytic sphere reference field."""
    finalize_grid_args(args)
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    geometry = workdir / "sphere.stl"
    h5 = workdir / "sphere_stl_coeff.h5"
    center = np.asarray(args.center, dtype=np.float64)

    make_args = argparse.Namespace(
        output=geometry, center=center, radius=args.radius,
        subdivisions=args.subdivisions, open_cap_angle=args.open_cap_angle,
    )
    make_sphere_stl(make_args)

    fluid_points_path = workdir / "sphere_fluid_points.txt"
    coeff_args = argparse.Namespace(**vars(args))
    coeff_args.geometry = geometry
    coeff_args.output = h5
    prepare_fluid_points(coeff_args, center=center, radius=args.radius, write_path=fluid_points_path)
    stl = coeff_from_stl(coeff_args)
    analytic = analytic_sphere_coeff(args, center, args.radius)

    diff = stl - analytic
    max_abs = float(np.max(np.abs(diff)))
    l2 = float(np.sqrt(np.mean(diff * diff)))
    active = np.abs(analytic) < 1.0e20
    local_rel = np.zeros_like(diff)
    local_rel[active] = np.abs(diff[active]) / np.maximum(1.0, np.abs(analytic[active]))
    max_rel_active = float(np.max(local_rel[active])) if np.any(active) else 0.0
    solid_mismatch = int(np.count_nonzero((np.abs(stl) > 1.0e20) != (np.abs(analytic) > 1.0e20)))

    mesh, _, _ = load_stl_mesh(geometry, repair=not args.no_repair)
    radial = np.linalg.norm(mesh.triangles_center - center, axis=1)
    facet_radial_error = float(np.max(np.abs(radial - args.radius))) if len(radial) else 0.0

    print(f"Analytic reference radius: {args.radius:.17e}")
    print(f"max facet-center radial error: {facet_radial_error:.17e}")
    print(f"max_abs_vs_analytic_sphere: {max_abs:.17e}")
    print(f"l2_vs_analytic_sphere:      {l2:.17e}")
    print(f"max_rel_active_vs_analytic: {max_rel_active:.17e}")
    print(f"solid_mismatch_count:       {solid_mismatch}")
    if solid_mismatch > args.max_solid_mismatch:
        raise SystemExit(
            f"STL sphere test failed: {solid_mismatch} solid/fluid mismatches > {args.max_solid_mismatch}"
        )
    if max_abs > args.abs_tolerance and max_rel_active > args.rel_tolerance:
        raise SystemExit(
            f"STL sphere test failed: max_abs {max_abs:.3e} > {args.abs_tolerance:.3e} "
            f"and max_rel_active {max_rel_active:.3e} > {args.rel_tolerance:.3e}"
        )



def make_stl_sphere_mesh(center: np.ndarray, radius: float, subdivisions: int):
    _, trimesh, _ = require_stl_tools()
    mesh = trimesh.creation.icosphere(subdivisions=subdivisions, radius=radius)
    mesh.apply_translation(np.asarray(center, dtype=np.float64))
    return mesh


def remove_stl_cap(mesh, center: np.ndarray, radius: float, angle_deg: float) -> None:
    if angle_deg <= 0.0:
        return
    cutoff = center[2] + radius * math.cos(math.radians(angle_deg))
    keep = mesh.triangles_center[:, 2] < cutoff
    mesh.update_faces(keep)
    mesh.remove_unreferenced_vertices()


def mesh_edge_counts(mesh) -> tuple[int, int]:
    if len(mesh.faces) == 0:
        return 0, 0
    edges = np.sort(mesh.faces[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    return int(np.count_nonzero(counts == 1)), int(np.count_nonzero(counts > 2))


def choose_faces(rng: np.random.Generator, nfaces: int, fraction: float) -> np.ndarray:
    count = max(1, int(round(fraction * nfaces)))
    count = min(count, nfaces)
    return rng.choice(nfaces, size=count, replace=False)


def create_stl_defect(case: str, center: np.ndarray, radius: float, subdivisions: int,
                      rng: np.random.Generator):
    """Create controlled mesh defects for robustness tests without relying on external files."""
    _, trimesh, _ = require_stl_tools()
    mesh = make_stl_sphere_mesh(center, radius, subdivisions)
    description = "clean watertight reference"

    if case == "clean":
        return mesh, description

    if case.startswith("cap_"):
        angle = float(case.split("_")[1])
        remove_stl_cap(mesh, center, radius, angle)
        return mesh, f"removed polar cap, half-angle {angle:g} degrees"

    if case.startswith("delete_"):
        fraction = float(case.split("_")[1])
        drop = choose_faces(rng, len(mesh.faces), fraction)
        keep = np.ones(len(mesh.faces), dtype=bool)
        keep[drop] = False
        mesh.update_faces(keep)
        mesh.remove_unreferenced_vertices()
        return mesh, f"deleted {fraction:.3g} of triangles randomly"

    if case.startswith("flip_"):
        fraction = float(case.split("_")[1])
        flip = choose_faces(rng, len(mesh.faces), fraction)
        faces = np.asarray(mesh.faces).copy()
        faces[flip] = faces[flip][:, [0, 2, 1]]
        mesh = trimesh.Trimesh(vertices=np.asarray(mesh.vertices).copy(), faces=faces, process=False)
        return mesh, f"flipped orientation of {fraction:.3g} of triangles"

    if case.startswith("duplicate_"):
        fraction = float(case.split("_")[1])
        dup = choose_faces(rng, len(mesh.faces), fraction)
        faces = np.vstack((np.asarray(mesh.faces), np.asarray(mesh.faces)[dup]))
        mesh = trimesh.Trimesh(vertices=np.asarray(mesh.vertices).copy(), faces=faces, process=False)
        return mesh, f"added duplicate copies of {fraction:.3g} of triangles"

    if case.startswith("reversed_duplicate_"):
        fraction = float(case.split("_")[2])
        dup = choose_faces(rng, len(mesh.faces), fraction)
        reversed_faces = np.asarray(mesh.faces)[dup][:, [0, 2, 1]]
        faces = np.vstack((np.asarray(mesh.faces), reversed_faces))
        mesh = trimesh.Trimesh(vertices=np.asarray(mesh.vertices).copy(), faces=faces, process=False)
        return mesh, f"added reversed duplicate copies of {fraction:.3g} of triangles"

    if case.startswith("crack_band_"):
        width = float(case.split("_")[2])
        vertices = np.asarray(mesh.vertices).copy()
        faces = np.asarray(mesh.faces).copy()
        centers = np.asarray(mesh.triangles_center)
        band = (np.abs(centers[:, 2] - center[2]) < 0.08 * radius) & (centers[:, 0] > center[0])
        selected = np.flatnonzero(band)
        if selected.size == 0:
            return mesh, "crack band requested, but no faces were selected"
        old_vertices = np.unique(faces[selected].ravel())
        remap = np.full(vertices.shape[0], -1, dtype=np.int64)
        remap[old_vertices] = np.arange(vertices.shape[0], vertices.shape[0] + old_vertices.size)
        duplicated = vertices[old_vertices].copy()
        radial = duplicated - center
        norm = np.linalg.norm(radial, axis=1)
        radial[norm > 0.0] /= norm[norm > 0.0, None]
        duplicated += width * radius * radial
        faces[selected] = remap[faces[selected]]
        mesh = trimesh.Trimesh(vertices=np.vstack((vertices, duplicated)), faces=faces, process=False)
        return mesh, f"detached an equatorial patch by {width:g} radii"

    if case.startswith("degenerate_"):
        fraction = float(case.split("_")[1])
        bad = choose_faces(rng, len(mesh.faces), fraction)
        faces = np.asarray(mesh.faces).copy()
        faces[bad, 1] = faces[bad, 0]
        mesh = trimesh.Trimesh(vertices=np.asarray(mesh.vertices).copy(), faces=faces, process=False)
        return mesh, f"collapsed {fraction:.3g} of triangles to zero area"

    if case == "floating_patch":
        vertices = np.asarray(mesh.vertices).copy()
        faces = np.asarray(mesh.faces).copy()
        patch = np.array([
            center + np.array([0.32, -0.015, -0.015]),
            center + np.array([0.32,  0.015, -0.015]),
            center + np.array([0.32,  0.000,  0.020]),
        ], dtype=np.float64)
        patch_face = np.array([[vertices.shape[0], vertices.shape[0] + 1, vertices.shape[0] + 2]], dtype=np.int64)
        mesh = trimesh.Trimesh(vertices=np.vstack((vertices, patch)), faces=np.vstack((faces, patch_face)), process=False)
        return mesh, "added one disconnected open triangle patch"

    raise ValueError(f"unknown STL stress case {case!r}")


def stl_error_metrics(coef: np.ndarray, analytic: np.ndarray) -> dict[str, float | int]:
    diff = coef - analytic
    active = np.abs(analytic) < 1.0e20
    local_rel = np.zeros_like(diff)
    local_rel[active] = np.abs(diff[active]) / np.maximum(1.0, np.abs(analytic[active]))
    return {
        "max_abs": float(np.max(np.abs(diff))),
        "l2": float(np.sqrt(np.mean(diff * diff))),
        "max_rel_active": float(np.max(local_rel[active])) if np.any(active) else 0.0,
        "solid_mismatch": int(np.count_nonzero((np.abs(coef) > 1.0e20) != (np.abs(analytic) > 1.0e20))),
    }


def stress_stl_watertightness(args: argparse.Namespace) -> None:
    """Exercise STL classification on clean and intentionally damaged sphere meshes."""
    finalize_grid_args(args)
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    center = np.asarray(args.center, dtype=np.float64)
    rng = np.random.default_rng(args.seed)
    analytic = analytic_sphere_coeff(args, center, args.radius)

    default_cases = [
        "clean",
        "cap_5", "cap_10", "cap_15",
        "delete_0.001", "delete_0.01",
        "flip_0.01", "flip_0.05",
        "duplicate_0.05", "reversed_duplicate_0.02",
        "crack_band_0.01",
        "degenerate_0.002",
        "floating_patch",
    ]
    cases = args.cases if args.cases else default_cases
    rows = []

    print(f"STL watertightness stress test: {len(cases)} cases")
    print(f"workdir: {workdir}")
    for case in cases:
        print(f"\n[{case}] generating mesh...")
        source_mesh, description = create_stl_defect(case, center, args.radius, args.subdivisions, rng)
        stl_path = workdir / f"{case}.stl"
        h5_path = workdir / f"{case}_coeff.h5"
        source_mesh.export(stl_path)

        mesh, vertices, faces = load_stl_mesh(stl_path, repair=not args.no_repair)
        boundary_edges, nonmanifold_edges = mesh_edge_counts(mesh)
        coeff_args = argparse.Namespace(**vars(args))
        coeff_args.geometry = stl_path
        coeff_args.output = h5_path
        prepare_fluid_points(coeff_args, center=center, radius=args.radius, write_path=workdir / "stress_fluid_points.txt")
        coef, stats = stl_coeff_from_mesh(mesh, vertices, faces, coeff_args)
        metrics = stl_error_metrics(coef, analytic)
        if args.write_coefficients:
            write_hdf5(h5_path, coef, coeff_args, stl_path, stl_extra_attrs(mesh, stats, coeff_args))

        status = "ok"
        if metrics["solid_mismatch"] > args.max_solid_mismatch:
            status = "fail"
        elif metrics["max_abs"] > args.abs_tolerance or metrics["max_rel_active"] > args.rel_tolerance:
            status = "warn"

        row = {
            "case": case,
            "status": status,
            "description": description,
            "stl_file": str(stl_path),
            "h5_file": str(h5_path) if args.write_coefficients else "",
            "vertices": int(len(mesh.vertices)),
            "faces": int(len(mesh.faces)),
            "watertight": bool(mesh.is_watertight),
            "euler_number": int(mesh.euler_number),
            "boundary_edges": boundary_edges,
            "nonmanifold_edges": nonmanifold_edges,
            "ambiguous_points": int(stats["ambiguous_points"]),
            "crossing_segments": int(stats["crossing_segments"]),
            "fallback_segments": int(stats["fallback_segments"]),
            "fluid_ray_checks": int(stats.get("fluid_ray_checks", 0)),
            "fluid_ray_ambiguous": int(stats.get("fluid_ray_ambiguous", 0)),
            "fluid_ray_disagreements": int(stats.get("fluid_ray_disagreements", 0)),
            "fluid_ray_overrides": int(stats.get("fluid_ray_overrides", 0)),
            **metrics,
        }
        rows.append(row)
        print(
            f"status={status} watertight={row['watertight']} boundary_edges={boundary_edges} "
            f"nonmanifold_edges={nonmanifold_edges} fallback={row['fallback_segments']} "
            f"fluid_overrides={row['fluid_ray_overrides']} solid_mismatch={row['solid_mismatch']} "
            f"max_rel={row['max_rel_active']:.3e}"
        )

    csv_path = workdir / "stl_watertightness_stress.csv"
    json_path = workdir / "stl_watertightness_stress.json"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with json_path.open("w") as f:
        json.dump(rows, f, indent=2)

    print("\nSummary")
    print("case,status,watertight,boundary_edges,nonmanifold_edges,fallback_segments,fluid_ray_disagreements,fluid_ray_overrides,solid_mismatch,max_rel_active")
    for row in rows:
        print(
            f"{row['case']},{row['status']},{int(row['watertight'])},{row['boundary_edges']},"
            f"{row['nonmanifold_edges']},{row['fallback_segments']},{row['fluid_ray_disagreements']},"
            f"{row['fluid_ray_overrides']},{row['solid_mismatch']},{row['max_rel_active']:.6e}"
        )
    print(f"\nCSV report:  {csv_path}")
    print(f"JSON report: {json_path}")

    if args.strict and any(row["status"] != "ok" for row in rows):
        raise SystemExit("STL watertightness stress test failed in strict mode")



def bent_pipe_basis() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fixed oblique basis used to avoid axis-aligned special cases in pipe tests."""
    ex = np.array([0.83, 0.30, 0.47], dtype=np.float64)
    ex /= np.linalg.norm(ex)
    ey = np.array([-0.18, 0.94, 0.29], dtype=np.float64)
    ey -= np.dot(ey, ex) * ex
    ey /= np.linalg.norm(ey)
    ez = np.cross(ex, ey)
    ez /= np.linalg.norm(ez)
    return ex, ey, ez


def bent_pipe_path_samples(straight_length: float, bend_radius: float,
                           straight_sections: int, bend_sections: int) -> tuple[np.ndarray, np.ndarray]:
    """Sample the local centerline and tangents of the bent-pipe validation geometry."""
    centers = []
    tangents = []
    for s in np.linspace(0.0, 1.0, straight_sections + 1):
        centers.append([-straight_length + straight_length * s, 0.0, 0.0])
        tangents.append([1.0, 0.0, 0.0])
    for theta in np.linspace(-0.5 * math.pi, 0.0, bend_sections + 1)[1:]:
        centers.append([bend_radius * math.cos(theta), bend_radius + bend_radius * math.sin(theta), 0.0])
        tangents.append([-math.sin(theta), math.cos(theta), 0.0])
    for s in np.linspace(0.0, 1.0, straight_sections + 1)[1:]:
        centers.append([bend_radius, bend_radius + straight_length * s, 0.0])
        tangents.append([0.0, 1.0, 0.0])
    centers = np.asarray(centers, dtype=np.float64)
    tangents = np.asarray(tangents, dtype=np.float64)
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]
    return centers, tangents


def bent_pipe_total_length(args: argparse.Namespace) -> float:
    return 2.0 * args.straight_length + 0.5 * math.pi * args.bend_radius


def bent_pipe_sample_arclengths(args: argparse.Namespace) -> np.ndarray:
    values = []
    values.extend(np.linspace(0.0, args.straight_length, args.straight_sections + 1))
    theta = np.linspace(-0.5 * math.pi, 0.0, args.bend_sections + 1)[1:]
    values.extend(args.straight_length + args.bend_radius * (theta + 0.5 * math.pi))
    tail = np.linspace(0.0, args.straight_length, args.straight_sections + 1)[1:]
    values.extend(args.straight_length + 0.5 * math.pi * args.bend_radius + tail)
    return np.asarray(values, dtype=np.float64)


def bent_pipe_inner_radius_profile(s: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    """Return the pipe bore radius, optionally including a smooth physical restriction."""
    radius = np.full_like(np.asarray(s, dtype=np.float64), args.inner_radius, dtype=np.float64)
    throat = getattr(args, "restriction_radius", 0.0)
    width = getattr(args, "restriction_width", 0.0)
    if throat <= 0.0 or width <= 0.0:
        return radius
    center = getattr(args, "restriction_center", 0.5) * bent_pipe_total_length(args)
    half_width = 0.5 * width
    d = np.abs(np.asarray(s, dtype=np.float64) - center)
    active = d < half_width
    blend = np.ones_like(radius)
    blend[active] = 0.5 * (1.0 - np.cos(math.pi * d[active] / half_width))
    radius = throat + (args.inner_radius - throat) * blend
    return radius


def bent_pipe_center_at_s(s: float, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    if s <= args.straight_length:
        center = np.array([-args.straight_length + s, 0.0, 0.0], dtype=np.float64)
        tangent = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    elif s <= args.straight_length + 0.5 * math.pi * args.bend_radius:
        theta = (s - args.straight_length) / args.bend_radius - 0.5 * math.pi
        center = np.array([args.bend_radius * math.cos(theta),
                           args.bend_radius + args.bend_radius * math.sin(theta), 0.0], dtype=np.float64)
        tangent = np.array([-math.sin(theta), math.cos(theta), 0.0], dtype=np.float64)
    else:
        y = args.bend_radius + s - args.straight_length - 0.5 * math.pi * args.bend_radius
        center = np.array([args.bend_radius, y, 0.0], dtype=np.float64)
        tangent = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    tangent /= np.linalg.norm(tangent)
    return center, tangent


def bent_pipe_global_center_at_s(s: float, args: argparse.Namespace) -> np.ndarray:
    center, _ = bent_pipe_center_at_s(s, args)
    basis, translation = bent_pipe_transform(args)
    return center @ basis.T + translation


def bent_pipe_local_mesh_arrays(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    """Build a triangulated finite pipe wall in local coordinates."""
    centers, tangents = bent_pipe_path_samples(
        args.straight_length, args.bend_radius, args.straight_sections, args.bend_sections
    )
    n_theta = int(args.around)
    angles = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False)
    plane_normal = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    vertices = []
    inner_radii = bent_pipe_inner_radius_profile(bent_pipe_sample_arclengths(args), args)
    for radius_values in (np.full(centers.shape[0], args.outer_radius), inner_radii):
        for c, t, radius in zip(centers, tangents, radius_values):
            n = np.cross(plane_normal, t)
            n /= np.linalg.norm(n)
            b = plane_normal
            for a in angles:
                vertices.append(c + radius * (math.cos(a) * n + math.sin(a) * b))
    vertices = np.asarray(vertices, dtype=np.float64)
    n_ring = centers.shape[0]

    def outer(i: int, m: int) -> int:
        return i * n_theta + (m % n_theta)

    def inner(i: int, m: int) -> int:
        return n_ring * n_theta + i * n_theta + (m % n_theta)

    faces = []
    for i in range(n_ring - 1):
        for m in range(n_theta):
            m1 = m + 1
            faces.append([outer(i, m), outer(i + 1, m), outer(i + 1, m1)])
            faces.append([outer(i, m), outer(i + 1, m1), outer(i, m1)])
            faces.append([inner(i, m1), inner(i + 1, m1), inner(i + 1, m)])
            faces.append([inner(i, m1), inner(i + 1, m), inner(i, m)])
    for m in range(n_theta):
        m1 = m + 1
        faces.append([outer(0, m1), outer(0, m), inner(0, m)])
        faces.append([outer(0, m1), inner(0, m), inner(0, m1)])
        faces.append([outer(n_ring - 1, m), outer(n_ring - 1, m1), inner(n_ring - 1, m)])
        faces.append([outer(n_ring - 1, m1), inner(n_ring - 1, m1), inner(n_ring - 1, m)])
    return vertices, np.asarray(faces, dtype=np.int64)


def bent_pipe_transform(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    vertices, _ = bent_pipe_local_mesh_arrays(args)
    ex, ey, ez = bent_pipe_basis()
    basis = np.column_stack((ex, ey, ez))
    global_vertices = vertices @ basis.T
    target_center = np.asarray(args.pipe_center, dtype=np.float64)
    translation = target_center - 0.5 * (global_vertices.min(axis=0) + global_vertices.max(axis=0))
    return basis, translation


def make_bent_pipe_mesh(args: argparse.Namespace):
    _, trimesh, _ = require_stl_tools()
    vertices, faces = bent_pipe_local_mesh_arrays(args)
    basis, translation = bent_pipe_transform(args)
    vertices = vertices @ basis.T + translation
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
    try:
        trimesh.repair.fix_normals(mesh, multibody=False)
    except Exception:
        pass
    return mesh


def make_bent_pipe_stl(args: argparse.Namespace) -> None:
    mesh = make_bent_pipe_mesh(args)
    Path(args.output).resolve().parent.mkdir(parents=True, exist_ok=True)
    mesh.export(Path(args.output).resolve())


def bent_pipe_local_points(points: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    basis, translation = bent_pipe_transform(args)
    return (np.asarray(points, dtype=np.float64) - translation) @ basis


def bent_pipe_radial_distance_and_s(points: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project points onto the bent-pipe centerline and return radius and arclength."""
    q = bent_pipe_local_points(points, args)
    x = q[:, 0]
    y = q[:, 1]
    z = q[:, 2]
    best = np.full(points.shape[0], np.inf, dtype=np.float64)
    best_s = np.zeros(points.shape[0], dtype=np.float64)
    valid = np.zeros(points.shape[0], dtype=bool)

    mask = (x >= -args.straight_length) & (x <= 0.0)
    d = np.sqrt(y * y + z * z)
    update = mask & (d < best)
    best[update] = d[update]
    best_s[update] = x[update] + args.straight_length
    valid |= mask

    dx = x
    dy = y - args.bend_radius
    theta = np.arctan2(dy, dx)
    mask = (theta >= -0.5 * math.pi) & (theta <= 0.0)
    d = np.sqrt((np.sqrt(dx * dx + dy * dy) - args.bend_radius) ** 2 + z * z)
    update = mask & (d < best)
    best[update] = d[update]
    best_s[update] = args.straight_length + args.bend_radius * (theta[update] + 0.5 * math.pi)
    valid |= mask

    mask = (y >= args.bend_radius) & (y <= args.bend_radius + args.straight_length)
    d = np.sqrt((x - args.bend_radius) ** 2 + z * z)
    update = mask & (d < best)
    best[update] = d[update]
    best_s[update] = args.straight_length + 0.5 * math.pi * args.bend_radius + (y[update] - args.bend_radius)
    valid |= mask
    return best, valid, best_s


def bent_pipe_radial_distance(points: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    distance, valid, _ = bent_pipe_radial_distance_and_s(points, args)
    return distance, valid


def analytic_bent_pipe_wall_mask(points: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    """Analytic solid mask for the finite pipe-wall validation case."""
    distance, valid = bent_pipe_radial_distance(points, args)
    solid = valid & (distance > args.inner_radius) & (distance < args.outer_radius)
    tol = args.reference_tolerance
    certain = ~valid | ((np.abs(distance - args.inner_radius) > tol) & (np.abs(distance - args.outer_radius) > tol))

    q = bent_pipe_local_points(points, args)
    start_r = np.sqrt(q[:, 1] * q[:, 1] + q[:, 2] * q[:, 2])
    end_r = np.sqrt((q[:, 0] - args.bend_radius) ** 2 + q[:, 2] * q[:, 2])
    near_start_annulus = (np.abs(q[:, 0] + args.straight_length) <= tol) & \
        (start_r >= args.inner_radius - tol) & (start_r <= args.outer_radius + tol)
    near_end_annulus = (np.abs(q[:, 1] - (args.bend_radius + args.straight_length)) <= tol) & \
        (end_r >= args.inner_radius - tol) & (end_r <= args.outer_radius + tol)
    certain &= ~(near_start_annulus | near_end_annulus)
    return solid, certain


def bent_pipe_reference_mask(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    masks = []
    certains = []
    for var in (VAR_U, VAR_V, VAR_W):
        points, shape = stl_extended_grid_points(var, args)
        solid, certain = analytic_bent_pipe_wall_mask(points, args)
        masks.append(solid.reshape(shape)[1:args.nx + 3, 1:args.ny + 3, 1:args.nz + 3])
        certains.append(certain.reshape(shape)[1:args.nx + 3, 1:args.ny + 3, 1:args.nz + 3])
    return np.stack(masks, axis=-1), np.stack(certains, axis=-1)


def bent_pipe_validation_fluid_points(args: argparse.Namespace) -> np.ndarray:
    """Create generic fluid probes in the pipe bore and on the box boundary."""
    centers, tangents = bent_pipe_path_samples(args.straight_length, args.bend_radius, 8, 12)
    ex, ey, ez = bent_pipe_basis()
    basis, translation = bent_pipe_transform(args)
    bore = centers[::max(1, centers.shape[0] // 12)] @ basis.T + translation
    box = []
    for x in (0.0, args.lx):
        for y in (0.0, args.ly):
            for z in (0.0, args.lz):
                box.append([x, y, z])
    return np.ascontiguousarray(np.vstack((bore, np.asarray(box, dtype=np.float64))), dtype=np.float64)


def bent_pipe_inner_surface_arrays(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    """Build only the pipe inner surface used for closed-cavity STL tests."""
    centers, tangents = bent_pipe_path_samples(
        args.straight_length, args.bend_radius, args.straight_sections, args.bend_sections
    )
    n_theta = int(args.around)
    angles = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False)
    plane_normal = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    vertices = []
    radii = bent_pipe_inner_radius_profile(bent_pipe_sample_arclengths(args), args)
    for c, t, radius in zip(centers, tangents, radii):
        n = np.cross(plane_normal, t)
        n /= np.linalg.norm(n)
        b = plane_normal
        for a in angles:
            vertices.append(c + radius * (math.cos(a) * n + math.sin(a) * b))
    vertices = np.asarray(vertices, dtype=np.float64)

    def idx(i: int, m: int) -> int:
        return i * n_theta + (m % n_theta)

    faces = []
    for i in range(centers.shape[0] - 1):
        for m in range(n_theta):
            m1 = m + 1
            faces.append([idx(i, m1), idx(i + 1, m1), idx(i + 1, m)])
            faces.append([idx(i, m1), idx(i + 1, m), idx(i, m)])
    return vertices, np.asarray(faces, dtype=np.int64)


def bent_pipe_lid_arrays(args: argparse.Namespace, side: str) -> tuple[np.ndarray, np.ndarray]:
    """Build one disk lid that can exactly or approximately close the pipe bore."""
    if side == "start":
        center = np.array([-args.straight_length, 0.0, 0.0], dtype=np.float64)
        tangent = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        normal_sign = -1
    elif side == "end":
        center = np.array([args.bend_radius, args.bend_radius + args.straight_length, 0.0], dtype=np.float64)
        tangent = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        normal_sign = 1
    else:
        raise ValueError(f"unknown pipe lid side {side!r}")
    if args.lid_gap != 0.0:
        center = center + args.lid_gap * tangent * normal_sign
    n_theta = int(args.around)
    angles = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False)
    plane_normal = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    n = np.cross(plane_normal, tangent)
    n /= np.linalg.norm(n)
    b = plane_normal
    radius = args.inner_radius * args.lid_radius_scale
    vertices = [center]
    for a in angles:
        vertices.append(center + radius * (math.cos(a) * n + math.sin(a) * b))
    faces = []
    for m in range(n_theta):
        m1 = 1 + ((m + 1) % n_theta)
        m0 = 1 + m
        if side == "start":
            faces.append([0, m1, m0])
        else:
            faces.append([0, m0, m1])
    return np.asarray(vertices, dtype=np.float64), np.asarray(faces, dtype=np.int64)


def transformed_mesh_from_arrays(vertices: np.ndarray, faces: np.ndarray, args: argparse.Namespace):
    _, trimesh, _ = require_stl_tools()
    basis, translation = bent_pipe_transform(args)
    vertices = vertices @ basis.T + translation
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
    try:
        trimesh.repair.fix_normals(mesh, multibody=False)
    except Exception:
        pass
    return mesh


def make_bent_pipe_inner_mesh(args: argparse.Namespace):
    vertices, faces = bent_pipe_inner_surface_arrays(args)
    return transformed_mesh_from_arrays(vertices, faces, args)


def make_bent_pipe_lid_mesh(args: argparse.Namespace, side: str):
    vertices, faces = bent_pipe_lid_arrays(args, side)
    return transformed_mesh_from_arrays(vertices, faces, args)


def make_lidded_bent_pipe_stls(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    """Write the inner surface and two lid STLs as separate generic input meshes."""
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    inner = output_dir / "bent_pipe_inner_surface.stl"
    start_lid = output_dir / "bent_pipe_start_lid.stl"
    end_lid = output_dir / "bent_pipe_end_lid.stl"
    make_bent_pipe_inner_mesh(args).export(inner)
    make_bent_pipe_lid_mesh(args, "start").export(start_lid)
    make_bent_pipe_lid_mesh(args, "end").export(end_lid)
    return inner, start_lid, end_lid


def make_lidded_bent_pipe_stl(args: argparse.Namespace) -> None:
    inner, start_lid, end_lid = make_lidded_bent_pipe_stls(args)
    print(f"inner surface STL: {inner}")
    print(f"start lid STL:     {start_lid}")
    print(f"end lid STL:       {end_lid}")


def analytic_bent_pipe_cavity_mask(points: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    """Analytic fluid mask for the lidded-pipe cavity validation case."""
    distance, valid, path_s = bent_pipe_radial_distance_and_s(points, args)
    local_inner_radius = bent_pipe_inner_radius_profile(path_s, args)
    fluid = valid & (distance < local_inner_radius)
    tol = args.reference_tolerance
    q = bent_pipe_local_points(points, args)
    near_inner = valid & (np.abs(distance - local_inner_radius) <= tol)
    start_r = np.sqrt(q[:, 1] * q[:, 1] + q[:, 2] * q[:, 2])
    end_r = np.sqrt((q[:, 0] - args.bend_radius) ** 2 + q[:, 2] * q[:, 2])
    near_start_lid = (np.abs(q[:, 0] + args.straight_length) <= tol) & (start_r <= args.inner_radius + tol)
    near_end_lid = (np.abs(q[:, 1] - (args.bend_radius + args.straight_length)) <= tol) & (end_r <= args.inner_radius + tol)
    certain = ~(near_inner | near_start_lid | near_end_lid)
    return fluid, certain


def bent_pipe_cavity_reference_mask(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray]:
    masks = []
    certains = []
    for var in (VAR_U, VAR_V, VAR_W):
        points, shape = stl_extended_grid_points(var, args)
        fluid, certain = analytic_bent_pipe_cavity_mask(points, args)
        masks.append((~fluid).reshape(shape)[1:args.nx + 3, 1:args.ny + 3, 1:args.nz + 3])
        certains.append(certain.reshape(shape)[1:args.nx + 3, 1:args.ny + 3, 1:args.nz + 3])
    return np.stack(masks, axis=-1), np.stack(certains, axis=-1)


def test_stl_lidded_bent_pipe(args: argparse.Namespace) -> None:
    """Validate inside-is-fluid STL classification for a closed lidded pipe cavity."""
    finalize_grid_args(args)
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    args.output_dir = str(workdir)
    inner, start_lid, end_lid = make_lidded_bent_pipe_stls(args)
    h5 = workdir / "lidded_bent_pipe_coeff.h5"

    mesh, vertices, faces = load_stl_meshes([inner, start_lid, end_lid], repair=not args.no_repair)
    boundary_edges, nonmanifold_edges = mesh_edge_counts(mesh)
    coeff_args = argparse.Namespace(**vars(args))
    coeff_args.geometry = [str(inner), str(start_lid), str(end_lid)]
    coeff_args.output = h5
    coeff_args.inside_is_fluid = True
    coef, stats = stl_coeff_from_mesh(mesh, vertices, faces, coeff_args)
    write_hdf5(h5, coef, coeff_args, ";".join(str(p) for p in (inner, start_lid, end_lid)), stl_extra_attrs(mesh, stats, coeff_args))

    analytic_solid, certain = bent_pipe_cavity_reference_mask(args)
    stl_solid = np.abs(coef) > 1.0e20
    mismatch = stl_solid != analytic_solid
    solid_mismatch = int(np.count_nonzero(mismatch & certain))
    ignored_mismatch = int(np.count_nonzero(mismatch & ~certain))
    probe_points = bent_pipe_validation_fluid_points(args)
    probe_fluid, _ = analytic_bent_pipe_cavity_mask(probe_points, args)

    print(f"Lidded bent pipe STLs: {inner}; {start_lid}; {end_lid}")
    print(f"HDF5 coefficients: {h5}")
    print(f"XDMF coefficients: {coefficient_xdmf_path(h5)}")
    print(f"vertices/faces: {len(mesh.vertices)} / {len(mesh.faces)}")
    print(f"watertight: {mesh.is_watertight}")
    print(f"boundary_edges: {boundary_edges}")
    print(f"nonmanifold_edges: {nonmanifold_edges}")
    print(f"inside_is_fluid: True")
    print(f"fluid-point ray disagreements: {stats.get('fluid_ray_disagreements', 0)}")
    print(f"fluid-point ray overrides: {stats.get('fluid_ray_overrides', 0)}")
    print(f"solid_mismatch_count_certain: {solid_mismatch}")
    print(f"ignored_near_surface_mismatches: {ignored_mismatch}")
    print(f"reference fluid probe points outside cavity: {int(np.count_nonzero(~probe_fluid[:max(1, probe_fluid.size - 8)]))}")
    if args.lid_gap == 0.0 and args.lid_radius_scale == 1.0 and (not mesh.is_watertight or boundary_edges != 0):
        raise SystemExit("lidded bent-pipe test failed: exact generated cavity mesh is not watertight")
    if solid_mismatch > args.max_solid_mismatch:
        raise SystemExit(
            f"lidded bent-pipe test failed: {solid_mismatch} certain solid/fluid mismatches > {args.max_solid_mismatch}"
        )


def test_stl_bent_pipe(args: argparse.Namespace) -> None:
    """Validate the open pipe-wall STL where both bore and exterior are fluid."""
    finalize_grid_args(args)
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    geometry = workdir / "bent_pipe.stl"
    h5 = workdir / "bent_pipe_coeff.h5"
    mesh = make_bent_pipe_mesh(args)
    mesh.export(geometry)

    mesh, vertices, faces = load_stl_mesh(geometry, repair=not args.no_repair)
    boundary_edges, nonmanifold_edges = mesh_edge_counts(mesh)
    coeff_args = argparse.Namespace(**vars(args))
    coeff_args.geometry = geometry
    coeff_args.output = h5
    if args.auto_fluid_points:
        coeff_args._fluid_points = bent_pipe_validation_fluid_points(args)
    coef, stats = stl_coeff_from_mesh(mesh, vertices, faces, coeff_args)
    write_hdf5(h5, coef, coeff_args, geometry, stl_extra_attrs(mesh, stats, coeff_args))

    analytic_solid, certain = bent_pipe_reference_mask(args)
    stl_solid = np.abs(coef) > 1.0e20
    mismatch = stl_solid != analytic_solid
    certain_mismatch = mismatch & certain
    ignored = mismatch & ~certain
    solid_mismatch = int(np.count_nonzero(certain_mismatch))
    ignored_mismatch = int(np.count_nonzero(ignored))
    bore_points = bent_pipe_validation_fluid_points(args)
    bore_solid, _ = analytic_bent_pipe_wall_mask(bore_points, args)

    print(f"Bent pipe STL: {geometry}")
    print(f"HDF5 coefficients: {h5}")
    print(f"XDMF coefficients: {coefficient_xdmf_path(h5)}")
    print(f"vertices/faces: {len(mesh.vertices)} / {len(mesh.faces)}")
    print(f"watertight: {mesh.is_watertight}")
    print(f"boundary_edges: {boundary_edges}")
    print(f"nonmanifold_edges: {nonmanifold_edges}")
    print(f"fluid-point ray disagreements: {stats.get('fluid_ray_disagreements', 0)}")
    print(f"fluid-point ray overrides: {stats.get('fluid_ray_overrides', 0)}")
    print(f"solid_mismatch_count_certain: {solid_mismatch}")
    print(f"ignored_near_surface_mismatches: {ignored_mismatch}")
    print(f"reference fluid probe points inside wall: {int(np.count_nonzero(bore_solid))}")
    if not mesh.is_watertight or boundary_edges != 0:
        raise SystemExit("bent-pipe STL test failed: generated wall mesh is not watertight")
    if solid_mismatch > args.max_solid_mismatch:
        raise SystemExit(
            f"bent-pipe STL test failed: {solid_mismatch} certain solid/fluid mismatches > {args.max_solid_mismatch}"
        )
    if np.count_nonzero(bore_solid) != 0:
        raise SystemExit("bent-pipe STL test failed: a nominal fluid probe lies inside the pipe wall")


def test_stl_restricted_lidded_pipe(args: argparse.Namespace) -> None:
    """Validate that a real narrow throat is preserved as fluid, not treated as a defect."""
    finalize_grid_args(args)
    if args.restriction_radius <= 0.0:
        args.restriction_radius = max(abs(args.lid_gap), 0.001)
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    args.output_dir = str(workdir)
    inner, start_lid, end_lid = make_lidded_bent_pipe_stls(args)
    h5 = workdir / "restricted_lidded_bent_pipe_coeff.h5"

    mesh, vertices, faces = load_stl_meshes([inner, start_lid, end_lid], repair=not args.no_repair)
    boundary_edges, nonmanifold_edges = mesh_edge_counts(mesh)
    coeff_args = argparse.Namespace(**vars(args))
    coeff_args.geometry = [str(inner), str(start_lid), str(end_lid)]
    coeff_args.output = h5
    coeff_args.inside_is_fluid = True
    coef, stats = stl_coeff_from_mesh(mesh, vertices, faces, coeff_args)
    write_hdf5(h5, coef, coeff_args, ";".join(str(p) for p in (inner, start_lid, end_lid)), stl_extra_attrs(mesh, stats, coeff_args))

    analytic_solid, certain = bent_pipe_cavity_reference_mask(args)
    stl_solid = np.abs(coef) > 1.0e20
    mismatch = stl_solid != analytic_solid
    solid_mismatch = int(np.count_nonzero(mismatch & certain))
    ignored_mismatch = int(np.count_nonzero(mismatch & ~certain))

    throat_s = args.restriction_center * bent_pipe_total_length(args)
    probe_s = np.linspace(throat_s - 0.45 * args.restriction_width, throat_s + 0.45 * args.restriction_width, 9)
    probe_s = np.clip(probe_s, 0.0, bent_pipe_total_length(args))
    coeff_args._check_fluid_points = np.vstack([bent_pipe_global_center_at_s(float(s), args) for s in probe_s])
    check_report = generic_stl_checks(mesh, vertices, faces, coeff_args)

    print(f"Restricted lidded bent pipe STLs: {inner}; {start_lid}; {end_lid}")
    print(f"HDF5 coefficients: {h5}")
    print(f"XDMF coefficients: {coefficient_xdmf_path(h5)}")
    print(f"restriction_radius: {args.restriction_radius}")
    print(f"restriction_width: {args.restriction_width}")
    print(f"vertices/faces: {len(mesh.vertices)} / {len(mesh.faces)}")
    print(f"watertight: {mesh.is_watertight}")
    print(f"boundary_edges: {boundary_edges}")
    print(f"nonmanifold_edges: {nonmanifold_edges}")
    print(f"inside_is_fluid: True")
    print(f"solid_mismatch_count_certain: {solid_mismatch}")
    print(f"ignored_near_surface_mismatches: {ignored_mismatch}")
    print_generic_stl_check_report(check_report)
    if args.lid_gap == 0.0 and args.lid_radius_scale == 1.0 and (not mesh.is_watertight or boundary_edges != 0):
        raise SystemExit("restricted lidded-pipe test failed: exact generated cavity mesh is not watertight")
    enforce_generic_stl_checks(check_report, coeff_args, label="restricted lidded-pipe probe")
    if solid_mismatch > args.max_solid_mismatch:
        raise SystemExit(
            f"restricted lidded-pipe test failed: {solid_mismatch} certain solid/fluid mismatches > {args.max_solid_mismatch}"
        )


def add_bent_pipe_restriction_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--restriction-radius", type=float, default=0.001,
                        help="minimum bore radius at the throat; set comparable to a lid gap")
    parser.add_argument("--restriction-width", type=float, default=0.08)
    parser.add_argument("--restriction-center", type=float, default=0.5,
                        help="throat center as fraction of pipe centerline length")


def add_generic_stl_check_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--check-fluid-points", default=None,
                        help="points that must be fluid; text/CSV, .npy, or HDF5")
    parser.add_argument("--check-solid-points", default=None,
                        help="points that must be solid; text/CSV, .npy, or HDF5")
    parser.add_argument("--max-fluid-check-failures", type=int, default=0)
    parser.add_argument("--max-solid-check-failures", type=int, default=0)
    parser.add_argument("--expect-watertight", action="store_true")
    parser.add_argument("--max-boundary-edges", type=int, default=10**18)
    parser.add_argument("--max-nonmanifold-edges", type=int, default=10**18)


def add_stl_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--classification", choices=["fast-winding", "winding"], default="fast-winding")
    parser.add_argument("--winding-threshold", type=float, default=0.5)
    parser.add_argument("--winding-margin", type=float, default=1.0e-3)
    parser.add_argument("--scale", type=float, nargs="+", default=[1.0],
                        help="scale STL coordinates before classification; pass one value or x y z")
    parser.add_argument("--translate", type=float, nargs=3, default=[0.0, 0.0, 0.0],
                        help="translation added after STL scaling")
    parser.add_argument("--bbox-padding-cells", type=float, default=2.0,
                        help="padding, in grid cells, added around the STL bounding box before culling")
    parser.add_argument("--no-bbox-cull", action="store_true",
                        help="classify the full grid instead of only the padded STL bounding box")
    parser.add_argument("--fluid-points", default=None,
                        help="optional known-fluid point cloud that reinforces winding classification; text/CSV, .npy, or HDF5")
    parser.add_argument("--auto-fluid-points", action="store_true",
                        help="for sphere validation/stress tests, generate known-fluid points on the box boundary")
    parser.add_argument("--fluid-ray-policy", choices=["check", "ambiguous", "override"], default="ambiguous",
                        help="how sure-fluid ray votes affect winding: report only, fix ambiguous points, or override all confident disagreements")
    parser.add_argument("--fluid-ray-scope", choices=["ambiguous", "all"], default="ambiguous",
                        help="which grid points get expensive fluid-probe ray votes; default only checks winding-ambiguous points")
    parser.add_argument("--fluid-ray-max-seeds", type=int, default=17,
                        help="maximum number of fluid seeds used for parity votes; <=0 uses all points")
    parser.add_argument("--fluid-ray-majority", type=float, default=0.5)
    parser.add_argument("--fluid-ray-margin", type=float, default=0.0)
    parser.add_argument("--fluid-ray-hit-tol", type=float, default=1.0e-10)
    parser.add_argument("--chunk-size", type=int, default=2000000)
    parser.add_argument("--jobs", type=int, default=1,
                        help="parallel worker processes for tiled coefficient generation")
    parser.add_argument("--tile-size", type=int, nargs=3, default=None,
                        help="write one chunked HDF5 file from independent coefficient tiles")
    parser.add_argument("--h5-compression", choices=["none", "lzf", "gzip"], default="none",
                        help="compression used for chunked coefficient datasets")
    parser.add_argument("--gzip-level", type=int, default=1,
                        help="gzip compression level when --h5-compression=gzip")
    parser.add_argument("--no-tiled-output", action="store_true",
                        help="for stl-ibm-coeff, use the legacy dense in-memory HDF5 writer")
    parser.add_argument("--distance-mode", choices=["hybrid", "winding-root"], default="hybrid")
    parser.add_argument("--inside-is-fluid", action="store_true",
                        help="interpret the inside of the STL shell as fluid and the outside as solid")
    add_generic_stl_check_args(parser)
    parser.add_argument("--no-repair", action="store_true")

def add_bent_pipe_lid_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--lid-gap", type=float, default=0.0,
                        help="offset lids along the pipe axis; nonzero creates an approximate/non-watertight closure")
    parser.add_argument("--lid-radius-scale", type=float, default=1.0,
                        help="scale lid radius relative to the pipe inner radius")


def add_bent_pipe_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--pipe-center", type=float, nargs=3, default=[0.5, 0.5, 0.5])
    parser.add_argument("--inner-radius", type=float, default=0.035)
    parser.add_argument("--outer-radius", type=float, default=0.055)
    parser.add_argument("--bend-radius", type=float, default=0.18)
    parser.add_argument("--straight-length", type=float, default=0.22)
    parser.add_argument("--straight-sections", type=int, default=18)
    parser.add_argument("--bend-sections", type=int, default=28)
    parser.add_argument("--around", type=int, default=48)
    parser.add_argument("--reference-tolerance", type=float, default=0.012)


def add_common_grid_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--grid-file", required=True,
                        help="HDF5 grid exported by the mobygrid executable")
    parser.add_argument("--re", type=float, default=100.0)


def main(argv: list[str] | None = None) -> int:
    """Build the command-line interface and dispatch to the selected geometry utility."""
    parser = argparse.ArgumentParser(description="mobyDiff STL-to-IBM coefficient utility")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("make-sphere-stl", help="write an STL sphere test body")
    p.add_argument("--output", required=True)
    p.add_argument("--center", type=float, nargs=3, default=[0.5, 0.5, 0.5])
    p.add_argument("--radius", type=float, default=0.2)
    p.add_argument("--subdivisions", type=int, default=5)
    p.add_argument("--open-cap-angle", type=float, default=0.0,
                   help="remove a polar cap with this half-angle, useful for non-watertight tests")
    p.set_defaults(func=make_sphere_stl)

    p = sub.add_parser("make-bent-pipe-stl", help="write an oblique open-ended 90-degree bent-pipe wall STL")
    p.add_argument("--output", required=True)
    add_bent_pipe_args(p)
    p.set_defaults(func=make_bent_pipe_stl)

    p = sub.add_parser("make-lidded-bent-pipe-stl", help="write inner pipe surface plus two lid STLs")
    p.add_argument("--output-dir", required=True)
    add_bent_pipe_args(p)
    add_bent_pipe_lid_args(p)
    p.set_defaults(func=make_lidded_bent_pipe_stl)

    p = sub.add_parser("stl-ibm-coeff", help="compute IBM coefficients from one or more STL triangle meshes")
    p.add_argument("--geometry", required=True, nargs="+")
    p.add_argument("--output", required=True)
    add_common_grid_args(p)
    add_stl_args(p)
    p.set_defaults(func=coeff_from_stl)

    p = sub.add_parser("block-active", help="write per-block solid-removal flags into a coefficient file")
    p.add_argument("--geometry", required=True, nargs="+")
    p.add_argument("--output", required=True,
                   help="existing coefficient HDF5 to amend with block_active")
    p.add_argument("--block-nb", type=int, required=True,
                   help="cubic block edge in cells ([blocks] nb in the solver input)")
    add_common_grid_args(p)
    add_stl_args(p)
    p.set_defaults(func=block_active_from_stl)

    p = sub.add_parser("block-table", help="write block-table IBM coefficients for refined runs")
    p.add_argument("--geometry", required=True, nargs="+")
    p.add_argument("--output", required=True)
    p.add_argument("--block-nb", type=int, required=True)
    p.add_argument("--levels", type=int, default=2,
                   help="number of levels (refine_levels + 1)")
    p.add_argument("--refine-box", type=float, nargs="+", action="append", default=None,
                   help="x0 x1 y0 y1 z0 z1 [level]: box refinement ADDED to the "
                        "body-driven one; repeatable; the optional 7th value is "
                        "the box's target level (default the finest)")
    p.add_argument("--refine-dims", choices=("xyz", "xz"), default="xyz",
                   help="refined directions ([blocks] refine_dims): xyz octree "
                        "(default) or xz quadtree (y keeps the global line)")
    p.add_argument("--no-dwall", action="store_true",
                   help="skip the per-leaf dwall_blocks wall-distance dataset (RANS)")
    p.add_argument("--keep-buried", action="store_true",
                   help="keep leaves buried inside the body (required when the "
                        "penalization-force statistic is used: a removed core's "
                        "closed faces absorb pressure loading outside the coef "
                        "bookkeeping)")
    add_common_grid_args(p)
    add_stl_args(p)
    p.set_defaults(func=block_table_from_stl)

    p = sub.add_parser("check-stl-geometry", help="run generic STL mesh/probe classification checks")
    p.add_argument("--geometry", required=True, nargs="+")
    add_stl_args(p)
    p.set_defaults(func=check_stl_geometry)

    p = sub.add_parser("test-stl-sphere", help="compare STL sphere coefficients with analytic sphere coefficients")
    p.add_argument("--workdir", default=str(ROOT / "tools" / "mobygeom_tests"))
    p.add_argument("--center", type=float, nargs=3, default=[0.5, 0.5, 0.5])
    p.add_argument("--radius", type=float, default=0.2)
    p.add_argument("--subdivisions", type=int, default=6)
    p.add_argument("--open-cap-angle", type=float, default=0.0,
                   help="remove a polar cap with this half-angle, useful for non-watertight tests")
    p.add_argument("--abs-tolerance", type=float, default=2.5e1)
    p.add_argument("--rel-tolerance", type=float, default=5.0e-2)
    p.add_argument("--max-solid-mismatch", type=int, default=0)
    add_common_grid_args(p)
    add_stl_args(p)
    p.set_defaults(func=test_stl_sphere)

    p = sub.add_parser("test-stl-bent-pipe", help="validate an oblique open-ended 90-degree bent-pipe wall STL")
    p.add_argument("--workdir", default=str(ROOT / "tools" / "mobygeom_tests" / "bent_pipe"))
    p.add_argument("--max-solid-mismatch", type=int, default=0)
    add_common_grid_args(p)
    add_stl_args(p)
    add_bent_pipe_args(p)
    p.set_defaults(func=test_stl_bent_pipe)

    p = sub.add_parser("test-stl-lidded-bent-pipe", help="validate a lidded bent pipe where only the bore is fluid")
    p.add_argument("--workdir", default=str(ROOT / "tools" / "mobygeom_tests" / "lidded_bent_pipe"))
    p.add_argument("--max-solid-mismatch", type=int, default=0)
    add_common_grid_args(p)
    add_stl_args(p)
    add_bent_pipe_args(p)
    add_bent_pipe_lid_args(p)
    p.set_defaults(func=test_stl_lidded_bent_pipe)

    p = sub.add_parser("test-stl-restricted-lidded-pipe",
                       help="validate a lidded bent pipe with a very small physical throat")
    p.add_argument("--workdir", default=str(ROOT / "tools" / "mobygeom_tests" / "restricted_lidded_pipe"))
    p.add_argument("--max-solid-mismatch", type=int, default=0)
    add_common_grid_args(p)
    add_stl_args(p)
    add_bent_pipe_args(p)
    add_bent_pipe_lid_args(p)
    add_bent_pipe_restriction_args(p)
    p.set_defaults(func=test_stl_restricted_lidded_pipe)

    p = sub.add_parser("stress-stl-watertightness", help="run damaged-STL robustness stress tests")
    p.add_argument("--workdir", default=str(ROOT / "tools" / "mobygeom_tests" / "stl_stress"))
    p.add_argument("--center", type=float, nargs=3, default=[0.5, 0.5, 0.5])
    p.add_argument("--radius", type=float, default=0.2)
    p.add_argument("--subdivisions", type=int, default=5)
    p.add_argument("--cases", nargs="*", default=None,
                   help="optional list such as clean cap_10 delete_0.01 crack_band_0.01")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--abs-tolerance", type=float, default=2.5e1)
    p.add_argument("--rel-tolerance", type=float, default=5.0e-1)
    p.add_argument("--max-solid-mismatch", type=int, default=0)
    p.add_argument("--write-coefficients", action="store_true")
    p.add_argument("--strict", action="store_true",
                   help="exit nonzero if any case exceeds the solid-mismatch threshold")
    add_common_grid_args(p)
    add_stl_args(p)
    p.set_defaults(func=stress_stl_watertightness)

    args = parser.parse_args(argv)
    result = args.func(args)
    return 0 if result is None else 0


if __name__ == "__main__":
    raise SystemExit(main())
