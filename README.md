# mobydiff
<img align="left" src="https://github.com/davecats/mobydiff/blob/main/.mobydiff.png" width="400"> <br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/>
A WIP simple and (reasonably) efficient solver of the incompressible Navier-Stokes equations.

<br clear="left"/>

## Overview

mobydiff integrates the incompressible Navier–Stokes equations with a
second-order finite-difference discretisation on a staggered (MAC) Cartesian
grid. The solver state lives on a set of equal-size structured blocks (BCM
style) that tile the domain, which enables 2:1 local mesh refinement and the
removal of blocks buried inside an immersed body.

- **Spatial discretisation:** 2nd-order finite differences on a staggered
  grid; each direction can be uniform or stretched (`uniform`, `cosine`,
  `tanh`, or the `natural` near-wall channel stretching).
- **Time integration:** explicit 3-stage Runge–Kutta, with optional CFL /
  diffusive time-step limiting.
- **Pressure–velocity coupling:** coupled velocity–pressure projection solved
  with a red-black SOR sweep; on refined grids the 2:1 interfaces are treated
  by a composite projection.
- **Immersed boundaries:** volume-penalization IBM. Geometry comes either from
  an analytic body or from a precomputed coefficient file generated from an
  STL mesh.
- **Turbulence modelling:** optional Smagorinsky LES.
- **Mesh refinement:** BCM-style equal-size blocks with 2:1 refinement (by box
  or driven by the immersed geometry); blocks fully inside a solid body are
  removed.
- **Parallelism:** MPI with a 3D Cartesian decomposition and 26-neighbour halo
  exchange, plus OpenMP `target` offload for GPUs (GPU-aware MPI on the GPU
  path).
- **I/O:** parallel HDF5 fields, usable as restarts.

## Dependencies

- _cmake_ ≥ 3.20
- A Fortran compiler: _gcc/gfortran_ ≥ 13, or NVHPC, or Intel
- _libhdf5_ ≥ 1.14 with **parallel** (MPI) support
- An MPI library (always required, even for single-rank runs)

For GPU support:
- NVIDIA NVHPC SDK + CUDA 12+ (GPU builds need a GPU-aware MPI stack)

## Build

`compile.sh` configures and builds into `build_cpu/` or `build_gpu/`:

```bash
./compile.sh cpu     # CPU build  -> build_cpu/main,  build_cpu/mobygrid
./compile.sh gpu     # GPU build  -> build_gpu/main,  build_gpu/mobygrid
```

Build the CPU path as well as the GPU path: the CPU build is the reference for
debugging and can be compiled with `-Mbounds` to catch out-of-bounds accesses
that show up only as opaque illegal-address faults on the GPU.

Set `HDF5_ROOT=/path/to/parallel-hdf5` if CMake cannot locate parallel HDF5.

## Usage

Always launch through `mpirun`, even on a single rank:

```bash
mpirun -n 1 ./build_gpu/main path/to/input.ini      # single GPU
mpirun -n 8 ./build_cpu/main path/to/input.ini      # 8 MPI ranks (CPU)
```

The simulation is fully described by an `.ini` input file (see below). Output
fields are written as HDF5 and can be fed back in via `[restart]`.

## Input file

The input is an INI-style file grouped into sections. A minimal channel or
generic case needs only a handful of them; the cases under `tutorials/` are
complete, runnable examples. The main sections:

| Section        | Purpose |
|----------------|---------|
| `[case]`       | `name = generic` or `name = channel` (case-specific setup/forcing); `[case.channel]` holds channel options. |
| `[grid]`       | `nx,ny,nz` cells and `lx,ly,lz` domain lengths. |
| `[grid.x/y/z]` | Per-direction `distribution` (`uniform`/`cosine`/`tanh`/`natural`) and stretching parameters; `subdivided = true` builds a line as the midpoint subdivision of its half-resolution line. |
| `[mpi]`        | `dims = 0 0 0` lets MPI choose the Cartesian decomposition. |
| `[blocks]`     | Block refinement: `nb` (block size in cells; even, ≥4, must divide the grid), `refine = x0 x1 y0 y1 z0 z1` boxes, `refine_levels`, `refine_body = true` (geometry-driven), `remove_solid`. Omit for one block per rank. |
| `[flow]`       | `re` (Reynolds number), body forcing `forcing_x/y/z`, optional `initial_*` / `initial_noise`. |
| `[time]`       | `dt`, `nsteps`, `t_final`, and limits `cflmax`, `pecletmax`, `dtmax`. |
| `[pressure]`   | `niter` (SOR iterations per projection) and `sor` (relaxation factor). |
| `[ibm]`        | `enabled`, and `coeff_file` for file-based geometry. |
| `[les]`        | `model = smagorinsky`, `cs` (Smagorinsky constant). |
| `[boundary]`   | Per-face periodicity (`periodic_x/y/z`) and Dirichlet/Neumann conditions with values (e.g. `x_min_u_value`, `x_max_p_type = dirichlet`). |
| `[output]`     | `field_interval` (steps between HDF5 dumps; 0 disables) and `field_prefix`. |
| `[restart]`    | `file = field_NNN.h5` to continue from a saved field. |

## Tutorials

Ready-to-run cases in `tutorials/`:

- `min_channel/` — small turbulent channel (single- and block-refined inputs).
- `channel_kmm180/` — Re_τ = 180 channel restart/reference data.
- `wavychannel/` — wavy-wall channel with an analytic immersed boundary.
- `sailplane/` — external aerodynamics with a file-based STL geometry.
- `interface_decay/` — white-noise decay on a refined patch (2:1-interface
  stability check).
- `interface_shear_mode/` — interface shear-mode check.

## Tools

In `tools/` (and the `mobygrid` companion executable):

- `mobygrid` — exports the staggered grid (node lines + metadata) from an
  `input.ini` to HDF5; the geometry preprocessor uses it to stay consistent
  with the solver grid.
- `mobygeom.py` — generates static IBM coefficient files from STL meshes (see
  `tools/README_mobygeom.md`).
- `compare_fields.py` — compares two HDF5 field files (used for bit-exact
  verification; `--export-global` reassembles block-table output).
- `check_parabolic_channel.py`, `check_interface_decay.py`,
  `channel_interface_validation.py`, `make_channel_restart.py` — validation and
  setup helpers.

## Repository layout

```
src/                solver sources
  main.f90          solver entry point
  mobygrid.f90      serial grid / preprocessing tool
  modules/          numerics, blocks, comm, IBM, LES, I/O, config, flow cases
tools/              pre/post-processing and validation scripts
tutorials/          runnable example cases
validation/         turbulence-validation setups
docs/               design notes and strategy documents
```
