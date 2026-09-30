# Running the solver

## Invocation

The solver reads a single `.ini` configuration file and is **always** launched through
`mpirun`, even on a single rank:

```bash
mpirun -n <ranks> ./build_gpu/moby_solve path/to/input.ini
```

- Use `build_cpu/moby_solve` for the CPU reference build, `build_gpu/moby_solve` for the
  GPU build (`main` in each build directory is a compatibility symlink to `moby_solve`).
- On the GPU build, one MPI rank drives one GPU; launch as many ranks as you have devices.
  On multi-node runs the rank-to-device order can be overridden with the `MOBY_GPU_ORDER`
  environment variable; `tools/moby_tune.sh` finds a good mapping by measurement.

Example — a single-GPU minimal channel and an 8-rank CPU channel:

```bash
mpirun -n 1 ./build_gpu/moby_solve tutorials/min_channel/input_gpu.ini
mpirun -n 8 ./build_cpu/moby_solve tutorials/channel_kmm180/input.ini
```

## Domain decomposition

The domain is split by a 3D MPI Cartesian decomposition. The process grid is set in the
`[mpi]` section:

```ini
[mpi]
dims = 0 0 0     ; 0 = let MPI choose the factorization for this direction
```

A `0` lets MPI pick the number of ranks along that axis; fixed non-zero values pin the
process grid. Blocks are distributed over the ranks along a Z-order (Morton) space-filling
curve, and each rank owns a contiguous run of block ids. Results are independent of the
number of ranks and of the block count (`[blocks] nb`) — see
[Numerical methods](numerical-methods.md#block-structured-grid-and-21-refinement).

## Grid, blocks and refinement

The base grid resolution and extent come from `[grid]`, with per-direction stretching in
`[grid.x]`, `[grid.y]`, `[grid.z]`. Block decomposition and 2:1 refinement are configured
in `[blocks]`. The full key list is in the [configuration reference](configuration.md).

## Preprocessing immersed bodies (`moby_prepare`)

Analytic immersed walls (`[ibm] wall_shape`) can be set up inline by the solver. STL
geometry — and, optionally, analytic geometry too — goes through the MPI-parallel
preprocessor, which reuses the solver's own grid, block and classification code:

```bash
mpirun -n 4 ./build_cpu/moby_prepare input.ini          # writes [case] file, else <field_prefix>.case.h5
mpirun -n 4 ./build_cpu/moby_solve   input.ini          # reads it; stops if missing or stale
```

`input.ini` is the run's `.ini` plus, for an STL body, the geometry declaration
(`[ibm] stl_file`, optional `stl_scale` / `stl_translate`). The output is the block-table
**case file** — node lines, leaf block table, per-block IBM coefficients, solid-removal
masks and wall distance — named by `[case] file` (an explicit second argument to
`moby_prepare` overrides it). EVERY run reads one, body or not: an unset `[blocks] nb` is one
block per rank of the PREPARE run, so prepare on the rank count you solve with. Prepare with
the CPU build (the canonical one). Regenerate the case file whenever the grid, the block layout or `[flow] re`
changes (the coefficients carry the 1/Re scaling). The design is in
[`prepare_solve_strategy.md`](prepare_solve_strategy.md); the retired Python `mobygeom`
pipeline is kept only as a cross-implementation reference.

## Time stepping

Time advancement is controlled by `[time]`: a nominal `dt`, a step or time limit
(`nsteps` / `t_final`), and stability caps (`cflmax`, `pecletmax`, `dtmax`). The actual
step is the largest value satisfying the CFL and Péclet limits, capped by `dtmax`.

## Output

Field snapshots are written as HDF5 by `[output]`:

```ini
[output]
field_interval = 50000            ; steps between field dumps
field_prefix   = channel_field    ; → channel_field_<step>.h5
```

Each snapshot holds the velocity components and pressure (`un vn wn pn`), plus whatever
the active models carry (`nut`, the RANS scalars `k`/`omega`/…, passive scalars by name).
Snapshots use the block-table layout (one dataset row per block, plus the `blocks` table);
no XDMF is written. To reassemble a field onto a single global grid for visualization or
comparison, use `tools/compare_fields.py --export-global` (for very deep refinement the
finest-lattice reassembly can exhaust memory).

Channel and boundary-layer cases accumulate turbulence statistics (the `[case.channel]` /
`[case.boundarylayer]` `stats_*` keys), passive scalars their own (`[scalar] stats_*`), and
the airfoil case writes C_L/C_D to its `runtime_file`. `[output] profile = true` prints a
per-phase timing breakdown after the run.

## Restart

To continue from a snapshot, point `[restart]` at a previous field file:

```ini
[restart]
file = channel_field_50000.h5
```

Restart works on any number of ranks — the block table in the file is redistributed over
the current process grid. Legacy single-level (global-3D) restart files are still read.
Solver parameters (e.g. the pressure `sor` / `niter`) come from the current `.ini`, not
from the restart file, so you can change them on continuation.
