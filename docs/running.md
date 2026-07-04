# Running the solver

## Invocation

The solver reads a single `.ini` configuration file and is **always** launched through
`mpirun`, even on a single rank:

```bash
mpirun -n <ranks> ./build_gpu/main path/to/input.ini
```

- Use `build_cpu/main` for the CPU reference build, `build_gpu/main` for the GPU build.
- On the GPU build, one MPI rank drives one GPU; launch as many ranks as you have devices.

Example — a single-GPU minimal channel and an 8-rank CPU channel:

```bash
mpirun -n 1 ./build_gpu/main tutorials/min_channel/input_gpu.ini
mpirun -n 8 ./build_cpu/main tutorials/channel_kmm180/input.ini
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

For IBM (immersed-body) cases, the grid must first be exported with `mobygrid` so the
Python preprocessor can classify the STL geometry against the solver's exact node lines —
see [Tutorials](tutorials.md) and the [tools reference](tools.md).

## Time stepping

Time advancement is controlled by `[time]`: a nominal `dt`, a step or time limit
(`nsteps` / `t_final`), and stability caps (`cflmax`, `pecletmax`, `dtmax`). The actual
step is the largest value satisfying the CFL and Péclet limits, capped by `dtmax`.

## Output

Field snapshots are written as HDF5 by `[output]`:

```ini
[output]
field_interval = 50000            ; steps between field dumps
field_prefix   = channel_field    ; → channel_field_1.h5, channel_field_2.h5, ...
```

Each snapshot holds the velocity components and pressure. For refined (multi-level) runs
the file uses the block-table layout (one dataset row-range per block); a companion `.xdmf`
is written for cases where it applies. To reassemble a refined field onto a single global
grid for visualization or comparison, use `tools/compare_fields.py --export-global`.

Channel cases additionally accumulate turbulence statistics (see the `[case.channel]`
keys), written to the configured `stats_file`.

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
