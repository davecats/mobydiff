# mobyGeom STL-to-IBM Utility

**RETIRED FOR PRODUCTION (prepare/solve split P1b, 2026-07-16,
`docs/prepare_solve_strategy.md`).** The geometry subcommands below
(`stl-ibm-coeff`, `block-active`, `block-table`) are superseded by the
MPI-parallel Fortran preprocessor: `moby_prepare case.ini case.h5` with
`[ibm] stl_file` (one line per STL, `stl_scale`/`stl_translate` for
transforms, `[blocks] keep_buried` replacing `--keep-buried`) — no venv,
no separate mobygrid handshake, identical file contract. mobygeom is KEPT
working as the independent cross-implementation reference: the
`validation/prepare/` gates compare moby_prepare's output against it.
The STL generators/checkers (`make-*-stl`, watertightness tests) stay in
normal use.

`mobygeom.py` generates static IBM coefficient fields from STL triangle meshes.
The grid definition always comes from a Fortran-generated `mobygrid` HDF5 file,
so the preprocessor uses the same node lines, stretching, periodicity, and
staggered coordinates as the solver.

## Grid Export

The `mobygrid` executable was absorbed by the prepare/solve split (P3):
`moby_prepare` writes the node lines and the mobygrid-format attributes
into every case file, so a prepared case file serves directly as
`--grid-file`:

```bash
mpirun -n 4 ./build_cpu/moby_prepare case.ini case.h5
# mobygeom cross-checks read the grid straight from it:
#   --grid-file case.h5
```

Legacy `mobygrid`-written grid files (e.g. the committed
`validation/channel_interface/les_ibm/grid.h5`) remain readable; they
contain `/x_nodes`, `/y_nodes`, `/z_nodes` and the grid metadata.

## Coefficient Generation

```bash
/home/davide/ibmc/bin/python tools/mobygeom.py stl-ibm-coeff \
  --geometry body.stl \
  --output ibm_coeff.h5 \
  --grid-file grid.h5 \
  --re 100
```

Use `--scale` and `--translate tx ty tz` when the STL coordinates need to be
converted or positioned in the solver domain before coefficient generation.

The HDF5 file contains `coef` with shape `(nx+2, ny+2, nz+2, 3)` plus component
virtual datasets `coef_u`, `coef_v`, and `coef_w`. A sidecar `.xdmf` file is
written for ParaView and uses the staggered coordinates derived from the
`mobygrid` node lines.

`stl-ibm-coeff` writes tiled/chunked HDF5 by default. Use `--tile-size ix iy iz`
to tune tile memory and `--jobs N` to process independent tiles in parallel.
`--no-tiled-output` restores the dense in-memory writer for small debugging
cases.

## Classification Checks

For external STL bodies, the utility classifies only a padded geometry bounding
box by default and treats points outside that box as fluid. This keeps large
domains cheap when the body occupies only a small fraction of the volume. Use
`--bbox-padding-cells N` to change the padding, or `--no-bbox-cull` for a full
grid classification debug run.

Use point probes as fail-fast checks:

```bash
/home/davide/ibmc/bin/python tools/mobygeom.py stl-ibm-coeff \
  --geometry body.stl \
  --output ibm_coeff.h5 \
  --grid-file grid.h5 \
  --check-fluid-points fluid_points.txt \
  --check-solid-points solid_points.txt
```

The probe files can be whitespace text, CSV, `.npy`, or HDF5. HDF5 probe files
may use dataset `points`; fluid probe files may also use `fluid_points`, and
solid probe files may also use `solid_points`.

Known-fluid ray reinforcement is still available with `--fluid-points`, but it
is intentionally separate from the cheaper `--check-fluid-points` validation.
By default, ray checks are restricted to winding-ambiguous grid points. Add
`--fluid-ray-scope all` only for deliberate full ray-vote debugging.

For closed cavity STL sets where the inside is the fluid region, add
`--inside-is-fluid`.

## STL Tests

All grid-based tests require a `mobygrid` file:

```bash
./build_cpu/mobygrid smoke/cpu_vs_gpu/tiny.ini tools/mobygeom_tests/tiny_grid.h5
/home/davide/ibmc/bin/python tools/mobygeom.py test-stl-sphere \
  --grid-file tools/mobygeom_tests/tiny_grid.h5 \
  --re 100 \
  --subdivisions 6
```

A small non-watertight smoke test can be run by removing a cap from the
generated sphere:

```bash
/home/davide/ibmc/bin/python tools/mobygeom.py test-stl-sphere \
  --grid-file tools/mobygeom_tests/tiny_grid.h5 \
  --subdivisions 5 \
  --open-cap-angle 12 \
  --abs-tolerance 100 \
  --rel-tolerance 1 \
  --max-solid-mismatch 20
```

The STL robustness stress test generates a clean sphere and controlled damaged
variants:

```bash
/home/davide/ibmc/bin/python tools/mobygeom.py stress-stl-watertightness \
  --grid-file tools/mobygeom_tests/tiny_grid.h5 \
  --re 100
```

The report includes watertightness, Euler number, boundary/non-manifold edge
counts, ambiguous classification points, fluid-ray disagreements/overrides,
fallback segment count, solid/fluid mismatches, and coefficient errors.
