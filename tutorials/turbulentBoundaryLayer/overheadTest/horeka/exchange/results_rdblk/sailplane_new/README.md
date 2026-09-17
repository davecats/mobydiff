# Sailplane Generic-Flow Tutorial

The source STL is in millimetres and is symmetric about its original `y = 0`
plane. The tutorial uses metres, keeps the symmetry plane at computational
`y = 0`, and simulates only the `y >= 0` half-domain.

The original STL bounding box is:

```text
min  = (3.226455e-10, -9.744634, -0.250183) m
max  = (8.644325,      9.744634,  1.700000) m
size = (8.644325,     19.489268,  1.950183) m
```

The computational domain in `input.ini` is:

```text
Lx = 5.0 * length = 43.22162499838677
Ly = 2.5 * span   = 48.72317
Lz = 5.0 * height = 9.750915
```

The coefficient generation command scales the STL by `0.001` and translates it
by `(17.288649999032064, 0.0, 4.150549)`. This places the full mirrored STL at
the centre of the corresponding full domain, while the solver only uses the
positive-`y` half.

Generate the IBM coefficients with moby_prepare (prepare/solve split,
`docs/prepare_solve_strategy.md` — a `[blocks] nb = 10` block layout plus
the STL declaration replace the retired mobygrid+mobygeom pipeline):

```bash
cd tutorials/sailplane
# prep_blocks.ini = input.ini with:
#   [blocks] nb = 10
#   [ibm] stl_file = "FRUE V0 ohneRundung.stl"
#         stl_scale = 0.001
#         stl_translate = 17.288649999032064 0.0 4.150549
mpirun -n 2 ../../build_cpu/moby_prepare prep_blocks.ini sailplane_case.h5
# solve with [blocks] nb = 10 and [ibm] coeff_file = sailplane_case.h5
```

The committed `sailplane_ibm_coeff.h5` (legacy mobygeom global layout,
usable without `[blocks] nb`) remains for the original one-step tutorial;
a 1-step solve from the prepared case file is bit-exact against it
(validation/prepare/run_gates_big.sh). If the grid changes in `input.ini`,
rerun moby_prepare. The retired mobygeom cross-check can read the grid
straight from the case file (`--grid-file sailplane_case.h5`).

The tutorial is currently sized as a large one-step smoke case for the 6 GB
Quadro RTX 3000 in this workstation. It is intentionally kept a little below
the estimated device-memory ceiling so both the solver and the Python STL
preprocessing remain stable under WSL.
