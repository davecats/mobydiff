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

## Boundary conditions (changed 2026-09-27, numerics review F4)

`x_min` is a uniform inflow (`u = 1`, Neumann p); `x_max` is a declared
**outlet** (`x_max_patch = outlet`: zero-gradient normal velocity in the
predictor, Dirichlet `p = 0` in the projection, Neumann tangential
velocities); the y and z faces are free-slip (Dirichlet 0 normal, Neumann
tangential). The ini used to spell the outlet as `x_max_{u,v,w}_type =
neumann` + `x_max_p_type = dirichlet`. A Neumann condition on the NORMAL
velocity component is now a config error: `apply_bc` rewrote that face from
the corrected interior inside every projection iteration, which broke the
SPD operator at that row. The change alters the outflow face's arithmetic on
purpose, so the 1-step `validation/prepare/run_gates_big.sh` sailplane leg
(prepared case file vs the committed legacy file, same ini on both sides)
was re-run rather than compared against older snapshots. With the outlet
patch the case also runs on the default damped-Jacobi / Chebyshev-Jacobi
projection (see the step-3 MEASURED entry in
`docs/numerics_review_2026-09-26.md` section 10); `solver = redblack` is kept
in the ini as the faster choice for this outflow case.

The tutorial is currently sized as a large one-step smoke case for the 6 GB
Quadro RTX 3000 in this workstation. It is intentionally kept a little below
the estimated device-memory ceiling so both the solver and the Python STL
preprocessing remain stable under WSL.
