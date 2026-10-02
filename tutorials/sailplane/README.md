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

The IBM coefficients live in the CASE FILE `sailplane_case.h5` (`[case] file`,
prepare/solve split, `docs/prepare_solve_strategy.md` + numerics review step
7): grid, `[blocks] nb = 10` block layout (40 x 45 x 10 = 18000 leaves), the
STL declaration (`stl_file`, `stl_scale = 0.001`, `stl_translate`) and the
coefficients themselves. `moby_prepare` writes it (minutes on 2 CPU ranks),
`moby_solve` only reads it and refuses a stale one -- grid, nb, transform or
`re` changed -- with the key named (re-run moby_prepare):

```bash
cd tutorials/sailplane
mpirun -n 2 ../../build_cpu/moby_prepare input.ini    # -> sailplane_case.h5 ([case] file)
mpirun -n 2 ../../build_gpu/moby_solve input.ini
```

The committed legacy `sailplane_ibm_coeff.h5` (the retired mobygeom's global
layout) was retired with the legacy reader at step 7-4; it stays in git
history (last at the step 7-3 commit). A 1-step solve from the prepared case
file was bit-exact against it (validation/prepare/run_gates_big.sh, P1b), and
the retired mobygeom cross-check still reads the grid straight from the case
file (`--grid-file sailplane_case.h5`).

## Boundary conditions (changed 2026-09-27, numerics review F4)

`x_min` is a uniform inflow (`u = 1`, Neumann p); `x_max` is a declared
**outlet** (`x_max_patch = outlet`: the face-normal velocity is an unknown,
advanced by the predictor and corrected by the projection against the held
`p = 0` -- `docs/next_session_outlet.md`; Neumann tangential velocities); the y and z faces are free-slip (Dirichlet 0 normal, Neumann
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
