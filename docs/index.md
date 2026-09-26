# mobydiff documentation

`mobydiff` is an incompressible Navier–Stokes solver: second-order finite differences on a
staggered Cartesian grid, RK3 time stepping, a damped-Jacobi / Chebyshev–Jacobi pressure
projection, block-structured 2:1 local refinement, a volume-penalization immersed boundary
method, optional LES / RANS / IDDES, and passive scalars with conjugate heat transfer. It
runs on CPUs (MPI) and NVIDIA GPUs (OpenMP target offload). The executables are
`moby_prepare` (preprocessor) and `moby_solve` (solver).

## User documentation

- **[Installation](installation.md)** — dependencies, CPU vs GPU toolchains, building.
- **[Running the solver](running.md)** — the run workflow, `moby_prepare`, MPI
  decomposition, output and restart.
- **[Configuration reference](configuration.md)** — every `.ini` section and key.
- **[Numerical methods](numerical-methods.md)** — discretization, time stepping, the
  pressure projection, block refinement and the 2:1 interface, the IBM, turbulence models
  and scalars.
- **[Tutorials](tutorials.md)** — the shipped cases, with walk-throughs of the
  `channel_kmm180` turbulent channel and the `sailplane` immersed-boundary case.
- **[Tools reference](tools.md)** — geometry helpers, verification checks, and
  post-processing/plotting scripts.
- **[Validation & verification](validation.md)** — how correctness is established; the case
  index is [`validation/README.md`](../validation/README.md).
- **[Developer guide](developer-guide.md)** — source layout, the block data model, the GPU
  programming model, and coding conventions.

## Design notes & internal history

The following documents are the design records and engineering notes accumulated while
building the block-refinement and 2:1-interface machinery. They are lower-level and more
historical than the user documentation above, but they are the authoritative record of
*why* the numerics are the way they are.

- [`block_refinement_strategy.md`](block_refinement_strategy.md) — the master design for
  BCM-style equal-size blocks and 2:1 refinement.
- [`interface_projection_derivation.md`](interface_projection_derivation.md) — derivation
  of the SPD 2:1-interface pressure projection.
- [`interface_review.md`](interface_review.md) — review of the interface treatment.
- [`corner_reconstruction_strategy.md`](corner_reconstruction_strategy.md) — edge/corner
  halo reconstruction across refinement interfaces.
- [`prepare_solve_strategy.md`](prepare_solve_strategy.md) — the `moby_prepare` /
  `moby_solve` split and the case-file contract.
- [`conjugate/`](conjugate/) — derivation of the conjugate heat-transfer scheme.
- The `momentum_interface_handout.md`, `jacobi_interface_handout.md` and the
  `next_session_*.md` files are working handouts for individual development sessions.
