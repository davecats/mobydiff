# Next session: an immersed body that reaches an outlet face

STATUS: **B0-B3 DONE 2026-10-05 -- and the gate they did not measure
FAILED the same evening: the body on the outlet plane is UNSTABLE in a
turbulent trough-cut geometry.** The rough-wall boundary layer on the
ELECTED CaNS-grid case (`tutorials/turbulentBoundaryLayer/
production_stats_rough.ini`, MacDonald et al. 2016 k+ 10 / lambda+ 113
egg-carton, crests 0.77 = 50 cells high) blew up within 2000 steps from the
developed field with the roughness running through the outlet, and a
4 M-cell cut-down of it at t = 52 with skew AND divergence convection; the
same cut-down with the roughness ramped out before the outlet runs 300
t.u. The onset is in the last column at the troughs: outlet face fluid or
graded, the face behind it solid (the `r ~ 1e27` limit of B1, gated only
in laminar Poiseuille). Write-up with numbers:
`tutorials/turbulentBoundaryLayer/tests_record.md` section 14; reproducer
`~/outlet_runs/rough_cans/cut_skew.ini` (10 min on one A6000, field dumps
around the onset beside it). **THE NEXT SESSION ON THIS TOPIC STARTS
THERE: stage-by-stage, the last fluid cell whose upstream face is solid,
with the stored pressure and the projection's dead-end diagonal in the
picture.** The production rough case meanwhile ends its roughness at
x = 580 (`wall_x_end`); the HoreKa pair runs in `rough_tbl_run`.** B0 is the section "B0.
Investigation write-up" below; what was built and every gate number is in
"Implementation (2026-10-05)". The handout as written on 2026-10-02 follows,
unchanged. The order was the one the outlet work used and the user asked for
again: **investigate first, write down what each stage does, choose the
simplest design that is correct on a low and on a high face, then implement
and gate.**

CONTEXT FROM THE USER (2026-10-05): a NEW `turbulentBoundaryLayer` tutorial is
planned that compares mobydiff directly with other codes run at the same
parameters. It is not part of this session, but the rough-wall boundary layer
is the production case behind the body-at-outlet work, and `rough_jacobi.ini`
is the configuration to keep running.

## How to start

1. Read this file, then `docs/next_session_outlet.md` ("Implementation" and
   "Follow-up, 2026-10-02": the outlet face as it is, every gate and where its
   numbers are), then CLAUDE.md's bullets "The outlet face is predicted" and
   "Outlet follow-up". The investigation write-up of the outlet
   (`next_session_outlet.md`, "One substage, stage by stage") is the model
   for the write-up this session owes.
2. Build all four variants (`./compile.sh cpu gpu cpu_nofma gpu_nofma`,
   nvhpc 25.9 loaded in the same shell as every `mpirun`). The reference set
   is `~/cfl_ref_binaries` (`cb47d35`; the head differs from it only by
   docs, data and this handout).
3. Reproduce the refusal on the cut-down rough-wall case (recipe in
   "Gates") before changing anything: that case is the smoke test of the
   whole session.
4. Increments, each gated before the next:
   - **B0, the investigation**: the stage-by-stage table for a face with a
     coefficient on it and on its neighbour, low and high, written into this
     file; the nine questions below answered; the design confirmed or
     changed. No code.
   - **B1, the kernel**: `predict_outlet_faces` in the one-statement form,
     `check_outlet_faces_fluid` deleted, `init_outlet_faces` reconsidered.
     Gate: every outlet-free and outlet case `max_abs 0` against
     `~/cfl_ref_binaries` (nofma AND production flags: the argument is
     exactness, so measure it both ways), then the body gates.
   - **B2, the inlet inside a body** (question 7): decide, implement, gate
     with the rough-wall case.
   - **B3**: docs (`docs/configuration.md` patch types, the tutorial README
     of the sailplane if it is affected), CLAUDE.md, `validation/README.md`,
     this file's STATUS, a reference set cut from the commit.

## Why this is the next thing

Since the outlet face is predicted (2026-10-01, `949148f`) the solver STOPS
when an immersed body touches an outlet face:

    error: the immersed body reaches the outlet face 2
    ERROR STOP [boundary] outlet: an immersed body touches the outlet face

(`boundary.f90 check_outlet_faces_fluid`, called from `moby_solve.f90`.) That
was a deliberate stop, not an oversight: the outlet predictor borrows its
neighbour's increment as it stands, which is right only where neither face is
penalized. But it took a capability away. **Every rough-wall or immersed-wall
boundary layer has its wall crossing the outlet**: the campaign config
`overheadTest/horeka/configs/rough_jacobi.ini` (egg-carton roughness on an
immersed wall, inlet at x_min, outlets at x_max and the top) no longer starts,
and ran before (measured 2026-10-02 on a 256 x 44 x 48 cut-down of it:
`154e48f` runs, the present head stops with the message above). The same holds
for a plate or sting leaving through the outlet, and for `les_ibm`-type wall
slabs in an inflow/outflow channel.

## What is known (read in the code, 2026-10-02; confirm before relying on it)

- **The projection is already written for it.** The outlet-face correction
  multiplies by the face's own `ibm%mu` on both sides
  (`pressure_solver.f90 jacobi_apply`: a low outlet face through the ordinary
  low-face statement, `cf*ibm%mu(i,...)` with `cfLow(1)`; a high one through
  the plane branch, `cf*ibm%mu(ip,...)` with `cfHigh`; the red-black sweep
  reads the same two factors); the diagonal is `(dnLow*mu + dnHigh*mu)*d1P`. `ibm%mu` and `ibm%coef` are
  ghost-inclusive (`update_ibm_mu` loops `lbound..ubound`), so the HIGH
  outlet face, index `nb+1`, has its own coefficient and factor. The normal
  velocity at `nb+1` sits ON the boundary plane (`slice_grid_direction`), so
  its coefficient needs no geometry beyond the domain.
- **The predictor is the one stage that assumes a fluid face.**
  `step.f90 predict_outlet_faces`:
  `q(f) = qs(n) + ((q(f) - q(n)) + dt_gamma (dp(n) - dp(f)))`,
  with `(f, n) = (1, 2)` or `(nb+1, nb)`. `qs(n) - q(n)` is the neighbour's
  increment only when the neighbour is not penalized.
- **What `qs(n)` holds at that point.** The fused kernel writes
  `qs = (q + I)*mu` with `I = dt_alpha N + dt_beta N_old - dt_gamma dp` and
  `mu = incr(x)`, `x = coef*dt_gamma` (`ibm.f90 penal_incr_factor`);
  `add_penalization_state_correction` adds `(state - incr) q` on the body
  blocks; the SGS and body-force passes add `incr`-scaled increments; the
  band filter rewrites `qs` in the near-body band. So, band filter aside,
  `qs(n) = state_n q(n) + incr_n I_n` when the outlet kernel runs, and
  `oldrhs(n)` already holds the NEW `N`, the old one is gone: the
  neighbour's unpenalized increment can only be recovered from `qs`.
- **`init_outlet_faces`** gives a face the run was not handed its neighbour's
  value, once. Restart planes (`un_xmax` ...) carry whatever the face held.
- **Not looked at:** an INLET face inside a body. A Dirichlet inlet is a
  pinned face with a prescribed value; if that value is not zero inside the
  solid (a Blasius profile measured from the domain wall, an immersed wall
  standing above it) the first solid cell has a mass flux the penalized
  faces cannot pass on. The rough-wall case has this at x_min: the mean
  plane of the egg-carton wall is y = 0.01 (`ibm.f90 wavy_wall_height`,
  `y_offset`), the Blasius inlet profile is measured from y = 0, so the inlet
  faces of the first row(s) are inside the solid with u > 0. It ran before,
  which says nothing about whether it was right.

## The proposed strategy (to be confirmed by the investigation, not assumed)

**The outlet face keeps its equation and gains its own penalization:**

    du_f/dt = N_n - (dp/dn)_f - lambda_f u_f

`N_n` is the non-pressure, NON-PENALIZATION right-hand side of the
neighbouring face (the zero-normal-gradient closure the outlet already
uses); the pressure gradient and the penalization are the face's own. In the
two-factor form of the scheme, with `x = lambda dt_gamma`:

    I_n   = (qs(n) - state_n q(n)) / incr_n            the neighbour's unpenalized increment
    q(f) := state_f q(f) + incr_f (I_n + dt_gamma (dp(n) - dp(f)))

Written as ONE statement, with `r = incr_f/incr_n`:

    q(f) = r qs(n) + ((state_f q(f) - r state_n q(n)) + incr_f dt_gamma (dp(n) - dp(f)))

- **No branch, and bit-identical without a body.** At `lambda = 0` the
  factors are exactly 1 (`penal_*` return 1.0 at x = 0) and the statement is
  today's, operation for operation; a multiplication by 1.0 is exact, with or
  without FMA contraction. Every outlet gate of today must come out
  `max_abs 0`, production flags included. That is the first gate.
- **Low and high alike**: the side still enters through `(f, n)` only.
- **The factors are there**: `incr = ibm%mu(f)`, `ibm%mu(n)` (refreshed each
  substage, ghost-inclusive); `state = penal_state_factor(coef*dt_gamma)`.
  The kernel additionally maps `ibm%coef` and `ibm%mu`.
- **Limits to check on paper and in a test**: `n` solid and `f` fluid
  (`r ~ x_n ~ 1e25` times `qs(n) ~ I_n/x_n`: finite, but look at it); both
  solid (the face goes to zero through `state_f`, as an interior solid face
  does); `f` solid and `n` fluid (a body starting exactly at the plane).
  `incr > 0` always (P3 is positive for x >= 0), so the division is safe.
- **`check_outlet_faces_fluid` is deleted**, not relaxed.

What the investigation has to settle before this is built:

1. Walk one substage and init / restart for a face with `coef(f) /= 0` and
   for one with `coef(n) /= 0`, low and high: who reads and writes what
   (the table of `docs/next_session_outlet.md`, "One substage, stage by
   stage", with a body in it). In particular the band filter (`qs` is not
   `state q + incr I` in the band: is the recovered `I_n` still the right
   thing to borrow?) and the SGS / body-force passes.
2. Is "zero normal gradient of `N`" acceptable when the body differs between
   `n` and `f` (roughness that varies along x)? It is the closure of the
   outlet as it stands; say what it costs, do not replace it silently.
3. The projection: confirm, with a test and not by reading, that the
   denominator and the correction use the same `mu` of the outlet face on
   both sides and in all three paths (Jacobi / Chebyshev, `interface_correct`,
   red-black).
4. The pressure of a solid cell on the outlet plane: the Dirichlet ghost is
   written there too. Harmless if the solid pressure is decoupled (`mu ~ 0`
   on all faces); show it.
5. Removed blocks (`remove_solid`, buried under `refine_body`) next to an
   outlet: a removed block has no boundary points and its neighbours' faces
   toward it are `FACE_CLOSED`. Is anything else needed?
6. Geometry beyond the plane. The normal face needs none. Tangential
   velocities, pressure and `dwall` ghosts sit at mirrored coordinates
   OUTSIDE the domain: an analytic indicator is defined there; an STL must
   extend past the plane or the body appears to end. Decide what the solver
   requires and say so in the config check or the docs.
7. The inlet inside a body (above): what value does a pinned inlet face
   inside the solid hold, and what should it? If the answer is "zero where
   the face coefficient says solid", decide which stage owns that (the
   boundary row's constant, at init) and make it the same on a low and a
   high inlet.
8. Scalars (`apply_scalar_bc_q` copy rows, the Dirichlet-penalized and the
   conjugate body at an outlet), the RANS wall-cell classification and
   `dwall` at the last cell, the statistics of the boundary-layer case.
9. A level jump on the outlet next to the body (`refine_body` reaching the
   outlet): the tangential-velocity edge halo there is read by the momentum
   stencil, and the coarse side of that halo is first order
   (`validation/channel_interface/les/README.md`, last section). Measure
   whether it matters here before deciding anything about it.

Alternatives, and why they are not the proposal:

- *Skip the borrowed increment where the neighbour is penalized* (copy the
  face, or freeze it): a second rule for a subset of faces, and it brings the
  reset back exactly where the body is.
- *Evaluate the face's own right-hand side*: needs the momentum stencil for
  one plane outside the fused kernel's range on the high side; rejected for
  the outlet itself for the same reason.
- *Keep the stop and ask for a fluid gap before the outlet*: not possible for
  a wall.

## Gates

- Body-free: the 7- and 9-case suites, `run_bitexact_outlet.sh` and
  `validation/freestream/run_gates.sh` unchanged, `max_abs 0` against
  `~/cfl_ref_binaries` (nofma), CPU and GPU. The one-statement form makes
  this exact by construction; measure it anyway.
- A zero-coefficient twin: a body whose case file has `coef = 0` crossing the
  outlet equals the body-free run at `max_abs 0`.
- **Poiseuille between immersed walls with inflow and outflow**: the
  `les_ibm` wall slabs (or the analytic plane wall) in a `pois_io`-type
  channel, the walls crossing BOTH x faces. The profile between the walls
  must be the one of the periodic twin with the same immersed walls, the
  pressure linear and equal to the held outlet value at the face (last-cell
  p = G dx/2 in the fluid), the velocity inside the solid at round-off.
- The mirrored pair of it (outlet at x_min): equal to the x_max image to
  round-off; the y and z directions under permutation.
- Restart at mid-length equals the continuous run (`run_restart_outlet.sh`
  with the body case added).
- A body that ENDS one cell before the plane and one that STARTS at it (the
  three limits above): stable, the mass balance closed to the projection
  residual, the solid faces at round-off.
- **The rough-wall boundary layer**: the cut-down `rough_jacobi` (256 x 44 x
  48, recipe below) runs; then the momentum-integral balance to the outlet
  (`assets/postpro/momentum_integral.py`) on a developed run, against the
  smooth-wall result of 2026-10-02 (balance within ~20 % up to 2 units from
  the outlet).
- 1 rank == 4 ranks, CPU == GPU.

The cut-down case used today (it stops at init on the present head):

    sed -e 's/^nx = .*/nx = 256/' -e 's/^ny = .*/ny = 44/' -e 's/^nz = .*/nz = 48/' \
        -e 's/^lx = .*/lx = 46.875/' -e 's/^lz = .*/lz = 8.0/' \
        -e 's/^n_wave_x = .*/n_wave_x = 3/' -e 's/^n_wave_z = .*/n_wave_z = 1/' \
        tutorials/turbulentBoundaryLayer/overheadTest/horeka/configs/rough_jacobi.ini > rough_small.ini

## Landmines

- The fused predictor kernel is not to be edited (register schedule). The
  outlet kernel is a surface kernel: its cost is not the issue, its extra
  maps are (`ibm%coef`, `ibm%mu` are resident; do not map by value).
- nvfortran 25.9 drops the assignment of a variable used only as a map
  bound (CLAUDE.md, coding conventions): keep looping over all boundary
  points.
- Since 2026-10-02 `cflmax` and `pecletmax` bound the SUM over directions.
  A reference binary older than `~/cfl_ref_binaries` differs BY DESIGN on
  every case whose step is set by the Courant or the diffusion limit; gate
  against that set, or with `dtmax` binding on both sides.
- The conjugate gate drivers rewrite tracked record files
  (`validation/conjugate/*.dat`, `*.heat.txt`): `git status` before
  `git add -A`.
- A solver-minted template carries the outlet plane of its own one-step
  run; an IC tool must write the plane or delete it.
- Everything in `docs/next_session_outlet.md` "Landmines" holds.

## State of the tree (2026-10-05)

- `main` is pushed through `11760d4`. On top of `949148f` (the outlet face
  predicted): `61a6ca4` the outlet follow-up, `f39c139` `pecletmax` on the
  sum over directions, `ced8e91` `cflmax` on the sum, `cb47d35` this handout,
  `11760d4` the regenerated boundary-layer tutorial data.
- Reference sets: **`~/cfl_ref_binaries`** (`cb47d35`, all four builds,
  PROVENANCE inside) is the current one. `~/limiter_ref_binaries` is before
  the `cflmax` change (`beltrami_yslab` differs), `~/outlet_runs/ref_949148f`
  before both limiter changes (six diffusion-bound cases differ).
- The boundary-layer tutorial: `data.nc` and the figures are the 10000-t.u.
  statistics of HoreKa job 5175523; `tutorials/turbulentBoundaryLayer/`
  holds `production_stats.h5` and `restart_field.h5` (the step-1587500
  field, 5.6 GB, untracked). On HoreKa, `tbl_stats_run` (11 snapshots) and
  `outlet_tbl_run` (6) under `optimiseBlockRefinement/` can be deleted;
  the worktrees `moby-outlet-new`, `moby-outlet-ref`, `moby-tbl-stats` too.
- `~/Codes/mobydiff.bl/.../production_stats.h5` is no longer the tutorial's
  old statistics (overwritten by another case on 2026-10-02); the old data
  live in git history only.

## B0. Investigation write-up (2026-10-05, main `f54cce8`; no code changed)

Read in the code and checked on the cut-down rough-wall case (prepared with
`build_cpu`, 4 ranks: 4 leaves of 64 x 44 x 48, every block a body block).
The refusal reproduces on the head (`error: the immersed body reaches the
outlet face 2`). What the cut-down's geometry actually is: `dy_1 = 0.19`
with the roughness `0.01 +- 0.05`, so on BOTH x planes only the ghost row
`j = 0` is solid-centred (`coef = SOLID/Re = 2.2e27`) and row `j = 1` is
graded; no interior face of the plane is inside the solid. It exercises the
graded limit, not the solid one -- a second variant with `amp_x = 0.5`
(crest-to-trough ~5 cells) is used below for the solid faces.

### One substage, stage by stage, with a body on the plane

`f` the outlet face, `n` its interior neighbour, `(f, n) = (1, 2)` low /
`(nb+1, nb)` high; `x_i = dt_gamma coef_i`, `incr_i = ibm%mu(i)`,
`state_i = penal_state_factor(x_i)`. Only what DIFFERS from the body-free
table of `docs/next_session_outlet.md` is listed; everything else is as
there.

| stage | LOW outlet, `f = 1` | HIGH outlet, `f = nb+1` |
|---|---|---|
| `update_ibm_mu` | `mu` refreshed over the ghost-inclusive range `0..nb+1` of every BODY block (a block whose only coefficient sits on the `nb+1` plane IS a body block: `ibm_body_blocks` scans the ghost-inclusive range). Non-body blocks keep `mu = 1` from allocation and `coef = 0`. So `incr_f`, `incr_n`, `coef_f`, `coef_n` are all available at the face on both sides | same |
| fused predictor | face 2: `qs(2) = (q(2) + I_2) incr_2`, reading `q(1)` (old outlet value) in its stencil | face `nb`: the same with `q(nb+1)` |
| `add_penalization_state_correction` | body blocks, under the start masks: face 2 gets `(state_2 - incr_2) q(2)`; face 1 is outside the mask | face `nb` the same; `nb+1` out of range |
| SGS, body force | `qs += dt_alpha term * incr` at the predicted faces | same |
| band filter | `qs += theta/4 del * incr` over the band list, which is interior-only (`1..nb`); face 1 can be in it (harmless: the commit skips it and the outlet kernel reads `qs(2)`, not `qs(1)`) | `nb+1` never in the list |
| **so, when the outlet kernel runs** | **`qs(n) = state_n q(n) + incr_n I_n` EXACTLY, with `I_n` the complete non-pressure increment minus `dt_gamma dp(n)`: every correction pass scales by `incr` (the band filter too -- the handout's worry is moot)** | same |
| `predict_outlet_faces` today | `q(1) = qs(2) + ((q(1) - q(2)) + dt_gamma (dp(2) - dp(1)))`: borrows `qs(2) - q(2)`, which is `(state_2 - 1) q(2) + incr_2 I_2`, NOT the increment, once `coef_2 /= 0`; and `q(1)` is never penalized | same with `(nb+1, nb)` |
| commit `qs -> q` | face 1 skipped (pinned mask) | `nb+1` out of range |
| exchange | unchanged | unchanged |
| projection: diagonal | `(dnLow(1) mu(1) + dnHigh(1) mu(2)) d1P`, `dnLow(1) = 2 d1f` on an outlet: the face's OWN `mu(1)` | `dnHigh(nb) mu(nb+1)`: the face's own `mu(nb+1)` (read at `ip = nb+1`, ghost-inclusive) |
| projection: correction | `q(1) += (phi(0) - phi(1)) cfLow(1) mu(1)` | `q(nb+1) += (phi(nb) - phi(nb+1)) cfHigh mu(nb+1)` (the `i == nx` plane branch) |
| red-black sweep | `gLo1 = face_grad_denom(.., d1x(1), outLow)`, `cLo1` with `mu_u_i`: the same pair in place | `gHi1 = face_grad_denom(.., d1x(ip), outHigh)` with `mu_u_ip` |
| final `apply_bc` | pressure ghost `p(0) = 2 p_out - p(1)` also in a solid cell; tangential ghosts copy the (penalized, ~0) interior | same |
| init | `init_outlet_faces`: `q(f) = q(n)`; a cold start has `u_inf` INSIDE the body too, and the first substage's `state` factor removes it from every solid DOF, the outlet face included once it is penalized | same |
| restart | low face in `un`; high in `un_xmax` | same |

### Answers to the nine questions

1. **The table above.** The predictor's stage is the only one that assumes
   a fluid face; every correction pass preserves `qs(n) = state_n q(n) +
   incr_n I_n`, so the neighbour's unpenalized increment IS recoverable from
   `qs` by one division, and the band filter is inside `I_n` like the SGS
   term (it is applied `x mu`). Nothing else on the list reads the outlet
   face differently with a body.
2. **The zero-normal-gradient closure with a body that differs between `n`
   and `f`.** The face inherits `N_n` (convection + diffusion of its
   neighbour, built from the neighbour's stencil) -- the neighbour's
   PRESSURE gradient cancels algebraically (`I_n` carries `-dt_gamma dp(n)`,
   the face adds `+dt_gamma dp(n)`), so a solid neighbour's decoupled
   pressure never reaches a fluid face. What is left when `n` is solid and
   `f` fluid: `N_n` has `q(n) ~ 0` in its diffusion stencil, i.e. the face
   sees `nu (q(f) - 0)/dx^2`, a wall-like sink -- the cost of a body that
   ends one cell before the plane, O(1) in that one face column. When `f`
   is solid, `N_n` is multiplied by `incr_f ~ 1/x_f` and is irrelevant.
   Roughness varying along x: the face's `N` is the neighbour's, an O(dx)
   error on the plane -- exactly the outlet's existing closure. Kept; not
   replaced.
3. **The projection uses the face's own `mu` on both sides in all three
   paths** -- read in `compute_rdenom` (`mu_u_ip` at `ip = nb+1`),
   `jacobi_apply` (`ibm%mu(ip,..)` under the `i == nx` plane branch,
   `ibm%mu(i,..)` with `cfLow(1)`), `redblack_sweep` (`mu_u_ip` with `gHi1`
   / `cHi1`), and `interface_correct` carries no outlet branch since O5.
   Tested below, not only read: the Poiseuille-between-immersed-walls gate
   (solid faces on the plane at round-off after the projection, fluid
   divergence at the projection residual) and the mirrored pair.
4. **The pressure ghost in a solid cell on the plane.** Harmless: every
   face of a solid cell has `mu ~ 1/x ~ 1e-27`; the diagonal is `~1e-27 d^2`
   and the divergence source scales the same way (the solid face velocities
   are `incr * I ~ 1e-27`), so `phi` stays O(dt dx) and the correction
   `cf mu dphi ~ 1e-27`. The outlet ghost enters only through
   `dnHigh mu(f)`, i.e. not at all. The solid-cell pressure on the plane is
   the same decoupled Brinkman pressure as anywhere in the body. Measured
   below (max |p| in solid cells bounded).
5. **Removed blocks next to an outlet.** A removed block has no boundary
   points, no exchange entries toward it, and its neighbours' faces toward
   it are `FACE_CLOSED`. The outlet face of a SURVIVING block is `nb+1` in
   the NORMAL direction, which has nothing to do with a tangential
   neighbour's absence; the tangential edge halos of the outlet row come
   from the same-level extension that keys on "combined neighbour absent"
   (Phase 2) and is unchanged. Nothing needed.
6. **Geometry beyond the plane.** The normal face sits ON the plane
   (`slice_grid_direction`): no geometry outside the domain. Tangential
   velocity, pressure and `dwall` ghosts sit half a cell outside; the
   analytic indicator is defined there (the egg-carton is a function of
   `(x, z)`); an STL must extend past the plane -- `moby_prepare` samples
   the ghost tiles from the mesh, so a mesh that ends at the plane gives a
   fluid ghost and the body appears to end half a cell before the boundary
   (a tangential ghost mirrored against a fluid value; the normal face
   itself is unaffected). Required of the user, stated in the docs
   (`docs/configuration.md`); not checkable by the solver (a ghost tile
   with `coef = 0` is a legal fluid ghost).
7. **The inlet inside a body.** A Dirichlet face is pinned to its datum
   whatever the coefficient says: no stage penalizes it (not in any
   predicted range; `cfLow = dnLow = 0` on a non-outlet physical face). So
   today a Blasius inlet face in the solid injects `u_in(y) dy` into a
   solid cell whose other faces have `mu ~ 1e-27`; the projection absorbs
   it with `phi ~ rdenom div ~ 1e27 div` and a correction `mu cf dphi =
   O(div dx)`, i.e. the solid passes the flux on through penalized faces at
   O(u_in), with a stored pressure that grows by `phi/dt_gamma` -- a conduit
   through the body. (The cut-down case has no such face: its inlet faces
   in row 1 are graded, fluid-centred; the production grid, `dyw_plus`
   0.15 on 176 cells, very likely has one or two.) DECISION: the datum of a
   Dirichlet VELOCITY row is ZERO where the row's own staggered location is
   solid-centred (`coef >= 1e20`, the marker threshold rans/scalar already
   use): the normal row at the face itself, a tangential row at its ghost
   position (the mirror `2v - q` with `v = 0` makes the plane no-slip
   there). Graded faces keep the datum (fluid-centred; the same staircase
   every cut cell has). Owner: the boundary row's constant `bcC`, masked
   ONCE at init from the host coefficients, in the slot
   `check_outlet_faces_fluid` leaves (between `read_ibm_coeff_file` and
   `enter_ibm_data`), with `enter_boundary_data` moved below it -- nothing
   between its present position and the first `apply_bc` touches `bc` on
   the device. Low and high alike (the row loop already resolves `gi` per
   side); every Dirichlet velocity face, not only declared inlets (a wall's
   datum is 0 already, a moving wall's would be wrong inside a body).
   Restart: the rows are rebuilt from the ini each run. Not: a time-
   dependent relaxation of the datum (invents a scale), a solver check that
   refuses it (it is a legitimate configuration -- the rough wall).
8. **Scalars, RANS, statistics.** The scalar outlet row is a Neumann copy:
   in a solid cell it copies the penalized body value (Dirichlet mode) or
   the conjugate solid temperature (zero-gradient = adiabatic at the domain
   boundary), nothing to add; `apply_scalar_bc_q` no longer exists (step 6:
   one `apply_bc` with scalar columns). RANS `wallcell`/`dwall` at the last
   cell come from the ghost-inclusive tiles that prepare writes from
   geometry: unchanged by the outlet. The momentum-integral statistics of
   the boundary-layer case are a post-processing matter (`momentum_integral.py`
   measures the balance to the outlet; the roughness is a form-drag term it
   does not contain, so for the rough case only the TREND to the outlet is
   comparable with the smooth result).
9. **A level jump on the outlet next to the body.** Not measured in this
   session: the production case is single-level, and the first-order
   restriction ghost row (6.8 % in an LES wall cell) is a known, documented
   residual of the exchange, not of the outlet. The `refined` freestream
   group (level jump on the outlet, from rest) stays the gate for the
   exchange side. Left as it was.

### Design, confirmed

The proposed one-statement form stands. With `r = incr_f/incr_n`, written
in the kernel as four locals so that every FMA contraction the compiler can
choose reproduces today's arithmetic when the factors are exactly 1:

    r   = mu_f/mu_n                      (ibm%mu: refreshed this substage)
    sf  = penal_state_factor(x_f)
    rsn = r*penal_state_factor(x_n)
    dtf = mu_f*dt_gamma
    q(f) = r*qs(n) + ((sf*q(f) - rsn*q(n)) + dtf*(dp(n) - dp(f)))

Body-free: `r = sf = rsn = 1`, `dtf = dt_gamma`, each product by 1.0 exact,
and every contraction candidate (`fma(r, qs(n), T)`, `fma(sf, q(f), -rsn
q(n))`, `fma(dtf, dpn - dpf, S)`) rounds exactly as today's
`fma(dt_gamma, dpn - dpf, S)` / `qs(n) + T` -- hence `max_abs 0` at
PRODUCTION flags is expected, not only at nofma. Limits: `f` solid, `n`
fluid: `r ~ 1e-27`, `q(f) -> state_f q(f) + 1e-27 (..)`, the face dies like
any solid DOF; `n` solid, `f` fluid: `r ~ 1e27`, `r qs(n) = I_n` to round-
off, `rsn q(n) ~ 6e-54 q(n)`; both solid: `r ~ 1`, `sf ~ 1e-81`. `incr > 0`
for `x >= 0` (P3 positive), `x^3 ~ 1e81` far from overflow. The kernel
gains `ibm` (maps `ibm%coef`, `ibm%mu`, both resident).
`init_outlet_faces` stays as it is (the face is state; in a solid it is
penalized away in the first substage like every other solid DOF);
`check_outlet_faces_fluid` is deleted.

## Implementation (2026-10-05, on top of main `f54cce8`)

Reference binaries: `~/cfl_ref_binaries` (`cb47d35`). Gate drivers: the 7-
and 9-case suites, `validation/freestream/run_bitexact_outlet.sh` (1 and 4
ranks), and the new `validation/body_outlet/run_gates.sh` (README there has
every number).

### B1. The outlet face carries its own penalization

`step.f90 predict_outlet_faces(blk, bc, ibm, dt_gamma)`, the one-statement
form of the design (`r = mu_f/mu_n`, `sf`, `rsn = r state_n`, `dtf = mu_f
dt_gamma`, all exactly 1 / dt_gamma without a body); `boundary.f90
check_outlet_faces_fluid` deleted with its call; `init_outlet_faces`
unchanged (the face is state; a solid one is penalized away in the first
substage). The kernel maps `ibm%coef` and `ibm%mu` (resident).

- Body-free, nofma: the 7-case suite (4 ranks), the 9-case suite, the
  outlet suite (1 and 4 ranks) `max_abs 0` against the reference, CPU AND
  GPU (30 + 23 comparisons).
- Body-free, PRODUCTION flags: CPU 30/30 `max_abs 0`; GPU 22/23 -- the 7-
  and 9-case suites and six outlet cases at 0, `blasius2d` at 1.4e-15.
  Run down: the first difference is ONE ulp (8.7e-19 on 4.8e-3) at 2 of
  384 columns of the TOP outlet (the stretched y line) at step 37, the x
  outlet of the same run exact through step 50; a form with the scaled
  operands in statements of their own gives bit-identical output. So the
  GPU compiler contracts the modified statement differently for rare
  operands -- the FMA class, and the reason CLAUDE.md gates an expression
  change at nofma. The B0 claim "exact at production flags too" holds on
  the CPU and not on the GPU; the kernel's header says so.
- The cut-down rough-wall case runs (4 ranks CPU, 400 steps, L2 divergence
  6e-8); a tall-roughness variant (`amp_x = 0.5`, phase pi/2, so rows 1..2
  are solid on both planes) runs 200 steps stable.

### B2. Dirichlet velocity data inside the body

`boundary.f90 mask_dirichlet_velocity_in_solid`: every Dirichlet velocity
row (the normal row at its face, a tangential row at its ghost) whose
location is solid-centred (`abs(coef) > SOLID_FACE_THRESHOLD`) gets `bcC =
0`, once at init on the host coefficients, between `read_ibm_coeff_file` and
the maps; `enter_boundary_data` moved below it (nothing in between touches
`bc` on the device). The solver prints the number of NON-ZERO data it zeroed
(a wall row inside the body is 0 already and would swamp the count: the
flat cut-down read 32828 of 43200 rows "inside" before the count was
restricted). `SOLID_FACE_THRESHOLD` now has ONE definition (`init.f90`),
replacing five identical local copies (les, rans, scalar x2, scalar_stats).
Body-free and bodies without a Dirichlet face: nothing changes (the B2
binaries reproduce every B1 gate above).

### Gates with a body (`validation/body_outlet/`, all PASS)

Poiseuille between two immersed slabs (0.3 dy off a face, gap 1.0) in an
inflow/outflow channel, uniform inlet through the slabs, `nb = 4` with 32 of
192 blocks removed (some at the outlet): outlet-face profile = the periodic
twin's to **1.35e-3** of 1.5 (0.09 %, the bulk ratio 0.9991), last-cell p =
G dx/2 to **0.08 %**, gradient 1.19987 vs 1.2, the 32 solid faces ON the
outlet plane **1.6e-28**, interior solid DOFs <= 9.4e-28, flux in = out to
2e-16; mirrored (outlet on the LOW face) = the image to **3.3e-15**;
stream along y = the x case at **0.0**, along z to 2.4e-15; the body
ENDING one cell before the plane (`r ~ 1e27`) and STARTING at it
(`r ~ 1e-27`) both stable 20000 steps with the solid faces at 3.4e-28 /
9.8e-28 and the flux closed; 4 == 1 ranks, restart == continuous, CPU ==
GPU (nofma pair, same case file), zeroed-coefficient twin == body-free,
each at **max_abs 0** including the `un_xmax` plane.

### The rough-wall boundary layer

A production-shaped cut-down (the production y line; 320 x 176 x 24 on 60 x
100 x 8) and its smooth twin, 300 t.u. each on istmcetus: both run, dtheta/dx
holds to the last half unit on the smooth twin (1.20 -> 1.15e-3; the old
outlet collapsed it to 5 %) and on the rough one stays at its upstream level
through the roughness phase (last unit +17..26 % over the same phase one and
two wavelengths upstream, against the smooth twin's +4..12 % in the same
band -- the layer's own modulation), theta growing uniformly through the
last six cells. **NOT measured: the turbulent balance** -- both layers are
still laminar at t = 300 (H 2.6; the trip does not transition them in 60
units). That needs the production domain on HoreKa, with
`rough_jacobi.ini` itself; it is the one open item of this handout.

### Not done, by decision

Question 9 (a level jump on the outlet next to the body): not measured; the
production case is single-level and the first-order restriction ghost is a
documented exchange residual (`validation/channel_interface/les/README.md`).

## After this: where the code and numerics review goes next

The outlet was item 3 of `docs/next_session_after_step9.md`. With it and the
body-at-outlet case closed, the review (`docs/numerics_review_2026-09-26.md`)
has these left, in the order suggested:

1. **The remaining simplifications of its section 6** (performance-neutral,
   each bit-exact, one session together). Still open, checked in the code
   today:
   - item 6, *split the scalar kernel by mode*: `scalar_transport_kernel`
     (`scalar.f90`, 600 lines, ~60 privates) carries the plain, turbulent,
     penalized, wall-function, conjugate and banded branches in one body. A
     plain kernel and a conjugate kernel, expressions moved verbatim: easier
     to read, and the register lesson of `jacobi_apply` says the dormant
     branches cost occupancy. Gate: the 9-case suite and the conjugate gates
     at nofma, then ncu on both kernels.
   - item 5, the three `interface_correct` kernels into one launch over a
     face list (~1 ms/step at 16 ranks).
   - item 8, `cheb_combine` fused into `jacobi_compute_phi`.
   - item 7 (double-buffering `q` to drop the `qs -> q` copy) stays "measure
     first".
2. **Step 10, implicit solid conduction for conjugate scalars** (section 10;
   two to three sessions; depends on nothing open). The one numerics step
   with a measured customer: the CHT pipe at `solid_rhocp = 0.0625` steps at
   the solid's limit. NOTE what changed today: the conjugate rate is now half
   the Gershgorin diagonal at every cell, so the bulk solid is bounded as
   strictly as the cut cells were. A conjugate run is slower than before
   unless `pecletmax` is raised toward 0.6, which makes step 10 worth more.
3. **The one-by-one re-measurement pass** (`validation/README.md` is the
   checklist; the long tier is still at niter-6 values). It now has a second
   reason: every recorded number of a run bound by the Courant or diffusion
   limit was taken with a larger step. Discharged so far: turb180, Blasius,
   the cylinder Re 40 drag, the `les_ibm` developed statistics, the
   boundary-layer tutorial (regenerated in full, 2026-10-04).
4. Step 11 (line-implicit wall-normal diffusion) only if a viscous-limited
   production case appears; step 12 (second halo layer, bounded scalar
   convection) last.

Smaller items left by the outlet sessions, none blocking:

- **`cflmax = 0.8` in every ini is now 46 % of the limit** (it was the
  per-component number: effectively 1.0 .. 1.2 on the sum in a turbulent
  flow, 2.4 in the worst case, against sqrt(3)). A Courant-bound run takes
  1.2 .. 1.5 times more steps than before at the same key. The user's call
  per case: 1.2 restores the old step with a 30 % margin. NOT changed in the
  inis.
- The two limits are bounded separately, and RK3's stability region is not a
  rectangle: at `pecletmax` 0.5 the convective margin shrinks (z = -2 + 1.2 i
  has |R| = 0.99). Worth a line in the docs or a combined bound; not
  measured.
- The first-order ghost row of a restriction at a level jump on a physical
  face (6.8 % in the LES wall cell of the coarse side): unequal gather
  weights (3/2, -1/2), or the boundary row applied to the halo column after
  the exchange.
- The control rows of the C3 transient gate use a binary of August that no
  longer reads today's inputs (`validation/conjugate/run_gates_c3.sh`,
  `~/c2_ref_binaries`): repair or retire the control.
- `rans.f90` warns that "the Peclet limiter is molecular-only": check it
  against `get_timestep_rates`, which adds the eddy viscosity whenever a
  turbulence model is on.
- `tutorials/naca/rans` from scratch (days of GPU time).
