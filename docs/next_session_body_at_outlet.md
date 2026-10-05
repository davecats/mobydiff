# Next session: an immersed body that reaches an outlet face

STATUS: **OPEN -- the next session's task (user, 2026-10-05).** Written
2026-10-02 at the end of the outlet work (`docs/next_session_outlet.md`, DONE
including its follow-up; the boundary-layer tutorial data were regenerated
with the predicted outlet on 2026-10-04 and the user accepted that as the
outlet's validation). Nothing below is implemented. The order is the one the
outlet work used and the user asked for again: **investigate first, write
down what each stage does, choose the simplest design that is correct on a
low and on a high face, then implement and gate.** This is the code
improvement of the next session; section "After this" says where the review
goes afterwards.

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
