# Next session: the outlet face, done properly

NEXT: `docs/next_session_body_at_outlet.md` (an immersed body reaching an
outlet face is refused since this work; rough-wall boundary layers need it).

STATUS: **DONE 2026-10-01, increments O0 .. O5 implemented and gated; the
design was ratified by the user the same day** (restart: store the planes;
the level-jump hole: fixed in this series). What was built, every gate
number and what is still owed are in the last section, "Implementation".
The open items of its "Not done" list were followed up on 2026-10-02 (the
boundary-layer restart pair, the diffusion limit behind `pecletmax`, LES at
a level jump on a wall): "Follow-up, 2026-10-02", near the end.
The sections before it are the handout as written, then the investigation
that preceded the code.

## The decision, in the user's words

> Take 1, make the next session investigate the whole timestep first, to
> find the simplest, most elegant solution to outlet faces at high or low
> indices, without working around errors!

"1" is the root fix: **the outlet face's normal velocity is owned by the
momentum predictor**, not written by a boundary row that is then
compensated. The prototype on the local branch `proto/outflow-incremental`
(`4908e7c`, worktree `~/moby_outlet_proto2`) is NOT to be merged: it keeps
the defective row and pre-loads its constant. It is the REFERENCE the new
implementation is measured against, and a statement of what the right
arithmetic is; it is not the code.

So the order of work is fixed:

1. **Investigate the whole time step first** (section "What to investigate").
   No code before the write-up of what each stage does to an outlet face, at
   a LOW index and at a HIGH index, and where that is wrong or accidental.
   The write-up also answers the two further checks the user asked for:
   whether the low/high indexing itself is the right one, and whether
   boundary-normal velocities should be placed on the boundaries.
2. **Choose the simplest design that is correct on both sides** and say why
   the alternatives lose. Record it here before implementing.
3. Implement, gate, re-measure every outlet case (section "Gates").

A fix that needs another fix to hide it is not acceptable: if a stage is
wrong, that stage changes.

## What is known (2026-10-01)

The defect, from the code. Take an outlet at `x_max`; `u_N` is the last
interior normal-velocity face, `u_out` the outlet face, `p_N` the pressure
of the last cell. Per RK substage, on main:

1. `momentum()` advances `u_N` (old pressure gradient included) and never
   touches `u_out`.
2. `apply_bc(..., outflow_copy = .true.)` RESETS `u_out := u_N`
   (`BCK_FACE_OUTFLOW`, `dst = 1*src + 0`).
3. The projection corrects `u_out` by `2 phi_N d1f` (mirrored `phi` ghost,
   i.e. `phi_face = 0`) and `p += phi/dt_gamma`.

Continuity in the last cell needs `u_out − u_N = −dx (dv/dy + dw/dz)`. The
reset discards it every substage, the projection rebuilds it every substage,
and the INCREMENTAL stored pressure integrates those `phi` in the last cell
column with nothing tying `p_N` to the held outlet pressure. The discrete
outlet is "the outflow must be parallel", stiffer the smaller dt.

History, which is why this is an error and not a design: the original A0
outlet had the predictor NEVER write the outlet face (the comment at
`pressure_solver.f90` ~line 699 still says so). That is consistent with the
incremental pressure, but the face then has no convection and no diffusion,
and the Poiseuille gate converged to a plug profile
(`docs/next_session_airfoil.md`, "AMENDMENT to §THE-gap item 5"). The
amendment added the zero-gradient copy. Neither state gives the face a
momentum equation; the copy merely borrows the neighbour's VALUE.

**The projection's outlet correction is NOT the error**: the half-cell
Dirichlet gradient (`face_grad_corr`, `d1f` against the mirrored ghost)
matches its denominator (`face_grad_denom`, `2*d1f`), and `p += phi/dt_gamma`
is the standard incremental update. Only its comment is stale.

The arithmetic that works (the prototype):

    u_out* = u_out + (u_N* − u_N) + dt_gamma [ dp/dn|face N − dp/dn|outlet face ]

i.e. the face's own momentum predictor with convection / diffusion /
forcing / penalization closed at zero gradient (the neighbour's increment)
and its OWN pressure gradient against the held `p_face = 0`. Without the
pressure terms the last column's pressure is unanchored (measured: a frozen
relic and a drifting outflow imbalance). Numbers, main against prototype:

| | main | prototype |
|---|---|---|
| cylinder Re 100, stored p rms after one time unit, dt 5e-3 / 6.25e-4 | 0.89 / 1.65 | 0.080 / 0.073 |
| velocity, dt 5e-3 against 6.25e-4 at T: max | 0.15 | 6.0e-5 |
| rms dv/dy, last column over x = 15.0 | 0.23 / 0.08 | 1.12 / 1.12 |
| 100 time units, production settings: St / mean C_D / CV C_L amplitude | 0.1670 / 1.4474 / 0.655 | 0.1744 / 1.4383 / 0.393 |
| freestream oblique | exact | exact |
| Poiseuille last-cell p (exact 1.25e-3) / profile dev at x/lx 0.9 | 2.41e-3 / 1.539e-3 | 1.247e-3 / 1.521e-3 |
| Lamb-Oseen reflected fraction | 2.2e-2, regrowing | 9.8e-5, monotone |
| Blasius worst theta / H / du/Ue / dv/v_edge | 1.33 % / 0.45 % / 2.40e-3 / 0.124 | 1.65 % / 0.82 % / 1.48e-3 / 0.013 |

## What to investigate (before any code)

Walk ONE substage and the init / restart path, stage by stage, and write
down for an outlet face at a LOW index and at a HIGH index what is read,
what is written, and by whom. Facts already collected, each to be confirmed
and followed to its consequences:

- **The two sides are not symmetric in storage.** The normal velocity of
  cell `i` lives on its LOW face. A low-side outlet face is index 1: inside
  the predictor's loop range, with an `oldrhs` slot, and IN the field file.
  A high-side outlet face is index `nx+1`: outside the loop `1..nx`, with no
  `oldrhs` slot (`oldrhs(1:nx, ...)`, `blocks.f90`), and NOT in the field
  file. Every test so far (cylinder, Blasius, freestream, the boundary
  layer) has outlets on HIGH faces only. The low side has never run.
- **`momentum_face_start` pins every `FACE_PHYS` low face** (start index 2),
  an outlet included: the predictor skips face 1 of an outlet at `x_min`.
- **The predictor's final `qs → q` copy runs over `1..nx` unconditionally**,
  so it overwrites a pinned low face with whatever `qs(1)` holds (never
  written by the predictor). Today `apply_bc` repairs that afterwards, for
  walls, inlets and outlets alike. That is a hidden dependence of one stage
  on a later one; decide whether it is acceptable or the same class of error.
- **`apply_bc` writes the outflow row** post-momentum and at init / restart
  (`moby_solve.f90`, three call sites with `outflow_copy = .true.`), and NOT
  inside the projection (two call sites in `pressure_solver.f90`).
- **A Neumann normal velocity is refused by the row resolver** (numerics
  review F4). Relevant if the clean design is "the face is an ordinary
  predicted DOF and its closure comes from ghost values": on the low side
  the fused predictor would then predict face 1 as it stands, reading
  `q(0)` and the pressure ghost `p(0)` — and with the Dirichlet pressure
  ghost mirrored, `dpx` at face 1 IS the face's own gradient against
  `p_face = 0`. Check whether that is the whole low-side answer, and what
  the high-side equivalent costs (a plane outside the fused kernel's range;
  the fused predictor must stay untouched, see CLAUDE.md on its register
  schedule).
- **The projection**: `cfLow(1,d,b)` / `cfHigh(d,b)` (`init_face_metrics`),
  the mirrored `phi` ghost re-applied every iteration, `dnLow/dnHigh` in the
  denominator; `interface_correct` and the red-black path
  (`pressure_solver.f90` ~line 990) carry their own copies of the outlet
  logic. Confirm low and high are treated alike in all three.
- **The halo exchange** between the predictor and the projection
  (`syncface = .true.`): what it does to `q(nx+1)` on a block whose high
  face is physical, and to the edge / corner halos next to an outlet (the
  tangential extension copies a neighbour block's ghosts).
- **Init and restart**: the high-side face is reconstructed by the copy
  today (a one-time kick once the reset is gone). Continuity in the last
  cell gives it exactly; decide where that lives.
- **Who else reads the outlet face**: scalar convection (called BEFORE
  momentum, reads `q`), the RANS scalars, the control-volume forces, the
  CFL limiter, the 2:1 interface machinery if an outlet face is refined,
  `ibm%mu` if a body reaches the outlet.
- **The other outlet rows**: tangential velocities (Neumann ghosts), the
  stored pressure ghost (Dirichlet 0), `phi` (mirror), scalars. Are they
  consistent with a face that now has a momentum equation?

### Two further checks, asked for by the user (2026-10-01)

These widen the investigation from "the outlet face" to "how boundary faces
are indexed and stored". They belong in the same write-up, BEFORE a design
is chosen, because the outlet answer may fall out of them.

> * if the low/high indices make sense as presently done, or better solution
>   that can simply handle periodicity and other bc with the least number of
>   if conditions exist
> * if it is worth to change convention and have boundary-normal velocities
>   be placed on boundaries. Would this simply everything?

**(A) Does the present low/high indexing make sense, or is there a scheme
that handles periodicity and every other boundary condition with the fewest
conditionals?** The present convention, as found in the code:

- The normal velocity of cell `i` lives on its LOW face; a block stores
  `q(0:nb+1)` and owns faces `1..nb`. Its HIGH face `nb+1` is a halo slot:
  the neighbour block's face 1, the periodic image, or — on a physical high
  boundary — a face that NO block owns.
- So one and the same physical object, a boundary-normal face, is an
  interior index on the low side and a halo index on the high side. The
  special cases this produces are worth listing in full; known so far:
  `momentum_face_start` (start index 2 on `FACE_PHYS` / `FACE_CLOSED` low
  faces), the final `qs → q` copy over `1..nb` regardless of pinning,
  `cfLow(idx,d,b)` per normal index against `cfHigh(d,b)` per block,
  `dnLow/dnHigh`, the `if (i == nx)` plane branches of `jacobi_apply`, the
  `side == SIDE_MIN` branches in
  `apply_bc` (face 1 / neighbour 2 against face `nb+1` / neighbour `nb`),
  the exchange's "+axis face" divergence subset, the field file holding the
  low boundary face but not the high one, `oldrhs(1:nb)`.
- Count the conditionals each boundary kind costs today (periodic, wall,
  inlet, outlet, closed, 2:1 coarse / fine) per stage, low and high
  separately, and ask of each whether it is essential or an artefact of the
  indexing. The target is the treatment with the least number of `if`s that
  is still exact for periodicity (where low and high are the same face).

**(B) Is it worth changing the convention so that boundary-normal velocities
are PLACED ON the boundaries — owned, stored unknowns on both sides — and
would that simplify everything?** To be assessed, not assumed:

- What it would mean concretely: the high boundary face becomes a real
  degree of freedom of the last block (a predictor range, an RK memory
  slot, a row in the field file), exactly like the low one; a Dirichlet
  face is then a pinned unknown, an outlet face a predicted one, on either
  side, by the same code.
- What it would remove (candidates: the outflow row and `outflow_copy`, the
  restart hole of the high face, the low / high asymmetry of the projection
  tables, the pinned-face repair by `apply_bc`) and what it would cost: the
  array extents and the "block owns `1..nb`" rule that every kernel, the
  exchange (a high face is today ALWAYS the neighbour's), the 2:1
  low-block-owns-face rule, the io layout and existing restart files rest
  on; periodic directions, where an extra stored face would be a duplicate;
  bit-exactness of every outlet-free case.
- A middle form to weigh against both: keep the storage, but give physical
  high faces an explicit owner (per-block plane tables as `cfHigh` already
  is), so that "who predicts and who stores this face" has one answer.
- State the verdict with the count from (A): how many conditionals and
  special tables each option leaves, and what each breaks.

Questions the write-up must answer:

1. What is the momentum equation of an outlet face, stated once, for both
   sides? (Candidates: the neighbour's increment plus the own pressure
   gradient, as in the prototype; or a genuine one-sided predictor with
   ghost closure. They differ in the RK memory and in what `oldrhs` holds.)
2. Which existing stage is the natural owner on each side, such that no
   later stage has to repair an earlier one?
3. Can low and high be ONE code path? If storage forbids it, what is the
   smallest asymmetric piece?
4. What becomes of `BCK_FACE_OUTFLOW` and the `outflow_copy` argument? If
   nothing uses them after the change, they go.
5. Is the pinned-low-face overwrite by the final copy part of this, or a
   separate finding to record?
6. (A) Is the present low/high indexing the one with the fewest
   conditionals across periodic and all other boundary kinds? If not, which
   is, and what does moving to it cost?
7. (B) Should boundary-normal velocities be stored unknowns on both
   boundaries? Verdict with what it removes, what it breaks, and whether
   the outlet fix should wait for it or precede it.

## Gates

- Cases WITHOUT an outlet: untouched by construction; the 7-case and 9-case
  nofma suites `max_abs 0` against `~/step9b_ref_binaries`, CPU and GPU.
- The new implementation against the prototype binary
  (`~/moby_outlet_proto2/build_{cpu,gpu}`) on every high-side outlet case:
  agreement to round-off where the arithmetic is the same, or an explained
  difference where the design departs from it.
- `validation/freestream/run_gates.sh` (oblique exact, Poiseuille last-cell
  p = G dx/2, Lamb-Oseen ≤ 1e-4, 1 == 4 ranks, CPU == GPU), plus NEW legs
  with the outlet on a LOW face (mirror the Poiseuille and the vortex case:
  flow in −x, outlet at `x_min`) — the mirrored case must reproduce the
  high-side numbers.
- Restart: a run restarted at mid-length equals the continuous run (it
  cannot today on main either; state the tolerance reached).
- `validation/cylinder`: `run_outlet_mode.sh` legs (stored p rms dt-independent,
  the two step sizes within 1e-4), the Re 40 steady drag, the Re 100
  Strouhal gate on the control-volume series (decide which of St 0.167 and
  0.174 is this domain's number: vary the lateral far field).
- `validation/blasius`, then the boundary-layer tutorial, `tutorials/naca/rans`
  and the sailplane re-measured; `validation/README.md` is the checklist.
- Profile: no new kernel launch on an outlet-free case; on an outlet case
  the cost stated.

## Landmines

- The fused predictor kernel is not to be edited (register schedule,
  CLAUDE.md "predictor register guard"). Its final copy kernel and anything
  outside it are fair game, measured.
- The clean-p protocol (zero `pn`, 300 steps) is not clean on main: start
  cylinder legs from a state produced by the new scheme, run long enough.
- A pre-registered absolute time or force must come from the head it is
  compared with (the 77.7 ms miss of 2026-10-01).
- Everything in `docs/next_session_after_step9.md` "Landmines" holds.

## Investigation write-up (2026-10-01, main `154e48f`; no code changed)

Everything below is from reading the code and from eight small CPU runs
(`build_cpu`, and `~/moby_outlet_proto2/build_cpu` for the prototype). The
runs are listed in "Measurements"; their inputs are `sed` variants of
`validation/freestream/oblique.ini` and of a 16x32x4 inlet/outlet channel,
recipes given there. This section is the proposal as it stood before any
code; what was built from it is in "Implementation" below.

### One substage, stage by stage

`f` is the outlet face, `n` its interior neighbour: `(f, n) = (1, 2)` on a
low face, `(nb+1, nb)` on a high face. "U" stands for the component normal
to the face.

| stage (file) | LOW outlet, `f = 1` | HIGH outlet, `f = nb+1` |
|---|---|---|
| readers at the start of the substage: scalar and RANS convection, SGS, the Courant reduction | read `q(f)` as left by the last projection | same |
| fused predictor (`step.f90 momentum`) | `momentum_face_start` returns 2 for `FACE_PHYS`: the face is skipped; `qs(1)`, `oldrhs(1)` are never written; face 2 reads `q(1)` | the face is outside `1..nb`; there is no `oldrhs(nb+1)`, no `lap*(nb+1)`, no `q(nb+2)`; face `nb` reads `q(nb+1)` |
| the correction kernels (penalization state, SGS, body force, band filter) | same start mask | not in range |
| final copy `qs -> q` over `1..nb` | **`q(1) := qs(1)`, which is 0 since allocation** | untouched |
| `apply_bc(outflow_copy = .true.)` | `q(1) := q(2)`: repairs the copy AND resets the face | `q(nb+1) := q(nb)`: resets the face |
| `exchange_halos(syncface)` | no entry targets a physical face; same-level entries copy the neighbour block's face into the edge halos (tangential extension); a CROSS-LEVEL tangential entry spans `1..nb`, which includes face 1 | same, except that the cross-level span `1..nb` does NOT include `nb+1` |
| projection tables (`init_face_metrics`) | `dnLow(1) = 2 d1f`, `cfLow(1) = d1f` | `dnHigh(nb) = 2 d1f`, `cfHigh = d1f` |
| each iteration | `phi(0) := -phi(1)`; corrected by the ordinary low-face statement of `jacobi_apply`; `p += phi/dt_gamma` | `phi(nb+1) := -phi(nb)`; corrected by the `if (i == nx)` plane branch |
| red-black sweep | in place, metric `2 d1f` | same |
| last iteration: `apply_bc` (no outflow), full exchange | face untouched; tangential and pressure ghosts refreshed | same |
| snapshot | the face IS in the file | the face is NOT in the file |
| init / restart: `apply_bc(outflow_copy = .true.)`, twice | **the value read from the file is overwritten by `q(2)`** | filled with `q(nb)` |

Geometry, confirmed in `slice_grid_direction`: `coord(i,U) = face_at(first +
i - 2)`, so `U(1)` of the first block sits ON the low boundary and
`U(nb+1)` of the last block ON the high boundary; ghost node lines are
mirrored, so `d1(1,U)` and `d1(nb+1,U)` are both `1/dx` of the boundary
cell and `(p(f) - p(f-1)) d1(f,U)` IS the half-cell gradient against the
held outlet pressure on either side (the pressure ghost is `-p + 2 value`).

### What is wrong, and what is accidental

1. **The reset** (the known defect): `apply_bc` writes the outlet face
   after the predictor. The face has no momentum equation and the stored
   pressure no anchor.
2. **The pressure LEVEL is not tied to the outlet value at all** (new,
   measured): a constant added to the whole pressure field changes nothing
   on main, because the only face that would feel it is the outlet face,
   whose "predictor" carries no pressure gradient. Run `touchLow` below ends
   with uniform velocity to 7e-10 and `p = 0.2254937` everywhere (min = max
   to seven digits) where the outlet pressure is 0; run `touchS` drifts at a
   constant rate to a level of 31 after 2400 steps.
3. **The final copy writes faces the predictor did not predict.** At every
   pinned low face (wall, inlet, outlet, closed) `q(1) := qs(1) = 0`, and
   the next `apply_bc` puts the value back. Same class as 1: a stage writes
   what it does not own and a later stage repairs it. On a high face nothing
   of the kind happens.
4. **Init / restart**: on the high side the outlet face is missing from
   the file; on the low side it is in the file and is thrown away by the
   init copy. Both sides lose the state, by different routes.
5. **A level jump touching a HIGH inlet or outlet face** (new, measured):
   the edge halo `q(nb+1, 0 | nb+1)` of the normal component is written by
   same-level entries only. Next to a 2:1 neighbour nothing ever writes it,
   and the predictor of the tangential component reads it (`vu_p`, `wu_p`).
   It keeps its initial value for the whole run. The uniform-flow gates
   cannot see this: they initialise every index to the right value.
6. Stale or dead: the comment in `jacobi_apply` ("the predictor never
   writes the outlet face"), the `BC_OUTFLOW` comment in `boundary.f90`
   ("apply_bc skips the face write"), and the `outLow/outHigh` arguments of
   `interface_correct` (an interface face never takes the outlet branch).

**Not wrong**, confirmed: the projection treats low and high alike in the
Jacobi / Chebyshev path, in `interface_correct` and in the red-black sweep
(denominator `2 d1f`, correction `d1f` against the mirrored ghost or `2 d1f`
in place, `ibm%mu` of the face on both). Measured below to round-off.

### Measurements

Mirror pair: 16x32x4, `lx 0.5`, `nb 4`, Re 100, fixed `dt 1e-3`, 400 steps
from rest, parabolic inlet `u = +-1` with `v = 0.3` (the outflow is not
parallel), walls in y; `high` = inlet `x_min` / outlet `x_max`, `low` = the
mirror image. Compared under `u(i) <-> -u(nx+2-i)`, `v, p (i) <-> (nx+1-i)`.

| run | result |
|---|---|
| main, Chebyshev 12: `high` against mirrored `low` | u 2.8e-15, v 1.6e-15, p 5.9e-14 |
| prototype, same | u 3.4e-15, v 1.5e-15, p 3.1e-14 |
| main, red-black 12 (`sor 1.5`) | u 1.3e-7, p 1.7e-6: not round-off, consistent with the colour order not being mirror-invariant (even `nx` swaps red and black); not pursued |
| main against prototype, same case | u 0.10, v 0.12, p 0.11: the two schemes are different flows when the outflow is not parallel |
| restart, main: 200 + 1 steps against 201 | u 1.9e-9, p 1.3e-7 (small because main forces the gap `u_out - u_n` to zero anyway) |
| restart, prototype: 200 + 1 / + 20 / + 200 | u 5.6e-4 / 5.7e-4 / 2.5e-4, p 6.0e-2 / 2.2e-4 / 5.6e-5: the one-time kick of the init copy |

**A low-side outlet runs, for the first time, and mirrors the high side to
round-off, on main and on the prototype.** So the `(f, n)` index pair is all
that distinguishes the sides, for the reset as for the prototype's row.

Level jump: `oblique.ini` (inlets `x_min`, `y_min`, `y_max` at
`(0.940, 0.342)`, outlet `x_max`, 32^3, `nb 8`), Chebyshev 12,
**`initial_u = initial_v = 0`**, one refinement box, `dtmax 2.5e-3`.

| run | box | velocity against uniform | stored pressure |
|---|---|---|---|
| `innerS` | `0.25 0.75 0.25 0.75 0 1`, touching nothing | 1.4e-11 after 2400 steps, falling a decade per 200 | 1e-10 |
| `touchS` | `0.5 1.0 0.25 0.75 0 1`, touching the outlet | **2.28e-2, frozen** (fine block at the outlet, on the level jump) | **level growing linearly: 7.9, 15.7, 23.4, 31.1** at steps 600..2400 |
| `touchS90` | same, fields initialised at 0.9 of the freestream | **2.20e-3**: ten times smaller, as a halo frozen at its initial value must give | growth ten times slower |
| `touchLow` | mirror image: outlet at `x_min`, box `0.0 0.5 ...` | 6.9e-10 after 1800 steps, falling | frozen uniform 0.2254937 |

So item 5 is a high-side-only hole (the low face is index 1, inside every
entry's range), and it is independent of the outlet row: it would hit a
high-side inlet with a non-zero normal velocity just as well. No committed
case has a level jump on a physical face other than a wall, where the
frozen value 0 happens to be right.

### The proposed design: the outlet row moves into the predictor

The momentum equation of an outlet face, stated once for both sides:

    du_f/dt = R_n - (dp/dn)_f

`R_n` is everything but the pressure gradient in the predictor of the
neighbour face (convection, diffusion, forcing, SGS, body force,
penalization): the closure is "zero normal gradient of the right-hand
side". `(dp/dn)_f` is the face's own gradient against the held outlet
pressure. In the predictor's own spelling, with
`dp(i) = (p(i) - p(i-1)) d1(i,U)`:

    q(f) = q(f) + [qs(n) - q(n)] + dt_gamma [dp(n) - dp(f)]

One kernel over the outlet points (`predict_outlet_faces`, in `step.f90`,
called inside `momentum()` after the last correction of `qs` and before the
final copy, so that `q(n)` is still the old value). `(f, n)` is the only
thing that depends on the side, exactly as in `apply_bc`. It needs nothing
the high side lacks: no `oldrhs`, no `q(nb+2)`, no Laplacian row. The RK
memory of the face is its neighbour's. This is the prototype's arithmetic
(its `C + u_n*`), owned by the predictor instead of being pre-loaded into a
boundary row. Differences from the prototype: no call before `momentum()`,
no row, no flag; the pressure ghost is READ, so a non-zero outlet pressure
value is honoured (the prototype hard-codes 0).

What goes: `BCK_FACE_OUTFLOW`, the `outflow_copy` argument and its branch in
`apply_bc`, the three `outflow_copy = .true.` call sites. `BC_OUTFLOW` stays
as the TYPE of the row (it is in restart metadata) and resolves to "no
boundary row"; the outlet points become a compact list in `boundary_type`.

Why the alternatives lose:

- **The face as an ordinary DOF with ghost closure.** On the LOW side it is
  nearly free: give outlet faces a kind that `momentum_face_start` does not
  pin and the fused kernel predicts face 1 with its full stencil, reading
  `p(0)` (the right gradient already) and a new ghost row for `q(0,U)`. On
  the HIGH side the same equation needs a second copy of the whole momentum
  stencil for one plane, in three directions, plus SGS and body force, or
  the storage change of (B). And the two sides would then solve different
  discrete equations: the mirrored case could not reproduce the high-side
  numbers. Rejected: duplicated code or an asymmetric scheme.
- **The prototype as it is.** Keeps the reset row and compensates it;
  an extra call. Rejected by the user.
- **Eliminating the face through continuity.** The projection needs the
  face flux as an unknown under a Dirichlet pressure. Not an option.

### Answers to the seven questions

1. `du_f/dt = R_n - (dp/dn)_f`, above. The RK memory is the neighbour's;
   `oldrhs` holds nothing for the face.
2. The predictor (`momentum()`) owns the predicted value, the projection
   the correction (`cfLow(1)` / `cfHigh`), on both sides. `apply_bc` owns
   nothing of an outlet-normal face.
3. Yes, one code path, through the `(f, n)` pair. The asymmetric pieces
   that remain are storage facts, not outlet logic: the start mask (low),
   `cfHigh` and its plane branch (high), the file (high).
4. `BCK_FACE_OUTFLOW` and `outflow_copy` go. Nothing uses them.
5. It is part of this. The copy must commit only the faces the predictor
   predicted (the start masks it already has); then a pinned face is
   written at init and never again, low or high, and the outlet kernel can
   write `q(f)` directly on both sides. Without that fix the low side would
   need the kernel to write `qs(1)` and the high side `q(nb+1)`.
6. (A) The present indexing is the one with the fewest conditionals. For a
   periodic or open face it costs none: face 1 is owned, face `nb+1` is the
   neighbour's, every face is stored once. A non-periodic line has `n+1`
   faces on `n` cells, so with equal-size blocks one boundary face must sit
   in a halo slot; "cell owns its high face" is the mirror image, not an
   improvement. The census of what the asymmetry costs today:
   - low only: the start mask (5 sites: fused predictor, three correction
     kernels, the band list), `zero_closed_halos` zeroing `q(1,U)`, the
     red-black start index;
   - high only: `cfHigh` and three plane branches in `jacobi_apply` (outlet
     only), the `+axis` entries and the divergence subset of the exchange,
     the cross-level reach (item 5), the file, `oldrhs` / Laplacian extents;
   - symmetric through one index pair or table: `apply_bc`'s `side` branch,
     `dnLow/dnHigh`, `interface_correct`, the RANS and scalar wall tests;
   - accidental, removed by this plan: the final copy, `outflow_copy`.
   The only scheme with fewer `if`s is "predict every face, then pin": no
   start mask at all. It edits the fused kernel and computes with unset
   ghosts. Not proposed.
7. (B) No. The boundary-normal velocities already SIT on the boundaries on
   both sides; what differs is who stores them. Making the high face an
   owned unknown needs a predictor range of `nb+1` for one component (an
   edit of the fused kernel), a second halo layer (`q(nb+2)`), `oldrhs` and
   Laplacian rows at `nb+1`, a file layout change and a duplicate face in
   every periodic direction. It would remove the plane branch, the file
   hole and item 5. The proposed row needs none of that storage, so the
   outlet fix does not wait for it; the file hole and item 5 are closed
   separately below, each where it lives.

### Increments (each gated before the next)

- **O0. The copy commits what the predictor predicted.** `qs -> q` under the
  start masks. Pure refactor: `max_abs 0` at nofma on the 7- and 9-case
  suites AND on the outlet cases (the value `apply_bc` restored is the value
  now kept). Measure the copy kernel on the GPU; it gains three predictable
  comparisons per cell.
- **O1. `predict_outlet_faces`**, the removal of the row and the flag, the
  outlet point list, a hard init error when `ibm%coef` is non-zero on an
  outlet face or its neighbour (the row assumes `mu = 1`). Init keeps the
  zero-gradient value as the documented INITIAL condition of a face the file
  does not carry. Pre-registered: outlet-free suites `max_abs 0`; **equal to
  the prototype binary at `max_abs 0` (nofma) on every high-side case** if
  the sum is written in the prototype's order, `qs(n) + ((q(f) - q(n)) +
  dt_gamma (dp(n) - dp(f)))`, so the prototype's table carries over
  (Poiseuille last cell 1.247e-3, Lamb-Oseen 9.8e-5, the cylinder and
  Blasius rows); mirrored low-side legs within 1e-13 of the high-side ones
  (measured today for both existing schemes: 3e-15 velocity, 6e-14
  pressure); `touchLow` ends at `p = 0`, not 0.2255.
- **O2. Restart.** Two ways, the user's choice. (a) Store the high outlet
  planes in the snapshot as one optional dataset per direction and read the
  low one that is already there instead of overwriting it: restart equals
  the continuous run exactly. (b) Rebuild the face from continuity in its
  cell: no format change, equal only to the projection residual, and
  undetermined in a cell with two outlet faces (the Blasius corner).
  Recommended: (a). The kick to beat, measured on the prototype: 5.6e-4
  velocity, 6.0e-2 pressure.
- **O3. Drop the post-predictor `apply_bc`.** After O0 and O1 it writes
  only values that the projection's own `apply_bc` rewrites before anything
  reads them (tangential and pressure ghosts; pins that no longer move). A
  hypothesis from reading, to be decided by the gate: `max_abs 0` on every
  suite, both smoothers. Worth one kernel launch per substage.
- **O4. Level jump on a physical high face.** Let the cross-level
  tangential entries reach `nb+1` in a direction whose combined neighbour is
  absent, as the same-level entries do. Gate: `touchS` recovers uniform flow
  to round-off like `touchLow`; the interface suites `max_abs 0`. Until it
  is in, a level jump touching a high inlet or outlet face should be a
  config error.
- **O5.** Dead arguments and stale comments (item 6), the new low-side legs
  in `validation/freestream/run_gates.sh` (mirrored Poiseuille and vortex,
  plus the zero-initialised refined legs above), then every outlet case of
  `validation/README.md` re-measured, the cylinder Strouhal question
  included.

### Side findings (not outlet matters)

- **`pecletmax` is a per-direction number.** RK3's real-axis limit is
  `4 nu dt sum(1/h^2) <= 2.51`, i.e. 0.63 / 0.31 / 0.21 per direction for
  one / two / three equally fine directions. The first refined runs here
  used the gate's `pecletmax = 0.5` and `dtmax = 1e-2` (0.41 on the fine
  level) and settled into a sustained O(0.65) disturbance inside the patch:
  the diffusive instability held at the stability boundary by the Courant
  limiter (`0.8/(64 dt_lim)` gives `|u|` 1.63 = freestream + 0.69). At a
  quarter of the step the same case converges to 1e-11. Attributed, not
  proven. 829 of the 901 tracked inis carry 0.5; it is safe where one direction dominates
  (stretched wall-normal lines) or another limit binds first.
- `q(0,U)` on a physical low face is never written and never read by the
  solver proper; the Courant reduction, which spans `0..nb+1`, sees its
  initial value.

## Implementation (2026-10-01, on top of main `154e48f`)

Reference binaries: `~/outlet_ref_src` (a pinned worktree of `154e48f`,
all four builds) and the prototype's nofma builds in `~/moby_outlet_proto2`.
Gate drivers: the 7- and 9-case suites, the new
`validation/freestream/run_bitexact_outlet.sh` (seven OUTLET cases: oblique,
pois_io, lamboseen, the mirrored pair, blasius2d, cyl_re100) and a scratch
list of six refined cases (`refine2d` rans_xz / uniform_xz / bp_xz_32,
`redblack_interface`, `interface_decay`, `multilevel_body`). nofma on both
sides throughout.

### O0. The copy commits what the predictor predicted

`step.f90 momentum`: the final `qs -> q` runs under the start masks. A
pinned low face is written at init and never again.
`max_abs 0` against `154e48f` on the 7-case suite (1 and 4 ranks), the
9-case suite and the seven outlet cases (1 and 4 ranks), CPU and GPU.
Cost, RTX 3060, 128x64x128 channel, three runs each: momentum bucket
1.8758 -> 1.8713 ms per call, step 38.33 -> 38.26 ms. Not measurable.

### O1. The outlet face is predicted

`step.f90 predict_outlet_faces`, called inside `momentum()` between the last
correction of `qs` and the commit:

    q(f) = qs(n) + ((q(f) - q(n)) + dt_gamma (dp(n) - dp(f)))

over the boundary points of the outlet faces, `(f, n) = (1, 2)` or
`(nb+1, nb)`. `BCK_FACE_OUTFLOW` and the `outflow_copy` argument are gone;
`BC_OUTFLOW` remains as the type and resolves to no boundary row;
`boundary_type` carries `faceOutlet` and `nOutlet` (no outlet point on the
rank: no launch). Init: `init_outlet_faces` gives a face the run was not
handed its zero-gradient value, once. An immersed body touching an outlet
face or its neighbour is a hard error (`check_outlet_faces_fluid`).

- Outlet-free suites `max_abs 0` against `154e48f`, CPU (4 ranks) and GPU.
- **Equal to the prototype binary at `max_abs 0` on all seven outlet
  cases**, GPU and CPU, 1 and 4 ranks, the low-side case included. The
  prototype's table is therefore this scheme's table.
- One code path, checked in every direction: the mirrored pair rerun with
  the flow along y (inlet/outlet on the y faces, walls in x) and along z
  equals the x case under the coordinate permutation at **0.0** in u, v, w
  and p, for the outlet on the high side and on the low side (production
  CPU build, 400 steps).
- FOUND, nvfortran 25.9: a variable used ONLY as the bound of mapped array
  sections has its assignment dropped. The first form of the kernel looped
  over a compact outlet list, `npts` appeared only in
  `bc%pointFace(1:npts)`, and two runs in three stopped with "partially
  present on the device" (a garbage size). The kernel loops over all
  boundary points and skips, as `apply_bc` does; the list is gone.

### O2. Restart

`io.f90 write_outlet_planes / read_outlet_planes`, `field_hdf5.c
fdm_h5_append_block_rows / fdm_h5_read_block_rows`: one dataset per
direction with a HIGH outlet (`un_xmax`, `vn_ymax`, `wn_zmax`), one row of
the block's `q(nb+1)` plane per block. A low outlet face is index 1 of the
velocity dataset and is now USED on restart instead of overwritten. Files
of cases without a high outlet are unchanged; a file without the plane is
announced and falls back to the zero-gradient value.

`validation/freestream/run_restart_outlet.sh`: restart at mid-length equals
the continuous run at **`max_abs 0`** on outlet_high, outlet_low, pois_io,
lamboseen and blasius2d (two outlets), CPU (1 and 4 ranks; blasius2d has
no `[blocks] nb` and stays on the one rank its IC was minted on) and GPU.
With the plane stripped
the same restarts differ by 1.9e-4 / 1.5e-7 / 4.4e-5 / 7.3e-6 in velocity.

FOUND while gating: **a solver-minted template carries the plane of its own
one-step run**, so an IC tool that overwrites only `un/vn/wn/pn` restarts
the outlet face from an unrelated field (pois_io against the prototype:
2.3e-3 until the plane was removed). `make_freestream_ics.py` and
`make_blasius_ic.py` now write the analytic plane.
`tutorials/naca/rans/interp_restart.py` builds a fresh file without one and
gets the announced fallback.

### O3. No boundary write after the predictor

The `apply_bc` between `momentum()` and the exchange is removed, with its
profiling bucket. `max_abs 0` against the tree that still had it: 7-case,
9-case, seven outlet cases, six refined cases, and the mirrored pair under
red-black, GPU; the final tree against `154e48f` and the prototype on CPU.
`apply_bc` now runs once per substage (it was 0.34 % of the step on the
3060 channel).

### O4. Level jump on a physical face

`comm.f90 candidate_boxes`: cross-level entries extend tangentially into a
halo row whose combined neighbour is absent, as same-level entries always
did. An unrefined dim and a PROLONG dim extend inside the entry (their
affine gather maps land on the source's row 0 / nb+1). A refined RESTRICT
dim averages two fine rows and cannot; its halo row is a separate one-row
piece of the candidate (`N_PIECES`, `plane(d)` in `entry_gather_map`). No
kernel changed.

- `zr_high` (zero start, level jump on the x_max outlet): **2.9e-13** from
  uniform, stored p within 8e-14; it was 2.28e-2 frozen with p growing to
  31. `zr_low` 2.2e-13 with p within 8e-14 (it was p = 0.2255 everywhere).
  `zr_inner` 1.8e-13. 1 == 4 ranks EXACT on `zr_high`.
- The six refined cases `max_abs 0` against `154e48f`, GPU and CPU 4 ranks
  (`rans_xz` has level jumps along both walls and compares k, omega, nut).

### O5. Clean-up and re-measurement

`interface_correct` lost its two dead outlet arguments; stale comments
fixed. `validation/freestream/run_gates.sh` gained the groups `low`,
`refined`, `restart`.

| case | main `154e48f` | now |
|---|---|---|
| Poiseuille last-cell p (exact 1.25e-3) | 2.41e-3 | **1.247e-3** |
| Poiseuille profile dev at x/lx 0.5 / 0.9 | 1.639e-3 / 1.539e-3 | 1.637e-3 / 1.521e-3 |
| Lamb-Oseen reflected fraction | 2.2e-2, regrowing | **9.80e-5**, monotone |
| the same two cases with the outlet at x_min | never run | identical numbers; fields equal to 5.6e-16 / 4.7e-15 |
| mirrored pair, high against low | 2.8e-15 | 3.4e-15 |
| Blasius worst theta / H / du/Ue / dv | 1.33 % / 0.45 % / 2.40e-3 / 0.124 | 1.65 % / 0.82 % / 1.48e-3 / **0.0129** |
| Blasius, spread over four projections | theta 1.33 .. 1.74 %, dv 0.096 .. 0.124 | theta 1.65 .. 1.66 %, dv 0.0127 .. 0.0130 |
| cylinder Re 100, stored p rms after one t.u., dt 5e-3 / 6.25e-4 | 0.89 / 1.65 | **0.080 / 0.073** |
| cylinder, velocity dt 5e-3 against 6.25e-4 | 0.15 | **6.0e-5** |
| cylinder, 100 t.u.: St / mean C_D / C_L amplitude | 0.1670 / 1.4474 / 0.655 | 0.1744 / 1.4383 / 0.393 |
| cylinder Re 40, C_D (60 t.u. from main's converged state) | 1.69235 | 1.69682 ± 1e-5 (+0.26 %: the outflow is no longer forced parallel) |
| `rans_inlet` | pass | pass (k 0.14 %, omega 2.15 %, 1 == 4 exact) |
| `multilevel_body` | exact | exact, 1 == 4 ranks, GPU |
| `tutorials/naca/rans`, 0.5 t.u. from the converged state, both binaries | C_L 0.52007, C_D 0.01275 | C_L 0.52008, C_D 0.01275 |
| `tutorials/sailplane`, 200 steps from rest, red-black | stable | stable; within 5.9e-4 in velocity |

**The Strouhal question is settled: 0.174 is the 16 x 16 domain's number.**
Cold starts to t = 200 on a domain twice as long and on one twice as tall,
both binaries (`validation/cylinder/README.md`, `shedding_fit.py`):

| domain lx x ly | scheme | St | second line in the lift | mean C_D |
|---|---|---|---|---|
| 16 x 16 | main | 0.1677 | f 0.503, amplitude 0.415 | 1.4474 |
| 16 x 16 | now | 0.1737 | none | 1.4383 |
| 32 x 16 | main | 0.1739 | none | 1.4398 |
| 32 x 16 | now | 0.1743 | none | 1.4412 |
| 16 x 32 | main | 0.1663 | f 0.240, amplitude 0.842 | 1.5007 |
| 16 x 32 | now | 0.1723 | none | 1.4277 |

Main with the outlet at 26 D lands on what the new scheme gives at 10 D
and at 26 D. Its extra lift line is not a harmonic (3 St on one domain,
1.44 St on the other) and is the "3rd harmonic of the confined lift" of
the A2 notes.

The NACA row needs reading with care: the outlet is 77 chords behind the
airfoil and the pair ran only 20000 steps from the 2026-08 converged
state (case file re-prepared, both binaries of today), so it shows that
the outlet change does not disturb the converged solution (fields within
1.7e-5 in velocity), not what a from-scratch convergence would give. Both
binaries read C_D 0.01275 where the tutorial records 0.01294 ± 0.00017:
whatever that difference is, it is not the outlet.

Cost: outlet-free channel on the RTX 3060, 39.44 -> 39.32 ms per step (the
removed `apply_bc`); cylinder 4.2 M cells on an A6000, 87.8 -> 87.6 ms;
NACA on the RTX 5090, 0.171 -> 0.170 s.

### Follow-up, 2026-10-02: the open items of "Not done"

Reference for this section: `~/outlet_runs/ref_949148f` (the committed tree,
nofma CPU/GPU and production GPU). The three source changes below are a
report in the time-step limiter (fields untouched), two comments and the
gates; the 7-case, 9-case and outlet suites are `max_abs 0` against the
committed tree (nofma, CPU 1/4 ranks and GPU), and
`validation/freestream/run_gates.sh` passes in full with the rebuilt binaries.

**1. The turbulent boundary layer (`tutorials/turbulentBoundaryLayer`): two
restart pairs, and the old outlet was stalling the layer.** The production
case (4096 x 224 x 192, red-black 6) restarted from the developed field of
the tutorial (step 1050000) with `154e48f` and with `949148f`. Both outlets
of the case are high faces (x_max, and the top, where the entrainment
enters). Lengths in inlet displacement thicknesses; delta99 at the outlet
is 18.

- (a) 100 time units on one RTX 5090 (0.58 s per step, 26 GB),
  `~/outlet_runs/tbl/` (`compare_tbl.py`).
- (b) 750 time units = one flow-through on HoreKa, 4 A100 per side, both
  sides in the same allocation (0.194 / 0.198 s per step), three chained
  one-hour jobs 5175219-21
  (`overheadTest/horeka/exchange/submit_outlet_tbl.sh`, results in
  `~/outlet_runs/tbl/horeka/`). The statistics below are the window
  250 .. 750, after the first 250 in which the new run leaves the old
  scheme's state behind.

Both stable throughout. Health and stored pressure, (b), last chunk:

| | `154e48f` | `949148f` |
|---|---|---|
| L2 divergence, mean | 1.6e-5 | 3.2e-6 |
| net mass imbalance (runtime column), mean / max | 0.50 / 1.38 | 0.040 / 0.11 |
| stored p, last column, rms (final field) | 6.8e-2 | 6.9e-4 |
| stored p, column at 0.9 lx, rms | 8.0e-3 | 2.7e-3 |
| mean p on the outlet plane, max over y | 5.4e-3 | 8.6e-5 |

The momentum integral, `d theta/dx` against `c_f/2 - (H+2)(theta/U_e)
dU_e/dx` (`assets/postpro/momentum_integral.py`), and the mean wall pressure
(the freestream pressure stays within 7e-4 of zero in both):

| x band | old: d theta/dx | balance | p_wall | new: d theta/dx | balance | p_wall |
|---|---|---|---|---|---|---|
| 450 .. 600 | 2.0 .. 2.3e-3 | 2.2 .. 2.3e-3 | +2e-4 | 2.0 .. 2.3e-3 | 2.2 .. 2.3e-3 | -1e-4 |
| 650 .. 680 | 1.8e-3 | 2.1e-3 | -4.0e-4 | 2.2e-3 | 2.1e-3 | +1.7e-4 |
| 700 .. 715 | 1.0e-3 | 2.0e-3 | -1.5e-3 | 2.1e-3 | 2.0e-3 | +1.7e-4 |
| 730 .. 740 | 4.0e-4 | 2.0e-3 | -2.2e-3 | 2.3e-3 | 2.0e-3 | +6.8e-4 |
| 745 .. 748 | 2.2e-4 | 2.0e-3 | -2.8e-3 | 2.5e-3 | 2.0e-3 | +8.8e-4 |
| 749.5 .. 750 | 6e-5 | 2.0e-3 | -3.0e-3 | 4.1e-3 | 2.2e-3 | +1.3e-4 |

**With the old outlet the momentum thickness stops growing over the last
~100 units (5 delta99)**: the wall pressure falls by 3e-3 toward the outlet
while the freestream pressure does not, a favourable pressure gradient
inside the layer that no zero-pressure-gradient layer has (`p_wall = p_e`
there, `<v'v'>` vanishing at the wall). It is the outlet pressure mode in its
time-averaged form. The tutorial's own 4000-time-unit statistics, produced
by the old code, show the same stall to the digit that matters (d theta/dx
1.2e-3 at 700..715, 5e-5 in the last cell, p_wall -4.2e-3): it was in the
shipped data, in a region nobody compared. **With the predicted face the
balance holds, within ~20 %, up to 2 units from the outlet**; what remains is
a wall pressure 9e-4 high over the last ~20 units and the adjustment of the
last cells to the uniform face pressure (the mean pressure inside the layer,
`-<v'v'>`, is brought to 0 over the last ~5 units: near-wall acceleration,
c_f +13 % over the last 10 units, most of it in the last cell).

New against old, statistics of the window: theta within 0.06 % up to x = 500,
+0.3 % at 600, +1.0 % at 650, +2.3 % at 700, +6.0 % at the outlet; c_f
within the sampling scatter of two decorrelated 500-unit windows (+-2 %) up
to x = 650, then 1.5 .. 2 % LOWER to the outlet (-3 .. -5 % against the
tutorial's long statistics: the old favourable gradient raised it) and
+11 % in the last cell. At the comparison station of the tutorial (Re_theta
677, x = 400) nothing moves beyond that scatter, so the tutorial's table
stands; its fields downstream of x ~ 600 are the old outlet's.

In (a), where the two runs still share the realization, the same picture at
its start: statistics within 7e-4 for x in 100 .. 600, theta +0.3 % at 30
from the outlet and +1.1 % at it, stored p of the last column 8.1e-2 -> 4.1e-4.

**2. `tutorials/naca/rans` from scratch: still not run.** Days of GPU time,
and it belongs to the one-by-one re-validation of `validation/README.md`.

**3. `pecletmax`: proven, then fixed (the user's decision, the same day).**
The attribution of the side finding is a measurement,
`validation/diffusion_limit/`: a decaying Beltrami flow at Re = 1 in a
periodic box, fixed step just below and just above
`dt sum_d nu/h_d^2 = 0.628` on grids with three, two and one fine direction.
Decay below, blow-up above, on all three (limits on the single-direction
number 0.209 / 0.305 / 0.558, bracketed to 2.5 %). With the adaptive step and
the old `pecletmax = 0.5` on the isotropic grid the run does not end in NaN:
it sits in a bounded state held by the Courant limiter, which is the O(0.65)
disturbance of the refined patch above, without a patch.
**`pecletmax` now bounds the sum** (`precompute_peclet_rate`, the
eddy-viscosity reduction of `get_timestep_rates`), so one value is stable or
not on any grid: 0.60 decays and 0.66 does not on all three (gate 12/12).
The conjugate rate becomes half the Gershgorin diagonal at EVERY cell, which
is the sum for a regular cell and what cut cells already used (`share = 2`):
the cut-cell special case is deleted, its step unchanged. A `pecletmax` above
0.628 is announced at init and a step beyond the limit warned about once
(fixed step too).
What moves, against `949148f` (nofma): 17 of 23 suite comparisons
`max_abs 0` on CPU and 11 of 14 on GPU (bound by `dtmax` or the Courant
limit); turb180, lam30t,
conduction, prsweep, turbsst and pois_io take a smaller step, and all six
are `max_abs 0` again with `dtmax` binding on both sides, so they move
through the step alone. The S1 scalar gates, the conjugate gates C1 / C2 /
C3 (61 / 75 / 15 checks), the freestream and the penalization gates pass
unchanged. `validation/refine2d/bp_xz_{32,64}.ini` carried 0.8
(never binding, `dtmax` does) and are set to 0.5.
**`cflmax` followed the same evening** (the user's decision): it bounds
`dt sum_d |u_d|/h_d`, whose RK3 limit with central convection is sqrt(3) for
a flow in any direction (`validation/courant_limit/`, 12/12: 1.65 stable and
1.80 NaN along one, two and three directions). Against the `pecletmax`
commit, 22 of 23 suite comparisons `max_abs 0` on CPU and 13 of 14 on GPU;
`beltrami_yslab` is Courant-bound and `max_abs 0` with `dtmax` binding. The
same gates pass. The inis keep `cflmax = 0.8`, which is now 46 % of the
limit and costs a Courant-bound run 1.2 .. 1.5 times more steps; 1.2
restores the old step. Reference set after both: `~/cfl_ref_binaries`.

**4. LES at a level jump on a physical face: measured, one residual.** The
SGS kernels do not read edge ghosts of cell-centred quantities that matter
(the eddy viscosity's physical ghosts are never written and stay 0, before
and after O4). What they read is the VELOCITY edge halo: the gradient tensor
of the wall cell `(nb, 1)` takes `du/dy` from `u(nb+1, 0)`. So O4 changes an
LES with a level jump running into a wall, the first case in the tree it
changes: `validation/channel_interface/les/run_wall_jump.sh` (x-invariant
profile, Smagorinsky, one step). Wall-row eddy viscosity in the column next
to the jump against its level's median: fine side 25.8 % before, 0.94 % now;
coarse side 13.4 % before, **6.8 % now**. The coarse-side remainder is the
one-row restriction of O4: the coarse ghost receives the single fine ghost
row, `-u_f(1)`, where the mirror of its own halo column is
`-(u_f(1) + u_f(2))/2`. Exact for uniform flow, first order otherwise. Closing
it needs unequal weights `(3/2, -1/2)` on the fine rows `(0, 1)` in the
gather kernels, or the boundary row applied to the halo column after the
exchange; neither is done (WALE's wall-cell eddy viscosity is ~1e-4 of the
molecular one, resolved RANS multiplies that strain by zero, the momentum
stencil does not read the halo at a wall). It would matter for a Smagorinsky
run or for a tangential velocity at an inlet/outlet face next to a jump.

**Also fixed here:** the sailplane excerpt of `docs/tutorials.md` still
showed the Neumann rows that are a config error since F4; two comments
described the outlet as a zero-gradient velocity; the gate scratch of
`validation/redblack_interface/run_gates.sh` is ignored by git.

### Not done

- **The tutorial's data are REGENERATED (user's decision 2026-10-02; HoreKa
  job 5175523 ran 2026-10-03/04, 26.5 h on 4 A100, `f39c139`, from
  `outlet_tbl_run/tbl_new_1087500.h5`; data.nc and the figures recommitted
  2026-10-05).** At Re_theta 677: c_f 0.00465 (was 0.00462), H 1.499 (same),
  u'_rms 2.738 (2.707), -u'v' 0.866 (0.867) -- sampling, the station is
  untouched. The outlet zone of the 10000-t.u. statistics
  (`momentum_integral.py`): the balance within 5 % to x = 650, 15 % to 715,
  20 % to 740, the wall pressure rising 1e-3 over the last 60 units, the last
  2 units adjusting (c_f +14 %); README "The outlet zone". The final field is
  the tutorial's local `restart_field.h5` now (step 1587500). The run
  directory on HoreKa (`tbl_stats_run`, 11 snapshots of 5.6 GB) can go.
- `tutorials/naca/rans` from scratch (days of GPU time): only the restart
  pair above.
- The first-order ghost row of a restriction at a level jump on a physical
  face (follow-up item 4).
- FOUND while gating the limiter, not caused by it: the CONTROL rows of the
  C3 transient gate (`validation/conjugate/run_gates_c3.sh`, the archived C2
  binary `~/c2_ref_binaries`, commit 7aa1c7b of 2026-08-28) no longer
  reproduce. The gate itself passes and its C3 rows are unchanged to every
  digit, but the control reads errors of 0.55 where the committed
  `c3_transient.dat` records 7e-3 .. 1e-3: a binary of August run on inis
  and case files of today (unknown keys are silently ignored). The committed
  records were left as they are; the control needs a binary that reads the
  present inputs, or retiring.
- Long output runs are in `~/outlet_runs/` (cylinder far-field pair `st/`,
  Re 40 `re40/`, NACA pair `naca/`, sailplane pair `sail/`); the pinned
  reference worktree is `~/outlet_ref_src`.
