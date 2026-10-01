# Next session: the outlet face, done properly

STATUS: **NOT STARTED. Decided by the user 2026-10-01.** Read this first,
then `validation/cylinder/README.md` (last section: the evidence) and
`CLAUDE.md` ("THE OUTLET PRESSURE MODE IS RUN DOWN").

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
