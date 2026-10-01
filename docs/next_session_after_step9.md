# Next session: after numerics review step 9

STATUS: **ITEMS 1 AND 2 DONE (2026-10-01), item 2 measured on HoreKa at
4, 8 and 16 ranks; item 3 RUN DOWN with a
measured prototype and DECIDED by the user the same day (the root fix,
investigation first: `docs/next_session_outlet.md`); items 4a and 4b DONE.**
Steps 0–9 of `docs/numerics_review_2026-09-26.md` section 10 are done. What
is open for the next session: (i) **`docs/next_session_outlet.md`, which
comes FIRST**; (ii) item 5, the plan. The commits of 2026-10-01 after `5c53a96` are LOCAL (not pushed).

## What holds now

- **Reference set `~/step9b_ref_binaries`** (PROVENANCE inside): F1, the
  predictor guard, the RATIONAL penalization factor (item 1). Body-free cases
  are `max_abs 0` (nofma) against `~/step8_ref_binaries` and
  `~/step9_ref_binaries` too; body cases are not (`~/step9_ref_binaries`
  carries the exponential).
- **Penalization.** `ibm%mu = (1 + x/2 + x²/6)/P3`, x = λ dt_γ,
  P3 = 1 + x + x²/2 + x³/6, is the factor on every velocity INCREMENT
  (predictor, projection correction, SGS / body-force / band-filter passes).
  The state factor `1/P3` is completed by
  `add_penalization_state_correction`, a separate kernel over
  `ibmm bodyBlocks`. The fused predictor is untouched; do not fold the
  correction into it (traffic on the solver's occupancy-limited kernel, and
  the body-free bit-exactness argument). Gate: `validation/penalization/`.
- **The predictor's `if (skew)` guard is deliberate** (register schedule,
  −1.4 to −2.0 % of the step). `skew` must stay a runtime value.

## Items

**1. The penalization factor: adopt the AMPHIBIOUS rational form, for
momentum AND the Dirichlet scalar (decided 2026-10-01). DONE 2026-10-01.**

*Result* (numbers in the review's step-9 MEASURED block and the READMEs it
names). Implemented as specified below, with two additions: the scalar
branch is entered where `coef_p ≠ 0` only, and the cubic coefficient is the
constant multiply `PENAL_SIXTH`. Gates: unit test PASS; the decay gate third
order for velocity and scalar (2.77 / 2.88 / 2.94; 1.01e-4 at dt 0.04);
body-free 28/28 `max_abs 0` vs `~/step9_ref_binaries` (CPU + GPU); `les_ibm`
moves 2.5e-9 in 20 steps; cylinder Re 40 `C_D` 1.69234525 (1.6923452515 /
…19); CPU == GPU `max_abs 0` on the wavy-wall body case and on `les_ibm`
± refine; the S3/S4 scalar body gates and both `ibmwf` budgets hold.
Pre-registered cost, on the device the numbers below were taken on (the
RTX 3060 — the handout did not say): `ibm_mu` 0.40 ms (≤ 0.6), step
24.01 vs 24.01 ms for the implicit-Euler binary in the same sitting (within
1 %); A6000 +0.9 %. **The first implementation missed both** (`ibm_mu`
0.72 ms, step +2.5 %): `x/6.0d0` is a second fp64 divide per DOF. And the
`momentum` bucket's +0.3 ms is NOT the `exp`: it is the separate
state-correction pass, unchanged by the factor's form (the text below
attributed the whole +3.6 % to the exponential).

*Source.* `~/Codes/AMPHIBIOUS/src/modules/solver-timestep.cpl`, lines 4–11
and 58: `limiter(a) = 1/[(a/6 + 1/2) a + 1]`, `new = (B u + Δ)/(B + a)`,
and `limden` (the same factor in the projection); the temperature uses the
same line with `a/Pr`. That `B` is `a/(e^a − 1)` with `e^a` replaced by its
THIRD-ORDER Taylor polynomial. In our two-factor notation, x = λ dt_γ:

    P3    = 1 + x (1 + x (1/2 + x/6))
    state = 1/P3                       (replaces e^{-x})
    incr  = (1 + x (1/2 + x/6))/P3     (replaces (1 - e^{-x})/x; this is ibm%mu)
    state - incr = -x (1/2 + x/6)/P3   (the correction kernel's factor, no cancellation)

A rational function: positive and monotone for x ≥ 0, exactly (1, 1) at
x = 0, (≈0, 1/x) for a solid cell, and `incr = (1 − state)/x` holds
algebraically, so the steady fixed point `λ q = R` is unchanged.

*Why (measured 2026-10-01).* (a) Accuracy: `du/dt = f(t) − λu` integrated
through our RK3 substages, λ = 20, error at T:

| forcing | factor | dt 0.04 | 0.02 | 0.01 | 0.005 | order |
|---|---|---|---|---|---|---|
| f = 1 | implicit Euler | 1.2e-2 | 6.0e-3 | 2.9e-3 | 1.5e-3 | 1 |
| f = 1 | rational P3 | 1.0e-4 | 1.5e-5 | 2.0e-6 | 2.6e-7 | 2.8–3.0 |
| f = 1 | exact exp | round-off | | | | — |
| f = cos 2πt | implicit Euler | 1.1e-3 | 6.1e-4 | 3.2e-4 | 1.7e-4 | 1 |
| f = cos 2πt | rational P3 | 1.1e-4 | 2.5e-5 | 5.9e-6 | 1.4e-6 | 2 |
| f = cos 2πt | exact exp | 1.2e-4 | 2.7e-5 | 6.1e-6 | 1.5e-6 | 2 |

With a time-varying right-hand side the exponential and the rational factor
are indistinguishable: the RK3 coupling, not the factor, sets the error.
(b) Cost of the exp version, `les_ibm` (256/640 body blocks), GPU, 400
steps, `[output] profile`: `ibm_mu` 0.52 → 1.04 ms/step, `momentum` 3.82 →
4.16, step 22.69 → 23.51 ms (**+3.6 %**). The rational form is one divide
and a few multiplies. (c) No libm call on the device: CPU == GPU by
construction (the exp version happened to read max_abs 0 on one case; the
T3 `log()` class did not), no series branch, no `STACK:8` in the kernel.

*What to change.* `ibm.f90`: the bodies of `penal_incr_factor` and
`penal_state_factor`, plus a third pure function for `state − incr` (then
`add_penalization_state_correction` no longer reads `ibm%mu`). Nothing
else in the momentum path: the separate body-blocks pass and the untouched
predictor STAY. AMPHIBIOUS clamps its argument at `minImbArg = −0.616`
because its coefficients can be negative and P3 has a root at −1.596; ours
are ≥ 0 (checked: cylinder, les_ibm, wavy case files) — add a one-time
`min(coef) ≥ 0` check at init rather than a per-cell clamp.
`scalar.f90` (`scalar_transport_kernel`, the `SC_IBM_DIRICHLET` branch,
today `mus = 1/(1 + dt_γ coef_p/Pr)`): with x = dt_γ coef_p/Pr and `s0`
the start-of-substage value,

    ss = (s0 + Δ)*incr + (state - incr)*s0 + (1 - state)*s_body

which at x = 0 is today's arithmetic bit for bit (no branch needed).

*Gates.* `penalization_test`: new reference values for the rational form,
the identity `state + x incr = 1`, exact (1, 1) at 0. `validation/penalization/`
becomes a THIRD-ORDER convergence gate (expect the P3 row above; gate
order ≥ 2.7 and error ≤ 2e-4 at dt 0.04) and gets a scalar twin (uniform
scalar, uniform `source`, uniform `coef_p`). Body-free: 7-case minus
les_ibm and the 9-case suite `max_abs 0` (nofma, CPU + GPU) against
`~/step9_ref_binaries`. Body cases move at O(x⁴) — quantify on les_ibm
(20 steps) and re-read the cylinder Re 40 `C_D` (`run_f7_gates.sh steady`,
must stay 1.692345). CPU == GPU `max_abs 0` on the wavy-wall body case.
Profile `les_ibm` again: pre-registered, `ibm_mu` back to ≤ 0.6 ms and the
step within 1 % of 22.7 ms. Scalar: the S3 body gates of
`validation/scalar/README.md` re-measured (solid cell == body value to the
last bit, the heat diagnostics, `ibmwf180`/`ibmwf1000` closed-form budgets —
all steady statements, expected to hold). Then cut `~/step9b_ref_binaries`
and update the review's step-9 MEASURED block, `validation/penalization/README.md`
and CLAUDE.md.

**2. The block order: a minimum-surface bit order for the Morton key
(decided 2026-10-01). DONE 2026-10-01, measured on HoreKa (job 5174076,
`results_horeka_2026-10-01.md`): `base_jacobi` at 8 ranks, default
placement, 130.67 → 75.93 ms, `mpi_wait` 59.07 → 3.63 ms; init line
2262016 / 33792 as pre-registered; `rect_jacobi` and the refined config
unchanged (−0.4 %, +0.1 % at 8 ranks); fields `max_abs 0` between the two
orders at production flags (138 M points, rows matched). The pre-registered
77.7 ms ± 1 % was MISSED on the fast side: that number predates the
predictor guard. Block tax at 8 ranks from a default-placed `base_jacobi`:
0.997. 16 ranks (job 5174077): 67.40 → 41.52 ms, `mpi_wait` 27.90 → 3.61 ms,
101,376 cells across nodes, block tax 1.024 — prediction 7 hit.**

*Result.* As designed below, point by point:
(1) `blocks.f90`: `blk%keyOrder/keyPos`, `min_surface_key_order` (the rule),
`set_key_order`, `legacy_key_order`; `leaf_key` reads the table in xyz mode;
`morton_key` is gone. (2) Case file: attribute `block_key_order` (directions
1..3, most to least significant bit); absent = legacy order, and the table
reader checks that the recorded bits cover every leaf. (3) Restart:
`field_hdf5.c block_row_map` matches rows on (origin, level) for
`fdm_h5_read_field` / `fdm_h5_read_scalar` (so the body-force file and the
RANS scalars too); consecutive rows are still one hyperslab; the solver
prints a note when a restart file is in another order. (4) Mirrors:
`tools/partition_analysis.py --order file|legacy|minsurface --per-node N`
(it reproduces the table below from `lattice` input, no file needed);
`make_channel_restart.py` keeps the legacy order, which the row map makes
immaterial; mobygeom not ported (noted in its README). (5) xz mode keeps
its key and writes no record. (6) Init line
`partition: face cells shared across ranks N (R ranks), across nodes M (K nodes)`
(`blocks.f90 partition_surface`, `comm.f90 report_partition`; node id = the
smallest world rank sharing the node's memory).
Also needed, not foreseen: `tools/h5maxdiff`, `validation/prepare/compare_case.py`
and `compare_snapshots.py` compared rows BY POSITION and now match them on
(origin, level) — every comparison of a new-order file with a legacy-order
or mobygeom one would otherwise read as different fields (the P1 `ransgeom`
check did, until fixed).
Gates, all PASS (`validation/prepare/README.md`, last section):
`keyorder_test`; `run_gates_order.sh` 43/43 on CPU and on GPU (legacy vs
re-prepared file: identical fields at 1/4/8 ranks, identical tiles, restart
across the orders, the init line == the Python mirror, incl. two refined
cases); the 7-case + 9-case nofma suites vs `~/step9b_ref_binaries`
`max_abs 0`, CPU and GPU (32/32); P0 26/26, P1 15/15; S3 scalar body gates.
**One statement below is not what the rule gives:** "for cubic blocks on any
lattice this IS the ordinary Morton curve" holds on a NON-periodic lattice
(x y z x y z, x on top); a periodic direction's doubled first cut groups its
top two bits (periodic cube: x x y y z z; the 256 x 128 x 256 channel:
x x y z z x y z). The cut sizes are the same in the cases of the table (66 k
both ways); the rule was implemented as specified.
HoreKa: `horeka/exchange/submit_order.sh`, `PREREGISTERED_order.md`
(ref `5bdc5eb`, new `d2ac839`). SUBMITTED 2026-10-01 on `accelerated`
(`dev_accelerated` refused with `QOSMaxSubmitJobPerUserLimit`: four other
jobs of the user were already there): job **5174076** (2 nodes: the field
gate + 4 and 8 ranks, run directory `order_run`) and job **5174077** (4
nodes, 16 ranks, `order16_run`, starts after the first). Worktrees on
HoreKa: `moby-2to1-order` (`d2ac839`) and `moby-2to1-order-ref` (`5bdc5eb`),
transferred by `git push` into the HoreKa clone (branches `gate-order`,
`gate-order-ref`), NOT through GitHub — `origin/main` is still at `5c53a96`,
the commits of this session are local. Job 5174076 ran the same evening
(the scheduler's estimate had been 2026-10-05), 5174077 forty minutes
later; both are written up in `results_horeka_2026-10-01.md`.

*The defect* (`results_horeka_2026-09-30.md` section 3). The leaf order is
the Morton key of the finest-lattice block coordinates with x in the lowest
bit and z on top (`blocks.f90 morton_key`). Where the lattice has the same
number of bits per direction, the coarsest cut is therefore a z plane,
whatever the geometry: `base_jacobi` (nb unset, 2x2x2 blocks of 2048x88x96)
at 8 ranks on two nodes puts 1,442 k face cells on the node boundary instead
of 34 k, `mpi_wait` 3.3 → 58.8 ms, step +70 %; `--map-by node` on the same
binary and case file restores 77.6 ms.

*The design.* Keep a pure bit-permutation key (so the closed-form linear
split `zorder_owner`, the canonical row order and the exchange enumeration
are untouched), but choose the SIGNIFICANCE of every coordinate bit once,
from the geometry, by recursive bisection:

    box = the global domain in cells; bits_d = bit length of the finest-lattice tile count
    repeat until no bits are left:
        among directions with bits left, take the one whose bisecting plane
        carries the fewest cells = product of the other two box extents,
        counted TWICE if the direction is periodic and this is its first cut
        (a periodic bisection makes two interfaces); ties in x, y, z order
        give its highest remaining bit the next most significant key position
        halve the box in that direction

For cubic blocks on any lattice this IS the ordinary Morton curve (the long
direction's extra bits come first, then the interleave); it differs only
where today's order is arbitrary. It depends on grid size, block size and
periodicity alone — not on the rank count — so it stays canonical.

*Evidence* (face cells crossing a node boundary, 4 ranks per node, scored
with a 60-line model that should move into `tools/partition_analysis.py`
as an `--order` option):

| layout | current | best fixed x/y/z slot order | minimum-surface bit order |
|---|---|---|---|
| `base_jacobi` 8 ranks (2x2x2) | 1442 k | 34 k | 34 k (x0 y0 z0, high → low) |
| `base_jacobi` 16 ranks (4x2x2) | 1476 k | 101 k | 101 k (x1 x0 y0 z0) |
| `rect_jacobi` 8 / 16 ranks (64x4x4) | 34 k / 101 k | same | same |
| channel 256x128x256, nb 32, 8 ranks | 66 k | same | same |

A fixed slot order ties with the bit rule on all of these; the bit rule is
chosen because it also covers strongly non-cubic blocks (64x1024 domain in
16x256 blocks, 4 ranks: 192 cells against 1088 by hand count) and needs no
second mechanism.

*What it touches.* (1) `blocks.f90`: a per-run bit-position table in `blk`
built from grid, nb, refine mask and periodicity; `leaf_key` reads it. It
has exactly two users: the leaf-table builder's sort (line ~897) and the
reader's order check (line ~359). (2) Case file: record the order (an
attribute listing the directions from most to least significant bit); a
file WITHOUT it is read with the legacy order, so every existing case file
keeps working unchanged and identically. (3) Restart snapshots carry their
own `blocks` table in the writer's order: replace the row-by-row
cross-check (`fdm_h5_check_block_table`) by a row MAP matched on (origin,
level), so a snapshot written under one order restarts under another.
(4) Python mirrors of the key: `tools/partition_analysis.py`,
`tools/make_channel_restart.py` (mobygeom is retired; note it, do not
port). (5) `refine_dims = xz` keeps its own key (x,z Morton above y,
measured 2026-08-28) unless the rule is shown to reproduce it on the
production boundary layer. (6) Print at init, next to the body-block
fraction: face cells shared across ranks and across NODES (node ids as
`comm.f90 select_target_device` already derives them) — a bad distribution
must be visible, not silent.

*Gates.* Fields are independent of ownership by construction: the 7-case
and 9-case suites `max_abs 0` (nofma, CPU + GPU) against the then-current
reference set, min_channel 1 == 4 ranks, a legacy-order case file solving
identically to a re-prepared one, a restart across the two orders. The
init print must reproduce the table above (34 k at 8 ranks). HoreKa,
pre-registered: `base_jacobi` n=8 with the DEFAULT placement back to
77.7 ms ± 1 % with `mpi_wait` ≈ 3.3 ms; `rect_jacobi` and the refined
configs unchanged within 0.5 %; then 16 ranks, which was never measured
post-step-7.

*What it is not.* A greedy minimum of exchanged cells at every
power-of-two level, not a proof of minimum run time: load balance is the
other half (increment 7-5, the per-leaf weight column + weighted split,
still a plan entry). Until item 2 lands: no block tax at 8+ ranks from a
default-placed post-step-7 `base_jacobi`.

**3. The outlet pressure mode on the cylinder: RUN DOWN 2026-10-01, a
prototype fix measured, the decision to ship it is OPEN.**
`validation/cylinder/README.md`, last section, has the tables.
(a) It is not the projection: Chebyshev, plain Jacobi and red-black at 60
and 240 iterations give the same pressure to 1 % and the same velocity to
1e-3 on the same leg. (b) It is the outlet row: `apply_bc(outflow_copy)`
resets the outlet face to its interior neighbour every substage, the
projection re-supplies `u_out − u_nx = −dx dv/dy` each time, and the
incremental stored pressure integrates those `phi` in the last cell column,
with no `p_face = 0` to anchor it. The discrete outlet is "the outflow must
be parallel", stiffer the smaller dt (rms dv/dy in the last column over its
value at x = 15: 0.23 at dt 5e-3, 0.08 at 6.25e-4). Not the corner. (c)
Prototype on the local branch `proto/outflow-incremental` (`4908e7c`,
worktree `~/moby_outlet_proto2`): the outlet face takes its neighbour's
increment plus its own pressure gradient against the held outlet pressure
(one kernel over the boundary points, `save_outflow_gap`). Cylinder: stored
p rms 0.89 / 1.65 → 0.080 / 0.073 at dt 5e-3 / 6.25e-4, the two step sizes
agree to 6e-5 in velocity (0.15 before); freestream gates: oblique exact,
Poiseuille last-cell p 2.41e-3 → 1.247e-3 (exact 1.25e-3), Lamb-Oseen
reflected fraction 2.2e-2 → 9.8e-5.
**DECIDED 2026-10-01 (user): do NOT merge the prototype.** The predictor
is to own the outlet face (it can, without an extra call: old and new values
coexist inside `momentum()`), and the next session first investigates the
whole time step for the simplest correct treatment of outlet faces at low
and high indices, without working around errors —
`docs/next_session_outlet.md`. The prototype compensates for the outflow row
of `apply_bc` (a reset bolted onto the original "predictor never writes the
face" design) instead of correcting it; the projection's outlet correction
is sound. What the prototype would have needed, for the record: it changes
every outlet case. Work it needs first:
the outlet face at restart (reconstruct from continuity; it is not in the
field file), `ibm%mu` at a face next to a body, 2:1 interfaces on an outlet
face, then re-measure the cylinder (C_D, St, C_L amplitude: 0.40 with the
prototype against the old 0.51), `validation/blasius`, the boundary-layer
tutorial, `tutorials/naca/rans`, the sailplane. Cases without an outlet are
untouched by construction. More evidence taken the same day: the Blasius
gate passes with it and its top-entrainment deficit disappears (dv/v_edge
0.124 → 0.013; theta error 1.33 % → 1.65 %, gate 2 %); 100 time units of
shedding at production settings read St 0.1670 / C_D 1.4474 / CV C_L
amplitude 0.655 on main and 0.1744 / 1.4383 / 0.393 with the prototype — it
is a different flow at the 4 % level in St, not only a cleaner pressure
(this also closes item 4b's Re 100 leg: main's CV lift is inflated 28 % over
the old penalization series by the mode).
The clean-p protocol (zero `pn`, 300 steps) is NOT clean on main either: the
CV lift is still ringing at ±0.8 when it ends.

*The text of the handout, for the record:* (`validation/cylinder/README.md`,
last section). From a clean pressure, one time unit at Chebyshev niter 60
takes p rms from 0.30 to 0.9–1.6 and the CV lift to ±3.9 at dt 6.25e-4;
velocity differences between step sizes peak in the OUTLET CORNERS. It is
what limits any time-convergence statement on that case, it is identical
with both penalization factors, and it is the same family as the clean-p
protocol and the B0 outlet finding. A real investigation: is it the
Chebyshev bounds on N = 512 (`lmin` 2.5e-5: 60 iterations do not touch the
long waves), the outlet row, or the corner where two Dirichlet-p faces
meet? Start with plain Jacobi and red-black on the same 1-t.u. leg.

**4. Small, each under an hour.** (a) DONE 2026-10-01: the banner
`commit:` is the build's commit (`src/modules/build_info.c`; CMake defines
`MOBY_COMMIT` from `git describe --always --dirty --exclude *` on that file
alone, so a new commit relinks rather than rebuilds). (b) Re 40 DONE
2026-10-01: the `steady` gate of `validation/cylinder/README.md` is
re-measured on the control-volume series (`f7/forces_st_new.txt`, C_D 1.6924
± 6.7e-4 over the last 20 %). The Re 100 `strouhal` gate still quotes the
old penalization series; it needs a clean-p leg long enough for a spectrum,
which item 3's outlet pressure mode limits.

**5. The plan.** Review section 10: step 10 (implicit solid conduction for
conjugate scalars, two to three sessions) does not depend on anything open;
step 11 only if a viscous-limited production case shows up; step 12 last.
The one-by-one re-measurement pass (`validation/README.md` is the checklist;
the long tier is still at niter-6 values) resumes whenever a session has a
GPU and no plan step: discharged so far are turb180, Blasius, the cylinder
Re 40 steady drag and the `les_ibm` developed statistics (case a_wale).

## Landmines, new this session

- **A suite MODE label reused across sessions used to reuse the case file**
  (the drivers deleted `.${pfx}.case.h5`, the solver reads `${pfx}.case.h5`);
  fixed in the three `run_bitexact*.sh`, but any other driver that
  re-prepares should be checked for the same typo before a "stale case
  file" stop is believed.
- **A gate that was not run since step 7 may stop on the input echo**: the
  solve ini must carry every key the case file was prepared with (the
  cylinder inis lacked `stl_file` and `remove_solid = false`). Fix the ini,
  do not derive a second prepare ini.
- **HoreKa**: code `…/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-code`
  (main), run trees beside it; `dev_accelerated` runs ONE job per user at a
  time, so a queued moby job waits behind any other job of yours there.
  `h5dump`/h5py are not available on the login node: dump case files locally.
- **`validation/prepare/run_gates_step7.sh` is not a `max_abs 0` gate any
  more** (found 2026-10-01): its reference is the 7-0 inline binary, which
  predates F1 and the penalization factors. 18 pass / 26 fail with ANY
  current binary (control: `~/step9b_ref_binaries` reads the same pattern).
  Do not chase those failures; use the suites and `run_gates_order.sh`.
- **A field or case-file comparison by ROW POSITION is wrong across block
  orders** (2026-10-01): `h5maxdiff`, `compare_case.py`,
  `compare_snapshots.py` and the `run_gates_stl.sh` ransgeom check now match
  rows on (origin, level); `tools/compare_fields.py` always reassembled.
  Any NEW script that reads two block-layout files must do the same
  (`validation/prepare/h5rows.py`). An `h5maxdiff` binary built before
  2026-10-01 compares by position: rebuild it (the HoreKa submit scripts
  that build it only when missing would keep a stale one).
- **Local `mpirun -n 1` launches all bind to core 0** (the istmcetus
  landmine holds on the workstation): three suite legs in parallel ran at
  20 % CPU each until `OMPI_MCA_hwloc_base_binding_policy=none` was exported.
- **`[ibm] enabled` defaults to TRUE**: a `generic` case with no `[ibm]`
  section gets the default flat analytic wall and coefficient tiles (the
  order gates' "box" case is therefore a body case).
- `nvhpc 25.9` (workstation) folds the unguarded predictor to 78 registers,
  `25.3` (HoreKa) to 100: register counts are per compiler, read them from
  the machine that runs.
- Everything in `docs/next_session_after_step8.md` "Landmines" still holds.
