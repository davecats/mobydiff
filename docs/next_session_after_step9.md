# Next session: after numerics review step 9

STATUS: **NOT STARTED (handout written 2026-10-01, items 1 and 2 DECIDED
the same day).** Steps 0–9 of `docs/numerics_review_2026-09-26.md` section 10
are done. Read `CLAUDE.md` (the "After step 8" bullet), then work items 1
and 2 in that order: both are decided designs with their evidence below,
each is a session's worth including gates. Item 3 is an open investigation,
item 4 small fixes, item 5 the plan.

## What holds now

- **Reference set `~/step9_ref_binaries`** (PROVENANCE inside): F1, the
  predictor guard, the exact penalization factor. Body-free cases are
  `max_abs 0` (nofma) against `~/step8_ref_binaries` too; body cases are not.
- **Penalization.** `ibm%mu = (1 − e^{−x})/x`, x = λ dt_γ, is the factor on
  every velocity INCREMENT (predictor, projection correction, SGS /
  body-force / band-filter passes). The state factor `e^{−x}` is completed
  by `add_penalization_state_correction`, a separate kernel over
  `ibmm bodyBlocks`. The fused predictor is untouched; do not fold the
  correction into it (traffic on the solver's occupancy-limited kernel, and
  the body-free bit-exactness argument). Gate: `validation/penalization/`.
  Item 1 replaces the two factor FUNCTIONS (not this structure).
- **The predictor's `if (skew)` guard is deliberate** (register schedule,
  −1.4 to −2.0 % of the step). `skew` must stay a runtime value.

## Items

**1. The penalization factor: adopt the AMPHIBIOUS rational form, for
momentum AND the Dirichlet scalar (decided 2026-10-01).**

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
(decided 2026-10-01).**

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

**3. The outlet pressure mode on the cylinder** (`validation/cylinder/README.md`,
last section). From a clean pressure, one time unit at Chebyshev niter 60
takes p rms from 0.30 to 0.9–1.6 and the CV lift to ±3.9 at dt 6.25e-4;
velocity differences between step sizes peak in the OUTLET CORNERS. It is
what limits any time-convergence statement on that case, it is identical
with both penalization factors, and it is the same family as the clean-p
protocol and the B0 outlet finding. A real investigation: is it the
Chebyshev bounds on N = 512 (`lmin` 2.5e-5: 60 iterations do not touch the
long waves), the outlet row, or the corner where two Dirichlet-p faces
meet? Start with plain Jacobi and red-black on the same 1-t.u. leg.

**4. Small, each under an hour.** (a) The banner `commit:` in `init.f90`
is the literal `7aa1c7b`: make it a build-time define or delete it.
(b) `validation/cylinder/forces_re40.txt` / `forces_re100.txt` are still
the old penalization series (README says so); the Re 40 one can be
regenerated from `run_f7_gates.sh steady`.

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
- `nvhpc 25.9` (workstation) folds the unguarded predictor to 78 registers,
  `25.3` (HoreKa) to 100: register counts are per compiler, read them from
  the machine that runs.
- Everything in `docs/next_session_after_step8.md` "Landmines" still holds.
