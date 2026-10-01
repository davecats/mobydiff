# Next session: after numerics review step 9

STATUS: **NOT STARTED (handout written 2026-10-01).** Steps 0–9 of
`docs/numerics_review_2026-09-26.md` section 10 are done. Read `CLAUDE.md`
(the "After step 8" bullet), then choose: items 1–3 are decisions and small
fixes that came out of the last session; item 4 is the next plan step.

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
- **The predictor's `if (skew)` guard is deliberate** (register schedule,
  −1.4 to −2.0 % of the step). `skew` must stay a runtime value.

## Items

**1. DECIDE: the block order of an nb-less case** (`results_horeka_2026-09-30.md`
section 3). One block per rank is numbered along the Morton curve with x in
the lowest bit; on an elongated domain with ranks filled node by node the
node boundary then cuts the LARGEST faces (`base_jacobi` 8 ranks: +70 %,
`mpi_wait` 3.3 → 58.8 ms; `--map-by node` restores it). Options: (a)
Cartesian rank order for `block_nb_auto` layouts (the leaf-table readers
check Morton order — relax for auto); (b) Morton bit significance by
block-face area (changes canonical ids wherever the lattice has equal bit
counts per direction; mobygeom mirror); (c) keep the order and document a
placement for the benchmark. Until decided: no block tax at 8+ ranks from a
default-placed post-step-7 `base_jacobi`. 16 ranks was not run.

**2. The outlet pressure mode on the cylinder** (`validation/cylinder/README.md`,
last section). From a clean pressure, one time unit at Chebyshev niter 60
takes p rms from 0.30 to 0.9–1.6 and the CV lift to ±3.9 at dt 6.25e-4;
velocity differences between step sizes peak in the OUTLET CORNERS. It is
what limits any time-convergence statement on that case, it is identical
with both penalization factors, and it is the same family as the clean-p
protocol and the B0 outlet finding. A real investigation: is it the
Chebyshev bounds on N = 512 (`lmin` 2.5e-5: 60 iterations do not touch the
long waves), the outlet row, or the corner where two Dirichlet-p faces
meet? Start with plain Jacobi and red-black on the same 1-t.u. leg.

**3. Small, each under an hour.** (a) The Dirichlet scalar penalization
(`scalar.f90`, `mus = 1/(1 + dt_γ coef_p/Pr)`) is the same first-order
factor F7 replaced for momentum: `penal_incr_factor` + a state correction,
gated by a scalar twin of `validation/penalization/`; the S3 gates compare
solid cells "to the last bit", so re-measure them. (b) The banner
`commit:` in `init.f90` is the literal `7aa1c7b`: make it a build-time
define or delete it. (c) `validation/cylinder/forces_re40.txt` /
`forces_re100.txt` are still the old penalization series (README says so);
the Re 40 one can be regenerated from `run_f7_gates.sh steady`.

**4. The plan.** Review section 10: step 10 (implicit solid conduction for
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
