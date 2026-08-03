# Next session: runtime C_L/C_D from the control-volume momentum budget

STATUS 2026-08-03: **DONE and COMMITTED.** The p_inf bug is fixed, the budget
is validated on `validation/cylinder` (Re 40 steady + Re 100 unsteady), the
penalization integral is REMOVED (`cv_box` now required), the keys are
documented in `docs/configuration.md`, and the 8-case bit-exactness suite
passes on CPU and GPU. Instrumentation stripped. Open items are the two
regeneration follow-ups listed under item 7/8 below.

THE FIX: `cv_pinf()` is a first reduction pass over the upstream border
alone (`NPI = 2`: p*dA and dA, through the same deterministic reduction, so
p_inf is rank-count independent too); the flux pass then subtracts it PER
FACE. The assembled-signed-area form and its `sAc/sAl/pA/aW` accumulators
are gone; `cv_reduce()` is the shared scatter/allreduce/ordered-sum helper
(`penalization_forces` reuses it).

RESULTS (`.cv*` scratch inis/logs left in validation/cylinder):

- Re 40, clean-p restart `cvpz_20301.h5` at niter = 60, box 4-8 x 6.5-9.5,
  step 20310: runtime **C_D = 1.70387475, C_L = -1.1986e-3** (penalization
  reference 1.6924 / 4e-4; the box scatter below covers the difference).
- CROSS-CHECK, exact: an independent reassembled-global-plane evaluation
  (scratchpad `cv_global.py`, third implementation) reproduces the runtime
  to **9 significant digits per border** (FN/FT, after the x4 span factor)
  and 1.703875 / -0.001199 overall.
- The offline `cv_forces.py` reads 1.70690 / -0.00307 on the SAME snapshot
  and box. **That 0.17 % gap is a defect in the OFFLINE TOOL, not the
  runtime**: its `ut_face` interpolates the tangential velocity only along
  the border NORMAL (`0.5*(A(T,t,i-1) + A(T,t,i))`), leaving a half-cell
  collocation offset in the TANGENTIAL direction; the runtime uses the
  correct 4-point average. Replicating exactly that one choice in the
  independent script reproduces the tool to every printed digit
  (+1.706898 / -0.003070). FIXED (user decision): the tool now collocates
  u_t at the tangential cell centre (`Tc()`) and reads 1.70405 / -0.00141,
  i.e. 0.011 % from the runtime — the rest is its first-order block-edge
  gradients, unavoidable without halos. The tutorial/validation numbers
  predate the fix and are annotated in their READMEs.
- Uniform flow (empty.ini + cv_box): 8.7e-16 with the box on block
  boundaries, EXACTLY 0.0 with it one cell inside (interior faces).
- Determinism: 1 vs 4 CPU ranks BYTE-IDENTICAL (empty and cylinder);
  CPU vs GPU identical to every printed digit.
- Box scatter (steady Re 40, runtime stencils): 1.6982 (5-7.5 x 7-9),
  1.7039 (4-8 x 6.5-9.5), 1.7256 (3-9 x 5-11), 1.7646 (2-12 x 4-12) —
  longer borders accumulate more collocation/gradient error, same trend as
  the README's offline numbers. Use a TIGHT box.
- d/dt term (item 4, was untested): Re 100 shedding, two runs from the same
  clean-p restart at niter = 60 (identical trajectory), one with cv_box and
  one with the penalization statistic. With the term: mean C_D 1.4469 vs
  1.4486 (**0.12 %**), C_L range [-0.572, +0.486] vs [-0.536, +0.483],
  pointwise rms deviation 0.101 (~20 % of the 0.51 amplitude; 0.050 at a
  2-sample shift, so much of it is a phase offset in the differenced term).
  WITHOUT the term the statistic is meaningless: C_L rms error 1.49 and the
  mean flips sign (+0.503 vs -0.031). The unsteady term is 1.15 in C_L
  units here — the budget is a difference of two O(1) terms.
- niter (item 3): from a CLEAN p, 2000 steps at production niter = 6 hold
  C_D = 1.692 +- 0.001 with NO growth. The catastrophe is the ACCUMULATED
  drift of a very long run, not niter itself: restarting the committed
  `cyl_re40_20001.h5` (20k steps at niter = 6, pn std 422 with a large
  SPATIALLY VARYING mode) still gives C_D = 261 AFTER the fix — no constant
  subtraction can rescue it. **Any CV-budget run must start from, or
  periodically re-converge, a clean pressure** (README recipe: zero pn,
  restart at niter = 60).

## Goal

Replace the airfoil case's runtime force statistic. It is currently the
penalization integral `F = int coef*u dV`, which is exact bookkeeping ONLY
when the solid interior is present. Production runs set `[blocks]
remove_solid` (the default, and unavoidable at cost), so the buried core is
gone and its share of the pressure-dominated loading is outside the coef
bookkeeping: the statistic UNDER-READS. The tutorial works around this by
post-processing snapshots with `tools`-style `cv_forces.py`. Move that
budget into the flow case so it is available at runtime.

    F = - d/dt int_V u dV - oint_S [ u (u.n) + (p - p_inf) n - tau.n ] dS

over a box around the body, spanning the full periodic extent (the two span
faces cancel, so only four lateral borders are integrated). The unsteady
term is differenced between consecutive SAMPLES, so all machinery stays in
the flow case and no extra state is needed.

## What is implemented

- `after_step` gained a `turb` argument (the CV viscous stress needs
  `turb%nut`: the wake crossing the downstream border carries a large eddy
  viscosity). Touched `flow_case_base.f90`, `generic_flow.f90`,
  `channel_flow.f90`, `moby_solve.f90`. `channel_stats` has its OWN
  `after_step` interface (no `ibm`) and is deliberately untouched.
- `airfoil_flow.f90`: `cv_forces()` (the budget), `cv_pinf()` (the p_inf
  pass), `cv_reduce()` (the deterministic reduction), `snap_border()`.
  `penalization_forces()` is GONE (item 7).
- `[case.airfoil] cv_box = c0 c1 l0 l1` is REQUIRED whenever forces are
  sampled; `force_sample_interval = 0` is the only way to run without it.
- Borders are snapped at setup to a face of the COARSEST level they cross
  (a coarse node is a node at every finer level), so a border never lands on
  a fine-only face and the box stays closed. Two `max`-reductions: the
  minimum level via its negative, then the coordinate.
- Only INTERIOR faces (i = 1..nb, the west face of cell i) are integrated,
  so each physical face is counted exactly once by its east-side block at
  that block's own level.
- Rank-count independence preserved: per-block partials scattered into the
  global block table, exact allreduce, ordered final sum.

## THE BUG (FIXED 2026-08-03 — kept for the record)

`p_inf` is subtracted from the ASSEMBLED signed areas (`fc = acc(3) +
p_inf*acc(7)`) instead of per face. The reasoning was that `oint n dS = 0`
for a closed box so the constant cancels analytically — which it does, IN
EXACT ARITHMETIC. In floating point it is catastrophic cancellation.

Measured on `validation/cylinder` Re 40 (`niter = 6`, the documented pn-drift
family): stored `pn` runs **-1200 .. +1207, mean +130**, while the force is
O(0.2). Per-border FN/area came out at 174.7. Runtime output: C_D = 261,
C_L = -1398, i.e. garbage. The offline `cv_forces.py` subtracts `p_inf` PER
FACE, so every term it sums is O(p - p_inf) ~ O(1) and the cancellation is
benign.

Same bug, milder, on the NACA nose case: there p is O(1e-2) against forces
of O(1e-3), so only ~1 digit is lost — C_L 0.483 vs 0.5155 (6 % low) and
C_D -0.0154 vs +0.0116 (sign flipped).

**Fix**: subtract `p_inf` inside the face loop. That needs it before the
flux sum, so a two-pass reduction — pass 1 accumulates `p*dA` and `dA` on
the upstream border and allreduces to get `p_inf`; pass 2 does the flux
integral with `(p - p_inf)`. One extra small allreduce per sample.
Consider subtracting the reference dynamic pressure / momentum flux too,
for the same conditioning reason; `p_inf` per face is the necessary part.

## What is already verified correct (do not re-litigate)

- **Formula**: checked term by term against `cv_forces.py` — `fn = u_n^2 + p
  - tau_nn`, `ft = u_n u_t - tau_nt`, `Fc -= sgn*fn` on chord borders,
  `Fc -= sgn*ft` on lift borders, `sgn = -1` on low borders.
- **Geometry/bookkeeping**: every border's integrated area is EXACT
  (NACA: 0.5625 / 0.75; cylinder: 0.75 / 1.0). Face ownership, areas and the
  2:1 handling are sound.
- **Shear and tangential terms**: FT agrees with the reference to 0.1-1 %
  on all four borders — the interpolations, gradients and `nut` are right.
- **Halos are current** at `after_step`: `pressure_solver.f90:162-164`
  exchanges `[VAR_U, VAR_V, VAR_W, VAR_P]` on the final Jacobi iteration,
  after `jacobi_apply` writes the pn increment (line 327). Stale halos were
  the first hypothesis and are EXCLUDED.
- **Blended 2:1 ghosts are not the problem**: switching `p_face` to a
  one-sided extrapolation from the owning side (never touching a halo)
  changed FN by ~2e-6. Almost no border face sits at a block edge.

## Test harness: use the CYLINDER, not the NACA case

`validation/cylinder`, Re 40: single level (no 2:1 interfaces), STEADY
(so d/dt ~ 0 and the surface integral alone must match), 512x512x8, minutes
per run, converged state `cyl_re40_20001.h5` (t = 100.0, step 20001) is
present. It is also `span = z` (liftDim = y), so it exercises the OTHER span
orientation from the NACA case.

    cv_box = 4.0 8.0 6.5 9.5          # in a copy of cyl_re40.ini
    [restart] file = cyl_re40_20001.h5
    t_final = 100.05                  # the restart is at t = 100, not 20
    field_interval = 10               # snapshot at the same step as a sample
    runtime_file = /tmp/...            # do NOT clobber the committed forces_re40.txt

reference for the identical box (the tool gained `--span-z` for this):

    python3 ../../tutorials/naca/rans/postProcess/cv_forces.py \
        cyl_re40_20010.h5 --boxes 1.5 --nose 5.5 8.0 --span-z

Cylinder geometry: D = 1 centred at (6.0, 8.02); expected C_D ~ 1.69 at
Re 40 (validation/cylinder README).

**Compare at the SAME STEP.** An early NACA comparison was runtime@640020 vs
reference@640000 — 20 steps apart, and stored `pn` drifts over those steps
while velocities do not. Write a snapshot at a force-sample step.

## After the fix

1. DONE — cylinder re-validated (see STATUS; the residual vs the offline
   tool turned out to be the tool's tangential-velocity collocation offset,
   not a block-edge stencil difference).
2. DONE — uniform flow: 8.7e-16 / exactly 0.0.
3. DONE — from a clean p, 2000 steps at niter = 6 hold C_D to +-0.001.
4. DONE — d/dt exercised on Re 100; essential and correct to 0.12 % in the
   mean C_D.
5. DONE — 1 vs 4 ranks byte-identical, CPU == GPU.
6. DONE — stripped to `NCV = 4`; the Re 40 gate number is unchanged
   (1.70387475) after stripping and rebuilding CPU+GPU.
7. DONE (user decision 2026-08-03: do NOT keep it) — `penalization_forces`
   is DELETED and `cv_box` is REQUIRED whenever forces are sampled;
   `force_sample_interval = 0` is the explicit opt-out and the only way to
   run the airfoil case without a box. The error is raised in
   `setup_after_grid`, which `moby_prepare` never calls, so prepare-only
   inis are unaffected. All 16 committed solver inis gained a `cv_box`
   (margin 1.5c around the profile = what `cv_forces.py --boxes 1.5`
   integrates; cylinder 4-8 x 6.5-9.5; `empty.ini` deliberately one cell
   inside the block boundaries so the uniform-flow closure is EXACTLY 0).
   Recover the removed integral from history if ever needed — do not
   reinstate it as a fallback.
   CONSEQUENCE: `validation/cylinder`'s committed `forces_re40.txt` /
   `forces_re100.txt` and the naca0012 / sd7003 published C_L/C_D are still
   penalization series. They are annotated as such in their READMEs and
   need regenerating before being quoted as control-volume numbers (the
   Re 100 one needs the clean-p protocol). `--keep-buried` was load-bearing
   only BECAUSE of the penalization integral and is no longer needed for
   forces, but it STAYS (user, 2026-08-03): the buried interior is wanted
   for future moving boundaries and for a heat equation in the solid.
8. DONE — bit-exactness suite, **8/8 PASS on CPU AND GPU** (max_abs
   EXACTLY 0 on every dataset present in the snapshots, `-Mnofma` /
   `-gpu=nofma` both sides). The standard 7 cases plus an 8th,
   `validation/cylinder/empty.ini`: the airfoil case is where the code
   actually changed, and it also byte-compares the forces file, which
   gates the `penalization_forces` reduction refactor (byte-identical).
   MECHANICS, since the working tree was already dirty: the reference
   binary comes from a detached `git worktree add --detach <dir> HEAD`
   (HEAD = 1bb4f10), NOT from the current tree; scratchpad
   `build_nofma.sh cpu|gpu <src> <build>` mirrors compile.sh's cmake call
   with `-DCMAKE_Fortran_FLAGS=-Mnofma` (+ `-DOPENMP_OFFLOAD_FLAGS="-mp=gpu
   -gpu=nofma"`), copying the MPI wrapper flags out of an existing
   CMakeCache; `gate_bitexact.sh cpu|gpu [case ...]` runs both sides and
   compares every dataset common to the two snapshots at tolerance 0.
   NOTE: `compare_fields.py` needs its dataset names as POSITIONALS BEFORE
   `--tolerance` (argparse rejects them after it).

Remaining before commit: item 7, and README/docs for the `cv_box` key.

## Landmines

- **The NACA converged state was destroyed** (`c11_nose_640000.h5`, and the
  L11 baseline before it) by a careless `rm -f c11_nose_6400[0-9]*.h5` that
  was meant to remove one stray test snapshot. There is no restart point for
  `tutorials/naca/rans`; regenerating it is `run_case.sh scratch`, ~60-70
  GPU-h on an RTX 3060. The committed tutorial and all its published numbers
  are unaffected (extracted before the deletion). Its prepared case file
  `assets/geometry/ibm_coeff_c11_nose.h5` survives.
- `module load toolkits/nvhpc/25.9` does NOT work on this machine (only
  <= 25.5 exists) and silently leaves `mpirun` as the system Open MPI, which
  dies at MPI_INIT. Set PATH explicitly to
  `/opt/Nvidia/nvhpc/Linux_x86_64/25.9/{compilers,comm_libs/12.9/hpcx/latest/ompi}/bin`.
- `moby_prepare` is MPI-parallel ONLY — the CPU build has no OpenMP, so the
  `!$omp parallel do` in `classify_block_geometry` is inert and
  `OMP_NUM_THREADS` does nothing. Give it MANY RANKS (`-n 20 --bind-to none`);
  the NACA nose case is ~40 min that way and hours at `-n 4`.
- `pkill -f <pattern>` matches the invoking shell's own command line when the
  pattern appears in its arguments — it will kill the very command issuing
  it. Kill by PID (`pgrep -x`).
- Mapping a component of the polymorphic `this` into a `target` region is not
  portable; `cv_forces` copies `cv_box` into a plain local first.
- The snapshot `pn` dataset is `q(:,:,:,VAR_P)` and accumulates `phi/dt` over
  the whole run — it carries a large velocity-neutral drift at low `niter`.
  Any pressure-based statistic must be referenced, never used raw.
