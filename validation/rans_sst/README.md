# RANS T2+T3+T4 gates — k-omega SST transport (resolved walls, wall functions, gamma-Re_thetat transition)

Physics gates for IDDES phases T2, T3 and T4 (docs/history/next_session_iddes.md).
The channels run in minutes on 4 CPU ranks; the two IBM channels (`ibm180`,
`ibm180wf`, 200 000 steps each) take ~1.5 h each on the RTX 3060:

```bash
cd validation/rans_sst
RANKS=4 BIN=/abs/build_cpu/moby_solve ./run_gates.sh            # all cases; or one name:
RANKS=1 BIN=/abs/build_gpu/moby_solve ./run_gates.sh ibm180
./run_gates.sh t3                 # just the T3 wall-function set
./check_gates.sh                  # python3 + h5py + numpy
```

Cases (each writes `<prefix>_<laststep>.h5` + `<prefix>.log`; `turb180` and
`ibm180wf` also write cadence snapshots, every 20 000 / 100 000 steps, for the
restart gate below).

RE-VALIDATED 2026-10-21, commit 532ef54 (banner), production flags, the
channels on 4 workstation CPU ranks (seconds to 17 min each; refine180 85 min
at 440 286 steps), `ibm180` / `ibm180wf` on the RTX 3060 (1 rank, 200 000
steps, 87 min / 26 ms per step): `check_gates.sh` ALL PASS. Settings changed
for the rerun: the twelve channels that declared plain Jacobi `niter = 6`
moved to `accel = chebyshev`, `niter = 12`; all sixteen inis `cflmax = 0.5 ->
1.0` and lost the ignored `sor` line (the step is Peclet- or dtmax-bound on
every case: turb180 132 655 steps against 132 565 before). Old value -> new
value where it moved:

- `laminar.ini`   gate (a): Re_tau 10 channel, SST on, tu = 1 %. Parabola to
                  2.444e-4 of u_max (unchanged); k decayed to 4.7e-46 (was
                  3e-47). (At Re_tau 30 / tu 5 % the no-transition SST finds
                  its weakly-turbulent branch, `lam30` below, so the gate
                  probes the genuinely subcritical regime.)
- `turb180.ini`   gate (b): developed channel Re_tau 180 (natural y, y+_1 ~
                  0.5). u_tau = 1.0000 (1.0008 before F1); log line to 5.0 %
                  (4.9 %); U+ centreline 18.16 vs DNS 18.20 (0.2 %); bulk U+
                  15.36.
- `turb395.ini`   gate (b): Re_tau 395. u_tau 1.0000 (1.0016); log line to
                  6.6 % (6.5 %); U+ centreline 20.17 vs DNS 20.13 (0.2 %,
                  was 20.16).
- `ibm180.ini`    gate (c), THE key IBM gate: the les_ibm off-grid file-IBM
                  plane channel (walls mid-cell at y = 0.259375 / 2.259375,
                  uniform y, y+_1 ~ 2-3) through the IBM wall treatment
                  (dwall_blocks + wall-cell omega pinning). Log line to 4.3 %
                  (unchanged), u_tau 0.966 (unchanged), U+ centreline 16.91.
                  Restarts from `../channel_interface/les_ibm/IC.h5`, whose
                  WALE `nut` is now read back (2026-10-21) while k and omega
                  take the cold-start IC -- the first steps only; the fixed
                  point is the same. Needs ../rans_geometry/ibm_coeff_blocks_l1.h5.
- `refine180.ini` gate (d): uniform-y Re_tau 180 with both wall bands
                  refined to level 1 (2:1 interfaces at y = 0.25 / 1.75). NO
                  interface band: jump / local-variation ratios 0.58 (k),
                  1.11 (omega), 0.90 (nut), 0.73 (u) -- unchanged; core
                  profiles match the RESOLVED turb180 reference to 0.6 % (u,
                  was 0.5), 0.5 % (k), 2.8 % (omega), 1.7 % (nut).
- `base180u.ini`  informational coarse-wall control (uniform y+_1 ~ 2.8, no
                  refinement): its under-resolved sublayer feeds a spurious
                  core-k plateau -- the demonstration of WHY the wall bands
                  get refined. U+ centreline 18.78, log line 8.7 %. Not a
                  pass/fail gate.

T3 wall-function cases (`[rans] wall_treatment = wall_function`; all
checked with `--mode wallfn` against the RESOLVED turb180 field, which
carries the DNS centreline anchor to 0.2 %, so the anchor is transitive):

- `wf180_y30/45.ini` T3 gate (a): turb180 coarsened to uniform ny = 6 / 4
                  (y+_1 = 30 / 45; ny = 6 is not nb-divisible, so it runs
                  with `[blocks] nb` unset -- rank-box blocks). Implied U+
                  centreline vs DNS 18.20: 1.2 % / 0.7 % (unchanged); u_tau
                  from the delivered wall stress (nu + nut_1) U_1/y_1 =
                  1.0000. The mid-gap row (y+ 90) sits 3.02 % above the
                  resolved profile (2.97 % at Jacobi 6): the documented
                  profile-shape bias of the buffer-anchored low-log rows,
                  whose `--tolerance-center` had been set AT the measured
                  value and is 0.035 since 2026-10-21 (check_gates.sh). y45:
                  3.75 % log-region, 0.94 % near-centre.
- `wf180_y05/y15/y22.ini` T3 gate (b): graceful degradation across the
                  buffer range (y+_1 = 5 / 15 / 22.5). Log-region rows vs the
                  resolved reference 4.3 / 7.9 / 5.5 % (tol 8), near-centre
                  3.2 / 3.7 / 3.7 % (tol 4): a mild overshoot, NO
                  double-counting dip (that would be a deficit). First cells
                  below y+ 30 sit 6-19 % above the resolved profile
                  (informational, the log-line error at the anchor cell).
- `ibm180wf.ini`  T3 gate (c): ibm180 through the wall-function blend (200k
                  steps, RTX 3060). The k-based y+ keeps every wall cell on
                  the viscous branch, whose arithmetic is exactly the
                  resolved treatment: at the mid-run cadence snapshot (step
                  100 000) the profile already matches the T2 resolved
                  ibm180 field to 1e-4 at every row; the converged
                  comparison and the restart gate are in the GPU block below.

RESTART == CONTINUOUS AT THE DEVELOPED STATE (the gate of the 2026-10-21
`nut` read-back on a converged RANS field; the cold-start gate of that day
was ten steps): `turb180` restarted from its step-120 000 cadence snapshot to
t_final reproduces the continuous run's final field at **max_abs 0** on
`un vn wn pn k omega nut` (4 CPU ranks, `tools/h5maxdiff`). The same for
`ibm180wf` from step 100 000 on the GPU: see the GPU block.

THE GPU BLOCK (RTX 3060, 1 rank, commit 532ef54, production flags):
`ibm180` 200 000 steps in 87 min and `ibm180wf` in 88 min (26 ms per step);
at the converged state `ibm180wf_200000.h5` == `ibm180_200000.h5` at
**max_abs 0** on `un vn wn pn k omega nut` -- the wall-function blend on a
wall whose every cell sits on the viscous branch IS the resolved treatment,
bit for bit (the 2026-07-08 record read 4.5e-16); `ibm180wf` restarted from
its step-100 000 cadence snapshot to step 200 000 reproduces the continuous
run at **max_abs 0** on the same seven fields (the nut read-back on a
converged RANS field, on the GPU); the 9-case suite at nofma against
`~/pass2_ref_binaries` (d821341): max_abs 0 on all nine cases, CPU (4 ranks)
and GPU, after the ini moves above.

Gates (d)/(e) of the T3 record (resolved-mode bit-exactness, wall-function
determinism) live in the suites: `turb180`, `wf180_y30` and `lam30t` are
three of the nine cases of `validation/scalar/run_bitexact.sh`, run after
every change at nofma against the reference set (max_abs 0, CPU and GPU);
wall-function 1 rank == 4 ranks is exact, CPU vs GPU <= 2e-13 (the `log()`
intrinsic of the wall-function branch differs by an ulp between host and
device libm -- resolved mode remains exactly CPU == GPU).

RESOLVED 2026-08-05 (found 2026-08-04 while gating the passive scalars'
thermal wall function, validation/scalar increment S5a): **a cold-started
RANS run was not rank- or nb-independent when the decomposition split a
direction in which the velocity is nonzero.** Symptom as recorded:
`wf180_y30.ini` with no `[scalar]` section, 1 rank vs 4 ranks (the default
dims, which split x and z) read `un` max_abs 8.800384e-03 / `k`
9.035115e-02 / `omega` 1.238466e+00 after 20 steps, and already `un`
7.6e-04 after ONE step; a pure z split (`[mpi] dims = 1 1 4`) was EXACT.

ROOT CAUSE — the k initial condition read an unfilled halo.
`init_rans_transport` (rans.f90) sets `k = 1.5 (tu/100 |u|)^2` with the
cell velocity interpolated from the two staggered faces,
`uc = ½(q(i) + q(i+1))`, so at a block's LAST interior cell it reads
`q(nb+1)` — a halo. `moby_solve.f90` called `apply_bc` +
`exchange_halos` only AFTER the whole init block, so that halo was still
zero: `k` came out a factor **4** low (`|u|` halved) on the last plane of
every block, and `omega = k/(nut_ratio nu)` fell with it until the viscous
limb `6 nu/(beta1 y_eff^2)` took over. Measured on the 8x6x8 case with
`dt = 1e-14` (i.e. the initial state itself): `k` 0.172673 vs 0.690693 =
exactly 1/4, `omega` 16.0 (the viscous limb, exactly) vs 22.204923, at
x = 7 on one rank and at x = 3 AND x = 7 on two x-ranks — one bad plane
per block, which is what made the answer depend on the decomposition.
Only x showed it because this channel's `v` and `w` are zero in the IC, so
the corresponding y/z halo reads contributed nothing; a case with a
nonzero `v`/`w` initial field (a restart missing its k/omega datasets, an
inlet case) was wrong on the high face of every block in all three
directions. Neither of the two suspects recorded above was involved.

FIX (`src/moby_solve.f90`): `apply_bc(..., outflow_copy=.true.)` +
`exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])` now run BEFORE the
`[rans]` init block. The calls that follow it are unchanged and
idempotent — both write only ghosts and halos, from interior data that
nothing in between modifies — so every non-RANS case stays bit-exact.

SECOND HALF OF THE FIX, and the reason to run the GPU suite rather than
argue from the CPU one: **those two calls write the DEVICE copy of `blk%q`
(mapped by `enter_block_data`), while `init_rans_transport` is HOST code
reading the host copy** — so on the GPU the fix above was a pure no-op.
Measured 2026-08-05: the fixed GPU binary reproduced the pre-fix GPU result
BIT-FOR-BIT on turb180 / wf180_y30 / lam30t and kept the x-dependent `k`
(x-spread 1.492e-02 after 20 steps, against 0.000e+00 on the CPU). The
completed fix adds `!$omp target update from(blk%q)` after the exchange,
under `USE_OPENMP_OFFLOAD` — the same pattern as the
`target update from(ibm%coef)` a few lines below. After it the GPU `k` is
x-uniform (x-spread **0.000e+00**) and turb180's 20-step deviation from the
pre-fix binary is IDENTICAL to the CPU's to every digit (un 2.2814559502704945e-02,
k 2.8996297594464926e-01, omega 5.2618740388259901e+00), i.e. the two paths
now agree on the corrected initial condition.

WHAT IT CHANGES: only the k/omega INITIAL CONDITION of a cold-started RANS
run (a restart carrying k/omega overwrites the IC, and was never
affected). The converged answer is unchanged: `wf180_y30` re-run to
`t_final` with the fixed binary reproduced the gate to EVERY PRINTED DIGIT
— the RANS fixed point does not remember the IC. Cold-started RANS field
snapshots taken before 2026-08-05 are therefore not reproducible
bit-for-bit with the current binary; their physics is.

GATE: `wf180_y30.ini` at 20 steps, 1 rank vs `dims = 4 1 1`, `1 1 4` and
`2 1 2`, all **max_abs 0** on `un vn wn pn k omega nut` (nofma CPU).
Independently re-confirmed on an IBM + RANS + scalar cold start
(`../scalar/ibmwf180.ini`, 8 blocks, 20 steps): 1 rank == 4 ranks
**max_abs 0** on `un vn wn pn k omega nut theta`.

T4 gamma-Re_thetat transition cases (`[rans] transition = true`,
Langtry & Menter 2009 = OpenFOAM kOmegaSSTLM, resolved walls only; run
group `t4`). The canonical flat plate is DEFERRED: an inlet can be
composed from the existing faces (Dirichlet velocity + Neumann pressure,
zero-gradient outlet), but the RANS layer is not inlet-aware yet (a
Dirichlet inlet classifies as a no-slip wall in domain_face_is_wall /
dwall, and the scalars lack inlet ghost values) and no inflow/outflow
case is validated — so the gates are channels. Re-validated 2026-10-21
(the block above) — all PASS:

- `laminart.ini`  gate (a) subcritical control: Re_tau 10 / tu 1 % with
                  transition on stays laminar exactly like T2 (parabola
                  2.444e-4, k -> 9.5e-16, gamma pinned at its 0.02 floor).
- `lam30.ini`     transition-OFF control at Re_tau 30 / tu 5 %: the plain
                  SST weakly-turbulent branch (parabola off by 12.2 %, k
                  self-sustained at 0.37; unchanged). Informational —
                  documents WHY lam30t is discriminating.
- `lam30t.ini`    gate (a), THE discriminating gate: same conditions with
                  transition ON laminarize (parabola 1.5e-3, was 1.6e-3;
                  wall-layer gamma 0.024, mean-k peak 8.7e-3, was 9.1e-3 =
                  40x below the branch; the k residual is the gamma-floor
                  2 % of P_k, hence --k-max 0.02).
- `turb180t.ini`  gate (b): turb180 with transition on preserves developed
                  turbulence: gamma >= 0.999 for y+ >= 30, U+ centreline
                  18.45 vs DNS 18.20 (1.4 %, tol 2 %; was 18.44), u_tau
                  1.0000 (1.0009); the kappa/B log-line fit is 6.7 % (6.6 %;
                  LM's sublayer D_k x 0.1 coupling lifts the low-log rows
                  ~2 % over T2 — gated at the turb395 tolerance 0.08; the
                  DNS anchor is the hard criterion).
- `t4_front_check.py` STEP-0 evidence for keeping first-order upwind on
                  the transported scalars: the gamma front is wall-normal
                  and the cross-front upwind diffusivity max|v| dy/2 is
                  7.0e-5 (lam30t) / 4.3e-14 (turb180t) of the physical
                  nu + nut (were 7.3e-5 / 7.8e-15) — first-order cannot
                  smear these fronts. Revisit together with the flat-plate
                  increment (inlet-vs-wall classification + scalar inlet
                  values + outflow validation).

T4 short gates: correlations unit-tested by `src/test_transition.f90` (26
tabulated values vs an independent transcription, every piecewise branch);
`transition = false` bit-exact (the suites); transition-on 1-rank == 4-rank
EXACT, CPU vs GPU <= 2.8e-14 (exp/pow intrinsics, the T3 log() class);
gamma/rethetat restart round-trip verified (restart with a CHANGED tu keeps
the read values) and legacy restarts (no gamma/rethetat datasets) warn +
reinitialize. Discretization note: the Re_thetat~ diffusivity sigma_thetat
(nu + nut) is TWICE the momentum diffusivity the Peclet dt controller
budgets for, so its diffusion diagonal is point-implicit (rans.f90 kernel
comment) — fully explicit it checkerboards to 1e6 within ~40 steps on lam30t.

`rans_channel_check.py` holds the pass criteria (tolerances overridable
per invocation in check_gates.sh).
