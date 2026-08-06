# RANS T2+T3+T4 gates — k-omega SST transport (resolved walls, wall functions, gamma-Re_thetat transition)

Physics gates for IDDES phases T2, T3 and T4 (docs/next_session_iddes.md).
These are the LONG runs, packaged to run on a remote machine (the T3
wall-function sweep is tiny and also runs fine on a laptop; ibm180wf ran
on the local GPU):

```bash
# on the big machine (after rsync + ./compile.sh cpu at the repo root)
cd validation/rans_sst
RANKS=8 ./run_gates.sh            # all cases, sequentially; or one name:
RANKS=16 ./run_gates.sh ibm180
./run_gates.sh t3                 # just the T3 wall-function set

# back home (or anywhere with python3 + h5py + numpy)
./check_gates.sh
```

Cases (each writes `<prefix>_<laststep>.h5` + `<prefix>.log`). RESULTS
2026-07-08 (runs on the remote machine, checks local) — all PASS:

- `laminar.ini`   gate (a): Re_tau 10 channel, SST on, tu = 1%. Parabola
                  to 2.4e-4 of u_max; k decayed 4e-3 -> 3e-47. (At
                  Re_tau 30 / tu 5% the no-transition SST finds its
                  textbook weakly-turbulent branch — k self-sustains at
                  nut/nu ~ 0.5 and the parabola flattens 12% — so the
                  gate probes the genuinely subcritical regime.)
- `turb180.ini`   gate (b): developed channel Re_tau 180 (natural y,
                  y+_1 ~ 0.5). Log line to 4.9%; U+ centreline 18.16 vs
                  DNS 18.20 (0.2%); u_tau = 1.0008.
- `turb395.ini`   gate (b): Re_tau 395. Log line to 6.5% (the pure
                  kappa/B line deviates from real profiles by several %
                  over the overlap — hence the DNS centreline anchor):
                  U+ centreline 20.16 vs DNS 20.13 (0.15%);
                  u_tau = 1.0016.
- `ibm180.ini`    gate (c), THE key IBM gate: the les_ibm off-grid
                  file-IBM plane channel (walls mid-cell, uniform y,
                  y+_1 ~ 2-3) through the IBM wall treatment
                  (dwall_blocks + wall-cell omega pinning). Log line to
                  4.3% with the walls at y = 0.259375/2.259375, u_tau
                  0.966. Needs ../rans_geometry/ibm_coeff_blocks_l1.h5
                  (../rans_geometry/setup.sh regenerates it).
- `refine180.ini` gate (d): uniform-y Re_tau 180 with both wall bands
                  refined to level 1 (2:1 interfaces at y = 0.25/1.75).
                  NO interface band: jump/local-variation ratios 0.58
                  (k), 1.11 (omega), 0.90 (nut), 0.73 (u); core profiles
                  match the RESOLVED turb180 reference to 0.5% (u/k),
                  2.8% (omega), 1.7% (nut).
- `base180u.ini`  informational coarse-wall control (uniform y+_1 ~ 2.8,
                  no refinement): its under-resolved sublayer feeds a
                  spurious core-k plateau (~1.9 vs the resolved ~0.9) —
                  the demonstration of WHY the wall bands get refined.
                  Not a pass/fail gate.

T3 wall-function cases (`[rans] wall_treatment = wall_function`; all
checked with `--mode wallfn` against the RESOLVED turb180 field, which
carries the DNS centreline anchor to 0.2%, so the anchor is transitive).
RESULTS 2026-07-08 (runs local: sweep on 4-rank CPU, ibm180wf on the
local GPU) — all PASS:

- `wf180_y30/45.ini` T3 gate (a): turb180 coarsened to uniform ny = 6/4
                  (y+_1 = 30/45; ny = 6 is not nb-divisible, so it runs
                  with [blocks] nb unset — rank-box blocks). Implied U+
                  centreline vs DNS 18.20: 1.2% / 0.7%; u_tau from the
                  delivered wall stress (nu + nut_1) U_1/y_1 = 1.0000.
- `wf180_y05/y15/y22.ini` T3 gate (b): graceful degradation across the
                  buffer range (y+_1 = 5/15/22.5). Implied centrelines
                  -3.1% / +2.8% / +3.0% — a mild overshoot, NO
                  double-counting dip (that would be a deficit). First
                  cells below y+ 30 sit 10-19% above the resolved
                  profile: the log-line error at the anchor cell,
                  informational. base180u (resolved, y+_1 2.8) is the
                  marginal control.
- `ibm180wf.ini`  T3 gate (c): ibm180 through the wall-function blend
                  (200k steps, local GPU, 18.6 ms/step). The k-based y+
                  keeps every wall cell on the viscous branch, whose
                  arithmetic is exactly the resolved treatment: the
                  steady profile matches the T2 resolved ibm180 field to
                  ROUND-OFF (u 4.5e-16, nut 1.1e-16) — zero regression.

Gates (d)/(e) (resolved-mode bit-exactness vs T2 8991192 + wall-function
determinism, nofma) are NOT here — they are short and run locally (cap
nsteps ~20 per the bit-exact-gates-short rule). T3 results: resolved
bit-exact (max_abs 0, CPU AND GPU) on min_channel/les_ibm ± refine_body/
Beltrami y-slab/turb180; wall-function 1-rank == 4-rank exactly, CPU vs
GPU <= 2e-13 (the `log()` intrinsic in the wall-function branch differs
by an ulp between host and device libm — resolved mode remains exactly
CPU == GPU).

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
`t_final` with the fixed binary reproduces the gate above to EVERY PRINTED
DIGIT (log-region dev 0.0297, near-centre 0.0297, implied U+ centreline
18.41, k/om/nut ranges identical) — the RANS fixed point does not remember
the IC. Cold-started RANS field snapshots taken before 2026-08-05 are
therefore not reproducible bit-for-bit with the current binary; their
physics is.

GATE: `wf180_y30.ini` at 20 steps, 1 rank vs `dims = 4 1 1`, `1 1 4` and
`2 1 2`, all **max_abs 0** on `un vn wn pn k omega nut` (nofma CPU). The
S5a determinism gate's `[mpi] dims = 1 1 4` pin is no longer load-bearing
(kept: it still passes, and a z split remains the cheaper decomposition
for that grid). Independently re-confirmed on an IBM + RANS + scalar cold
start (`../scalar/ibmwf180.ini`, 8 blocks, 20 steps): 1 rank == 4 ranks
**max_abs 0** on `un vn wn pn k omega nut theta`.

RE-BASELINE REQUIRED. `validation/scalar/run_bitexact.sh` against the
pre-fix (S5a) binaries now reads exactly the intended signature — the four
non-RANS cases untouched, the three cold-started RANS cases moved:

```
min_channel  les_ibm  les_ibm_refine  beltrami_yslab   -> max_abs 0, PASS
turb180      un 2.28e-02  k 2.90e-01  omega 5.26e+00   -> differs (the fixed IC)
wf180_y30    un 1.55e-02  k 1.58e-01  omega 2.21e+00   -> differs
lam30t       un 1.77e-02  gamma 2.80e-01  rethetat 8.0 -> differs
```

Those three are 20-step transients from a corrected initial condition, and
their CONVERGED answers are unchanged (above). So `~/s5a_ref_binaries/` is
stale for RANS cases: the next increment must build its reference set from
a commit that CONTAINS this fix, or it will chase this difference.

T4 gamma-Re_thetat transition cases (`[rans] transition = true`,
Langtry & Menter 2009 = OpenFOAM kOmegaSSTLM, resolved walls only; run
group `t4`). The canonical flat plate is DEFERRED: an inlet can be
composed from the existing faces (Dirichlet velocity + Neumann pressure,
zero-gradient outlet), but the RANS layer is not inlet-aware yet (a
Dirichlet inlet classifies as a no-slip wall in domain_face_is_wall /
dwall, and the scalars lack inlet ghost values) and no inflow/outflow
case is validated — so the gates are channels.
RESULTS 2026-07-09 (local 4-rank CPU) — all PASS:

- `laminart.ini`  gate (a) subcritical control: Re_tau 10 / tu 1% with
                  transition on stays laminar exactly like T2 (parabola
                  2.4e-4, k -> 8.6e-16, gamma pinned at its 0.02 floor).
- `lam30.ini`     transition-OFF control at Re_tau 30 / tu 5%: the plain
                  SST weakly-turbulent branch (parabola off by 12.2%, k
                  self-sustained at 0.37). Informational — documents WHY
                  lam30t is discriminating.
- `lam30t.ini`    gate (a), THE discriminating gate: same conditions with
                  transition ON laminarize (parabola 1.6e-3, wall-layer
                  gamma 0.024, mean-k peak 9.1e-3 = 40x below the branch;
                  state stationary t=150 -> t=300; the k residual is the
                  gamma-floor 2% of P_k, hence --k-max 0.02).
- `turb180t.ini`  gate (b): turb180 with transition on preserves developed
                  turbulence: gamma >= 0.999 for y+ >= 30, U+ centreline
                  18.44 vs DNS 18.20 (1.3%, tol 2%), u_tau 1.0009; the
                  kappa/B log-line fit is 6.6% (LM's sublayer D_k x 0.1
                  coupling lifts the low-log rows ~2% over T2 — gated at
                  the turb395 tolerance 0.08; the DNS anchor is the hard
                  criterion).
- `t4_front_check.py` STEP-0 evidence for keeping first-order upwind on
                  the transported scalars: the gamma front is wall-normal
                  and the cross-front upwind diffusivity max|v| dy/2 is
                  7.3e-5 (lam30t) / 7.8e-15 (turb180t) of the physical
                  nu + nut — first-order cannot smear these fronts.
                  Revisit together with the flat-plate increment
                  (inlet-vs-wall classification + scalar inlet values +
                  outflow validation).

T4 short gates (local): correlations unit-tested by
`src/test_transition.f90` (26 tabulated values vs an independent
transcription, every piecewise branch); transition = false bit-exact vs
T3 25ef6ed (nofma, max_abs 0 incl. k/omega/nut, CPU AND GPU) on
min_channel / les_ibm ± refine_body / Beltrami y-slab / turb180 /
wf180_y30; transition-on 1-rank == 4-rank EXACT, CPU vs GPU <= 2.8e-14
(exp/pow intrinsics, the T3 log() class); gamma/rethetat restart
round-trip verified (restart with a CHANGED tu keeps the read values) and
legacy restarts (no gamma/rethetat datasets) warn + reinitialize.
Discretization note: the Re_thetat~ diffusivity sigma_thetat (nu + nut)
is TWICE the momentum diffusivity the Peclet dt controller budgets for,
so its diffusion diagonal is point-implicit (rans.f90 kernel comment) —
fully explicit it checkerboards to 1e6 within ~40 steps on lam30t.

`rans_channel_check.py` holds the pass criteria (tolerances overridable
per invocation in check_gates.sh).
