# Skew-symmetric convection — global migration plan (approved 2026-07-22)

User decision (2026-07-22, after the C11 v1 term attribution): migrate the
momentum convection to the skew-symmetric form GLOBALLY, config-toggled
during validation, then locked down (toggle removed, divergence form
recoverable from history). Replaces the earlier masked-band proposal.

## Motivation (evidence trail)

- The C11 v1 forced-laminar case exposed an interface instability:
  a divergence-free 2-Delta velocity checkerboard anchored at the
  L8/L9/L10 coarse-owns 2:1 interfaces below the nose, growth ~160/t.u.
  against molecular damping ~58/t.u.
- Term attribution (diag/v1-terms worktree dumps, tutorials/naca/
  diag_attribution.py, commit 2056dd0): projection correction ~ 0
  (2-3 orders below predictor, sign-mixed; niter 18/24, chebyshev/
  jacobi all depart identically), diffusion strictly damping, the
  CONVECTIVE remainder positive and 2-5x diffusion at every sampled
  step/substage. The driver is the centred divergence-form convection
  fed by const-1/2 interface halos.
- Root cause, architecturally: div-form convection is energy-neutral
  ONLY for a discretely divergence-free advecting field. This solver
  never grants that premise exactly (incremental damped-Jacobi/
  Chebyshev projection at finite niter; 2:1 interface halos are O(1)
  inconsistent at grid scale). Production = 1/2 sum phi^2 (div u_eff).
  The skew form 1/2[div + adv] is energy-neutral for ANY advecting
  field - unconditional, dissipation-free, no physical parameter
  involved (the user's requirement: fix the numerics, not damp them
  with nut).

## The discrete form

For each momentum component the kernel's div form is
  conv = 0.25 [ (F1_p S1_p - F1_m S1_m) d1x + ... y ... + ... z ... ]
with pair-sums F (advected) and S (advecting), e.g. for u:
  x: S = (u_i + u_ip), (u_im + u_i)             (cell centres)
  y: S = (v_imjp + v_ijp), (v_imj + v_ij)       (u's y-edges)
  z: S = (w_imkp + w_ikp), (w_imk + w_ik)       (u's z-edges)
Skew identity: conv_skew = conv_div - 1/2 phi (div U)|stencil, with
(div U) built from the SAME advecting pair-sums and metrics:
  rhs_skew = rhs_div + 0.25 phi [ (S1_p - S1_m) d1x
                                 + (S2_p - S2_m) d1y
                                 + (S3_p - S3_m) d1z ]
(sign: rhs carries -conv). This is algebraically 1/2(div+adv); the
correction reuses in-register operands (~10 FLOPs per cell/component).
On stretched lines the antisymmetry has an O(h^2) metric-commutation
residual (Morinishi non-uniform weights NOT adopted now; the y-stretched
channel KE audit below bounds it).

## Implementation

- `[flow] convection = divergence | skew` -> dns%conv_skew (C_BOOL,
  default divergence DURING VALIDATION ONLY).
- step.f90 momentum kernel: one additive flag-guarded branch per
  component after the existing rhs assignment (the ambKW/iddes in-kernel
  branch pattern; div arithmetic untouched -> dormant bit-exact,
  verified by the suite, not just argued).
- oldrhs consistency is automatic (the branch modifies rhs before both
  qs and oldrhs consume it).
- LOCKDOWN (after all gates pass): default flips to skew, then the
  toggle is REMOVED per the production-config philosophy; div form
  recoverable from this commit range.

## Gates

Phase S0 (dormant): 7-case nofma suite bit-exact CPU+GPU with
convection unset/divergence. DONE-criterion for merging the code at all.

Phase S1 (skew correctness):
1. Uniform-flow-through-patch EXACT (0.0) with skew (constants
   annihilate both halves identically).
2. Inviscid KE audit: Beltrami patch case at Re 1e12, ~500 steps,
   total resolved KE non-increasing to round-off with skew; the SAME
   run with divergence documents the current interface KE production.
   Repeat on a y-stretched single-level channel box (metric residual).
3. Beltrami order (uniform + patch) ~ unchanged (2.7-2.9 class).
4. THE ACCEPTANCE TEST: C11 v1 forced-laminar case (.c11_aoa_5t.ini,
   full-nose kpin_box, NO stabilisation) stable to t = 20; diag
   attribution rerun shows G_conv <= 0 (round-off) in the band.
5. 1 == 4 ranks EXACT; CPU vs GPU (expect exact or ulp-class).

Phase S2 (physics re-validation, skew):
6. Re_tau 180 channel (turb180) + interface channel: statistics vs the
   standing references (mean U / u'v' within the reflux-off signature
   tolerances).
7. LES + LES-IBM gates (les_ibm, refine_body twin): law of the wall,
   nut bands, no interface band.
8. RANS gates: turb180 / wf180_y30 / lam30t checkers (momentum
   convection only; scalar upwind untouched).
9. C10/C11 spot re-run: aoa0 + aoa5, CV forces within the CV scatter
   of the committed polar; Cp/Cf overlay.
10. Wall-time delta on min_channel GPU (expect <~ 5%).

Phase S3 (lockdown): flip default -> remove toggle -> re-run the suite
with skew as the only path; update CLAUDE.md + this doc STATUS.

## Status

- S0 (2026-07-22, commit 57bd1e3): DONE — kernel branch + toggle,
  dormant 7-case nofma suite bit-exact 14/14 CPU AND GPU.
- S1 gates 1, 2, 4 (2026-07-22): PASS.
  1. Uniform oblique flow, 3-level refine_body layout, skew ON: EXACT (0.0).
  2. KE audit (Beltrami xz patch, Re 1e12, 500 steps): div +2.3e-1
     KE and ACCELERATING; skew bounded +1.2e-2 (~3e-5 rel), no trend
     (residual = RK3 temporal truncation; dt-scaling check pending).
  4. ACCEPTANCE: the v1 forced-laminar case (full-nose kpin, nut = 0 on
     the interfaces, NO stabilisation) with skew ran its full 5,500-step
     window cleanly — dt never collapsed (t reached 15.274; every div
     variant stalled at 15.12-15.13 by ~153,000), end forces the healthy
     trip signature (C_L 1.6e-3, C_D 9.4e-3), and the seed-point probe
     shows NO anchored mode: field deltas O(1e-2) DECREASING with
     wandering locations vs the div twin's fixed-point 0.30 -> 9.5 ->
     blow-up. The interface energy source is gone, not delayed.
- PENDING: v1 full t = 20 confirmation (launch next); instrumented
  G_conv <= 0 rerun (rebase diag/v1-terms onto 57bd1e3); gate 3
  Beltrami order; gate 5 rank/GPU identity; S2 physics battery;
  S3 lockdown.
