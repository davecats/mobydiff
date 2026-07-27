# BoostConv steady-state accelerator — implementation plan (PROPOSED 2026-07-27)

Goal: accelerate convergence of steady RANS cases (the C11 airfoil polar
class: the OF-ambient aoa-5 run needed >10 t.u. = 200k+ steps to converge
its circulation) by residual recombination, WITHOUT touching the
validated stepping. References: Citro & Palitta, "Residual Recombination
Methods as Anderson-like Acceleration" (tutorials/naca/boostconv.pdf —
the ROBUST Algorithm 3.1 with QR downdating and a linear-dependency
threshold tau) and the CPL reference implementation
(cplcode.net/Applications/Numerical/BoostConv; plain variant: Gram
matrix + PLU, circulating slots; sources in ~/fri/).

## The algorithm (robust variant, Alg. 3.1 of the paper)

The solver's own advance is the fixed-point map: one ACTIVE iteration =
p timesteps, x_{k+1} = x_k + r_k with r_k = Phi_p(x_k) - x_k. BoostConv
replaces r_k by xi_k = r_k + W_k c_k where c_k solves the small least
squares min ||r_k - V_k c||_2 over the last N stored residual-difference
directions:
  V_k = [..., r_{k-1} - r_k]            (residual differences)
  W_k = [..., xi_{k-1} + r_k - r_{k-1}] (matching state responses)
kept as a skinny QR of V_k (Q_k R_k): one Gram-Schmidt pass per active
iteration (O(nN)), Givens downdating when the window slides, and the new
column DISCARDED if ||(I - QQ^T) v_new|| < tau ||v_new|| (rank
protection). c_k = R_k^{-1} Q_k^T r_k. Update: x <- x_k + xi_k, keep
(r_k, xi_k) for the next column (paper Remark 2.1). Interpretation:
matrix-free multisecant/Anderson — kills the slowly-decaying modes
(exactly our circulation/ambient-decay tail) while leaving the solver
untouched.

## Module layout (per user spec)

1. `src/modules/boostconv.f90` — GENERIC, physics-blind:
   - `boostconv_type`: capacity N, interval p, tau, nDof; device arrays
     W(nDof,N), Q(nDof,N); host R(N,N) + bookkeeping; prev-state and
     prev-residual vectors x_prev(nDof), r_prev(nDof), xi_prev(nDof).
   - `boostconv_init(bc, nDof, N, p, tau)` + device maps
     (`enter_/exit_boostconv_data`, the established pattern).
   - `boostconv_apply(bc, x, comm)`: takes the PACKED state vector
     (device-resident), forms r, runs Alg. 3.1, overwrites x with
     x_prev + xi. All O(nDof) work as flat OpenMP-target loops (dot,
     axpy, pack of the new Q column); the N x N algebra (R update,
     Givens, triangular solve) on the host. MPI dot products via
     comm_allreduce_sum.
2. `boostconv_rans(...)` in rans.f90 — the RANS-facing wrapper:
   - pack u,v,w,p (blk%q interior) + k,omega (sst%k/omg interior) into
     the work vector with PER-VARIABLE DIAGONAL SCALING (see below),
     call boostconv_apply, unpack, then state hygiene: velocity halo
     exchange (syncface), scalar halo exchanges, apply_bc,
     rans constrained-cell pinning (solid k, wall omega, kpin boxes) —
     all existing calls, so the boosted state re-enters the loop
     exactly like a restart-read state does.
   - called from the moby_solve main loop at end-of-step when
     mod(step, p) == 0 (the ACTIVE iterations; W/V built only there,
     paper section 2.2).

## Design decisions embedded (flagged for review)

- SCALING: raw concatenation would be omega-dominated (wall omega
  ~1e5-6 vs u ~1). Per-variable diagonal weights 1/s_var with
  s_var = max(rms(field at activation), floor). Recomputed at each
  activation from the CURRENT field (cheap: the dot products are being
  formed anyway); frozen alternatives (fixed s) kept as a config
  fallback for strict reproducibility studies.
- DIVERGENCE: the boosted velocity is not discretely divergence-free.
  Default plan: rely on the incremental projection of the following
  steps (the same way a restart or interpolated IC is absorbed);
  evaluate in V1 whether one extra projection call after boost is
  needed (config `boostconv_project = true` if so).
- WHAT COUNTS AS ONE RESIDUAL: r over p SOLVER STEPS (p = interval).
  Paper guidance: apply every p > 0 iterations, p larger when the local
  Jacobian varies slowly. For RANS airfoil tails (slow circulation
  mode) p ~ 20-100 steps is the expected sweet spot; V1/V2 sweep
  p in {10, 25, 50, 100} x N in {6, 10, 20}.
- STEADY-ONLY: config-gated under [rans] (boostconv = true,
  boostconv_interval, boostconv_capacity, boostconv_tau, default
  1e-8); requires [turbulence] model = rans (hard error otherwise).
  Genuinely unsteady cases are the user's responsibility; a cheap
  guard prints a WARNING if ||r|| grows over 3 consecutive
  activations (suspend boost, keep stepping).
- MPI: plain allreduce dot products. CONSEQUENCE: the boosted
  trajectory is not bit-identical across rank counts (sum order); the
  boost-OFF path stays in the bit-exact regime, and the boost-ON gate
  is steady-state EQUIVALENCE (converged fields agree to ~1e-10)
  rather than trajectory identity. If trajectory-exactness is wanted
  later, the A2 block-ordered exact-reduction pattern can be dropped
  in (once per activation, cost negligible).
- MEMORY: nDof (C11) = 6 x 8.21M = 49.3M doubles = 394 MB per column;
  W + Q = 2N columns: N = 10 -> 7.9 GB device (fine on A6000 48 GB /
  5090 32 GB on top of the ~4 GB case; NOT on the 3060 12 GB).
  Config `boostconv_host_basis = true` keeps W,Q in host memory and
  streams per activation (~8 GB over PCIe ~ 0.3-0.5 s, negligible at
  p >= 10 steps between activations) — the 3060/loaded-GPU option.

## What is NOT touched

The fused predictor kernel, the projection, the RANS substage kernels,
the exchange machinery: zero changes. BoostConv is a pure end-of-step
state overwrite at active iterations — dormant-off bit-exactness is by
construction (no call, no state), verified by the suite anyway.

## Gates

V0 dormant: 7-case nofma suite bit-exact CPU+GPU with boostconv absent/off.
V1 correctness + first speedup (cheap, known fixed point): turb180 —
   boost-on from the standard IC must converge to the SAME steady state
   as the standing converged field (max|diff| <= ~1e-9 class) in
   measurably fewer steps (target >= 3x on the residual-vs-step curve);
   p/N sweep; check the divergence question (with/without extra
   projection); 1 vs 4 ranks: same converged state to ~1e-10.
V2 the payoff case: C11 aoa5 OF-ambient (restart from t = 20):
   time-to-converged C_L (band +-0.002 around the asymptote) vs the
   unaccelerated history (~10 t.u.); target >= 3x fewer steps. Then the
   polar campaign restarts WITH boost (user go).
V3 GPU == CPU sanity on a short boosted run (ulp-class tolerances);
   wall-time overhead per activation measured (expect < 1 % at p >= 25).

## Cost estimate

boostconv.f90 ~ 300-400 lines + ~80 in rans.f90/moby_solve/config;
V0-V1 gates ~ half a day of runs (turb180 is minutes per try on the
A6000); V2 one cetus run. No solver-kernel risk.

## Open questions for the user

1. Basis residency default: device (fast, 8 GB at N = 10) or host
   (fits everywhere, +0.5 s per activation)?
2. Scaling: dynamic per-activation rms (adaptive, default) vs fixed
   scales (strictly reproducible)?
3. V2 acceptance threshold: is >= 3x steps-to-converged-C_L the right
   bar to justify running the polar with boost?

## STATUS — V1 NEGATIVE (2026-07-27): BoostConv loses to plain marching here

Five configurations on turb180 (baseline: plain marching reaches the
fixed point GEOMETRICALLY, dist 1e-3/1e-5/1e-10 at 16k/44k/116k steps;
per-25-step contraction ~0.9935):

| config | result |
|---|---|
| V1a p=25 N=10, 6-var pack, drift-flush metric | plateau 8.8e-5; metric thrash (35 flushes) |
| V1b p=25 N=10, NO p in pack, frozen metric | plateau 1.5e-4, 156 growth flushes |
| V1c + tau 1e-3 + ||xi|| cap | IDENTICAL to V1b (guards never trigger) |
| V1d 6-var pack (pn = projection MEMORY, must be packed), pinned dt | best: converges, but 2.7x SLOWER than plain |
| V1e paper-faithful p=1, N=20 | worst: 0.25-0.29x, 8244 guard events |

The ALGORITHM is verified correct: a rule-for-rule Python mirror on a
stiff linear fixed point (rho = 0.999) converges 1e-10 in ~100 its
where plain iteration stalls — textbook BoostConv (scratchpad
bc_unit.py, reproduced in this doc's history).

DIAGNOSIS: mismatch between BoostConv's premises and an EXPLICIT
time-marching map. Phi_p (p RK3 steps at CFL ~0.8) carries a large
cloud of marginally-damped, strongly non-normal fast modes; every
recombination kick excites the cloud, the secant pairs sample healing
transients rather than the smooth slow-mode response, and the
correction quality collapses (guards fire massively at p = 1).
Published BoostConv successes ride heavily-damped implicit /
pseudo-compressibility inner iterations. Fixes that helped but could
not flip the sign: pack pn (the incremental projection makes it state
MEMORY), freeze the scaling metric per basis, pin dt (map autonomy),
tau/cap safeguards. Also documented: k-floor at unpack corrupts secant
pairs at near-wall cells (bounded, secondary).

DECISION (user, pending): (A) keep the module dormant-gated as-is
(V0 bit-exact 14/14; zero solver risk) and run polars unaccelerated;
(B) research direction — wrap BoostConv around a damped INNER
iteration (pseudo-time smoother) instead of the raw explicit map;
(C) optionally still try V2 on the C11 airfoil (its circulation mode
is ~5x slower than turb180's: rho_25 ~ 0.9988 — the slow-mode gap is
larger, but the explicit-cloud objection stands); one 16 h cetus run.

## UPDATE (2026-07-28): Dicholkar package in; alpha sweep verdict

Implemented after reading boostconv2.pdf (Dicholkar, Zahle & Sorensen,
JWEIA 2022 — BoostConv + SIMPLE-like RANS + k-omega SST airfoils to
Re 1e7; their headline finding = our V1: the ORIGINAL alpha = 1 method
is not robust on turbulent RANS):
- [rans] boostconv_alpha (default 0.02): xi = r + alpha W c (their
  relaxation; convergence valley at 0.01-0.05).
- [rans] boostconv_start: delayed activation (their 3000-iteration
  warm-up; we used 4000 on turb180).
- SPECIFIC-RESIDUAL formulation (user suggestion): the basis stores
  r/(elapsed pseudo-time), the update rescales by the current elapsed
  time — identical under pinned dt, well-posed under adaptive dt.

turb180 alpha sweep (p = 25, N = 20, tau 1e-3, start 4000; plain
baseline 1e-3/1e-5/1e-10 at 16k/44k/116k):
  alpha 0.005: 0.96-0.98x, final 4.8e-12
  alpha 0.02 : 0.88-0.91x, final 4.7e-11  (first machine-class boosted run)
  alpha 0.05 : 0.76-0.80x, final 7.3e-10
Monotone: alpha -> 0 recovers plain; no alpha accelerates turb180.

CONCLUSION: the Dicholkar modification delivers exactly what the paper
claims — ROBUSTNESS (converging where alpha = 1 fails), NOT speedup of
already-converging cases; their gains are convergence where the plain
solver LIMIT-CYCLES (post-stall polars, shedding cylinders). For us
the method's value proposition is therefore: (i) stabilizing
steady-RANS points that limit-cycle (high-alpha polar angles), and
(ii) possibly the C11 circulation tail (rho_25 ~ 0.9988, a genuinely
isolated slow mode) — V2 remains the one untested question. turb180
was the wrong yardstick for acceleration; the 3x bar there is
unreachable by construction.

## V2 (2026-07-28): the C11 slow-tail test — usable, cheap, converged

V2 = the finalized aoa-5 case restarted from the t = 26 near-steady
state with boost (alpha 0.02, N 20, p 25, specific residuals), t 26-29
on the healed cetus A6000: overhead +2.4 % s/step (0.457 vs 0.446; the
earlier 0.93 was the dying machine, not boost); C_L plateaus at
0.516 +- 0.005 within ~1 t.u. of activation (plain-trend estimate:
~2-3 t.u. to the band, though the plain asymptote extrapolation was
itself low — no concurrent plain twin, so claim conservatively:
boost converged the tail at no measurable cost and no downside).
FINAL VALIDATION NUMBERS (vs OpenFOAM): C_L 0.516/0.5142 (+0.4 %),
C_D 0.0134/0.0134 (exact), Cp_min -1.787/-1.780 (+0.4 %).
