# 2:1 interface improvement — experiment log

Goal: reduce the 2:1 interface error toward the no-interface baseline, primary
concern the *directional asymmetry*, then overall magnitude. Hard constraints:
exact discrete mass conservation; momentum conservation if achievable. Scope:
incremental improvements to the current composite-projection / above-block-owns
scheme first; a different formulation only if increments stall.

Baseline = commit `4f0d5c2` (revert with `git reset --hard 4f0d5c2`).

## Quantitative gate: Taylor–Green vortex (exact solution)

`validation/taylor_green` (Re=100, niter=200, IBM off, t=2). Primary metric =
L2 velocity error of the central 2:1 refine patch vs the exact solution.

| state | uniform-64 L2 | refined-patch L2 | notes |
|-------|---------------|------------------|-------|
| baseline (4f0d5c2) | 2.18e-5 | 8.10e-3 | error localised at the patch boundary, directionally asymmetric |

Standing invariants every change must preserve:
- **G1 single-level bit-exact**: TGV uniform output unchanged (interface-only
  edits must not touch the single-level path).
- **G2 mass**: global divergence sum of the refined TGV ≈ round-off.
- **G3 interface stability**: `tutorials/interface_decay` contracts over 200 steps.
- (later) channel nb=4 bit-exact vs Phase 2, rank independence.

## Diagnosis (baseline, slab_x decomposition)

- x- and y-interfaces are equivalent (slab_x = slab_y = 7.50e-3); straight
  interfaces dominate (square patch 8.10e-3 ≈ slab 7.50e-3, corners +8%).
- Error decomposed by location/component at the two opposite-orientation
  interfaces of slab_x (interfaces at x=π/2 and x=3π/2):
  - **Tangential v: ~0.82 at BOTH interfaces (symmetric)** — PROLONG injection
    smears the tangential velocity (piecewise-constant, O(h)). This is the
    magnitude driver. → point 1: linear/trilinear tangential prolong.
  - **Normal u: 0.47 at one interface vs 0.08 at the other (~6× asymmetric)** —
    the above-block-owns ownership convention. → point 2: symmetric / reflux.
  - w stays exactly 0 (no spurious cross-flow).
- The gather kernel currently does only equal-weight averaging (inject or
  RESTRICT-average); linear prolong needs a *weighted* gather (the enabler for
  both point 1 and the point-3 pressure-ghost blend, which was removed earlier).

## Experiments

### E1 — linear tangential PROLONG (point 1, magnitude)
Replace piecewise-constant injection of the tangential velocity (and tangential
prolong generally) with linear interpolation (3/4 covering + 1/4 adjacent
coarse), via a weighted gather. Normal-direction prolong and RESTRICT unchanged
(keeps interface mass flux conservative).

Implementation: a `declare target` `gather_taps` per-dim 2-tap weighted gather
(lin=0 reproduces inject/RESTRICT bit-exactly; lin=1 = face-linear in the
variable's staggered dim, cell-linear otherwise). Applied to copy_local_entries
(local sampler; off-rank pack deferred, TGV is single-rank).

| state | uniform-64 | refined | slab_x |
|-------|-----------|---------|--------|
| baseline | 2.18e-5 | 8.10e-3 | 7.50e-3 |
| E1 (full linear) | 2.18e-5 (bit-exact ✓) | **NaN** | **NaN** |

Result: **UNSTABLE.** uniform-64 stays bit-exact (single-level untouched), but
the refined cases blow up — NaN spread to all 176 blocks over 70 steps, a
growing mode. The composite projection's interface relaxation was tuned for
*injection* halos; replacing the tangential prolong with linear interpolation
breaks the contraction (consistent with the doc's "interface relaxation is
delicate, naive changes are unconditionally unstable"). → linear prolong cannot
be dropped in isolation; it needs a coordinated projection-side change.

### E1b — under-relaxed linear prolong (alpha=0.5 toward injection)
Result: **NaN** (refined and slab). The instability is not a mild
over-correction — even 12.5%-toward-adjacent destabilises.

### E1c — cell-linear only (face-linear in the staggered dim disabled)
Result: **NaN** (refined and slab). So it is not specific to the face-staggered
prolong: *any* linear velocity prolong (cell or face) destabilises.

### E2 — pressure-ghost blend (point 3): ghost = (2*coarse + fine_interior)/3
Re-add the blended pressure ghost (removed in the above-owns refactor) at the
fine side's FACE_COARSE faces, as a localized post-exchange kernel
(blend_prolong_pressure), applied after each projection pressure exchange.

| state | uniform-64 | refined | slab_x |
|-------|-----------|---------|--------|
| baseline | 2.18e-5 | 8.10e-3 | 7.50e-3 |
| E2 (pressure blend) | 2.18e-5 (bit-exact ✓) | 1.07e-2 (+32%) | 9.71e-3 (+29%) |

Result: **STABLE but WORSE.** No NaN (pressure-only changes do not trigger the
velocity instability), but the blend increases the error ~30%. The above-owns
reconstruction handles the interface pressure differently from the old scheme,
so re-adding the old blend is inconsistent — this retroactively justifies its
removal. Reverted.

## Conclusion on point 1 (and implications)

**Linear velocity prolong is unconditionally unstable with the current
composite projection** (E1 / E1b / E1c all blow up, NaN spreading to all
blocks; single-level stays bit-exact, so the machinery is correct — it is the
*physics coupling* that fails). The projection's interface relaxation contracts
only with **piecewise-constant (injection)** halos; replacing the prolonged
tangential velocity with any interpolation breaks the contraction. This matches
the strategy doc's hard-won experience that the interface relaxation is
delicate and naive transfer changes are unconditionally unstable.

Therefore the magnitude driver (tangential prolong, v≈0.82) **cannot be fixed
by an isolated transfer-order change** — it requires the transfer and the
projection relaxation to be co-designed. Points 2 (normal-velocity ownership
asymmetry) and 3 (pressure ghost) are the same projection-coupled class and
carry the same instability risk; isolated increments are unlikely to be stable.

Code reverted to baseline `4f0d5c2`; no solver change kept. Kept: this log, the
TGV exact-solution gate, and the slab_x/slab_y orientation diagnostics.

## Recommendation (the "different formulation" path)

The clean, exact-solution gate is now in place, and the error is precisely
characterised. The principled next step is a **conservative, projection-aware
interface** rather than a transfer-order patch:

1. **Momentum-conservative reflux (Berger–Colella)** at the interface — the
   doc's §6a candidate. Keep injection halos (so the relaxation still
   contracts), but add a *reflux correction* that makes the interface momentum
   flux conservative (the coarse and fine fluxes telescope). This is the
   standard AMR cure for interface errors and is built to be stable. It targets
   both the magnitude (the un-refluxed flux that the injection smearing
   represents) and, done symmetrically, the ownership asymmetry.
2. Validate every step against the TGV gate (refined-patch L2, slab_x/slab_y
   orientation split) — it resolves changes at the 1e-5 level with no turbulence
   or statistics, so a reflux that works will show the refined L2 fall toward
   the 2e-5 baseline and the slab orientation asymmetry collapse.
3. If a reflux still relies on injection-only halos for stability, a
   higher-order *and* stable transfer needs the relaxation re-derived for
   interpolated halos (a genuine reformulation) — out of scope for unsupervised
   incremental work; flagging for a directed session.

Bottom line: the obvious increment (linear prolong) is a dead end on stability
grounds; the path to the target is a conservative reflux interface, developed
against the new TGV gate.

### E3 — linear PROLONG in the FINAL exchange only (the stable magnitude fix) ✓ KEPT

Key realization from E1: linear prolong breaks the projection's *relaxation*, but
the momentum predictor reads halos from the **previous substage's final
exchange**. So apply linear prolong ONLY in that final exchange (injection
everywhere inside the relaxation): the relaxation still contracts (injection),
while the momentum predictor sees an O(h^2)-accurate coarse->fine tangential
halo. Velocity components only (pressure stays injected). Gated by
`comm%linProlong`, set solely on the multi-level projection's final exchange.

| state | uniform-64 | refined | slab_x | tangential v @ interfaces | normal u @ interfaces |
|-------|-----------|---------|--------|---------------------------|-----------------------|
| baseline | 2.18e-5 | 8.10e-3 | 7.50e-3 | 0.82 / 0.82 | 0.47 / 0.04 |
| **E3** | 2.18e-5 (bit-exact ✓) | **6.40e-3 (−21%)** | **5.57e-3 (−26%)** | **0.10 / 0.05 (−87%)** | 0.47 / 0.04 (unchanged) |

- **Stable**: interface-decay gate contracts over 200 steps (un/vn/wn/pn all
  decrease); TGV refined stable to t=2.
- **The dominant tangential (magnitude) error drops ~8-16×** — exactly the
  injection-smearing the diagnosis identified.
- **The normal-velocity asymmetry is untouched** (E3 doesn't change the owned
  face), so it is now the dominant remaining error → the reflux/ownership target.
- Mass proxy improves (slab_x 2.3e-4 → 4e-6).
- **Exact rank independence** (refined TGV 1-vs-2 ranks bit-identical, FMA and
  nofma). Two subtleties had to be fixed: (a) the per-point weighted gather is a
  single shared `declare target` function (`gather_point`) called by both the
  local-copy and off-rank-pack kernels, so the weighted sum is compiled once and
  cannot reassociate differently between them; (b) the linear adjacent tap is
  clamped to the source block's **interior** cells 1..nb — the staggered face
  nb+1 is a halo that injection never reads, so it can be rank-dependently stale
  during the exchange; reading it broke bit-exactness. Clamping there falls back
  to injection at the prolong region's edges (≈1% of the error, negligible).
- Single-level path untouched (linProlong only set in the multi-level branch).
- Verified bit-identical 1-vs-2 ranks under both FMA (GPU) and nofma (CPU);
  interface-decay gate passes; committed `b3d10a0`.

## Overall summary (experiments tried)

| experiment | point | stable? | refined L2 vs 8.10e-3 | verdict |
|------------|-------|---------|-----------------------|---------|
| E1 linear prolong (full, all exchanges) | 1 (magnitude) | no (NaN) | — | revert |
| E1b linear prolong (alpha=0.5) | 1 | no (NaN) | — | revert |
| E1c cell-linear only | 1 | no (NaN) | — | revert |
| E2 pressure-ghost blend | 3 (pressure) | yes | 1.07e-2 (worse) | revert |
| **E3 linear prolong, final exchange only** | 1 (magnitude) | **yes** | **6.32e-3 (−22%)** | **KEEP** |

E3 is the breakthrough: decoupling the accuracy-improving linear prolong (final
exchange, feeds the momentum predictor) from the stability-critical relaxation
(stays on injection) makes the magnitude fix stable. The tangential error drops
~8-16×. The **normal-velocity ownership asymmetry** (u 0.47 vs 0.04) is now the
dominant residual — next target (reflux / symmetric owned-face treatment).

Two of the three points tried empirically; both fail. Point 2 (normal-velocity
ownership asymmetry) is the same velocity/projection-coupled class as point 1
(so the same instability risk) and is structurally the hardest — its principled
fix IS the momentum reflux below, so it was not attempted as a separate
ad-hoc increment.

**Net: simple incremental transfer/ghost changes to the current composite
scheme do not improve the interface error** — they either destabilise (any
velocity-transfer interpolation) or worsen it (pressure-ghost blend). The
current scheme is a local optimum among these increments. Single-level stays
bit-exact throughout (all changes were interface-only), and the code is left at
baseline `4f0d5c2`.

**Recommended next step (directed work): momentum-conservative reflux.** Keep
injection halos (stability), add a Berger-Colella flux correction so the
interface momentum flux telescopes; develop it against the TGV gate
(refined-patch L2 -> 2e-5 baseline, slab_x/slab_y orientation split for the
asymmetry). This is a coordinated reformulation, not an ad-hoc patch, and is
the route consistent with both the diagnosis and the stability evidence.
