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
| **E3 + two-phase** | 2.18e-5 (bit-exact ✓) | **6.32e-3 (−22%)** | **5.53e-3 (−26%)** | **0.10 / 0.05 (−87%)** | 0.47 / 0.04 (unchanged) |

- **Stable**: interface-decay gate contracts over 200 steps (un/vn/wn/pn all
  decrease); TGV refined stable to t=2.
- **The dominant tangential (magnitude) error drops ~8-16×** — exactly the
  injection-smearing the diagnosis identified.
- **The normal-velocity asymmetry is untouched** (E3 doesn't change the owned
  face), so it is now the dominant remaining error → the reflux/ownership target.
- Mass proxy improves (slab_x 2.3e-4 → 4e-6).
- **Exact rank independence** (refined TGV 1-vs-2 ranks bit-identical, FMA and
  nofma). The per-point weighted gather is a single shared `declare target`
  function (`gather_point`) called by both the local-copy and off-rank-pack
  kernels, so the weighted sum is compiled once and cannot reassociate
  differently between them.
- **Two-phase final exchange (exact second order, no clamp).** The linear
  stencil's adjacent tap `base±1` can land in the source block's halo (e.g. the
  staggered face nb+1, or the coarse cell across a block boundary). The final
  exchange runs in two passes so that halo is current: Pass A is a plain
  injection exchange (refreshes every coarse halo from the post-sweep interiors
  via same-level copies + restricts), Pass B does the linear PROLONG on the
  velocity, whose 2-point stencil reads those now-current halos. This is exact
  O(h²) to the patch edges — generic across orientations AND "wedged" coarse
  blocks (a coarse block bordering fine patches on two perpendicular sides: its
  halo-from-restrict is set in Pass A before the Pass-B prolong reads it; no
  special-casing). Verified: refined 6.32e-3 (the earlier edge clamp's ~1% loss
  recovered), 1-vs-2 ranks bit-identical, interface-decay stable, and an
  L-shaped two-patch wedged case bit-identical across ranks. Cost: one extra
  exchange pass on the final exchange (~13% s/step on a small refined channel,
  the second pass's fixed launch/MPI overhead). The earlier interim solution
  clamped the adjacent tap to the interior 1..nb (injection at edges); the
  two-phase supersedes it.
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

## Deferred work (after E3 — for a later session)

E3 fixed the tangential magnitude. The dominant residual is the **normal-velocity
ownership asymmetry** (fine-owns interfaces u≈0.47 vs coarse-owns 0.04). Two
complementary directions, neither yet built:

1. ~~Bespoke normal-velocity interface reconstruction~~ **TRIED — FAILS (mass).**
   Implemented a normal-direction reconstruction of the fine owner's coarse-side
   normal-velocity halo (average the two bracketing coarse faces, 1/2,1/2,
   instead of injecting the single covering face; faces only; Pass B only). It
   made the error WORSE (refined 6.32e-3 -> 1.16e-2) and blew up the mass
   residual ~10x (net-mass 2.5e-4 -> 2.5e-3). Root cause: the normal velocity is
   the interface FLUX. The projection makes the field divergence-free *with the
   injection halos*; overwriting u(0) in the final exchange destroys the
   divergence-free property -> mass error. So the normal halo is mass-locked and
   CANNOT be naively reconstructed (unlike the tangential halo, which is free).
   The fine-owns accuracy is in genuine tension with mass conservation: injection
   conserves mass but makes the interface diffusion/convection of the owned face
   read a piecewise-constant coarse halo. A fix must improve momentum accuracy
   WITHOUT changing the (mass-locked) normal halo — i.e. a momentum correction
   (reflux), not a halo change. Reverted.
2. **Momentum-conservative reflux** (targets conservation, the hard constraint).
   Flux register at the interface, restrict the fine convective+viscous fluxes to
   the coarse face, correct the coarse cell by the mismatch. Guarantees momentum
   conservation but corrects the coarse side, so it likely does NOT directly
   reduce the fine-side 0.47; pursue for conservation, measure its asymmetry effect.
3. ~~Two-phase exchange to drop the linear-prolong edge clamp~~ **DONE** (see
   the two-phase bullet under E3 above). Implemented as Pass A injection (which
   already refreshes the post-sweep coarse halos via copies AND restricts) +
   Pass B linear prolong. Putting restricts in Pass A (the injection pass), not
   alongside the prolong, is what makes the wedged case work with no phase
   hierarchy: the only halo-consumer is PROLONG, everything else reads interiors,
   so the dependency graph is acyclic and two phases suffice.

### Diagnosis of the fine-owns error (Reynolds sweep, slab_x)

| Re | fine-owns u | coarse-owns u | ratio |
|----|------------|---------------|-------|
| 30 | 2.69 | 1.27 | 2.1 |
| 100 | 0.47 | 0.037 | 12.5 |
| 300 | 0.61 | 0.012 | 50.6 |

- Coarse-owns error ~ 1/Re (vanishes at high Re): pure **diffusion**, and it is fine.
- Fine-owns error drops then RISES (0.47 -> 0.61, Re 100->300) while diffusion
  vanishes: a Re-INDEPENDENT component (~0.4-0.6) dominates at the relevant Re.
  Subtracting the diffusion part (≈ coarse-owns), the residual is ~constant ->
  **convection**. It is fine-owns-specific (coarse-owns uses accurate RESTRICT).
- Conclusion: the fine-owns asymmetry is the **convective momentum flux** reading
  the mass-locked injected normal halo u(0). It is a MOMENTUM error
  (Re-independent, convective), not a mass error -> correctable conservatively
  WITHOUT changing the stored (mass-critical) halo. The fix is a momentum-side
  correction: either a convective-flux reflux, or a reconstructed normal halo
  used ONLY in the momentum RHS (convection/diffusion) while the divergence keeps
  the injection halo. Both keep mass exact; both are real implementations.

Both/all develop against the TGV gate (refined L2 + slab_x/slab_y orientation split).

### Convective-flux reflux (fine-owns interface) — KEPT

The fine owner's interface convective flux uu_m=(u(0)+u(1))^2 reads the stored
coarse-side normal velocity u(0), which is injected one fine cell deeper than the
fine halo position. Its correct value at the halo position is 0.5*(u(0)+u(1)),
computable on the fly from data the fine cell already has (no extra exchange, no
dual halo). Use it ONLY in the interface convective flux (momentum RHS); the
stored halo, read by the divergence/projection, is untouched -> mass stays exact.
Applied at i/j/k=1 when physLow(d)==FACE_COARSE (the fine-owns - edges), for the
normal velocity (uu_m/vv_m/ww_m).

Results (GPU): refined 6.32e-3 -> 5.85e-3 (-7.5%; cumulative -28% from the 8.10e-3
baseline), slab_x 5.53e-3 -> 5.22e-3; mass residual 2.5e-4 -> 1.3e-5 (**18x
better** -- the corrected flux makes the interface genuinely more conservative);
uniform-64 bit-exact; rank-independent (1v2 = 0); interface-decay stable.

Honest assessment: this captures the **conservation** part of the fine-owns error
(hence the large mass gain) but only modestly moves the **asymmetry** peak
(fine-owns u 0.465 -> 0.453; ratio 12.5 -> 8.9).

### Convergence study (resolution sweep) — the residual is a SCHEME DEFECT, not inherent

Ran uniform + slab_x at nx = 32 / 64 / 128 (with the reflux):

| nx | uniform L2 | slab_x L2 |
|----|-----------|-----------|
| 32 | 8.72e-5 | 5.36e-3 |
| 64 | 2.18e-5 | 5.22e-3 |
| 128| 5.61e-6 | 5.38e-3 |

uniform converges at **order 2.0**; the slab_x interface error is **flat
(order ~0)** -- resolution-INDEPENDENT. So the residual interface error does NOT
converge: it is a **consistency defect (an O(1)-in-a-finite-region scheme error),
NOT inherent coarse-grid truncation** (this corrects the earlier "inherent"
guess). A consistent interface treatment would converge with the grid; order 0
means a term in the fine-owns interface treatment is simply wrong at leading
order and survives refinement. This is fixable in principle -> a reformulation
(deferred option 2: change the ownership/coupling so the fine owner does not
import coarse-grid data through the owned-face reconstruction) is worth pursuing,
not a fundamental limit. The convective reflux removed the conservation part; the
consistency part remains for a future directed effort.

### Truncation-error probe — DEFECT FOUND: tangential pressure injection

`MOBY_TRUNC` (step.f90:truncation_probe + main.f90 hook): with the exact TGV
field + scheme halos, prints the u-momentum operator terms RMS'd over y-z at each
x of a fine-owns and a coarse-owns interface; the inviscid balance conv+pres
should -> 0 in a consistent scheme. At the owned face (i=1):

| | \|conv\| | \|pres\| | \|conv+pres\| |
|---|---------|----------|---------------|
| fine-owns | ~0 | 7.24e-2 | 7.24e-2 |
| coarse-owns | 4.4e-3 | 1.9e-3 | 6.27e-3 |

The owned-face momentum imbalance is **12x larger on the fine-owns side**
(matching the 0.45/0.05 solution asymmetry) and is **entirely the pressure
term**. Mechanism: the projection reconstructs the owned face as
u(1) = qs - dt*ifGrad*(p_new(1) - p_new(0)); the fine owner's p_new(0) is the
**tangentially injected** (y-staircase) coarse pressure, which does not match the
fine p_new(1)'s y-variation -> a spurious O(h) gradient. The coarse owner reads
RESTRICT'd (accurate) pressure -> clean. CONFIRMED: giving the pressure the
tangential linear prolong dropped the fine-owns imbalance 7.24e-2 -> 1.67e-3
(43x), cleaner than coarse-owns.

CAVEAT (why the simple fix did not move the solution): doing this only in the
final exchange (the momentum predictor) is inert -- the predictor's start-of-
projection pressure gradient CANCELS in the reconstruction (the pStart terms
telescope, leaving u(1) = q + dt*rhs - dt*ifGrad*(p_new(1)-p_new(0))). The fix
must put the tangentially-smoothed coarse pressure into p_new(0) **inside the
projection** (the per-colour pressure exchange that feeds the sweep
reconstruction), not the predictor. Concerns to handle there: (1) stability --
it is in the relaxation (E1 showed velocity-linear there is unconditionally
unstable; pressure may differ but must be checked); (2) rank independence -- the
per-colour exchange is single-phase, so the linear adjacent tap needs the
two-phase (producers-then-prolong) treatment. This is the precise, minimal target
for the asymmetry fix.

## Original recommendation (superseded by E3 for the magnitude)

**Momentum-conservative reflux.** Keep
injection halos (stability), add a Berger-Colella flux correction so the
interface momentum flux telescopes; develop it against the TGV gate
(refined-patch L2 -> 2e-5 baseline, slab_x/slab_y orientation split for the
asymmetry). This is a coordinated reformulation, not an ad-hoc patch, and is
the route consistent with both the diagnosis and the stability evidence.
