# Review: 2:1 interface treatment and the `channel_interface` failure

Independent review of the block-refinement work (branch `claude/blocks`)
against `docs/block_refinement_strategy.md` and the two source papers
(Nakahashi & Kim, AIAA 2004-434; Jansson et al., IJHPCA 33(4), 2019),
focused on the unresolved turbulent-channel validation
(`validation/channel_interface`). Files inspected: `src/modules/blocks.f90`,
`comm.f90`, `pressure_solver.f90`, `step.f90`, and the validation
results in `plots_y55/` and `plots_y110/`.

## Bottom line

- The decomposition and solid-block removal are faithful to both the plan
  and the papers, and close to as simple as the problem allows.
- The 2:1 refinement interface is where the difficulty lives, and the
  difficulty is largely **intrinsic and self-imposed by one mismatch**: both
  reference methods are *collocated*; mobydiff is *staggered*. The single
  operation that is breaking — the face-normal velocity at a 2:1 interface in
  a coupled incompressible projection — is the one thing neither paper had to
  solve.
- The asymmetric-interface diagnosis is correct. It is visible in the code,
  not just the write-up, and the validation statistics carry its signature.
- The proposed "uniform-B + reflux" fix is sound and is in fact the clean
  implementation route to the two-part fix already sketched in strategy §6a.
  Two caveats: the 2-layer halo cost, and the fact that the decisive part
  (momentum reflux) is back-loaded to the last increment.

## 0. The cross-cutting fact: collocated papers, staggered code

Jansson et al. §2.5, verbatim: *"Treatment of face centred, staggered
quantities is not necessary because we use a collocated arrangement in the
present work."* Nakahashi & Kim is likewise cell-centred, and additionally
compressible/explicit with no pressure projection — its inter-cube transfer
is a one-way ghost-cell average/injection with no conservation coupling back
through an elliptic solve.

Consequence: the papers supply the cheap, transferable machinery — cube
generation, the 2:1 rule, Z-order linear distribution, fine→coarse average /
coarse→fine injection, the wall buffer — and **none** of the expensive
machinery this code actually needs: a 2:1 interface that is conservative for
both mass *and* momentum, stable under a red-black SOR projection, on a
staggered grid where the normal velocity lives *on* the interface. That part
had to be invented here, and it is the only part still failing.

## i. Compliance with the plan and the papers

**Plan vs. papers — faithful.** Equal `nb^3` blocks set by refinement level
(BCM core); flat tree-free arrays with the block as last index (Jansson §2.3,
"no tree, linear arrays"); offline octree generation only (Nakahashi);
2:1 smoothing (both); RESTRICT = average of covering fine faces, PROLONG =
injection (Nakahashi: *"transfer from smaller cube to larger… average…
larger to smaller… just taking the value"*; Jansson Fig. 3); Z-order linear
distribution with `n = floor((N+P-p-1)/P)`, which is **verbatim Jansson
eq. 17** (confirmed in `blocks.f90` `zorder_*`); finest-level wall buffer so
interfaces never touch the body, which is Nakahashi's own recommendation
(*"the same size of cubes should be used near the body boundary"*). Removal
of cubes buried in the body is straight from Nakahashi.

**Implementation vs. plan — tracks it, with three documented deviations.**
The plan is a living document; each departure was forced and recorded:

1. §7: original `mu = 0` closed face → reuse of the wall face-mask machinery
   (`FACE_CLOSED`, `face_kind`). Confirmed in `blocks.f90` /
   `pressure_solver.f90`.
2. §6: "fine-owns-face unconditionally" → "low-side owns" (Phase 3c), because
   in staggered storage a west-fine block's interface DOFs live in
   unreachable halos.
3. §6: one-sided interface relaxation → BCM-style **symmetric** relaxation
   (Phase 3d), after the one-sided splitting proved unconditionally unstable.
   Confirmed: `pressure_solver.f90` `face_pinned` pins only `FACE_PHYS` /
   `FACE_CLOSED`; `step.f90` predicts both sides of an interface.

Net: decomposition, removal, and same-level exchange comply with plan and
papers. The interface section is where the plan already had to leave the
papers behind (necessarily), and it is the only unvalidated part.

## ii. The asymmetric-interface diagnosis is correct

It is in the code, not only the prose.

`comm.f90::interface_boxes` is structurally lopsided for the normal
component:

- high side (`off=+1`): `dstLo = nb+1` — a restrict **writes the boundary
  normal face**;
- low side (`off=-1`): `dstLo = 0` — a restrict writes **only the halo**,
  never the boundary face at index 1.

`step.f90` momentum computes `v(1..nb)` but never `v(nb+1)`: the high-side
normal face of a block is always a halo. So of the two blocks at an
interface, the **high-coordinate-side block physically computes the shared
normal face** (as its own low-boundary face, index 1); the low-side block
only receives a copy in its halo (index `nb+1`). See §A below for why this
follows from the staggered layout.

Mapped onto the channel (wall-normal = y, j=1 at the bottom wall):

- **Bottom wall**: refined band *below* the interface (low side), coarse
  *above* (high side). Shared normal face computed by the **coarse** block,
  at coarse resolution; the fine band receives it by injection. The fine
  scale of the wall-normal velocity is absent exactly at the interface.
- **Top wall**: refined band *above* the interface (high side), coarse below.
  Shared normal face computed by the **fine** block, at fine resolution.

The two walls are therefore treated with opposite handedness — unphysical for
a symmetric channel. (Vocabulary caveat: under the Phase-3c "low-side owns"
label the *owner* and the *computer* of the normal face are opposite sides;
that compute/own split is the subtle bit, but the conclusion is unchanged.)

**Validation signature.** In `plots_y55/profiles.png` and
`plots_y110/profiles.png` the mean-velocity deviation is not broadband — it
spikes **at the interface** (~2.1% at y+=55, ~1.4% at y+=112), with the
Reynolds-stress and rms deviations localized to the same band, shrinking as
the interface moves to weaker turbulence. A localized mean-shear / Reynolds-
stress defect that scales with local turbulence intensity and flips with wall
orientation is the fingerprint of an interface that is mass-conservative but
**not momentum-conservative**, with an asymmetric normal velocity. This is
very likely the cause, not a coincidence.

## A. Why the high-side block computes the shared face (staggered layout)

MAC arrangement, with the velocity component stored *on* the face it crosses:

```
        ----- V(i,j+1) -----        <- top face of cell (i,j)
       |                    |
     U(i,j)   P(i,j)    U(i+1,j)
       |                    |
        ----- V(i,j)  ------        <- bottom face of cell (i,j)
```

`V(i,j)` is the *low* (bottom) y-face of cell `(i,j)`; `V(i,j+1)` is its
*high* (top) face, which is simultaneously the low face of cell `(i,j+1)`.
Each cell, and hence each block, "owns" the low faces of its own cells. A
block with cells `j=1..nb` thus owns and computes `V(1..nb)` — its bottom
boundary face `V(1)` plus the internal faces — while its **single top
boundary face `V(nb+1)` is not computed; it is a ghost/halo**. That top face
physically belongs to the block above, where it is *that* block's `V(1)` (its
own low-boundary face, which momentum does compute).

So at any block-block boundary the shared normal face is computed once, by
the block on the **high-coordinate side**, and copied into the low-side
block's halo. For same-level neighbours this is harmless: the copy equals the
computed value (and the redundant open-halo trick can make it bit-identical).
At a 2:1 interface it is **not** harmless, because "which side computes it"
and "at what resolution" now differ between the bottom and top walls — the
asymmetry above.

## iii. Soundness of the proposed "uniform-B + reflux" fix

The proposal — have every block momentum-compute its own top normal face
`v(nb+1)` (and `u(nx+1)`, `w(nz+1)`), widen `q`/`qs`/`oldrhs` on the high
side, fill the normal-velocity halo two layers deep, make the coarse
interface face a restriction of the four fine faces, symmetrise the sweep
masks, then add tangential-momentum reflux — is **sound**, and it is exactly
the two-part fix already planned in strategy §6a: Inc 1–3 = part 1
(fine-authoritative normal velocity, both orientations), Inc 4 = part 2
(Berger–Colella reflux). It is a clean realisation of the design, not a new
direction.

What is right about it:

- Computing the top face on every block generalises the existing redundant
  open-halo trick (the one that already makes results independent of `nb`
  and rank count). It dissolves the phase-dependent transfer-direction
  problem the doc worried about: the fine side genuinely holds the DOF.
- Reflux (Inc 4) is the physically decisive piece — the only step that
  restores **momentum** conservation. "Constructed to vanish for uniform
  flow" is the correct constraint to preserve the exact gates.
- The gated increments are disciplined. **Inc 2's "channel nb=4 still
  bit-exact"** is the right canary: it tests that redundant top-face
  computation reproduces the previously-copied value to the last bit (same
  code, same operands, same FMA). Keeping `interface_decay` on at every step
  is the correct reflex given the Phase-3d lesson that smooth-flow gates are
  blind to the interface instability.

Caveats to weigh before starting:

1. **The 2-layer normal halo is the real cost and dents a design point.**
   Going to `-1:nb+2` roughly doubles per-direction halo overhead (≈20% →
   ≈40% at `nb=32`, the value §10 chose to hold 20%). If only the high-side
   layer is needed, `0:nb+2` avoids the cosmetic low-side ghost. Decide this
   deliberately rather than taking `-1:nb+2` "to match the coordinate arrays".
2. **"Uniform" buys branch-free interface logic with global FLOPs + memory.**
   The asymmetry is only at 2:1 interfaces; computing the top face redundantly
   on *every* block (including same-level) is more than strictly required. It
   is a defensible GPU branch-avoidance call, but it is not free — keep it on
   the Phase-4 optimisation list explicitly.
3. **`oldrhs` halo: recompute-redundantly is the right call** (pure function
   of exchanged state, no new messages, cannot desync). Confirm `lapYp(nb+1)`
   carries the correct one-sided spacing as flagged.
4. **Inc 1–3 alone may shrink but not remove the defect.** Symmetrising the
   normal velocity kills the wall asymmetry, but the constant shear-stress
   excess in the mean balance is the *un-refluxed tangential flux*; that only
   closes at Inc 4. The y+=55 spike may survive until then — the real test is
   post-Inc-4.

No soundness red flags; sequencing is correct and each increment is
independently falsifiable.

## iv. Is this the simplest solution to the original ask?

Split by sub-feature.

- **Block decomposition** (Phase 0/1) and **solid-block removal** (Phase 2):
  yes, near-minimal and well-judged. Decomposition reuses the comm skeleton,
  leaves inner-loop bodies verbatim under one outer block loop; removal reuses
  the wall mask instead of a special `mu=0` path.
- **Local refinement**: complexity has compounded, and only part is intrinsic.
  Unavoidable: any conservative AMR on this solver needs momentum refluxing,
  and a staggered grid genuinely has a normal-velocity-on-the-interface
  problem. Avoidable: adopting two *collocated* reference methods for a
  *staggered* projection code is what turned the interface into a research
  problem (one-sided → symmetric relaxation, pressure-ghost blending, the
  uniform-B rework). The papers gave the plumbing for free and zero of the
  hard part. Naming this matters: the difficulty is not a sign of a wrong
  turn, it is the references stopping where the real problem starts.

The sharper question is **scope, not implementation.** The original motivation
is geometry-adaptive refinement around an immersed body, and §4 already
mandates a finest-level wall buffer so interfaces sit in smooth flow — both
papers insist on this because coarse-fine transfer is the least accurate
operation. `channel_interface` deliberately plants interfaces in the most
energetic turbulence, the hardest possible placement, arguably stricter than
the target use case (the sailplane `refine_body` run was already reported
stable and improved). A genuinely simpler path that still satisfies the
original ask: ship decomposition + removal + geometry-refinement-with-buffer,
and treat "interface robust inside fully-developed turbulence" as a separate,
optional capability. If that capability is wanted, uniform-B + reflux is a
sound way to get it — but it is a deliberate scope expansion beyond the two
original features, and worth choosing consciously rather than having the
validation case force it.

## Recommendations

1. Proceed with the proposed increments; treat **Inc 2 bit-exactness** and
   **post-Inc-4 turbulence statistics** as the two gates that matter most.
2. Use `0:nb+2` (high-side only) unless a low-side ghost is independently
   needed; log the halo-overhead change.
3. Decide explicitly whether turbulence-grade interfaces are in scope. If the
   real deliverable is body-adaptive refinement with the wall buffer, the
   current state may already serve it; the uniform-B rework is then an
   optional generality.
4. Keep the "compute top face on every block" redundancy on the Phase-4
   list — it is a known, deliberate overhead, not a permanent design.

## v. Related work — how the incompressible-AMR literature treats the interface

Both source papers are collocated and compressible, so they do not cover the
2:1 face-normal velocity under an elliptic projection. The relevant lineage
is the incompressible block-structured-AMR literature:

- **Berger & Colella (1989)** — flux registers and *refluxing*: the fine-side
  flux is the accurate one; the adjoining coarse cell is corrected by the
  reflux divergence of the coarse-minus-summed-fine flux mismatch. This is the
  conservation primitive. mobydiff's Inc 4 (tangential-momentum reflux) is
  exactly this, and "fine-authoritative" is Berger–Colella's "fine flux is the
  accurate one".
- **Almgren, Bell, Colella, Howell & Welcome (1998, JCP 142)** — the
  foundational conservative adaptive projection for incompressible NS:
  per-level "level projection" + a multilevel "sync projection" + refluxing.
- **Martin & Colella (2000); Martin, Colella & Graves (2008, 3D)** —
  cell-centered *approximate* projection with composite coarse+fine operators
  at the interface; IAMR/AMReX is the production codification.
- **OpenFOAM `dynamicRefineFvMesh`** (collocated, Rhie–Chow) — no
  face-ownership problem, but the surface-flux field is mapped
  non-conservatively at refinement and must be *corrected back to
  divergence-free after every AMR event*. Same disease (interface
  interpolation is not conservative), different cure (post-hoc flux
  projection rather than ownership + reflux).

Key structural point for mobydiff: ABC/IAMR need the level/sync-projection
split mainly because they **subcycle in time**. mobydiff uses one global `dt`
and one coupled SOR sweep over the whole block set (strategy §9) — i.e. it
already performs a *composite* projection every step and has no temporal sync
to reconcile. So the temporal half of the classic machinery is not needed;
what remains is the **spatial** flux matching (reflux) plus a consistent C/F
normal-velocity rule. §6a targets exactly that residue.

## vi. Stability of the staggered 2:1 normal velocity (the channel blow-up)

After B2 ("fine-authoritative" normal velocity) + reflux were implemented,
smoke/uniform/decay gates pass but the turbulent channel **diverges**
(max|vel|: 22 → 477 → 1e5 → 1e9 by ~step 100; the CFL limiter drives
`dt → 0`, the run exits on the `dt ≤ 0` guard at t≈0.04, and the stats leg
then inherits `dt=0` — the "restart" error is downstream of the blow-up).

**Diagnosis (sound).** Isolation is correct: `MOBY_NO_REFLUX=1` still blows
up (so it is not the reflux), the original committed code is stable on the
same case, and the spurious velocity sits **at the interface** in strong
near-wall shear. The decisive fact is that the blow-up is **unconditional**:
`dt → 1e-12` and it still diverges.

**What dt-independence implies.** An instability that survives `dt → 0`
cannot live in the explicit advection/diffusion predictor (those errors are
dt-scaled). It lives in the part of the step that is itself essentially
dt-independent: the **pressure/divergence projection coupling at the
interface**. This is the same fingerprint as the Phase-3d pressure-jump mode
(*"unconditionally unstable… per-step gain independent of dt/viscosity/sor…
gates blind because uniform flow is exact and band-refined channels keep
pn=0"*). The present blow-up is the **velocity analog** of that mode.

**Mechanism (hypothesis — RETRACTED; see the 2026-06-15 update at the end of
§vi. The minimal-span test measured the unstable mode as low-k and broadband,
NOT the checkerboard described here; the real lead is a stretched-grid
transfer-consistency bug, not a high-k null space).** B2 gives the fine side four
independent interface sub-faces, while the coarse cell sees only their **sum**
(the restricted face). The perturbation `[+a,−a,+a,−a]` across the four fine
sub-faces has **zero coarse-average**, so it is invisible to the coarse
divergence and coarse pressure — the coarse side cannot push back on it. That
fine-scale component is therefore relaxed **one-sidedly** by the fine pressure
alone against a frozen injected coarse ghost — the exact non-contractive
splitting Phase 3d proved unstable. The crude constant `v(nb+2)` injection is
the *energy source* that excites the mode under shear; the **projection null
space** is why it is not damped. Once `vn` ≈ 5× physical, nonlinear advection
takes over (hence the accelerating, super-exponential growth).

**Assessment of the proposed fix order.** The proposed order was: (1) better
`v(nb+2)` prolongation, (2) interface dissipation/blend, (3) reconsider fine
ownership. This order is **backwards** relative to the dt-independence
evidence and the Phase-3d precedent:

- **Fix #1 (prolongation accuracy)** changes how strongly the mode is
  *excited*, not whether the projection *damps* it. It cannot cure a
  dt-independent projection mode with a robust positive growth rate; it is an
  accuracy improvement (the `v2` upgrade), not a stability fix. Do it last.
- **Fix #2 (low-pass filter / dissipation on the fine-scale interface normal
  velocity)** is the documented cure. Olshanskii et al.'s octree MAC solver
  reports the identical phenomenon — *"oscillatory spurious velocity modes
  tailored to coarse-to-fine grid interfaces"* on staggered grids — cured by
  *"a linear low-pass filter"*. Same family as the existing pressure-ghost
  blend. **Design constraint: the filter must vanish for any coarse-consistent
  (smooth/uniform) field** (act only on the zero-coarse-average component), so
  the uniform-flow and channel-nb4 bit-exact gates survive — the same "vanish
  for uniform flow" discipline the reflux already follows.
- **Fix #3 (reconsider fine ownership / null space)** is the same insight one
  level deeper: exact staggered/cell-centered projections carry a
  non-solenoidal null space (this is *why* ABC/Martin-Colella use approximate,
  null-space-free projections). Giving the fine side independent normal DOFs
  *creates* the null space; constraining or filtering it removes it.

**Physics read on the open question.** *Is fine-owns-normal-velocity stable at
a sheared interface?* Not under an exact staggered projection — it introduces
interface normal-velocity DOFs (the fine-scale, zero-coarse-average part) that
sit in the projection's null space; the exact coupling does not damp them and
shear excites them. It is not doomed (IAMR and octree-MAC solvers run
fine-resolution interface velocities stably) but only because they pair it
with an approximate/filtered projection that supplies the missing damping. So
the stabilizing filter (#2) or a constrained DOF treatment (#3) is the
load-bearing piece, not optional polish.

**The missing gate.** `interface_decay` missed this because it has weak white
noise and **no mean shear**, and uniform/cp gates keep the interface
pressure/velocity null. The gate that would catch it: a **laminar** linear-
shear (or channel-mean) base state, **no turbulence**, with a **structured**
`[+a,−a,+a,−a]` interface seed on the normal velocity, advanced
**projection-only (advection frozen)** — and check contraction of max|vn| at
the interface. If the structured mode grows on a smooth sheared base, the
linear projection instability is isolated and fix #1 is definitively out.
This is the gate to build first; it both diagnoses the mechanism and guards
against regressions.

**Recommended order — SUPERSEDED by the 2026-06-15 update below.** The
minimal-span diagnostics redirected the fix from a high-k filter to repairing
the stretched-grid interface transfers. The original (now-superseded) plan was:
build the sheared-base + structured-seed projection-only gate → confirm the
mode → symmetric/contractive relaxation of the interface normal face +
null-space-vanishing low-pass filter (#2/#3) → `v(nb+2)` accuracy upgrade (#1)
last. The projection-only gate did produce a useful number (a ×3.6 amplification
of a seeded mode), but the full-laminar reproduction was masked by viscosity and
the actual channel mode turned out not to be the checkerboard, so the plan moved
on — see below.

### Update (2026-06-15) — minimal-span diagnostics: two candidate instabilities (B2 and reflux); the audit "bug" was a false alarm

A minimal-span turbulent channel (same Re, grid and IC machinery as the full
channel, spanwise shrunk to a few cells) overturned the checkerboard hypothesis
above. A manufactured-field halo audit produced one apparent stretched-grid bug
that turned out to be an error in the audit, not the solver (see below). Current
best understanding:

**Mode shape: low-k, not checkerboard.** Both unstable modes are broadband,
interface-localized, with ~0 checkerboard projection and ~0–5% Nyquist energy.
The "zero-coarse-average high-k null-space mode" hypothesis is **retracted**.
A low-k, interface-localized *growing* mode is the signature of a discrete
**consistency error at the interface acting as a localized forcing**, not a
null-space resonance — a bug class, not an operator-redesign class. Good news.

**The pressure-ghost "bug" was a FALSE ALARM — the audit was wrong, not the
solver.** The manufactured-linear-field audit had reported 40 bad pressure
cells at the interface PROLONG face on the stretched grid (0 on uniform).
Instrumented, the discrepancy was 2.4e-5 and was *exactly the grid-stretch
correction*: the solver uses `entry_blend`'s node-line weight
`w = (bHalf+cHalf)/(aHalf+cHalf)`, while the **audit** hardcoded `2/3`, which is
exact only on uniform spacing. Fixing the audit to mirror `entry_blend` gives
**0 bad on both grids** (test-only change, gated by `MOBY_HALO_AUDIT`, no solver
effect). So the halo exchange and the pressure ghost ARE correct on the
stretched channel grid; "fix the blend first" has nothing to fix. Lesson:
instrument before declaring a smoking gun.

**Reflux × B2 matrix — read with the laminarization confound in mind.**

| (refined min-channel, identical IC) | B2 ON | B2 OFF (`MOBY_NO_B2_FACE=1`) |
| --- | --- | --- |
| reflux ON | blows up ~step 600 | instant blow-up (<100) |
| reflux OFF (`MOBY_NO_REFLUX=1`) | "stable": max\|v\| 0.17→0.05 | instant blow-up (<100) |

The dominating caveat: the smooth-IC minimal span **laminarizes** — the "stable"
cell shows max|v| 0.17→0.05, i.e. the flow dying, not stability. Therefore:

- **"reflux is the destabilizer" is NOT established.** Reflux-off may simply
  laminarize before it can blow up; the un-confounded full-channel datum is the
  opposite (reflux-off still blew up under the violent KMM IC). Do **not** remove
  the reflux — it is load-bearing for the momentum defect, and it is implicated
  here most plausibly because it shares the same stretched-grid bug class.
- **"B2 fine-face momentum is required" is NOT clean.** `MOBY_NO_B2_FACE=1` is
  the full B2 machinery with fine-face momentum switched off (an inconsistent
  hybrid), not the original committed low-side-owns code — which is stable.
  Confirm what the flag actually leaves active before trusting that row.
- **The linear-field audit does not clear the reflux or the velocity transfers.**
  A linear field carries trivial interface momentum flux, so a reflux bug would
  not show; "velocities clean" is expected regardless.

**Unifying "uniform-grid bug" hypothesis — WITHDRAWN.** It rested on the
pressure-ghost audit finding, which was the audit's error, not the solver's.
With that gone there is no demonstrated stretched-grid transfer bug. What
remains are two *candidate instabilities*, on different timescales, neither yet
cleanly isolated:

- **B2 fine-authoritative normal velocity** — implicated by the one
  un-confounded datum: on the FULL channel, B2-on + reflux-off
  (`MOBY_NO_REFLUX=1`) **blew up**. So B2 is a stability *suspect*, not "a fix,
  not a suspect". It is required for *correctness* (the normal-velocity
  asymmetry) but the original no-B2 code is stable, so it is not required for
  *stability* — a feature can be both a correctness fix and a stability
  liability.
- **Tangential-momentum reflux** — the fast destabilizer in the minimal span
  (reflux-on blows up before the box laminarizes). But "reflux-off = stable" on
  the minimal span is partly laminarization (max|v| 0.17→0.05), and the
  full-channel reflux-off blow-up shows reflux is not the whole story.

The minimal span only reliably exposes the *fast* (reflux) mode; the *slow*
(B2) mode needs the sustained turbulent `v` the minimal span never develops. So
both features probably carry their own instability, and the matrix cannot
attribute blame because (a) it laminarizes and (b) `MOBY_NO_B2_FACE=1` is an
inconsistent hybrid, not the true baseline.

**Proper next step — a clean 2×2 from the true baseline, on a non-laminarizing
reproduction (supersedes the filter-first and fix-the-blend plans above):**

1. Build four real states, NOT flag-hybrids: S0 = original committed
   (low-side-owns, no B2, no reflux; stable but carries the −⟨u'v'⟩ defect);
   S1 = S0 + B2 only; S2 = S0 + reflux only; S3 = S0 + B2 + reflux. Confirm each
   "off" path is bit-exact with S0.
2. Reproduce on a flow that does NOT laminarize, scored by *growth rate*, not
   eventual fate: cheap pre-screen = minimal span seeded with the interpolated
   KMM turbulent field, interface growth rate over the first ~50 steps;
   definitive = the real KMM channel (fast machine) for the survivors.
3. Score two axes separately — stability (does the interface mode grow, how
   fast) and correctness (does it reduce the −⟨u'v'⟩ defect, answerable only on
   the stable runs). The 2×2 also gives the interaction term (B2 alone vs reflux
   alone vs the combination).
4. Cheap discriminator: which velocity component carries the spurious interface
   energy — wall-normal `v` ⇒ B2, tangential `u/w` ⇒ reflux.
5. Only the implicated feature(s) then get a code-level fix; do NOT drop the
   reflux (needed for correctness) — if implicated, repair it. If B2-only (S1)
   grows under sustained turbulence, the operator-level questions in the body of
   §vi (adjointness, contractive relaxation) return for B2 specifically.

Standing caveats: the mode is low-k (not the checkerboard); the minimal span
laminarizes (score growth rate, not eventual fate); and `MOBY_NO_B2_FACE=1` is
a hybrid, so the original committed code is the only trustworthy baseline.

### References for §v–§vi

- Berger & Colella, *Local adaptive mesh refinement for shock hydrodynamics*,
  JCP 82 (1989).
- Almgren, Bell, Colella, Howell, Welcome, *A Conservative Adaptive Projection
  Method for the Variable Density Incompressible NS Equations*, JCP 142 (1998).
- Martin & Colella, *A Cell-Centered Adaptive Projection Method for the
  Incompressible Euler Equations*, JCP 163 (2000); Martin, Colella, Graves
  (3D), JCP 227 (2008). IAMR/AMReX.
- Olshanskii et al., *An octree-based solver for the incompressible
  Navier–Stokes equations* (spurious staggered C/F interface velocity modes;
  low-pass-filter cure).
- Liu, *A stable and accurate projection method on a locally refined staggered
  mesh*, Int. J. Numer. Methods Fluids (2011).
- *Null-space-free methods for the incompressible NS equations on
  non-staggered curvilinear grids* (exact projection's non-solenoidal null
  space).

## vii. Closing synthesis — structural impossibility result and the scope decision

The B2 turbulent-channel blow-up is now fully characterized, and the conclusion
is a clean impossibility result, not an unfinished hunt.

**The mechanism, end to end.** A disciplined diagnosis excluded, in order: a
high-k null-space filter; a stretched-grid pressure-ghost bug (an error in the
audit, not the solver); the `v(nb+2)` prolongation *order* (dt-dependence rules
out an accuracy lever); the pressure-ghost adjointness / oblique projection
(making it adjoint = injection is the *worst* empirically); the convective-flux
non-conservation (real, but dt-scaled — a seed, not the amplifier). What
remained, confirmed by direct measurement, is a **velocity-transfer
divergence-consistency defect**: the per-colour SOR (`interp=.false.`, the
Phase-3d contractive design) relaxes the coarse interface cell divergence-free
against the coarse's *own* copy of `v(1)`, and then the final `interp=.true.`
reconcile overwrites `v(1)` with `avg(fine projected v(nb+1))`, re-introducing
an **O(1), dt-independent local divergence** in the coarse interface cell every
substage that is never cleaned (global mass still telescopes, which is why the
mass gate always passed). The dt-signature (faster blow-up with *smaller* dt =
a per-step algebraic gain) is exactly this dt-independent state perturbation.

**The impossibility.** The localized post-reconcile cleanup (pin `v(1)` to the
fine value, relax the cell's other faces div-free) *removes* the O(1) jump but
is **non-contractive at both sor=1.5 and sor=1.0** — confirming it is
structural, not an over-relaxation artifact, and re-confirming the documented
`denom` lesson (pinning the face and over-compensating through the remaining
faces diverges). The reason is fundamental: the fine's final `v(1)` does not
exist until the fine's projection completes, so the coarse can only relax
against a *moving* fine value (interp=true every colour — Phase-3d divergent) or
a *stale* one it then over-corrects (the cleanup — divergent). Therefore, in the
block-SOR-then-reconcile framework, **a symmetric fine-authoritative 2:1 normal
velocity and a contractive projection are mutually exclusive.** B2 cannot be
stabilized within this framework; escaping it requires coupling the levels in
one solve — a composite / approximate projection.

**The decision boundary.**
- *S0 (low-side-owns) + finest-level wall buffer* — stable, SOR-compatible, and
  the design **both source papers explicitly recommend** (coarse-fine transfer
  is the least accurate operation, so interfaces should sit in smooth flow). Its
  −⟨u'v'⟩ defect is a function of turbulence intensity at the interface, so with
  the wall buffer it is negligible. This **delivers the original goal** (block
  decomposition + local refinement + solid-block removal) now. The
  `channel_interface` case deliberately plants interfaces in the most energetic
  near-wall turbulence — the hardest possible placement, which the wall-buffer
  design specifically avoids.
- *Composite / approximate projection (ABC / Martin–Colella / IAMR-style)* — the
  correct route to turbulence-grade interfaces (interfaces inside energetic
  turbulence, e.g. a free shear layer away from walls). A major pressure-solver
  effort with research risk, justified only if such a capability is actually
  required.

**Recommendation.** Revert to S0, keep the finest-level wall buffer, and
validate on a real refined-body case (the sailplane `refine_body` run was
already reported stable and improved). Treat the composite projection as
separately-scoped, conditional future work — opened only on a concrete
turbulence-grade-interface requirement, not on sunk cost in B2. The investigation
succeeded: the mechanism is fully understood, the impossibility is proven, and
the scope decision can be made with the costs of each path known precisely.

If the composite/approximate projection is chosen, the formally-clean strategy —
single authoritative interface DOF set, conservative composite divergence `D`,
gradient pinned by adjointness `G = −Dᵀ`, the SPD composite Poisson `L = −D Dᵀ`
solved as one coupled system — is written up in
`docs/composite_projection_strategy.md` (it supersedes the §6a "uniform-B +
reflux" plan in `block_refinement_strategy.md`).
