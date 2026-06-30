# 2:1 interface-NORMAL velocity treatment (no LES) — RESOLVED 2026-06-30

**OUTCOME: the ~9% coarse-owns v' asymmetry is ACCEPTED as the intrinsic price of
the energy-conserving const-1/2 interface. No code change. The production
interface is unchanged and stays validated.** The investigation below pinned the
mechanism conclusively and ruled out every lever; do not re-litigate without new
infrastructure (a 2-layer normal halo, and even that only partially helps).

## What was measured (no LES, `vface_asym.py`, committed)

Developed core-patch run, time-avg v'_rms over the box's two mirror-symmetric
y-faces (physics demands they be EQUAL):

| y-face | v_full | v_coarse (injectable) | v_fine (sub-coarse-cell) |
|---|---|---|---|
| LOWER = coarse-owns | 0.850 | 0.817 | **0.235** |
| UPPER = fine-owns   | 0.778 | 0.773 | **0.086** |
| interior (centreline) | 0.566 | 0.560 | 0.084 |

lower/upper v'_rms mean ratio ≈ 1.09–1.12. The excess is a **~2.7× spurious
sub-coarse-cell fine-structure spike at the coarse-owns face** (v_fine 0.235 vs
the ~0.086 interior / fine-owns baseline) — NOT a coarse-scale (mean-transport)
effect.

## Mechanism (conclusive)

At the **coarse-owns** face the FINE block is the high side, so it PREDICTS its own
interface-normal face `q(1)` — but its advective flux `vv_m = (q(0)+q(1))²` reads
the deep halo `q(0)`, which under const-1/2 is the **prolong-INJECTED coarse value**
(one coarse v replicated over the 4 fine cells). That pumps Jensen-type
(avg-of-squares ≫ square-of-avg) variance into the fine prediction → the v' spike.
At the **fine-owns** face the fine's interface face is its high HALO, not predicted
(only injected + projection-corrected), so it shows just the physical baseline.

## Why every lever is blocked (three converging probes)

1. **Deep-halo lever (the source) is violently unstable.** A diagnostic that
   gently replaced the injected normal deep halo with `q(0)=q(1)` (removing the
   coarse-injected variance, no cubic) blew the developed patch up exponentially —
   `max|v|` 4.3 → 4.2e3 → 5.7e7 … ~10⁴/snapshot, diverging within ~120 steps. The
   const-1/2 injection of the interface-normal deep halo is on a **stability
   knife-edge** and is LOAD-BEARING — the same failure mode the cubic
   reconstruction had. The asymmetry's source cannot be touched.
2. **The face reconciliation can't reach the source.** It changes the stored face
   value / divergence, not what the predictor reads, so a symmetric shared-face
   reconciliation would not remove a predictor-sourced sub-cell variance excess.
3. **The face reconciliation is independently storage-blocked.** `blk%q` carries a
   SINGLE halo layer `(0:nb+1)`. Making the fine-owns face fine-authoritative needs
   the fine block to PREDICT its high-halo interface-normal face `q(ny+1)`, whose
   flux/Laplacian read `q(ny+2)` — which does not exist. A 2-layer normal halo
   would be required, and even then the fluid ACROSS the fine-owns face is the
   COARSE region (no fine data exists there), so the normal velocity stays
   coarse-influenced. Partial payoff at best for a large infrastructure change.

Net: const-1/2's stability and its ~9% interface-normal v' asymmetry are two faces
of the same conservative low-order treatment. The deep halo is the source but is
untouchable; the face lever can't reach it and is storage-blocked. Accepted.

Tools (permanent record): `validation/channel_interface/core_patch/vface_asym.py`
(per-face v'_rms + peakedness) and the within-coarse-cell decomposition snippet in
the 2026-06-30 session transcript. Production settings unchanged:
`interface_constant_half = true`, `momentum_reflux = false`.

---

# ORIGINAL HANDOUT (superseded by the resolution above)

Branch `claude/jacobi-interface`. The 2:1 interface is validated in turbulence for
both the flat-face and edge/corner cases, and LES rides on it cleanly
(`validation/channel_interface/les/`, `../core_patch/`). The **one remaining
interface blemish** is an orientation asymmetry in the interface-NORMAL velocity
treatment. This session addresses it **without LES** (pure 2:1 mechanics).

## The problem (well localized)

The Phase-3c convention is **"the low-side block owns the 2:1 shared face."** For a
y-face this splits into two orientations:

- **fine-owns** (fine block is the low side): the fine block predicts and corrects
  its own shared face with its *resolved* normal velocity; the coarse neighbour
  gets the conservative RESTRICT (4-sub-face average). **Clean.**
- **coarse-owns** (coarse block is the low side): the coarse block owns the face;
  the fine cells across the interface receive their interface-normal velocity by
  **PROLONG = injection of the under-resolved coarse face value**. The fine grid
  then develops resolved cross-stream motion on top of a normal velocity pinned to
  the coarse value → a small **v'-only spike** in the fine cells at that face.

This is the residual that the momentum-reflux removal (`momentum_reflux=false`,
the production default) does NOT remove — reflux was the *large* band; this is the
ownership/prolong-injection residual underneath it.

### Evidence (measured, LES patch run but it is a NORMAL-treatment defect, not LES)

Embedded core patch (fine box in the channel core). The box's **lower** y-face is
coarse-owns, its **upper** y-face is fine-owns. Time-avg fluctuation rms, fine
cell at the interface, `max/mean` over the face plane (baseline interior = 1.45):

| component | lower face (coarse-owns) | upper face (fine-owns) |
|---|---|---|
| u' (tangential) | 1.36 | 1.46 |
| **v' (normal)**  | **1.69** | 1.45 |
| w' (tangential) | 1.45 | 1.40 |

Only v' (the component **normal** to the y-face) spikes, only at the coarse-owns
face. Reproduce: `validation/channel_interface/les/fig_interface_rms.py` (the
`max/mean` per interface plane is a few lines on top of its `accumulate`). The
same rule produced the earlier flat-wall-band asymmetry (there the core-side top
interface was coarse-owns). It is purely a 2:1-mechanics defect — drop LES for
this work (it only made the spike easy to see in a developed field).

## Goal

Make **both** orientations behave like the clean fine-owns one: the interface-
normal velocity on the fine side should be resolved/consistent, not an injected
coarse value, while staying **conservative** (single-valued shared flux) and
SPD-friendly for the Chebyshev-Jacobi projection. The most promising route (noted
across the validation): a **two-sided symmetric conservative reconciliation** of
the shared normal-velocity face — both blocks predict the face, then reconcile to
a single conservative value — instead of coarse-owns(inject)/fine-owns(restrict).

## Where the code is

- `src/modules/comm.f90`: the exchange entries (op COPY/RESTRICT/PROLONG,
  `src_samples`, the single weighted gather `entry_gather_map`; `interface_normal_dim`
  / the "skip prolong on the owned normal-velocity face" logic). This is where the
  ownership + transfer of the normal component lives.
- `src/main.f90`: the post-predictor face sync (`c%syncFace`, ~the line after the
  momentum predictor) and the momentum predictor's face pinning (only PHYS/CLOSED
  pinned; FACE_FINE/COARSE predicted on both sides).
- `src/modules/step.f90`: the sweep denominator/correction noflux masks
  (`noflux_low = {PHYS,CLOSED}`, `noflux_high = any non-open`) — the fine-owns-face
  split that 3c encodes; changing ownership touches these.
- CLAUDE.md "Phase 3c" + the `momentum-interface-todo`, `interface-validation-suite`,
  `interface-projection-adjoint-fix`, `jacobi-interface-restart` memories have the
  full lineage (incl. why unconditional fine-owns is storage-blocked: a fine-west
  face's DOFs sit in unpredictable halos; restriction would write the coarse
  INTERIOR plane against the prolong reading it).

## Gates (no LES, build `-Mnofma` / `-gpu=nofma`)

MANDATORY before/after (from `beltrami-interface-regression`): RK single-step dump
+ pressure-correction-alone (exact-overwrite) + whole-Beltrami t=8; baselines/tools
in that memory. Plus:
1. **Bit-exact no-interface** (single-level) vs Phase 2, CPU 8 ranks + GPU.
2. **`MOBY_HALO_AUDIT`** clean (the normal-velocity halo is the subtle one).
3. **Uniform-flow exactness** through a refined patch (max dev 0.0) — any new
   transfer MUST preserve constants.
4. **Global mass / divergence** to round-off (MOBY_STEPDIV); stable 1000 steps on
   the 3D patch and x-band channels.
5. **The asymmetry metric**: core-patch (reflux off), `v'` `max/mean` at the
   coarse-owns vs fine-owns face — target: coarse-owns drops from ~1.69 to the
   fine-owns/interior baseline (~1.45). (LES optional here; a no-LES developed
   patch works too — the spike is in the velocity field.)
6. **CPU == GPU** bit-identical.

## Risks

- The fix touches the most delicate part of the scheme (staggered normal velocity
  at a 2:1 face). Earlier attempts at "unconditional fine-owns" hit the storage
  block; a cubic/anti-diffusive normal reconstruction DESTABILIZED (see
  `corner-reconstruction-todo`). Favour a **conservative, symmetric, low-order**
  reconciliation over a higher-order one.
- Keep the projection SPD at the interface (`face_grad` composite stencil +
  conservative copy reconciliation) — that is what keeps Chebyshev-Jacobi stable.
- Validate `momentum_reflux=false` (production) AND on; the change must not
  reintroduce the reflux band.

---

## NEXT-SESSION PROMPT

> Read `docs/next_session_interface_normal.md` and the `momentum-interface-todo`,
> `interface-validation-suite`, `jacobi-interface-restart`,
> `beltrami-interface-regression` memories. Branch `claude/jacobi-interface`. The
> 2:1 interface + LES are validated; the one remaining blemish is an orientation
> asymmetry in the interface-NORMAL velocity treatment (Phase-3c low-block-owns-
> face): at COARSE-OWNS y-faces the fine cells get their normal velocity by
> prolong-injection of the under-resolved coarse value, giving a small v'-only
> spike (max/mean ~1.69 vs ~1.45 baseline); FINE-OWNS faces are clean. Work
> **without LES** (pure 2:1 mechanics). Goal: make both orientations behave like
> the clean fine-owns one via a two-sided symmetric CONSERVATIVE reconciliation of
> the shared normal-velocity face (not coarse-owns-inject / fine-owns-restrict),
> keeping the projection SPD and constants/mass exact. Run the mandatory Beltrami
> regression first, then the gates in the handout (bit-exact no-interface, halo
> audit, uniform-flow exactness, mass/divergence to round-off, the v' asymmetry
> metric coarse-owns→fine-owns baseline, CPU==GPU). Make a plan first; we execute
> after.
