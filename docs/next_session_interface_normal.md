# Next session — 2:1 interface-NORMAL velocity treatment (no LES)

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
