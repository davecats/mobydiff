# Next session — LES ↔ IBM coupling validation

Branch `claude/jacobi-interface`. LES (WALE) is validated across block refinement
and the 2:1 interface **for channel flow with grid-aligned walls** (see
`validation/channel_interface/les/`, CLAUDE.md "edge/corner + LES"). The standing
caveat is that **LES is IBM-UNAWARE in practice**: the `ibm_aware` solid-cell nut
masking was never exercised, and the LES↔IBM coupling is untested both on a single
grid and across the 2:1 interface. This session closes that gap.

## What the coupling actually is (read first)

`src/modules/les.f90`, `update_generic_les_viscosity` / `update_wale_viscosity`
(called from `main.f90` `update_les_viscosity`, which receives `ibm`):

- `les%ibm_aware` (config `[les] ibm_aware`, default **true**). When set AND
  `[ibm] enabled`, a cell is SOLID iff `|ibm%coef| > 1e20` (= `SOLID/Re`, the
  penalization coef of a fully-solid staggered face) at ANY of its 6 staggered
  velocity faces (`les.f90:353-361`); a solid cell gets `nut = 0` (`cycle`).
- Everything else (the IBM **band** cells, where `ibm%coef` is finite second-order
  interpolation weights, NOT `>1e20`) gets the NORMAL nut from the SGS model.
- The filter width is grid-based: `les%filter_{x,y,z}` from `1/d1{x,y,z}` (cube
  root of cell volume, `les.f90:215-243`), per level — it knows nothing about
  where the IBM wall sits inside a cell.

So the coupling is a binary mask, and the OPEN QUESTIONS are exactly where it is
NOT binary:
1. **Band cells.** At an off-grid wall the wall cuts THROUGH a cell. That cell is
   not solid (coef finite) so it gets full nut, but its velocity gradients are
   dominated by the IBM forcing dragging the velocity to ~0 across a fraction of
   the cell → the SGS model may read a huge spurious strain → a spurious nut SPIKE
   right at the wall (the opposite of the physical nut→0). WALE's `sd2` operator is
   designed to give nut→0 in pure shear, but the IBM stair-step is not pure shear.
2. **Filter width at the wall.** `delta` is the full cell size even in a cell the
   wall bisects, so the SGS length scale is too big near the IBM wall.
3. **2:1 interface + IBM together (the real untested combination).** `refine_body`
   refines at the wall, so the IBM band and a 2:1 interface coincide. nut must step
   by the filter-width ratio across the interface AND vanish into the wall, with no
   spurious band where the two meet. Never tested.

## Test case (the one Davide asked for)

A **plane-wall channel whose walls are described by the IBM and do NOT coincide
with any grid node.** Build it with the **file-based IBM** (`tools/mobygeom.py`),
NOT the analytic path: the analytic `isInBody` (`ibm.f90:173`) is a single
hardcoded bottom wavy wall (`amp_x=0.025`, `y_offset=0.01`), not a two-wall flat
channel. Use mobygeom to raster two horizontal planes at y = y_lo and y = Ly−y_lo
chosen so neither lands on a `grid.y` node (e.g. a uniform y-line and walls at a
non-multiple of dy), write the coefficient file, and run with `[ibm] enabled` +
`[ibm] coeff_file` + `[les] model = wale`. Re_tau ≈ 180 to reuse the existing
channel reference. The grid-aligned `tutorials/channel_kmm180` channel (no IBM)
and the `validation/channel_interface/les/` LES cases are the references.

Three runs, smallest first:
- (a) **single level, IBM wall, LES on** — isolates the IBM↔LES masking with no
  interface in play. The decisive first gate.
- (b) **single level, IBM wall, LES off** — the control (does LES change the mean
  flow correctly vs. no-SGS at this resolution).
- (c) **`refine_body` at the IBM wall, LES on** — the 2:1-interface ↔ IBM ↔ LES
  triple, the actual goal. Only attempt after (a) passes.

## Gates (build `-Mnofma` / `-gpu=nofma`)

1. **Solid-cell masking exact:** `nut == 0` in every cell with `|ibm%coef|>1e20`
   (dump nut + coef, check). Hard zero gate.
2. **No spurious wall nut:** the time-mean nut profile must →0 INTO the IBM wall
   and show no spike at the band cells (the open question #1). Compare to the
   grid-aligned LES channel's near-wall nut. If a band-cell spike appears, the fix
   is a band-aware nut damping (mask or scale nut by the cell's solid fraction /
   the IBM coef), or a van-Driest-style wall correction keyed to the IBM wall
   distance — design after measuring.
3. **Mean U log law:** with the IBM wall, the mean profile (shifted by the wall
   offset) matches the grid-aligned LES channel within the IBM wall-modelling
   error (`tools/channel_loglaw.py`, `channel_stats_profile.py`).
4. **nut across the 2:1 interface at the wall (run c):** steps by the physical
   filter-width ratio (δ²), no spurious band where the interface meets the IBM band
   (reuse `validation/channel_interface/les/fig_interface_rms.py` /
   `tools/patch_interface_stats.py`).
5. **Stability:** `MOBY_STEPDIV` bounded, mass round-off, 1000+ steps, no NaN.
6. **bit-exact no-LES / no-IBM:** `[les] model = none` and `[ibm] enabled = false`
   paths byte-identical to before (the change, if any, must be gated).
7. **CPU == GPU** bit-identical (the masking is a branch in the offloaded kernel).

## Risk / likely outcome

The most probable real defect is #1/#2: a **spurious nut spike at the IBM band
cells** because the SGS model sees the IBM-forced velocity drop as resolved strain.
If so the deliverable is a band-aware nut treatment (scale/mask nut by the IBM
solid fraction so it vanishes into the wall like it does on a grid-aligned wall),
validated by gates 2–4. If nut is already clean at the band (WALE's `sd2` may
already kill it), the task reduces to confirming gates 1–7 and documenting that
LES↔IBM is validated. Either way, START by MEASURING the nut field at the IBM wall
(run a, dump nut + coef) before writing any fix.

## Code pointers

- `src/modules/les.f90`: `update_generic_les_viscosity` / `update_wale_viscosity`
  (the `ibm_aware` mask, `solid_threshold=1e20`, the WALE `sd2` operator), filter
  widths `les%filter_{x,y,z}` / `inv_d{x,y,z}`.
- `src/modules/ibm.f90`: `set_ibm_coeff` (penalization coef, `SOLID/Re` in solid,
  finite second-order weights in the band), `isInBody` (analytic geometry — single
  wavy wall, NOT the test geometry).
- `src/main.f90`: `update_les_viscosity(les, blk, dns, ibm)` call site (predictor),
  `exchange_scalar_halos(c, les%nut, blk)`.
- Config: `[les] model = wale|smagorinsky|none`, `cw`, `cs`, `delta_scale`,
  `ibm_aware`; `[ibm] enabled`, `coeff_file`.
- `tools/mobygeom.py`: geometry → coefficient file (+ `block-active` / `block-table`
  for refine_body). `tools/channel_loglaw.py`, `channel_stats_profile.py`,
  `validation/channel_interface/les/` for the stats/figures.
- nut is written to snapshots via `fdm_h5_append_nut` (no-LES output byte-identical).

## Regression note

No interface-machinery change is expected, so the Beltrami interface regression
([[beltrami-interface-regression]]) is not the relevant guard here; the guards are
the no-LES / no-IBM bit-exactness (gate 6) and CPU==GPU (gate 7). If a fix touches
the shared LES kernel, keep the grid-aligned LES channel validation
(`validation/channel_interface/les/`) bit-exact or re-validated.

---

## NEXT-SESSION PROMPT

> Read `docs/next_session_les_ibm.md` and the `interface-normal-treatment`,
> `les-validation-plan` memories. Branch `claude/jacobi-interface`. LES (WALE) is
> validated for grid-aligned channel walls across block refinement + the 2:1
> interface, but the LES↔IBM coupling (`ibm_aware` solid-cell nut masking,
> `les.f90`) is **untested in practice**. Validate it. Test case = a plane-wall
> channel whose walls are described by the IBM and do NOT coincide with grid nodes
> (build via `tools/mobygeom.py`, file-based IBM — the analytic `isInBody` is a
> single hardcoded wavy wall, not a flat channel). Re_tau ≈ 180, reference the
> grid-aligned `channel_kmm180` + the `validation/channel_interface/les/` cases.
> Run (a) single-level IBM wall + LES, (b) the LES-off control, then (c)
> `refine_body` at the IBM wall (the 2:1-interface ↔ IBM ↔ LES triple) only after
> (a) passes. START by MEASURING the time-mean nut field at the IBM wall (dump nut
> + ibm%coef) — the prime suspect is a spurious nut spike at the band cells where
> the SGS model reads the IBM-forced velocity drop as resolved strain; the physical
> behaviour is nut→0 into the wall. Gates (in the handout): solid-cell nut==0 hard
> zero; no spurious wall-nut spike; mean-U log law vs the grid-aligned LES channel;
> nut steps cleanly across the 2:1 interface at the wall with no band; stability
> 1000+ steps (MOBY_STEPDIV); bit-exact no-LES/no-IBM; CPU==GPU. If a band-cell nut
> spike appears, the deliverable is a band-aware nut damping (scale/mask nut by the
> IBM solid fraction); design it after measuring. Make a plan first; execute after.
