# Next session — LES across the block decomposition and the 2:1 interface

Branch `claude/jacobi-interface`. This is the handout for **open item 2** of
`docs/next_session_edges_les.md` (LES validation). Open item 1 (edge/corner
turbulence) is DONE — see `validation/channel_interface/core_patch/` and the
`edge-corner-validation` memory. Read those first for the machinery this reuses.

## Goal (Davide)

Validate the LES path with block refinement, in order:
1. **LES ⊥ block decomposition** — confirm the LES path is unaffected by how the
   domain is split into blocks/ranks (no refinement yet).
2. **LES ⊥ 2:1 refinement interfaces** — confirm `nut` and the SGS stress are
   consistent across the resolution jump.
3. **Proper turbulence validation**, three channel cases:
   - **(i) uniform channel + LES** (no refinement) — the LES baseline.
   - **(ii) slab-refined channel + LES** — flat y-interfaces, NO edges (the
     existing wall-band geometry).
   - **(iii) block-refined channel + LES** — an embedded patch, WITH edges/corners
     (the `core_patch` geometry).

## How the LES path works (read before touching it)

`src/modules/les.f90`, `les_type`. Enabled by `[les] model = smagorinsky|wale`
(default `none`); `cs=0.10`, `cw=0.325`, `delta_scale=1.0`, `ibm_aware=.true.`
(config keys in `apply_les_value`, `src/modules/config.f90`).

- **`les%nut(0:nb+1, 0:nb+1, 0:nb+1, nBlocks)`** — eddy viscosity, per block with
  a trailing block index, halos included. Recomputed every substage.
- **Per step** (`src/main.f90:371`): `update_les_viscosity` (Smagorinsky/WALE
  strain from the velocity halos) → `exchange_scalar_halos(c, les%nut, blk)` →
  `momentum(..., les, les_prof)` which calls `add_les_momentum_correction`
  (adds the `2 nut S_ij` SGS stresses to the predictor RHS). `nut` also feeds the
  viscous timestep limit (`update_timestep_limits(blk, dns, c, les)`).
- **Filter width** `delta = delta_scale * (filter_x*filter_y*filter_z)`, where
  `filter_d = (1/d1_d)^(1/3)` is the local cell size from the per-block metric
  tables. These are sliced PER LEVEL, so `delta` automatically tracks the local
  cell size on each side of a 2:1 interface (coarse side ⇒ larger filter). Good
  — but VERIFY this is what actually happens at interface cells.
- **`nut` cross-level transfer**: `exchange_scalar_halos(c, les%nut, blk)` is
  called with DEFAULT args (no `interpProlong`, no `ifaceRow`), so at a 2:1
  interface `nut` RESTRICT = average of the fine samples, PROLONG = injection of
  the covering coarse value — it rides the SAME gather as the phi scalar but with
  plain injection prolong. The velocity `interpProlong` / blended-ghost machinery
  is NOT applied to `nut`. Whether injection is adequate (vs a `nut` jump at the
  interface) is an OPEN question to check (see risks).
- **IBM**: `ibm_aware` zeroes `nut` in solid cells (coef > 1e20 at any of the 6
  staggered faces). Irrelevant for the channel cases (no IBM) but keep in mind.

## Plan

### Phase A — LES ⊥ block decomposition (no refinement)

The block refactor's invariant is that results are EXACTLY independent of `nb`
and rank count (the redundant open-halo sweep). LES must preserve it. Gates,
build `-Mnofma`/`-gpu=nofma`, compare with `tools/compare_fields.py`:
- **nb-invariance**: a uniform channel + LES at `nb = default` vs a small `nb`
  (e.g. 4) — bit-exact (the `nut` field, exchange, and SGS stress are all
  per-block). This is the key new check (the LES kernels carry the trailing block
  index; confirm no block-boundary artefact in `nut` or the strain).
- **rank-invariance**: 1 vs 2 ranks, x- or z-split (NOT y) — bit-exact.
- **inert path** (regression for everyone else): with LES present but
  `model=none`, a run is bit-exact vs the pre-LES path (LES allocations/maps must
  not perturb the no-LES result).

### Phase B — LES ⊥ 2:1 interface (mechanics)

Before turbulence statistics, confirm the interface mechanics:
- **Inert-with-interface**: a SINGLE-LEVEL LES run with the interface code present
  is bit-exact vs Phase A (LES + const-1/2 both inert without a 2:1 interface).
- **`nut` transfer sanity**: on a refined patch, dump `nut` and check it is
  smooth/positive across the interface (RESTRICT=avg, PROLONG=inject); no
  negative or NaN `nut`, no order-1 discontinuity beyond the physical
  coarse/fine filter-width ratio. `MOBY_HALO_AUDIT` audits exchange-written
  halos against manufactured linear fields — run it with LES on to confirm the
  `nut` halo is filled correctly across levels.
- **`delta` per level**: verify `delta` at the coarse interface cell uses the
  coarse cell size and at the fine cell the fine size (the per-level metric
  slice) — print/inspect `nut` and `delta` for a manufactured strain.

### Phase C — turbulence validation (the three cases)

KEY SETUP ISSUE: **LES must be ACTIVE.** The current Re_tau 180, 128×64×128
channel is DNS-resolved, so `nut ≈ 0` and LES is inert (it would just reproduce
the DNS — a null test). Make the SGS term matter by running a deliberately
UNDER-RESOLVED grid: reuse the KMM180 restart interpolation
(`tools/make_channel_restart.py`) but onto a COARSE base (e.g. 64×48×64 or
96×48×96), where Re_tau 180 is under-resolved and `nut` is non-negligible. The
reference is the well-resolved no-LES DNS (the existing uniform-256 / -128
fields), filtered to the LES grid. (Alternative: raise Re_tau — but that needs a
new DNS restart source; coarsening reuses what we have. Decide first.)

- **(i) uniform + LES**: coarse base grid, `[les] model = smagorinsky` (then
  WALE). Gate: stable; mean U (log law) and resolved stresses match the
  filtered-DNS reference within the LES model's expected error; `nut` profile
  physically sane (→0 at the wall for WALE, finite for Smagorinsky without a
  wall damping — note which model needs van Driest). This establishes the LES
  baseline and the model constants on THIS code before adding refinement.
- **(ii) slab-refined + LES**: the existing wall-band geometry (two refine slabs
  spanning the full x–z plane, `validation/channel_interface/developed/`-style),
  `[les]` on. Flat y-interfaces only, NO edges. Gate: vs the uniform-fine LES
  reference — mean + stresses correct, and NO `nut` discontinuity/band at the
  y-interface (reuse the band tooling: `nut` is just another field for
  `patch_interface_stats.py` / `channel_band_profile.py`). reflux OFF (the
  open-item-1 conclusion).
- **(iii) block-refined + LES**: the `core_patch` geometry (embedded patch,
  edges+corners), `[les]` on. Gate: same as (ii) plus no `nut` band at the
  edges/corners. Reuse `run_gate5.py` (add `[les]` to the inis) +
  `patch_interface_stats.py` (it already aggregates any field by
  interior/face/edge/corner — point it at `nut` as well as u,v,w).

## Risks / open questions

- **`nut` injection prolong**: the coarse `nut` is injected onto the fine halo
  (not interpolated). At a sharp `nut` gradient this could step. If a `nut` band
  appears at the interface, try `exchange_scalar_halos(..., interpProlong=.true.)`
  for `nut` (the two-pass tangential interp already exists for phi) — but confirm
  it stays positive and conservative.
- **Strain at the interface cell** reads the const-1/2 velocity halos (mildly
  dissipative for small-scale tangential energy — the open-item-1 ~5% loss). The
  SGS model adds its own dissipation; check they don't compound into an over-
  dissipative interface band in `nut`.
- **`delta` discontinuity** is PHYSICAL (coarse cells have a larger filter), so a
  `nut` STEP across the interface is expected and correct — distinguish that
  physical step from a spurious band/oscillation.
- **Stability**: `nut` raises the viscous timestep limit cost; confirm the
  adaptive dt handles the coarse-side larger `nut*delta`.
- **Model choice**: Smagorinsky needs near-wall damping (van Driest) for a
  channel; WALE has the correct near-wall `nut~y^3` built in. Prefer WALE for the
  channel unless van Driest is wired in (check — it is not obviously in
  `les.f90`).

## Tooling (mostly reuse)

- `tools/make_channel_restart.py` — add a coarse-base target if coarsening (the
  `base`/`refined`/`patch` modes already exist; a coarser base is a small
  generalization or a new `--nx/--ny/--nz`).
- `validation/channel_interface/core_patch/run_gate5.py` + `base.ini` /
  `core_patch.ini` — add a `[les]` section; the driver already does two-leg
  transient+stats with snapshot dumps.
- `tools/patch_interface_stats.py` — already classifies face/edge/corner and
  aggregates ANY of un/vn/wn/pn; extend `VARS` to read `nut` (write `nut` to the
  field files, or add a `nut` dataset to the dumps) for the `nut`-band metric.
- `tools/plot_patch_slice.py` — cross-sections of `nut` across the interface.
- `MOBY_HALO_AUDIT=1`, `MOBY_STEPDIV=1` — halo/divergence sanity with LES on.
- Build: `module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3`,
  `./compile.sh cpu && ./compile.sh gpu`. One GPU case at a time; `pkill
  build_*/main` between runs.

---

## NEXT-SESSION PROMPT

> Read `docs/next_session_les.md` and the `edge-corner-validation` +
> `interface-validation-suite` memories. We are validating the LES path with
> block refinement on branch `claude/jacobi-interface`. Edge/corner (no-LES)
> validation is already DONE.
>
> Work in this order, gating each before the next:
> 1. **LES ⊥ block decomposition** (no refinement): nb-invariance (default vs
>    nb=4) and 1-vs-2-rank bit-exact for a uniform channel with `[les] model =
>    wale`; plus the `model=none` inert-path regression. Build `-Mnofma`.
> 2. **LES ⊥ 2:1 interface mechanics**: single-level-LES-with-interface-code
>    bit-exact; `MOBY_HALO_AUDIT` with LES on; confirm `nut` RESTRICT/PROLONG and
>    the per-level `delta` are sane on a refined patch (positive, no NaN, the
>    physical coarse/fine filter-width step only).
> 3. **Turbulence**, three cases, LES actually ACTIVE (coarsen the base grid so
>    Re_tau 180 is under-resolved — reuse the KMM180 interpolation; decide
>    coarsen-vs-raise-Re first): (i) uniform+LES baseline vs filtered DNS;
>    (ii) slab-refined+LES (flat, no edges); (iii) block-refined+LES
>    (`core_patch`, edges). reflux OFF. Reuse `run_gate5.py` + `base.ini` /
>    `core_patch.ini` with a `[les]` section, and `patch_interface_stats.py` /
>    `plot_patch_slice.py` pointed at `nut` as well as u,v,w. Gate: stable; mean
>    + stresses match the reference; NO `nut` band at the interface/edges/corners
>    beyond the physical filter-width step.
>
> Start by deciding the LES-active setup (coarse grid resolution + model:
> WALE preferred for the channel near-wall behaviour) and confirming the
> `model=none` inert path is bit-exact, then proceed through the phases.
