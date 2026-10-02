# LES validation with block refinement (2:1 interface)

**STATUS: VALIDATED (2026-06-30)** — the LES (`[les] model = wale`) path works
across the block decomposition and the 2:1 refinement interface (channel flow);
see "Results" below. Open item 2 of `docs/next_session_edges_les.md` (item 1,
edge/corner no-LES, is DONE in `../core_patch/`). Settings throughout: **WALE**,
`momentum_reflux = false`, `interface_constant_half = true` (both keys removed
2026-07-01; that configuration is now hardwired, so current inis omit them).

> **CAVEAT — LES is still IBM-UNAWARE in practice.** This validation is channel
> flow only (no immersed boundary). The LES↔IBM coupling — the `ibm_aware`
> solid-cell `nut` masking in `les.f90` — was never exercised here and is
> **untested with block refinement and across the 2:1 interface**. Validate the
> LES+IBM+refinement combination before trusting LES on IBM bodies.
>
> *(Lifted 2026-06-30: the LES↔IBM coupling was validated with and without
> refinement in `../les_ibm/`.)*

## Why WALE + a coarse grid

The standard 128×64×128 Re_tau-180 channel is DNS-resolved, so `nut≈0` and LES is
inert (a null test). These cases run a deliberately **coarse 64×48×64** grid
(dx⁺≈35, dz⁺≈18, dy⁺_wall≈2.5) where Re_tau 180 is under-resolved and the SGS
term is **active** (WALE `nut` peaks ~5× molecular). WALE is used for the correct
near-wall `nut~y³` (no van Driest needed) — verified `nut/ν = 1e-4` at the wall
cell, rising through the log layer, peaking in the core. ICs are interpolated from
the KMM180 DNS restart (`tools/make_channel_restart.py`, coarse base via the new
`--nx/--ny/--nz`).

## Status (2026-06-30)

### Mechanics — VALIDATED (this session, CPU + GPU)

- **LES ⊥ block decomposition** (no refinement) — all bit-exact (max diff 0.0;
  same binary, FMA baked in identically):
  - inert path: `[les] model=none` ≡ no `[les]` section.
  - nb-invariance: WALE, 1 block vs nb=8 (2048 blocks) vs nb=4 (16384 blocks).
  - rank-invariance: WALE, nb=8, 1 vs 2 ranks (x-split).
- **LES ⊥ 2:1 interface mechanics** (`../core_patch/` geometry, CPU):
  - `nut` cross-level exchange exact — `MOBY_HALO_AUDIT` (a diagnostic hook,
    removed in the 2026-06-30 cleanup; extended at the time to audit the
    `nut` scalar exchange): 0 bad / maxErr 0 / 0 unwritten at all levels.
  - `delta` per-level exact — same-footprint `delta_base/delta_fine = 2.000000`
    (midpoint subdivision); the fine filter width is exactly half the base.
  - `nut` positive, finite, no NaN; stable; no spurious interface band — the
    `nut` step across the interface is the **physical filter-width step**
    (`nut ∝ delta²`, ratio ≈ 4×; measured fine 0.060 → coarse 0.222 ≈ 3.7×),
    not a band.
- **Turbulence setup, all three case types run active+stable** (GPU, 500-step
  spot checks): uniform / slab / patch all stable (div decaying, mass ~1e-16),
  `nut` sane (≥0, no NaN, coarse-side > fine-side = physical step), LES active
  (`nut_max/ν` 4.4–5.1).

### Statistics campaign — READY TO RUN (developed-flow)

The quantitative gate (mean U log law + resolved stresses + `nut` profile,
time-averaged over t=5..25). Four cases via `run_les.py` (two legs each:
transient stats-off, then a stats leg with `channel_stats` ON + `nut` snapshots);
analyse with `les_stats.py`. The reference is the **128³ no-LES DNS**
(`reference.ini`, DNS-adequate at Re_tau 180, = the 64³ coarsened 2:1).

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
./compile.sh gpu          # build_gpu/main must carry the nut output (this branch)
cd validation/channel_interface/les
# on this host the system mpirun is wrong; pass the hpcx launcher:
MP=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/comm_libs/13.1/hpcx/hpcx-2.25.1/ompi/bin/mpirun
python3 run_les.py --arch gpu --case all --mpirun "$MP"   # reference,uniform,slab,patch
python3 les_stats.py                                       # -> les_profiles.png + ratios
```

`run_les.py` generates each IC (`--mode base` for reference/uniform; `--mode
patch` with the .ini's own refine boxes for slab/patch, so the IC leaf table is
bit-identical to the solver's). Defaults: `--t-transient 5 --t-average 20
--snap-interval 800`. `channel_stats` (sample every 50 steps) carries the
velocity statistics; snapshots carry `nut`.

Gates:
- **(i) uniform**: `les_profiles.png` — U⁺(y⁺) on the log law; resolved stresses
  track the DNS reference (LES sits slightly below the reference stress peaks by
  the SGS contribution, since the reference is *unfiltered* DNS — expected); the
  printed CORE ratio table LES/reference for U, u', v', w', −⟨u'v'⟩ near 1.
  `nut(y)` physically sane (→0 at the wall, WALE; verified `nut/ν=1e-4` at the
  wall cell, peak ~0.3 in the core).
- **(ii) slab**: same, plus NO `nut` band at the flat y-interface — `nut(y)`
  steps by ~the physical filter-width ratio (`delta²` ≈ 4× across a 2:1 face) and
  no more (the `nut(y)` panel marks it). Velocity stats from `channel_stats` are
  per-level (x,z homogeneous) so they cross the interface cleanly.
- **(iii) patch**: edge/corner band metric via the nut-aware
  `tools/patch_interface_stats.py` (`les_stats.py` runs it: patch run vs the
  uniform LES run as the matched base control) — gate: BAND ratios ≈ 1 for u,v,w
  AND `nut`, edge/corner not worse than face. Cross-sections:
  `tools/plot_patch_slice.py runs/patch/stats/patch_<N>.h5
  runs/uniform/stats/uniform_<N>.h5 --out slice.png` (renders a `nut` panel).

## Files

- `reference.ini` (128³ DNS), `uniform.ini`, `slab.ini`, `patch.ini` (64³ + WALE,
  reflux off, `channel_stats` on). `slab.ini` boxes give a **symmetric** fine band
  (wall block-rows {0,1,4,5}, interfaces at y⁺≈88 and mirror; verified via
  `make_channel_restart.box_leaf_table`).
- `run_les.py` — campaign driver (4 cases, two legs, channel_stats + nut snapshots).
- `les_stats.py` — analysis (channel_stats profiles + log law + nut(y) + patch band).
- `nut_interface_slice.png` — the Phase-B mechanics `nut` cross-section.
- `wall_jump.ini`, `check_wall_jump.py`, `run_wall_jump.sh` — the level jump that
  touches a wall (last section).

## Results (developed run, t=5..25, GPU)

Figures: `les_profiles.png` (profiles), `fig_{patch,slab}_slice_*.png`
(instantaneous cross-sections, all vars incl nut), `fig_{slab,patch}_rms_slice.png`
(time-avg fluctuation rms + mean nut cross-sections, `fig_interface_rms.py`).

**LES quality (uniform vs filtered-DNS reference).** CORE (0.3<y<1.7) ratio
LES/reference: U 1.01, −⟨u'v'⟩ 0.99 — mean velocity (log law recovered) and the
resolved Reynolds shear stress match the DNS. Normal stresses show the classic
coarse-LES/WALE bias at Re_tau 180: u' +5% (streak energy piles up), v'/w' −10%
(under-resolved cross-stream motions). Part of the v'/w' deficit is the legitimate
filtering gap (reference is *unfiltered* DNS). `nut` is physical: →0 at the wall
(WALE y³, no van Driest), peak ~0.26 ν in the core. Acceptable LES.

**Flat interface (slab).** `nut` steps down by the physical `delta²` (~4×) into the
fine wall bands — a clean sharp step, NO overshoot/band (see `fig_slab_rms_slice`
mean-nut panel and the `nut(y)` profile). The fine bands resolve MORE velocity
fluctuation (closer to DNS) — refinement improves the near-wall stresses; the
transition is smooth with no spurious rms ridge at the interface.

**Edge/corner interface (patch).** Band metric (patch vs uniform LES control):
velocity AND `nut` band ratios 0.98–1.03 at face/edge/corner — NO band,
edge/corner not worse than face. The mean-`nut` cross-section shows the fine core
box as a crisp low-`nut` square with no surrounding halo. The const-1/2 +
reflux-off interface treatment carries over cleanly to LES.

**Caveat — pressure null mode (NOT LES/interface related).** ALL runs, including
the no-LES 128³ reference, carry a large velocity-decoupled pressure null mode
(pn std 2e5–1.3e6, |max| up to 2e7) from the under-converged niter=6
Chebyshev-Jacobi projection (see the `pressure-volume-average-drift` note). It is
gradient-free: velocity, divergence and all turbulence statistics are healthy and
match DNS. It makes the pressure field meaningless and inflates the patch `p` band
ratio (1.2–1.5) — DISREGARD the `p` column; u/v/w/nut are the valid ones. For runs
where pressure matters, raise niter or pin/zero-mean the pressure.

## Notes

- The reference is *unfiltered* DNS, so LES resolved stresses sit slightly below
  the reference peak by the SGS contribution; mean U and −⟨u'v'⟩ (total momentum
  balance) should match closely. For a stricter test, box-filter the DNS field to
  64³ before computing reference stresses (not done here).
- `nut` is written to a field snapshot only when LES is active (a separate
  `fdm_h5_append_nut` HDF5 call); no-LES output is byte-identical.

## A level jump that touches a wall (2026-10-02)

`./run_wall_jump.sh` (`wall_jump.ini`, `check_wall_jump.py`; one step, CPU,
seconds). Every case above keeps its 2:1 interfaces away from the walls or
parallel to them. Here one refinement box spans the whole height and span, so
the two interfaces are x-normal planes that run INTO both walls. The field is
the x- and z-invariant mean profile and the model is Smagorinsky (non-zero in
plain shear), so within a level the eddy viscosity of a wall row has one value
and the column next to the jump shows what the SGS kernel reads there.

What it reads: the velocity-gradient tensor of the cell `(nb, 1)` takes
`du/dy` from `u(nb+1, 0)`, the x halo of the wall GHOST row. A same-level
neighbour has always delivered it (the tangential extension copies its ghost
row). Across a level jump nothing wrote it before 2026-10-01
(`docs/next_session_outlet.md`, increment O4) and it kept its initial value 0.

Deviation of the wall-row eddy viscosity in the column next to the jump from
the median of its level, both walls alike:

| | main `154e48f` (halo never written) | now |
|---|---|---|
| fine side (prolongation of the coarse ghost row) | 25.8 % | 0.94 % |
| coarse side (restriction of the fine ghost row) | 13.4 % | 6.8 % |
| every other column, either level | ≤ 1.3 % | ≤ 1.3 % |

The gate is ≤ 2 % fine, ≤ 8 % coarse. The 1.3 % is the first fine column on
the low-x side of the band, whose halo column is the injected coarse value
(the const-1/2 transfer); it is the same before and after.

**The coarse-side 6.8 % is a known first-order residual, not fixed.** A coarse
ghost cell covers two fine rows, of which only one exists beyond the wall, so
that halo row is restricted from the single fine ghost row: the coarse block
receives `-u_fine(1)` where the mirror of its own halo column is
`-(u_fine(1) + u_fine(2))/2`, about twice as large in the viscous sublayer.
The wall shear seen through that one face is 12.5 % low, the cell average
6.25 % low; the measured 6.8 % is that on the stretched line. Exact for a
uniform field (the freestream gates), first order otherwise. What it would
take: unequal sample weights `(3/2, -1/2)` on the fine rows `(0, 1)` in the
gather kernels (exact for a linear profile, equal weights are all they carry
today), or the boundary row applied to the halo column after the exchange.
Not done: with WALE the wall-cell eddy viscosity is ~1e-4 of the molecular
one, resolved RANS multiplies the wall-row strain by a zero eddy viscosity,
and the momentum stencil does not read that halo at a wall.

A 40-step WALE run from a disturbed profile on the same layout differs
between the two binaries by 1.3e-8 in velocity and 2.1e-10 in `nut`, all of it
starting in those columns: the first case in the tree in which O4 changes a
result, and the direction is the one above.
