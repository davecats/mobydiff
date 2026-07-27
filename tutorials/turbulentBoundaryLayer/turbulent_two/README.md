# turbulent_two — refined ZPG TBL DNS matched to the spectral reference

A second production run of the zero-pressure-gradient turbulent boundary
layer, set up to **reproduce the in-house pseudo-spectral TBL DNS**
(`../turbulent/tbl_uncontrolled.mat`, Schlatter & Örlü-class, Re_θ up to
~2500) up to **Re_θ = 900**, with two resolution/robustness improvements
over `../turbulent`.

## What the spectral data fixes

The reference `.mat` stores mean U⁺(y⁺) at 2400 streamwise stations plus
u_τ, Re_θ, Re_τ. It contains **only the turbulent region** — the first
station is x = 97.66 δ*₀ at Re_θ = 383.6 and u_τ decreases monotonically
from there, so the laminar inlet and the trip are upstream and cropped out.
Reconstructing the setup from the flow physics:

- **Length scale / Reynolds number.** From the ZPG momentum integral
  ν = u_τ²/(dRe_θ/dx); in the developed region this gives 1/ν = 456 ≈ 450,
  i.e. lengths are in the inlet displacement thickness δ*₀ with
  **Re_δ*,0 = 450**, U∞ = 1. → `[flow] re = 450`.
- **Blasius inlet at Re_θ = 173.7** (= Re_δ*/H_Blasius = 450/2.591). This is
  the boundaryLayer-case default (`theta_in = 1/H_Blasius`).
- **Trip at x₀ = 15 δ*₀** (boundaryLayer default). Transition completes by
  Re_θ ≈ 390, matching the spectral data's first station.

So the physics, nondimensionalization, inlet, and trip already match. The
only mismatch was streamwise extent and x-resolution.

## Changes vs `../turbulent`

| | `../turbulent` | `turbulent_two` | why |
|---|---|---|---|
| lx | 500 (→ Re_θ 734) | **750** | reach **Re_θ = 900** (at x ≈ 645, before the outflow zone) |
| nx | 1280 (uniform) | **2048, geometric** | hold Δx⁺ ≈ 8 (see below) |
| ny | 160 | **176** | keep Δy⁺_wall ≈ 0.2 over the taller resolved band |
| resolved_height | 30 | **36** | δ₉₉(Re_θ=900) ≈ 16, so ≈ 2.25 δ₉₉ |
| convection | divergence | **skew** | skew-symmetric form, energy-neutral for any advecting field |
| ly | 100 | 100 | domain height kept large |
| lz / nz | 32 / 192 | 32 / 192 | Δz⁺ ≈ 3.5 |

**Geometric x-grid.** A uniform Δx leaves the high-u_τ region
under-resolved: because Δx⁺ = Δx·u_τ/ν and u_τ peaks in the transition zone
(and falls downstream), a uniform grid gives Δx⁺ from 8.2 at the outlet up
to **15.8 in the transition region**. `[grid.x] distribution = geometric,
stretch = 1.5` grows the cells downstream (finest at the inlet), holding
Δx⁺ ≈ 8:

```
             Δx⁺(transition)   Δx⁺(core, Re_θ 390–900)   Δx⁺(outlet)
uniform            15.8              8.3–8.9                 8.2
geometric R=1.5    ~12               7.4–8.0                 ~9.4
```

`stretch = R` is the last/first cell-width ratio (dx₀ = 0.297 → dx_end =
0.445; R = 1 recovers a uniform grid). The core and the physically-compared
range (Re_θ ≳ 390) sit at Δx⁺ ≈ 7–8; the transition band (Re_θ < 390, not
compared) is roughly halved from the uniform case.

**Skew-symmetric convection.** `[flow] convection = skew`: the divergence
form is discretely energy-neutral only for exactly divergence-free
advection, which the incremental projection never grants; the skew form is
neutral for any advecting field (`docs/next_session_skew_convection.md`).
This is the first turbulent-BL exercise of the ported skew term.

Grid: 2048 × 176 × 192 = **69.2 M cells**.

## Run recipe (two phases, same as `../turbulent`)

```bash
module load toolkits/nvhpc/25.9
# Phase 1 — transition + transient wash-out (statistics OFF), ~100k steps:
mpirun -n 1 ./build_gpu/moby_solve production.ini
# Phase 2 — span+time statistics (niter=12), restart from production_100000.h5:
mpirun -n 1 ./build_gpu/moby_solve production_stats.ini
```

Estimated ~0.30–0.37 s/step on the RTX 5090 (istmcorax) → phase 1 ≈ 8–10 h,
phase 2 ≈ 19–24 h; ~1.2–1.4 days total. Chebyshev + niter 8/12 was stable
on the `../turbulent` run at this dt; keep niter ≥ 8 (the Chebyshev +
Dirichlet-p-outlet 2Δx resonance shows up only at niter = 6 / large dt).

## Analysis

```bash
python3 bl_stats.py production_stats.h5 --plot first_stats.png --retheta 677
                                        # 6-panel; SIMSON overlaid on every panel
python3 compare_passivewall.py    production_stats.h5   # c_f/H trends + U+ + Reynolds stresses vs SIMSON
python3 compare_spectral.py       production_stats.h5   # U+, c_f, H vs tbl_uncontrolled.mat (mean only)
python3 compare_schlatter_orlu.py production_stats.h5   # vs the KTH Re_θ=677 profile
python3 dpdx.py        production_stats.h5              # <p>_yz(x), d<p>_yz/dx (ZPG diagnostic)
python3 resolution.py  production_stats.h5 production_p2_200000.h5   # Δx+/Δz+/Δy+ vs Re_θ
python3 viz_flowfield.py production_p2_200000.h5        # instantaneous fields
```

Reference data:
- **`passivewall.hdf5`** — SIMSON spectral ZPG-TBL DNS (Schmitt/KIT; Re_δ*,0=450,
  same as this case; trip at x=10; passive-scalar wall). Carries mean profiles
  AND Reynolds stresses over the full x-development, so it is the primary
  reference (`bl_stats.py --ref`, `compare_passivewall.py`). **Not committed**
  (~99 MB) — keep it in this directory locally.
- `tbl_uncontrolled.mat`, `ref_schlatter_orlu_Re670.prof` — symlinked from
  `../turbulent` (mean-only spectral and the KTH Re_θ=677 profile).

The c_f offset vs the SIMSON reference is ~−4.6% at Re_θ=677 and shrinks
downstream (−8% near the trip → −3.6% at Re_θ=880): a tripping/development-history
effect (our trip is at x=15 vs SIMSON's x=10), not resolution or a Re_θ bias.
The Reynolds stresses match to ~2% (u'/v'/w'/−u'v' peaks 2.68/1.00/1.30/0.88
vs 2.63/1.02/1.30/0.88).
