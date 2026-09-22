# ZPG TBL — record of the tests performed

Chronological record of the study that produced the shipped case. The driving
question: **our DNS c_f sat ~5 % below the spectral reference at matched Re_θ —
why, and can appropriate setup close it?** The answer decomposed cleanly into
tripping + streamwise resolution + an intrinsic FD-vs-spectral floor.

Non-dimensionalisation throughout: length by inlet δ*₀, velocity by U∞, so
Re_δ*,0 = 1/ν = 450. "t.u." = time units (δ*₀/U∞).

## 0. Precursor — laminar Blasius + top boundary condition

A quasi-2D laminar Blasius layer validated the inflow/outflow setup and, crucially,
the **top BC**: only an `outlet` top (Dirichlet-p, which frees the entrainment v)
pins a true ZPG (U_e = 1.0000, θ within ~1 % of Blasius). Neumann-p variants
(displacement / prescribed-v) let the freestream pressure drift → a weak favourable
gradient and θ ~9 % high. Lesson: **the freestream is set by the top pressure.**

## 1. Turbulent setup

- Inlet: laminar Blasius at Re_θ = 173.7 (= Re_δ*=450 / H_Blasius); trip a bit
  downstream (x₀ = 15) via the **Schlatter–Örlü random wall-normal volume force**
  (implemented on-device: only the Fourier coefficients cross PCIe — 7.8× faster).
- Grid: a dedicated **blayer** wall-normal distribution (wall clustering over a
  resolved band, geometric coarsening above) so the tall ZPG domain (ly = 100 δ*₀)
  stays affordable. Standard DNS spacing.
- Reference the spectral data extracts: from the ZPG momentum integral
  ν = u_τ²/(dRe_θ/dx) → 1/ν ≈ 456 ≈ 450, confirming the shared nondimensionalisation.

## 2. Baseline production DNS (trip 0.15, standard grid)

2048 × 160/176 × 192, developed and averaged ~3400 t.u. Canonical TBL, but at
Re_θ = 677 vs SIMSON: **c_f −4.8 %, H +1.4 %**, u′_rms peak matching. The gap was
the puzzle.

## 3. Ruling out measurement artifacts

- **Wall-shear estimate** — 1st-order (U₀/y₀) vs 2nd-order 3-point stencil: identical
  to 0.00 %, because y₀⁺ ≈ 0.1 is deep in the linear sublayer. Not the cause.
- **Re_θ / abscissa shift** — the momentum integral dθ/dx = c_f/2 holds to ~2–3 %
  (a ~18 % Re_θ error would be needed to explain the gap as a shift; ruled out), and
  Re_δ* = 1032 vs the reference 998 (+3.4 %, all from H) — not +18 %. Re_θ is right.
- **Pressure gradient** — d⟨p⟩_yz/dx ≈ 3e-6 in the developed region (genuine ZPG);
  not a residual streamwise gradient.
- **von Kármán closure** — c_f/2 = dθ/dx + (δ*+2θ)/U_e·dU_e/dx closes to <0.5 % in the
  developed region (term2 ≈ 0, true ZPG); c_f from the integral balance matches the
  direct wall-shear c_f. So the deficit is **real mean-flow physics**, not a
  post-processing artifact.

## 4. The trip × resolution 2×2

At Re_θ = 677 vs SIMSON (~4000 t.u. windows):

| | standard grid | fine outer BL (Δy⁺_max 6→4) |
|---|---|---|
| **trip 0.15** | c_f −4.8 %, H +1.4 % | c_f −4.1 %, H +1.4 % |
| **trip 0.03** | c_f −3.5 %, H +0.1 % | c_f −3.2 %, H +0.2 % |

- **Tripping is the H driver and worth ~1.3 pts of c_f.** Reducing the trip 5×
  (0.15→0.03) puts H exactly on the reference and raises c_f; the transition
  overshoot collapses from c_f≈0.018 to ≈0.006 (matching the reference's gentle trip).
  Our original trip was simply more aggressive than the reference's.
- **Outer wall-normal resolution is minor** (~0.3–0.7 pts, no H change): Δy⁺_max≈6
  was already adequate. (Implemented by holding the wall spacing fixed —
  `dyw_plus` 0.15→0.28 as ny grows — so the diffusive/Peclet dt limit, set by the
  smallest wall cell, is unchanged: dt stays 0.02.)

## 5. Streamwise resolution — the dominant remaining lever

Best case (trip 0.03, std grid) had Δx⁺≈8 vs Δz⁺≈3.7 — anisotropic. Refining to
**Δx⁺≈4 ≈ Δz⁺** (nx 2048→4096) moved c_f from −3.5 % to **−1.6 %** at the full
window (H stays −0.4 %). So the residual was mostly **streamwise under-resolution**:
Δx⁺≈8 was too coarse for the 2nd-order FD to capture the streamwise gradients that
set the wall shear. The IC was the developed low-trip field interpolated in x
(`make_finex_restart.py`; u is x-face-staggered → face→face, v/w/p centre→centre).

Caveat learned the hard way: **provisional c_f over-reads at short windows** (this
case: −0.8 % @ 1000 t.u. → −1.1 % @ 2360 → −1.6 % @ 3920). Quote only near-full
windows.

## 6. Decomposition of the original ~5 % c_f gap

| contribution | Δc_f |
|---|---|
| over-aggressive tripping (0.15 → 0.03) | ~1.3 pts, and the entire H discrepancy |
| streamwise under-resolution (Δx⁺ 8 → 4) | ~1.9 pts |
| wall-normal outer resolution | ~0.5 pts (minor) |
| residual | **~1.6 %** — the intrinsic 2nd-order-FD-vs-spectral floor |

Not a solver deficiency: two fixable setup choices plus a small method floor.

## 7. Red-black SOR pressure solver (this case's smoother)

The original single-level red-black Gauss-Seidel projection was restored as a
**selectable** alternative (`[pressure] solver = jacobi | redblack`) to the
damped-Jacobi/Chebyshev path — kept mutually exclusive with the 2:1 block interface
(whose current consistency machinery is tied to the Jacobi phi buffer). Verified:
the Jacobi path stays **bit-exact**; red-black CPU==GPU; and on this case it is

- **stable** at niter=6, sor=1.5 over ~5000 t.u. — where **Chebyshev+niter=6+outlet
  was not** (the 2Δx pressure mode);
- **~1.8× faster** than Chebyshev-Jacobi niter=12 on the same 5090 (0.47 vs 0.85
  s/step), driven by needing half the iterations;
- **statistically identical** to the Jacobi run (c_f/H agree to <0.03 %) — both
  converge the projection to the same divergence-free flow.

This is the solver the shipped case uses.

## 8. Block-overhead / xz 2:1-refinement timing test

Asked whether coarsening x,z **above** the boundary layer via 2D
(`refine_dims = xz`) 2:1 block refinement would save wall clock on this 138 M-cell
production grid. Measured 2026-08-06 (cetus RTX A6000, 400-step cold starts, matched
dt): **no — a ~16 % net loss** at the only block size this grid permits. The 2:1
refinement delivers its cell saving, but the block-decomposition halo overhead it
requires (redundant halo-layer sweeps) exceeds the saving here. The block
decomposition is exactly result-invariant (`nb` unset vs `nb = 16` give identical
runtime lines), as designed. So the shipped case runs **single-level** (`nb` unset).

## 9. Finer wall-normal grid — the shipped "finey" case

Once Δx⁺ was adequate (§5), the wall-normal outer resolution turned out to be a
**bigger** lever than §4's coarse-Δx estimate suggested: refining the resolved band
from Δy⁺_max ≈ 6 to ≈ 4.3 (**ny 176 → 224**, wall spacing held fixed via
`dyw_plus` 0.15 → 0.221, so dt stays 0.02) moved c_f from **−1.6 % to −0.4 %** at
the full ~4000 t.u. window (H −0.6 %). The errors are **non-additive**: at Δx⁺≈8
the streamwise error dominated and masked the y effect (§4 measured it as ~0.3 pt);
at Δx⁺≈4 the y-refinement becomes the dominant remaining lever (~1.2 pt). This
`ny=224` grid is the **shipped case** (`make_finey_grid.py` ports the solver's
blayer line exactly; the IC is a developed field interpolated in y via
`make_finewall_restart.py`).

## 10. Trip-amplitude sweep — the u′_rms peak is not the trip

The one residual after §9 was the near-wall **u′_rms peak sitting ~3 % above
SIMSON** (finey 2.73 vs 2.63). Two runs on the ny=224 grid tested whether the trip
sets it (HoreKa, 4× A100, 6000 t.u. windows):

| trip_amp | 0.024 (−20 %) | 0.03 (finey) | 0.18854 (SIMSON params) |
|---|---|---|---|
| u′_rms peak | 2.709 | 2.727 | 2.724 |
| c_f @ 677 | −1.2 % | −0.4 % | −3.1 % |

An **8× span in trip amplitude leaves the peak flat at ~2.71–2.73.** The strong
SIMSON-parameter trip (applied verbatim; our `g(z)` is unit-rms so it is ~6× too
hot) also over-trips — a violent transition overshoot (c_f 17e-3 at Re_θ 304) and a
~3 % c_f deficit — without helping the peak. **Trip amplitude sets the transition
location/overshoot, not the developed near-wall peak.** finey (gentle trip) stays
the best match.

## 11. Code comparison — the u′_rms peak is an FV/FD signature

Reduced through one common post-processor (`compare_codes.py`) against **SIMSON
(spectral), CaNS and AMPHIBIOUS** at Re_θ ≈ 677:

| code | c_f | H | u′_rms peak |
|---|---|---|---|
| SIMSON (spectral) | 0.00471 | 1.507 | 2.631 |
| CaNS | 0.00463 | 1.503 | 2.706 |
| AMPHIBIOUS | 0.00460 | 1.508 | 2.723 |
| mobydiff (finey) | 0.00469 | 1.498 | **2.727** |

- **mobydiff's c_f is the closest of the non-spectral codes to SIMSON** (−0.4 %).
- **All three finite-volume/difference codes overshoot the spectral u′_rms peak by
  ~3 %** (2.71–2.73), collapsing onto each other. The peak excess is therefore an
  **intrinsic 2nd-order-FV/FD-vs-spectral signature**, reproduced by two independent
  established DNS codes — not a mobydiff defect. Combined with §10 (trip-independent)
  this fully closes the question.

The mobydiff data are exported to the reference NetCDF format at
`assets/mobydiff/xyz_4096_224_192/data.nc` (`make_mobydiff_nc.py`).

## Bottom line

The original 5 % c_f gap to the spectral reference was **not** a limitation of the
2nd-order FD solver. It decomposed into an over-aggressive trip (which also caused
the entire H discrepancy), streamwise under-resolution, and wall-normal outer
resolution — all fixable — leaving only the intrinsic FV/FD-vs-spectral floor. The
shipped **finey** configuration (gentle trip, Δx⁺≈Δz⁺≈4, Δy⁺_max≈4.3, red-black
niter=6) reproduces the SIMSON spectral reference in c_f (−0.4 %), H, mean profile
and Reynolds stresses to ~1–2 %, and sits squarely with CaNS and AMPHIBIOUS — the
residual ~3 % near-wall u′_rms peak is a signature shared by all non-spectral codes.
