# lowtrip_srdgridy_finex — best case + finer streamwise grid

Takes the campaign winner — **low trip (0.03) + standard (blayer) y-grid** =
`../turbulent_lowtrip` — and refines only the **streamwise** resolution:
**nx 2048 → 4096** so Δx⁺ ≈ Δz⁺ in the developed region (was Δx⁺≈8 vs Δz⁺≈3.7;
now Δx⁺≈4). Everything else identical (geometric x stretch 1.5, ny=176 blayer,
nz=192, skew convection, niter=12). 4096×176×192 = **138 M cells**.

**Why:** the 2×2 trip×resolution campaign left a residual ~3.5% c_f offset vs the
spectral reference under the best trip and finer *wall-normal* grid. This tests
whether it is **streamwise** under-resolution. dt is unaffected — finer x doesn't
touch the diffusive/Peclet limit (set by the wall cell) and the CFL limit stays
< 0.8 — so dt = 0.02, no penalty.

**IC:** case (i)'s developed field `../turbulent_lowtrip/production_p2_550000.h5`
(t=11000) **interpolated onto the finer x-grid** → `restart_finex.h5`, built by
`make_finex_restart.py` (u is x-face-staggered → interp face→face; v/w/p x-centred
→ centre→centre; geometric x rebuilds from nx+stretch, no attr to resync). A
laminar cold start is not an option — trip 0.03 only *sustains* turbulence, it
can't trip from Blasius. `restart_finex.h5` (4.4 GB) is generated, not committed.

**Run (two phases, corax):**
```bash
mpirun -n 1 build_gpu_corax/moby_solve production.ini        # phase 1: re-equilibration on the finer x
mpirun -n 1 build_gpu_corax/moby_solve production_stats.ini  # phase 2: statistics
```
