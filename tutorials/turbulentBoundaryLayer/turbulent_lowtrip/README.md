# turbulent_lowtrip — reduced-tripping variant of turbulent_two

Identical to `../turbulent_two` **except the trip intensity is reduced by a
factor 1/5** (`trip_amp` 0.15 → **0.03**). Same grid (2048×176×192, geometric
x, blayer y), same skew convection, same niter=12 Chebyshev.

**IC:** restarts from the developed `../turbulent_two/production_p2_250000.h5`
(t=5000), so the layer only re-equilibrates to the weaker trip rather than
transitioning from scratch.

**Purpose:** a direct test of the tripping/development-history hypothesis for
the c_f offset. Our trip overshoots to c_f≈0.018 at Re_θ≈300 while the SIMSON
reference (trip x=10) only bumps to ~0.006; a gentler trip here should move our
development toward the reference and shrink the ~4–5% c_f offset at matched Re_θ.

**Run (two phases, on corax):**
```bash
mpirun -n 1 build_gpu_corax/moby_solve production.ini        # phase 1: re-equilibration (stats off)
mpirun -n 1 build_gpu_corax/moby_solve production_stats.ini  # phase 2: statistics (restart from the phase-1 end field)
```
Phase 1 continues the step counter from 250000 → 350000; phase 2 restarts from
the latest `production_<step>.h5`. Analysis scripts and reference data are as in
`../turbulent_two` (passivewall.hdf5 symlinked).
