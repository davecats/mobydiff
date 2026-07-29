# turbulent_lowtrip_finewall — case (iii): reduced trip + fine outer BL

Combines the two variants:
- **reduced tripping** from `../turbulent_lowtrip` (`trip_amp = 0.03`, 1/5), and
- **fine outer-BL resolution** from `../turbulent_finewall` (ny=256, dyw_plus=0.28
  → Δy⁺_max≈3.85, wall spacing kept so dt=0.02 — no Peclet penalty).

Everything else matches turbulent_two (geometric x, skew convection, niter=12).
2048×256×192 = 100.7 M cells.

**IC:** the same interpolated field as case (ii),
`../turbulent_finewall/restart_interp.h5` (the developed turbulent_two field on
the ny=256 grid), restarted here with the weaker trip.

**Purpose:** isolate whether the residual c_f offset that survives the
reduced-trip case (case (i): −3.5% vs SIMSON, H matched) is a **resolution**
effect. If the fine outer BL closes the remaining gap, the offset was
resolution; if not, it's the 2nd-order-FD-vs-spectral floor.

**Run (two phases):**
```bash
mpirun -n 1 build_gpu/moby_solve production.ini        # phase 1: re-equilibration
mpirun -n 1 build_gpu/moby_solve production_stats.ini  # phase 2: statistics
```
Ran on istmcetus A6000 (GPU 1) in parallel with case (ii) on corax; the campaign
driver is `chain_cetus.log`. Analysis scripts / references as in turbulent_two.
