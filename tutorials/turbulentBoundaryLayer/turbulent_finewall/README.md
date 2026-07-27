# turbulent_finewall — Δy⁺_max=4 variant of turbulent_two

Identical to `../turbulent_two` **except the OUTER boundary layer is refined** so
the largest cell within the layer is **Δy⁺_max ≈ 3.85** (ny 176 → **256**). x and
z grids unchanged (2048×192, geometric x). 2048×256×192 = 100.7 M cells. Same trip
(0.15), skew convection, niter=12.

**Key: the wall spacing is held fixed, so there is no Peclet penalty.** The
diffusive/Peclet stability limit dt ≲ pecletmax·min(Δy²)/ν is set by the
*smallest* cell (the wall cell). The wall resolution was already good, so we keep
it: `dyw_plus` is raised 0.15 → **0.28** as ny grows, which holds Δy⁺_wall at
0.226 (turbulent_two: 0.218) and the wall cell at Δy=1.07e-2 (was 1.03e-2). The
Peclet dt limit stays 0.0257 ≥ dtmax, so **dt = 0.02, unchanged** — the extra
points only refine the outer BL (Δy⁺_max 6.2 → 3.85), which is the resolution we
actually need. (An earlier version refined the wall too, dropping dt to 0.0074 —
fixed.)

**IC:** the developed `../turbulent_two/production_p2_250000.h5` field
**interpolated onto the finer y-grid** → `restart_interp.h5`, built by
`make_finewall_restart.py`:
```bash
python3 make_finewall_restart.py ../turbulent_two/production_p2_250000.h5 new_y_nodes.npy restart_interp.h5
```
The staggered layout is respected (u/w/p interpolate cell-centre→cell-centre, v
lower-face→lower-face); `new_y_nodes.npy` is the exact ny=256 blayer node line
from the solver's own grid build, so the field lands on the grid the restart
rebuilds from config. `restart_interp.h5` (3.2 GB) and `new_y_nodes.npy` are
generated, not committed.

**Purpose:** check the c_f/stresses under stricter wall-normal resolution — is
the ~4–5% c_f offset sensitive to Δy⁺_max, or is it purely tripping/history?

**Run (two phases, on corax):**
```bash
mpirun -n 1 build_gpu_corax/moby_solve production.ini        # phase 1: re-equilibration on the finer grid
mpirun -n 1 build_gpu_corax/moby_solve production_stats.ini  # phase 2: statistics
```
Analysis scripts / reference data as in `../turbulent_two` (passivewall.hdf5 symlinked).
