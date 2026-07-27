# turbulent_finewall — Δy⁺_max=4 variant of turbulent_two

Identical to `../turbulent_two` **except the wall-normal grid is refined** so the
largest cell within the boundary layer is **Δy⁺_max ≈ 4** (ny 176 → **256**;
measured Δy⁺_max 6.2 → 3.9, Δy⁺_wall 0.22 → 0.12). x and z grids unchanged
(2048×192, geometric x). 2048×256×192 = 100.7 M cells. Same trip (0.15), skew
convection, niter=12.

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
