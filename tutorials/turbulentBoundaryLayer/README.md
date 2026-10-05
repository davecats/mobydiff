# Turbulent zero-pressure-gradient boundary layer (ZPG TBL)

A spatially-developing incompressible **ZPG turbulent boundary layer** DNS,
non-dimensionalised by the inlet displacement thickness (Re_δ*,0 = 450, U∞ = 1).
The best-resolved case from an extended trip / resolution / solver study (full
history in **[`tests_record.md`](tests_record.md)**) — the fine wall-normal grid
plus a **byte-faithful port of the CaNS/AMPHIBIOUS trip** — validated against three
independent reference DNS: the **SIMSON** pseudo-spectral code, **CaNS**, and
**AMPHIBIOUS**.

## The case

| ingredient | value | why |
|---|---|---|
| grid | 4096 × 224 × 192 (176 M) | Δx⁺≈4 ≈ Δz⁺≈3.7, Δy⁺_wall≈0.23, Δy⁺_max≈4.3 (in δ₉₉) |
| x-grid | geometric, stretch 1.5 | Δx⁺ ~ uniform as u_τ falls downstream |
| y-grid | blayer (wall-clustered + freestream coarsening) | resolve the BL, keep the tall domain affordable |
| domain | 750 × 100 × 32 δ*₀ | ly = 100 δ*₀ (very tall) so the top pins a true ZPG |
| trip | Schlatter–Örlü, **CaNS-exact port**, `trip_amp = 0.18854`, `nmodes = 16` | byte-faithful port of the CaNS/AMPHIBIOUS trip (`src/modules/bodyforce.f90`) |
| convection | skew-symmetric | energy-neutral under the incremental projection |
| **pressure solver** | **red-black SOR, niter = 6, sor = 1.5** | stable at low niter on this outlet case, ~1.8× faster than Chebyshev-Jacobi niter=12 |

## Code comparison — SIMSON / CaNS / AMPHIBIOUS

All four codes share the nondimensionalization (Re_δ*,in = 450), so the flow is
compared directly at a matched **Re_θ ≈ 677**. mobydiff now runs the **CaNS-exact
trip** (same parameters as CaNS/AMPHIBIOUS), so the **transition region matches too**
— a clean monotone rise to a single c_f overshoot at x ≈ 73, not the old weak trip's
spike-then-dip (see `tests_record.md`).

![Code comparison vs SIMSON, CaNS, AMPHIBIOUS](assets/figures/code_comparison.png)

| code | grid (nx·ny·nz) | c_f | H | u′_rms peak | −u′v′ peak |
|---|---|---|---|---|---|
| SIMSON (spectral) | 3072·301·— | 0.00471 | 1.507 | 2.631 | 0.882 |
| CaNS | 3200·384·135 | 0.00463 | 1.503 | 2.706 | 0.874 |
| AMPHIBIOUS | 3200·384·135 | 0.00460 | 1.508 | 2.723 | 0.899 |
| **mobydiff** | **4096·224·192** | **0.00465** | **1.499** | **2.738** | **0.866** |

Mean U⁺(y⁺) and the Reynolds stresses of all four codes overlay through the
sublayer, log region and wake. **Two takeaways:**

- **With the CaNS-exact trip, mobydiff reproduces CaNS** — c_f 0.00465 vs CaNS
  0.00463, −u′v′ peak 0.866 vs 0.874 — as it must, since the same trip drives the
  same developed layer. The three finite-volume/difference codes cluster tightly
  (c_f ≈ 0.0046, 1–2 % below SIMSON).
- **The near-wall u′_rms peak sits 3–4 % above SIMSON for _all three_ non-spectral
  codes** (mobydiff 2.74, CaNS 2.71, AMPHIBIOUS 2.72). This overshoot is a generic
  second-order finite-volume/difference signature relative to a spectral method,
  **not** a mobydiff artifact.

The table is the 10 000-time-unit statistics of 2026-10-04, produced with the
predicted outlet face (below). The previous data set (2026-09, the outlet face
reset every substage) read c_f 0.00462, u′_rms 2.707, −u′v′ 0.867 at the same
station: the outlet treatment moves nothing there, the rest is sampling.

The mobydiff comparison data, in the same NetCDF format as the reference files, is
committed at **`assets/mobydiff/xyz_4096_224_192/data.nc`**.

## The outlet zone

Do not use the last stretch of the plate. How long that stretch is depends on
the outlet treatment, and it was measured (2026-10-02 and 2026-10-04,
`assets/postpro/momentum_integral.py`: the growth of the momentum thickness
against `c_f/2 − (H+2)(θ/U_e) dU_e/dx`, and the mean wall pressure, per x
band). Two data sets: the old outlet (the face reset to its neighbour every
substage, 500 time units after a 250-unit settling, 4 GPUs) and the shipped
one (the face predicted, the committed 10 000 time units):

| x band | old outlet: dθ/dx | balance | p_wall | shipped: dθ/dx | balance | p_wall |
|---|---|---|---|---|---|---|
| 450 … 600 | 2.0 … 2.3e-3 | 2.2 … 2.3e-3 | +2e-4 | 2.2 … 2.3e-3 | 2.2 … 2.2e-3 | −1e-4 |
| 650 … 680 | 1.8e-3 | 2.1e-3 | −4.0e-4 | 2.2e-3 | 2.1e-3 | +1.5e-4 |
| 700 … 715 | 1.0e-3 | 2.0e-3 | −1.5e-3 | 2.4e-3 | 2.1e-3 | +4.4e-4 |
| 730 … 740 | 4.0e-4 | 2.0e-3 | −2.2e-3 | 2.4e-3 | 2.0e-3 | +8.6e-4 |
| 745 … 748 | 2.2e-4 | 2.0e-3 | −2.8e-3 | 2.7e-3 | 2.0e-3 | +9.5e-4 |
| 749.5 … 750 | 6e-5 | 2.0e-3 | −3.0e-3 | 4.3e-3 | 2.2e-3 | +1.4e-4 |

- **Old outlet** (data before 2026-10-04): the wall pressure falls toward the
  outlet and the momentum thickness stops growing over the last ~100 δ*₀
  (5 δ₉₉); the freestream pressure does not move. A zero-pressure-gradient
  layer has `p_wall = p_e`: this was the outlet pressure mode, time averaged.
- **Shipped** (the outlet face predicted): the balance holds within 5 % up to
  x = 650, within 15 % to 715 and within 20 % to 740, under a wall pressure
  that rises by 1e-3 over the last 60 units (a mild adverse gradient the
  freestream does not see); the last 2 δ*₀ adjust to the uniform outlet
  pressure (the mean pressure inside a turbulent layer is `−⟨v′v′⟩`, the face
  value 0): c_f rises 14 % there, most of it in the last cell. The stored
  pressure of the last column is 100 times quieter and the divergence
  residual of the 6-iteration red-black projection 5 times lower (3e-6).
- **Use x ≤ 650** for anything quantitative (Re_θ ≤ 940); the comparison
  station Re_θ ≈ 677 is at x ≈ 376.

## Layout

```
production.ini            phase 1: re-equilibration (stats off)
production_stats.ini      phase 2: statistics accumulation -> production_stats.h5
cold_start.ini           stage 0: from-scratch laminar->turbulent (regenerates restart_field.h5)
reproduce.py             regenerate assets/mobydiff/*.nc + all figures from the statistics
tests_record.md          the full trip/resolution/solver study behind the shipped case
assets/
  mobydiff/xyz_4096_224_192/data.nc   our reduced statistics (committed, res-study format)
  figures/                            the committed comparison figures
  postpro/                            post-processing + grid/IC tooling (below)
```

Large data are **kept locally** (all exceed GitHub's 100 MB limit), listed in
`.gitignore` and regenerable — see below:

- `restart_field.h5` — a developed instantaneous field (the phase-1 IC),
- `production_stats.h5` — the raw span+time statistics,
- `assets/postpro/passivewall.hdf5` — the SIMSON spectral reference.

The **CaNS / AMPHIBIOUS** reference data live on the group LSDF share
(`.../tbl-dns/res_study/data_processed/{cans,amphibious}/*/data.nc`); the committed
mobydiff `data.nc` is in the identical format so the comparison is symmetric.

## Reproduce from scratch

No large input file is needed — the grid, the Blasius inflow and the trip forcing
are all generated by the solver from the `.ini`.

```bash
# build (see the repo README): ./compile.sh gpu   (or cpu)
MPI="mpirun -n 4 ./build_gpu/moby_solve"      # one rank per GPU
PREP="mpirun -n 4 ./build_cpu/moby_prepare"   # the case file of each ini, once, SAME rank count

# stage 0 — cold start: laminar Blasius -> turbulent, ~2000 t.u. (~100k steps)
$PREP cold_start.ini; $MPI cold_start.ini
mv coldstart_100000.h5 restart_field.h5        # the developed field

# phase 1 — re-equilibrate on the shipped (grid, trip, solver), ~1000 t.u.
$PREP production.ini; $MPI production.ini                            # -> production_100000.h5

# phase 2 — accumulate statistics, 10000 t.u. (500000 steps, ~27 h on 4 A100) -> production_stats.h5
$PREP production_stats.ini; $MPI production_stats.ini
```

Then regenerate the comparison data and figures:

```bash
python3 reproduce.py        # -> assets/mobydiff/.../data.nc + assets/figures/*.png
```

`reproduce.py` exports the statistics to the res-study NetCDF format
(`assets/postpro/make_mobydiff_nc.py`), draws the code comparison
(`compare_codes.py`, falling back to SIMSON-only if the LSDF reference data are not
mounted) and the SIMSON single-code figures. To continue the run instead, restart
`production_stats.ini` from a later `production_p2_*.h5` — the statistics
accumulators continue seamlessly.

### Post-processing tooling (`assets/postpro/`)

| script | purpose |
|---|---|
| `make_mobydiff_nc.py` | export `bl_stats.h5` → res-study NetCDF (`data.nc`) |
| `momentum_integral.py` | dθ/dx against the von Kármán balance and the wall pressure, per x band (the outlet zone) |
| `compare_codes.py` | mobydiff vs SIMSON / CaNS / AMPHIBIOUS (the figure above) |
| `compare_passivewall.py` | mobydiff vs SIMSON (single-code, 4-panel) |
| `bl_stats.py` | boundary-layer profile plots from the statistics |
| `make_finey_grid.py` | build the blayer wall-normal node line (ny=224) |
| `make_finewall_restart.py` | interpolate a field onto a finer y-grid |
| `dpdx.py`, `vonkarman.py`, `resolution.py`, `dyplus_profiles.py`, `viz_flowfield.py` | diagnostics |
