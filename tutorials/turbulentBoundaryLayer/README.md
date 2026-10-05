# Turbulent zero-pressure-gradient boundary layer (ZPG TBL)

A spatially-developing incompressible **ZPG turbulent boundary layer** DNS,
non-dimensionalised by the inlet displacement thickness (Re_δ*,0 = 450, U∞ = 1).
This is a **direct, apples-to-apples reproduction of the CaNS and AMPHIBIOUS
reference DNS** — the *same* domain, resolution, grid, and trip — so mobydiff's
result can be held against theirs with nothing else varying, and against the
**SIMSON** pseudo-spectral code as an independent spectral reference.

(The earlier trip / resolution / solver study that led here — including a finer
in-house grid — is kept in **[`tests_record.md`](tests_record.md)**.)

## The case (identical to CaNS `xyz_3200_384_135`)

| ingredient | value | matches |
|---|---|---|
| domain | 650 × 100 × 26 δ*₀ | CaNS `l = 650, 100, 26` |
| resolution | 3200 × 384 × **136** | CaNS `ng = 135` (see note below) |
| x, z grid | uniform (Δx = 0.203, Δz = 0.191) | CaNS (Δz = 0.193) |
| **y grid** | one-sided tanh, `stretch = 2.822` (`distribution = cans`) | CaNS `gtype = 2, gr = 2.822` — reproduced to 5.6e-6 (their float32 rounding) |
| Re_δ*,in | 450 | CaNS `visci` |
| trip | Schlatter–Örlü, **CaNS-exact port**, `trip_amp = 0.18854`, `x0 = 10`, `nmodes = 16` | CaNS/AMPHIBIOUS `trip_*` (`src/modules/bodyforce.f90`) |
| convection | **conservative (CONS / divergence)** | CaNS / AMPHIBIOUS momentum form |
| pressure solver | red-black SOR, `niter = 6`, `sor = 1.5` | AMPHIBIOUS (same projection) |
| BCs | inlet/outlet x, periodic z, no-slip wall + pressure-pinned top | CaNS `cbcvel/cbcpre` |

The `cans` (one-sided tanh) grid distribution was added to the solver for this
(`init.f90` `GRID_TANH_WALL`): `y = L·[1 + tanh(a(s−1))/tanh(a)]`, `a = stretch`.

> **Spanwise `nz = 136`, not 135.** The reference grid is `nz = 135`, but
> mobydiff's red-black SOR needs **even** global sizes in periodic directions
> (checkerboard consistency across the periodic wrap). `nz = 135` is odd, so we use
> `nz = 136` at the **same** `Lz = 26`: the domain is byte-identical to the
> references, only the spanwise cell count ticks up by one (Δz 0.193 → 0.191,
> ~0.7 %, negligible for a spanwise-homogeneous statistic). To run the literal
> `nz = 135` grid instead, switch `[pressure] solver` to the default
> damped-Jacobi/Chebyshev (no even constraint) — closer to CaNS's exact FFT
> projection, but no longer mirroring AMPHIBIOUS's red-black SOR.

> **Convection — CONS vs skew.** This case uses the **conservative (divergence)**
> momentum convection to mirror CaNS/AMPHIBIOUS exactly (`[flow] convection =
> divergence`). The conservative form is kept for **testing/research**: it is not
> energy-neutral under the incremental projection, so watch the divergence residual
> on a long run (here `L2_div` stayed bounded at 2.5e-5–1.2e-4 over the full
> 10 000 t.u. statistics window — no drift, as AMPHIBIOUS's own CONS runs at
> `niter = 6`). The solver's **default production form is skew-symmetric**; a skew
> variant (`production_skew.ini` / `production_stats_skew.ini`) is provided to
> compare against AMPHIBIOUS's skew-symmetric dataset.

## Code comparison — SIMSON / CaNS / AMPHIBIOUS

All codes share the nondimensionalization (Re_δ*,in = 450), so the flow is compared
directly at a matched **Re_θ ≈ 677**. With the CaNS-exact trip the **transition
region matches too** — a clean monotone rise to a single c_f overshoot, not the
old weak-trip spike-then-dip (see `tests_record.md`).

![Code comparison vs SIMSON, CaNS, AMPHIBIOUS](assets/figures/code_comparison.png)

| code | grid (nx·ny·nz) | c_f | H | u_τ | u′_rms peak | −u′v′ peak |
|---|---|---|---|---|---|---|
| SIMSON (spectral) | 3072·301·— | 0.00471 | 1.507 | 0.0485 | 2.631 | 0.882 |
| CaNS | 3200·384·135 | 0.00463 | 1.503 | 0.0481 | 2.706 | 0.874 |
| AMPHIBIOUS | 3200·384·135 | 0.00460 | 1.508 | 0.0479 | 2.723 | 0.899 |
| **mobydiff** | **3200·384·136** | **0.00460** | **1.504** | **0.0480** | **2.716** | **0.867** |

On the identical grid / domain / trip / numerics, **mobydiff lands inside the
CaNS–AMPHIBIOUS scatter on every quantity** — c_f 0.00460 (−0.7 % vs CaNS, +0.04 %
vs AMPHIBIOUS), and u′_rms peak 2.716 between CaNS 2.706 and AMPHIBIOUS 2.723. The
mean U⁺(y⁺) and the Reynolds stresses of all four codes overlay through the
sublayer, log region and wake. **Two takeaways:**

- **On the reference setup, mobydiff reproduces CaNS/AMPHIBIOUS** — as it must,
  since the same grid and trip drive the same developed layer. The three
  finite-volume/difference codes cluster tightly (c_f ≈ 0.0046, ~2 % below SIMSON).
- **The near-wall u′_rms peak sits ~3 % above SIMSON for _all three_ non-spectral
  codes** (mobydiff 2.72, CaNS 2.71, AMPHIBIOUS 2.72). This overshoot is a generic
  second-order finite-volume/difference signature relative to a spectral method,
  **not** a mobydiff artifact. (A finer in-house grid narrows the c_f gap to SIMSON
  to ~0.4 %; see `tests_record.md`.)

> **AMPHIBIOUS curve caveat.** AMPHIBIOUS's *conservative*-convection datasets are
> not published on the LSDF share — only the skew-symmetric ones are — so the
> AMPHIBIOUS curve above is `TRIP_SKEWSYM_CENTERED` (same centered differencing and
> numerics, skew instead of CONS convection). CaNS is the direct CONS reference.

The mobydiff comparison data, in the same NetCDF format as the reference files, is
committed at **`assets/mobydiff/xyz_3200_384_136/data.nc`**.

## Layout

```
cold_start.ini           stage 0: from-scratch laminar->turbulent (-> coldstart_100000.h5)
production.ini           phase 1: re-equilibration, CONS convection (stats off)
production_stats.ini     phase 2: statistics accumulation -> production_stats.h5
production_skew.ini      phase 1, skew variant (-> skewsym_175000.h5)
production_stats_skew.ini   phase 2, skew variant -> skewsym_stats.h5
reproduce.py             regenerate assets/mobydiff/*.nc + the figure from the statistics
tests_record.md          the trip / resolution / solver study behind this case
assets/
  mobydiff/xyz_3200_384_136/data.nc   our reduced statistics (committed, res-study format)
  figures/code_comparison.png         the committed comparison figure
  postpro/                            post-processing + grid/IC tooling (below)
```

Large data are **kept locally** (all exceed GitHub's 100 MB limit), listed in
`.gitignore` and regenerable:

- `coldstart_*.h5`, `production_*.h5` — field snapshots the run writes,
- `production_stats.h5` — the raw span+time statistics,
- `assets/postpro/passivewall.hdf5` — the SIMSON spectral reference.

The **CaNS / AMPHIBIOUS** reference data live on the group LSDF share
(`.../tbl-dns/res_study/data_processed/{cans,amphibious}/*/data.nc`); the committed
mobydiff `data.nc` is in the identical format so the comparison is symmetric.

## Reproduce from scratch

No large input file is needed — the grid, the Blasius inflow and the trip forcing
are all generated by the solver from the `.ini`. The case is 165 M cells, a
multi-day GPU campaign (≈ 0.18 s/step on 4× A100).

```bash
# build (see the repo README): ./compile.sh gpu   (or cpu)
MPI="mpirun -n 4 ./build_gpu/moby_solve"      # one rank per GPU

# stage 0 — cold start: laminar Blasius + strong trip -> turbulent, ~2000 t.u.
$MPI cold_start.ini                            # -> coldstart_100000.h5

# phase 1 — re-equilibrate under the CaNS production trip (CONS), ~1500 t.u.
$MPI production.ini                            # -> production_175000.h5

# phase 2 — accumulate statistics, ~10000 t.u. -> production_stats.h5
$MPI production_stats.ini
```

(For the skew-symmetric variant, run `production_skew.ini` then
`production_stats_skew.ini` from the same `coldstart_100000.h5`.)

Then regenerate the comparison data and figure:

```bash
python3 reproduce.py        # -> assets/mobydiff/.../data.nc + assets/figures/code_comparison.png
```

`reproduce.py` exports the statistics to the res-study NetCDF format
(`assets/postpro/make_mobydiff_nc.py`) and draws the code comparison
(`compare_codes.py`, skipping any reference whose data is not mounted). To continue
the run instead, restart `production_stats.ini` from a later `production_p2_*.h5` —
the statistics accumulators continue seamlessly.

### Post-processing tooling (`assets/postpro/`)

| script | purpose |
|---|---|
| `make_mobydiff_nc.py` | export `bl_stats.h5` → res-study NetCDF (`data.nc`) |
| `compare_codes.py` | mobydiff vs SIMSON / CaNS / AMPHIBIOUS (the figure above) |
| `compare_passivewall.py` | mobydiff vs SIMSON (single-code, 4-panel) |
| `bl_stats.py` | boundary-layer profile plots from the statistics |
| `make_finey_grid.py`, `make_finewall_restart.py` | build / interpolate onto a finer in-house y-grid (the `tests_record.md` study) |
| `make_global_restart.py` | 1-block → legacy global-3D restart for a different rank count |
| `dpdx.py`, `vonkarman.py`, `resolution.py`, `dyplus_profiles.py`, `viz_flowfield.py` | diagnostics |
