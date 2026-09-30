# CaNS / AMPHIBIOUS exact-match ZPG-TBL

A direct, apples-to-apples reproduction of the CaNS and AMPHIBIOUS reference
boundary-layer cases: **same domain, resolution, grid, and trip**, so mobydiff's
result can be compared to theirs with nothing else varying.

## The setup (identical to CaNS `xyz_3200_384_135`)

| ingredient | value | matches |
|---|---|---|
| domain | 650 × 100 × 26 δ*₀ | CaNS `l = 650, 100, 26` |
| resolution | 3200 × 384 × **136** | CaNS `ng = 135` (see note) |
| x, z grid | uniform (Δx=0.203, Δz=0.191) | CaNS Δz=0.193 |
| **y grid** | **one-sided tanh, `stretch = 2.822`** (`distribution = cans`) | CaNS `gtype = 2, gr = 2.822` — reproduced to 5.6e-6 (their float32 rounding) |
| Re_δ*,in | 450 | CaNS `visci` |
| trip | CaNS-exact port, `trip_amp = 0.18854`, `x0 = 10`, `nmodes = 16` | CaNS/AMPHIBIOUS `trip_*` |
| BCs | inlet/outlet x, periodic z, no-slip wall + pressure-pinned top | CaNS `cbcvel/cbcpre` |

The `cans` (one-sided tanh) grid distribution was added to the solver for this
(`init.f90` `GRID_TANH_WALL`): `y = L·[1 + tanh(a(s−1))/tanh(a)]`, `a = stretch`.

> **Spanwise `nz = 136`, not 135.** The reference grid is `nz = 135`, but
> mobydiff's red-black SOR needs **even** global sizes in periodic directions
> (checkerboard consistency across the periodic wrap). `nz = 135` is odd, so we
> use `nz = 136` at the **same** `Lz = 26`: the domain is byte-identical to the
> references, and only the spanwise cell count ticks up by one (Δz 0.193 →
> 0.191, ~0.7%). For a spanwise-homogeneous statistic this is negligible. To run
> the literal `nz = 135` grid instead, switch `[pressure] solver` to the default
> damped-Jacobi/Chebyshev (no even constraint) — that path is closer to CaNS's
> exact FFT projection but no longer mirrors AMPHIBIOUS's red-black SOR.

## Convection — two variants

- **`production_skew.ini` / `_stats_skew.ini`** — skew-symmetric convection (the
  mobydiff production form). Compare against **AMPHIBIOUS_TRIP_SKEWSYM_CENTERED**,
  which uses the same numerics. This is the stable, recommended baseline.
- **`production_cons.ini` / `_stats_cons.ini`** — divergence (conservative)
  convection. Compare against **AMPHIBIOUS_TRIP_CONS_CENTERED** / CaNS. The
  conservative form was re-added for testing/research (`[flow] convection =
  divergence`); it is **not** energy-neutral under the incremental projection, so
  watch the divergence residual on the long run — it may drift.

## Run (cold start → phase 1 → phase 2)

```bash
MPI="mpirun -n 4 ./build_gpu/moby_solve"

# stage 0 — cold start (skew, strong trip): laminar Blasius -> turbulent, ~2000 t.u.
$MPI cold_start.ini            # -> coldstart_100000.h5

# phase 1 — re-equilibrate under the CaNS production trip (~1500 t.u.)
$MPI production_cons.ini       # (or production_skew.ini)

# phase 2 — statistics, ~10000 t.u.
$MPI production_stats_cons.ini # (or production_stats_skew.ini)
```

## Compare

```bash
python3 ../assets/postpro/make_mobydiff_nc.py cons_stats.h5 cons_p2_<last>.h5 cons.nc
python3 ../assets/postpro/compare_codes.py --mobydiff cons.nc \
        --amphibious .../amphibious/xyz_3200_384_135_TRIP_CONS_CENTERED/data.nc
```

165 M cells ≈ the shipped case, so a similar multi-day GPU campaign.
