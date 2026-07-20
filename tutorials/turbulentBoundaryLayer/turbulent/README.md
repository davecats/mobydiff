# turbulent — ZPG turbulent boundary layer DNS (Skote 2001, trip by Schlatter & Örlü)

Zero-pressure-gradient turbulent boundary layer: laminar Blasius inflow
tripped to turbulence by a random wall-normal volume force, then a
spatially developing turbulent BL. Setup follows the ZPG case of

  M. Skote (2001), *Studies of turbulent boundary layer flow through direct
  numerical simulation*, PhD thesis, KTH (Table 2.3),

with the trip forcing of Schlatter & Örlü (2012) implemented in the solver.

## Nondimensionalization and domain

Lengths by the **inflow displacement thickness** δ*₀ = 1, velocity by
U∞ = 1, so `[flow] re` = Re_δ*,₀ = **450** (⇒ Re_θ,₀ = 450/H = 174,
H = 2.5911, θ₀ = 1/H = 0.38594 — the `blasius_theta` used by the inlet).

| | Skote ZPG | here |
|---|---|---|
| reference length | δ*₀ | δ*₀ |
| Re_δ*,₀ | 450 | 450 |
| L_x × L_y × L_z | 600 × 30 × 34 (incl. fringe) | 500 × **100** × 32 |
| streamwise end | fringe (periodic) | physical **outlet** |
| top | — | validated ZPG **outlet** at y = 100 |

The top is raised to 100 δ*₀ (vs Skote's 30) since the free wall-normal
`natural` stretching keeps the extra height cheap; it uses the
`../topbc_outlet` boundary (Dirichlet-p pins ZPG, the best-validated top).

## Boundary conditions

- `x_min`: Blasius **inlet** (u + entrainment v), `blasius_theta = 0.38594`.
- `x_max`: **outlet** (zero-gradient velocity, Dirichlet p).
- `y_min`: no-slip **wall** (the plate); `y_max`: **outlet** (top).
- `z`: periodic (spanwise homogeneous).

## Resolution

Standard second-order-FD DNS at the turbulent end (Re_θ ≈ 700,
l⁺ = ν/u_τ ≈ 0.046 δ*): Δx⁺ ≈ 8.5, Δz⁺ ≈ 3.6, Δy⁺_wall ≈ 0.3 growing to
Δy_top ≈ 0.54 δ* by one-sided natural stretching. Production grid
(`zpg.ini`): 1280 × 256 × 192 = 63 M cells.

## The trip forcing (Schlatter & Örlü 2012)

Implemented as `[force] type = trip` (`src/modules/bodyforce.f90`). A random
volume force on the **wall-normal (v)** momentum, Gaussian-localized near the
inlet, random in span and time:

    f_v(x,y,z,t) = amp · exp[ -((x-x₀)/l_x)² - (y/l_y)² ] · g(z,t)
    g(z,t) = (1-b(t)) g_k(z) + b(t) g_{k+1}(z),   b = 3p² - 2p³,  p = t/t_s - k

- `g_k(z)` is a unit-rms random spanwise function (`trip_nmodes` Fourier
  modes, period L_z), redrawn every `trip_ts`; `b(t)` is the C¹ smooth step
  that makes the forcing continuous in time.
- The RNG is seeded deterministically (`trip_seed`) so a run is reproducible
  and rank-independent (the trip is a global spanwise function evaluated
  identically on every rank).

Config keys (`[force]`): `trip_x0`, `trip_lx`, `trip_ly`, `trip_amp`,
`trip_ts`, `trip_nmodes`, `trip_seed`. The force reuses the bodyforce device
array + correction kernel, so it is added exactly like any `f·μ` momentum
source (and is a no-op when `enabled = false`).

## Run and verify

```bash
# smoke (small + coarse, runs in-session): validates the trip + transition
mpirun -n 4 ../../../build_cpu/moby_solve template.ini   # mint template_1.h5
python3 ../make_blasius_ic.py --theta 0.38594            # -> IC_blasius.h5
mpirun -n 4 ../../../build_cpu/moby_solve smoke.ini
python3 check_turbulent.py smoke_4000.h5 --plot smoke.png

# production DNS (zpg.ini): same recipe at full resolution -- a long run
```

`check_turbulent.py` reports Re_θ(x), the shape factor H, the skin friction
c_f(x) vs the turbulent correlation 0.024 Re_θ^(-1/4), the spanwise-
fluctuation marker w_rms(x) (transition location), and the mean U⁺(y⁺)
profile against the log law. NOTE: a single snapshot is only z-averaged;
converged statistics need time-averaging over many snapshots (production).

## Smoke validation result

`smoke.ini` (256×96×64, `trip_amp = 0.15`) run from the Blasius IC. The trip
triggers transition and a turbulent boundary layer establishes and spreads
downstream at ≈U∞ (snapshots, z-averaged):

| t (step) | turbulent region (w_rms ≳ 3%) |
|----------|-------------------------------|
| 45 (1500)  | x ≈ 15–40 (trip disturbance holding, not decaying) |
| 90 (3000)  | x ≈ 20–70 |
| 135 (4500) | x ≈ 20–120, front near the outlet |

At t = 135 (`smoke.png`, `check_turbulent.py smoke_4500.h5`):

- **shape factor H drops from the laminar 2.59 to ≈1.5** in the transitioned
  region (x = 47–79) — the turbulent signature;
- **c_f tracks the turbulent correlation** 0.024 Re_θ^(−1/4): e.g. at
  x = 79, c_f = 5.6e-3 vs 5.5e-3;
- w_rms shows the trip spike (x ≈ 15) then a sustained ≈3.5% turbulent
  plateau (x = 20–105), decaying to 0 where the front has not yet arrived;
- the mean U⁺(y⁺) at x = 64 has the viscous sublayer (U⁺ = y⁺) turning
  toward the log region.

The run stays stable (u ∈ [−0.07, 1.36]) throughout. The mild under-log and
the U⁺ ≈ 17 plateau are the expected low-Re (Re_θ ≈ 360) + single-snapshot
(z-average only, no time averaging) behavior — they tighten with time/
spanwise averaging over many snapshots at the production resolution.

**This validates the setup end-to-end** (Blasius inlet, trip, transition,
turbulent BL, ZPG outlet top). Converged turbulence statistics — a smooth
c_f(Re_θ), a mean profile on the log law, second-order stats — require the
production `zpg.ini` run over several flow-throughs with time+spanwise
averaging.
