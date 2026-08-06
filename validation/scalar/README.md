# Passive scalars — S0/S1/S2/S3/S4/S5a gates

Gates for increments **S0** (layout, config, io, halos, BCs — the scalar is
carried but not advanced), **S1** (the transport kernel, molecular diffusivity
only), **S2** (the turbulent closure `D_face = 1/(Re Pr) + nut/Pr_t`,
constant `Pr_t` or Kays–Crawford), **S3** (the immersed body: cell-centred
coefficients, both wall modes, `coef_p_blocks`), **S4** (in-solver
statistics + tooling: `scalar_stats.f90`, `tools/scalar_stats.py`,
`compare_fields.py` dataset discovery) and **S5a** (the THERMAL WALL
FUNCTION: Kader/Jayatilleke under `[rans] wall_treatment = wall_function`,
which S2 rejected as a hard config error) of `docs/next_session_scalar.md`.

```bash
./compile_nofma.sh cpu && ./compile_nofma.sh gpu   # bit-exactness builds
./run_gates.sh    [uniform|conserve|conduction|wave|pr|det|restart]   # S1
./run_gates_s2.sh [kays|wferr|sst|les|band|det]                       # S2
./run_gates_s3.sh [solid|conserve|prep|missing|refine|det|cyl]        # S3
./run_gates_s4.sh [stats|accum|plane|levels|restart|det|noeffect|heat|adia|cyl|tools]
./run_gates_s5.sh [unit|ref|sweep|det]                                # S5a
REF=<pre-change nofma binary> MODE=cpu RANKS=4 ./run_bitexact.sh      # count = 0
REF=<S1 nofma binary>         MODE=cpu          ./run_bitexact_s1.sh  # turbulence off
REF=<S2 nofma binary>         MODE=cpu          ./run_bitexact_s3.sh  # no body
```

`run_gates.sh` takes `BIN` (default `../../build_cpu/moby_solve`), `GBIN`
(GPU binary, used by the `det` group) and `RANKS`. Helpers: `scalar_tools.py`
(per-level block geometry: cell centres, widths, volumes — the same midpoint
subdivision `blocks.f90` uses), `make_scalar_ic.py` (manufactured
`s = A sin(k·x)` initial condition written into a copy of a snapshot, since
the solver has no analytic non-uniform scalar IC — this also exercises the S0
restart path), `check_scalar.py` (all the checkers).

## Gate cases

| ini | gate |
|---|---|
| `uniform3.ini` | (a) uniform scalar through a **3-level** refined block layout — the cylinder-shaped `validation/multilevel_body` leaf table on its ZERO-FORCE twin, so every l0–l1 and l1–l2 interface, edge and corner is crossed while the body exerts no force. Two scalars, Pr 0.71 and 10, plus inlet (Dirichlet) and outlet (zero-gradient) physical faces. |
| `conserve.ini` | (b) global conservation: fully periodic Beltrami box, manufactured `sin(k·x)` scalar, single level (the flux form telescopes exactly only across matching faces). |
| `conduction.ini` | (c) pure conduction between Dirichlet walls on a **stretched** (tanh) grid, no flow. `cond_lin`: the linear profile is an exact fixed point of the discrete operator. `cond_16/32/64`: constant volumetric source ⇒ parabolic steady state ⇒ spatial order. |
| `wave.ini` | (d) advection–diffusion MMS: uniform oblique velocity in a periodic cube (an exact NS solution, so the advecting field stays constant to round-off) carrying `s = sin(k·x)`; exact solution `exp(-D\|k\|²t) sin(k·(x − ut))`. `cflmax = pecletmax = 0` freezes `dt`, so the (negligible) time error is identical on all three grids. |
| `prsweep.ini` | (e) Pr sweep 0.1/1/10 in a laminar Poiseuille channel: the scalar depends on y only and `u` does not vary in x, so the discrete convection vanishes identically and the transient is pure conduction at `D = 1/(Re·Pr)` **while a real flow is advanced** — it isolates the Pr scaling and checks the transport kernel does not corrupt it. Compared against the analytic series and the wall (Nusselt) flux. |
| `det.ini` | (f) determinism on a refined transport case (`conserve.ini` + a 2:1 y slab): 1 rank == 4 ranks, CPU == GPU; also reports the (non-zero) 2:1-interface conservation residual. |
| `smoke.ini` | (g) S0 restart round-trip: with the scalar dataset present the file must win over `[scalar.N] initial`; with it absent (renamed scalar) the solver must warn and reinitialise. |

## S2 gate cases

| ini | gate |
|---|---|
| (unit test) | (m) **Kays–Crawford correlation**: `src/test_scalar.f90` (CMake target `scalar_test`) checks `prt_kays` host-side against an independent mpmath transcription (50–1000 digits, so the large-`Pe_t` cancellation is resolved exactly) of the same formula. Every branch: the `Pe_t = 0` guard, the direct expression, the small-`x` series, both sides of the `x = 1/2` crossover plus its continuity, the `Pe_t → 0` limit (`Pr_t = 2 Prt_inf`) and the `Pe_t → ∞` limit at `Pe_t = 1e300`, where the direct expression overflows `a²` and cancels to nothing. Plus a monotonicity sweep over six decades. |
| `turbles.ini` | (h) **developed turbulent channel under LES**, `Re_tau` 180, Pr 0.71. Geometry, grid, WALE settings and velocity restart are `../channel_interface/les` (`uniform.ini`) — the coarse 64×48×64 grid where the SGS term is genuinely active. Scalar: antisymmetric isothermal walls (+1 / −1), no source, so in statistical steady state the TOTAL wall-normal flux `J = <v'θ'> − (D + nut/Pr_t) d<θ>/dy` is exactly CONSTANT. `θ_tau = J` (u_tau = 1) and `θ+ = (θ_w − θ)/θ_tau` are then measured against Kader's correlation, and the constancy of `J(y)` is a second, internal gate on both convergence and the transport operator. |
| `turbslab.ini` | (j) **2:1 wall-band-refined channel**, the same scalar, the `../channel_interface/les` `slab.ini` geometry (symmetric wall bands refined 2:1, flat y-interfaces at y+ ≈ 88). Band metric = the `tools/patch_interface_stats.py` method: the `θ'_rms` ratio of the refined run to the matched unrefined control on the COARSE rows adjacent to the interface, against the same ratio in the core. A spurious band is a localized excess — the u'/v' lesson of `../channel_interface/README.md`. |
| `turbsst.ini` | (i) **steady k-ω SST, resolved walls** — `../rans_sst/turb180.ini` with two scalars. Steady RANS collapses to a 1D fixed point with NO resolved fluctuations, so the scalar equation reduces to one ODE whose only coefficient is the `nut` the solver wrote: integrating `dθ/dy = −J/(1/(Re Pr) + nut/Pr_t(y))` from the snapshot's own `nut` and comparing with its own `theta` tests the S2 face diffusivity and the `Pr_t` model **with the turbulence model divided out** — the sharpest gate in the set, and free of time averaging. `theta` uses the constant `Pr_t`, `theta_kc` the SAME scalar with `prt_model = kays`, so the two differ ONLY by the correlation — the in-kernel exercise of `prt_kays`. (The log-layer slope is reported for context only: the scalar carries a CONSTANT flux while the channel's momentum flux falls linearly, so the reference is `Pr_t/(kappa(1 − y/h))`, which is `Pr_t/kappa` = 2.073 only as `y/h → 0`.) |
| `wferr.ini` | (k) `[scalar]` + `[rans] wall_treatment = wall_function` was a **hard config error** in S2 (a thermal wall function was deferred to S5). `../rans_sst/wf180_y30.ini` with one scalar bolted on. **S5a implemented it and FLIPPED this gate**: the same ini must now start and report the thermal wall function's constants — see the S5a table below. |
| `detles.ini` | (l) **determinism on a scalar + turbulence case**: `turbles.ini` for 20 steps without the `[mpi] dims` pin — 1 rank == 4 ranks and CPU == GPU, at tolerance 0, on nofma builds. |

`make_theta_ic.py` seeds the developed velocity restarts with the
Reynolds-analogy `θ(y)` built from the snapshot's OWN mean velocity (see its
header for why that is not circular with the Kader gate); `check_scalar_turb.py`
holds the three analyses (`channel`, `band`, `rans`).

## S3 gate cases (the immersed body)

| ini | gate |
|---|---|
| `ibmwavy.ini` | (n) + (o) the analytic wavy bottom wall with ONE SCALAR IN EACH WALL MODE. `theta` is `dirichlet` with `ibm_value = 1`: inside the body the penalization coefficient is `SOLID/Re`, so `mu_s` underflows and every solid cell must read **exactly** 1 — an equality, not a tolerance. `phi` is `adiabatic`, seeded through the restart path with a manufactured `2 + sin(k·x)` (a NON-ZERO mean, so `∫φ dV` is a real quantity): with x,z periodic, Neumann-0 y walls and every solid face masked there is no sink at all, so the integral must be conserved to round-off. Conservation alone would also hold if the mask never fired, so the checker adds the POSITIVE half — a cell whose six staggered faces are all inside the body is sealed on every side and must be **frozen exactly**. |
| `ibmwavy.ini` (prepared) | (p) `moby_prepare` on the SAME ini writes `coef_p_blocks` **and changes nothing else**: `h5same.py --ignore coef_p_blocks` against the case file prepared from the `[scalar]`-stripped twin. The tiles are then checked against an independent Python transcription of the graded sharp-interface formula `Σ((d0−d)/d)/d0²/Re` (including the ibm.f90 bisection), the prepare is repeated on 4 ranks, and the solve FROM the file is compared with the inline analytic solve at tolerance 0 — `theta`/`phi` included, which is what makes the file round-trip a bit-exactness statement rather than a plausibility one. (p2) the same case file WITHOUT `coef_p_blocks` must be a hard error naming the fix. |
| `ibmwavyr.ini` | (q) the same body under `refine_body` (2 levels): the cell-centred coefficients must be evaluated at EACH LEAF'S OWN LEVEL. Checked level by level against the same independent transcription, plus solve-from-file == inline analytic at tolerance 0. Conservation is deliberately NOT gated here — across a 2:1 interface the flux form telescopes only approximately (see the S1 note at the end). |
| `ibmwavy.ini` (nofma) | (r) determinism on a scalar + IBM case: 1 rank == 4 ranks and CPU == GPU at tolerance 0, through the ANALYTIC path and again through the FILE path (`coef_p_blocks` read on both sides). |
| `cylheat.ini` | (s) **heated cylinder, Re 40, Pr 0.71** — the quantitative body gate. Geometry, grid and far field are `../cylinder/cyl_re40.ini` verbatim (D = 1 at (6.0, 8.02) in a 16 × 16 × 0.25 box, uniform h = 1/32), restarted from that campaign's converged steady field, so only the thermal problem is new: an isothermal body (`theta = 1`) in a `theta = 0` freestream. The Nusselt number is measured TWICE from the same snapshot — from the penalization integral `Q = ∫coef_p (1 − θ)/Pr dV`, and from the Gauss/CV border flux of `(u θ − D_th ∇θ)·n` through a box in the fluid (the A2 cross-check method). The two share nothing but the `theta` field. Unlike the momentum cross-check this one never touches the pressure, so the `niter = 6` pn-drift caveat of `../cylinder/README.md` does not apply. |

## S4 gate cases (statistics and tooling)

| ini | gate |
|---|---|
| `s4stats.ini` | (t) (t2) (t3) (t4) (t5) (t6) the IN-SOLVER statistics. The case is `turbles.ini` (the S2 developed LES channel, 64×48×64 WALE at Re_τ 180, so `nut` is genuinely active) restarted from one of that campaign's mid-run snapshots, with TWO scalars differing ONLY in the `Pr_t` model (`theta` constant, `theta_kc` Kays–Crawford): one run gates the multi-scalar column layout AND the correlation branch of the statistics kernel's face diffusivity. THE GATE: the solver samples the END-OF-STEP field and the snapshot IS the end-of-step field, so `check_scalar_stats.py`, recomputing the same seven columns from the snapshots the same run wrote, must reproduce the solver's rows to ROUND-OFF — there is no statistical tolerance to hide in. Run with one sample (t), with four accumulated samples (t2), in the `plane` layout (t3), across a restart (t4, the accumulators must continue from the file), on 1 vs 4 ranks and CPU vs GPU (t5, a tolerance — atomics and the allreduce reorder the sums), and against a statistics-OFF twin whose fields must be bit-identical (t6). |
| `turbslab.ini` + statistics | (t3b) the 2:1 wall-band-refined case: the row tables are PER LEVEL (the `channel_stats` `lvlOff` layout), so the run writes one file per level. `check_scalar_stats.py rows` accumulates straight from the leaves — no global box — and therefore gates the level split that `profile` cannot reach; the face-flux columns stay gated at single level, where the halo values a snapshot does not carry are not needed. |
| `ibmwavy.ini` + `heat_interval` | (u) (u2) the immersed body's HEAT RELEASE. The runtime samples are compared with `check_scalar_ibm.py surface` — the cancellation-free pair (staircase interface flux + graded-cell penalization) that S3 validated against the full energy budget — on the very snapshots the same run wrote, and then the physics is closed with the SOLVER's own `Q`: with no boundary flux anywhere, `½[Q(t₁)+Q(t₂)] = d/dt ∫θ dV`. (u2) is the positive control on the other wall mode: the seeded `phi` (`2 + sin(k·x)`, manifestly non-zero) is `adiabatic`, so its heat columns must be EXACTLY zero — no penalization is applied and every body face is masked, and any non-zero number would be a flux the solver never applied. |
| `cylheat.ini` + `heat_interval` | (u3) the heated cylinder's RUNTIME Nusselt number: 20 steps from the S3 campaign's converged `t = 120` field, so the number is directly comparable with the S3 post-processed 3.3655 and with Churchill–Bernstein's 3.35. |
| (tooling) | (v) `compare_fields.py` with no dataset arguments discovers the datasets present in BOTH files (canonical `un vn wn pn` first, then the scalars / `nut` / the RANS variables alphabetically), and `tools/scalar_stats.py` reads the files this session wrote. |

## S5a gate cases (the thermal wall function)

| ini | gate |
|---|---|
| `wfsst.ini` | (w0) the RESOLVED REFERENCE: `turbsst.ini` (the S2/S4 resolved-wall SST scalar channel) restarted from its converged `t = 200` field with the S4 statistics on. `theta_tau` and the wall flux for BOTH sides of the comparison therefore come out of the same instrument — the solver's own statistics — instead of a hand-rolled post-processing, which is what S4 was built for. |
| `wfs180_y05/y15/y30/y45.ini` | (w) (x) the wall-function sweep: `../rans_sst/wf180_*.ini` (the T3 grids, `y+_1 = 5/15/30/45`, `u_tau = 1` by construction) with the `turbsst` scalar pair bolted on. Every run is step-bound, and the statistics are switched on with `interval = nsteps`, so exactly ONE sample is taken, at the last step, of the very field the snapshot carries — which makes the wall gate an IDENTITY rather than a comparison of two time levels. |
| (unit test) | (m2) the two correlations host-side (`src/test_scalar.f90`, the S2 `kays` group's driver): Jayatilleke's `P(Pr/Pr_t)`, the thermal sublayer thickness `y+_T` (a ROOT — so the solver's bisection is checked against mpmath's) and `wall_diffusivity` itself, including that it is exactly 0 on the conduction branch and continuous across the switch. |
| `wferr.ini` | (k, FLIPPED) the S2 gate was that `[scalar]` + `wall_treatment = wall_function` **error-stops**. S5a implements it, so the same ini is now the positive control: the solver must START and report the per-scalar `P` / `y+_T`. |
| `wfs180_y30.ini` (nofma) | (z) determinism: 1 rank == 4 ranks at tolerance 0 and CPU == GPU on the wall-function scalar path. |

## The thermal wall function at an IMMERSED wall (S5a's open item)

| ini | gate |
|---|---|
| `ibmwf180.ini`, `ibmwf1000.ini` | (y) the S5a closure where S5a never ran it: the wall is an immersed body, so the classified wall cell is a CUT cell and the Dirichlet penalization pins the same cell the wall function installs a diffusivity on. Geometry = the `../channel_interface/les_ibm/wall_{lo,hi}.stl` slabs (walls mid-cell at `y = 0.259375 / 2.259375`, fluid gap EXACTLY 2.0); `ly = 2.5` over `ny = 8`. The pair differs ONLY in `re` (180 / 1000), which is the `y+` lever at an immersed wall: the converged wall cells sit at `y+_k = 5.8` (conduction branch) and **39.7** (log branch) respectively, against the thermal switch `y+_T = 12.178` and the momentum switch `y+_lam = 11.530`. T3's only IBM wall-function case (`../rans_sst/ibm180wf.ini`) sits at `y+ ~ 2-3`, so the log branch at an immersed wall was never exercised before. `nx = nz = 8` (RANS here is 1-D; the STL slabs are padded well past the domain, so a small x/z box is covered by the same geometry). |
| `ibmwf180.ini` (`ibm_value` 0 vs 1) | (y2) the SAME physics run with the two body-value conventions — `theta` differs by a constant, so the staircase flux is bit-identical — as the decisive test of the S4 heat diagnostic's `ibm_value` invariance. This is what turned up the double count fixed 2026-08-05 (below). |

ONE body value serves both walls (`ibm_value` is per scalar, not per
surface), so S5a's antisymmetric +1/−1 walls are unavailable: the case
drives the scalar with isothermal walls plus a constant volumetric
`source` instead. That makes the steady budget CLOSED-FORM, which is the
gate's real instrument — summed over all interior cells the convective and
diffusive fluxes telescope to the domain boundary (periodic x/z, zero-flux
y), so at steady state the heat crossing into the body is `source *
V_fluid` exactly, with no reference run needed.

`check_scalar_ibmwf.py` holds the two analyses (`wall`, `budget`) with its
own transcription of the correlations and of the budget. The case file must
carry `coef_p_blocks`, so it is built from the ini itself:
`mpirun -n 1 ../../build_cpu/moby_prepare ibmwf180.ini ibmwf180_case.h5`.

`check_scalar_wf.py` holds the two analyses (`wall`, `compare`) and its own,
INDEPENDENT transcription of the wall-function correlations. Note that
`check_scalar_stats.py` (S4) is deliberately NOT used on a wall-function
case: its face diffusivity is the resolved `nut/Pr_t`, so it cannot
reproduce a wall row there — the `wall` analysis covers that instead, and
more sharply.

`check_scalar_stats.py` holds the five analyses (`profile`, `plane`, `rows`,
`heat`, `diff`); `tools/scalar_stats.py` is the production reader (mean profile, rms,
turbulent flux, wall flux → `theta_tau` / Nusselt; the body heat → Nusselt).

`check_scalar_ibm.py` holds the five analyses (`solid`, `conserve`, `coefp`,
`nusselt`, `cv`); `run_bitexact_s3.sh` is the "scalar WITHOUT a body is
bit-exact vs the S2 binaries" gate and doubles as the cheap, sharp form of
"every S1/S2 gate still reads the same number" (see its header).

`run_bitexact.sh` is the **`[scalar] count = 0` bit-exactness** gate: the
standard 7-case suite (min_channel, les_ibm ± refine_body, Beltrami y-slab,
turb180, wf180_y30, lam30t) run with a pre-change nofma binary and the new
one, compared at tolerance 0 on every field dataset (`un vn wn pn` plus
`nut`, `k`, `omega`, `gamma`, `rethetat` where the case has them).

## Results — S0/S1 (2026-08-03, branch `scalar`) — ALL PASS

Local runs (RTX 3060 GPU + nvhpc 25.9 CPU builds), branch `scalar`.

**(a) uniform scalar, 3-level refined layout** (`uniform3.ini`, 50 steps,
928 leaves = 224 l0 + 192 l1 + 512 l2, CPU 1 rank):

```
max|theta - 2.5|  = 0.000e+00        (Pr 0.71)
max|phi  + 1.25|  = 0.000e+00        (Pr 10)
max|un - u0| = max|vn - v0| = max|wn - w0| = 0.000e+00,  pn spread 0.000e+00
```

Every l0–l1 and l1–l2 interface, edge and corner is crossed, the physical
faces are a Dirichlet inlet and a zero-gradient outlet, and the transport
kernel runs every substage — the constants survive bit-for-bit.

**(b) global conservation** (`conserve.ini`, 32³ periodic Beltrami box,
manufactured `sin(k·x)`, 200 steps, t 0.002 → 0.402):

```
int s dV: -2.2759572004815709e-15 -> -2.4424906541753444e-15
relative drift (per max|s|·V, V = 248.05) = -6.838e-19        PASS
```

**(c) pure conduction on a stretched (tanh) grid, no flow**

- `cond_lin` — the LINEAR profile is an exact fixed point: after 324 steps
  `L2 = 1.96e-17`, `Linf = 1.11e-16` (GPU; the CPU build gives exactly
  `0.0`). Relaxing to it from `s = 0` instead reaches `2.4e-5` at `t = 4`,
  which is precisely the analytic slowest-mode residual
  `exp(-D π² t/L²) = 5e-5` — convergence, not discretisation.
- `cond_16/32/64` — constant source ⇒ parabolic steady state, `t = 12`
  (analytic transient residual 1.4e-13):

| ny | steps | L2 | order |
|---|---|---|---|
| 16 | 7756 | 2.814767e-03 | |
| 32 | 35154 | 7.054808e-04 | **2.00** |
| 64 | 149695 | 1.764823e-04 | **2.00** |

(`Q = 1`, `D = 1`, walls at 0, so the parabola peaks at 0.5 and the L2 above
is absolute. Error ratios 3.990 / 3.998. The step counts grow as `ny²` —
`dt` is Peclet-bound — so the `ny = 64` leg takes ~19 min on the GPU.)

**(d) advection–diffusion MMS** (`wave.ini`, uniform oblique flow
`u = (1, 0.5, 0.25)`, `s = sin(k·x)` with `k = 2π(1,1,1)/L`, `D = 0.01`,
fixed `dt = 1e-3`, `t = 0.5`, GPU):

| n | L2 | order |
|---|---|---|
| 16 | 1.554700e-02 | |
| 32 | 3.909091e-03 | **1.99** |
| 64 | 9.786703e-04 | **2.00** |

The advecting field is exact throughout: `max|un − 1| = max|vn − 0.5| =
max|wn − 0.25| = 0.0` and `pn` spread `0.0` after 500 steps, so the measured
error is the scalar discretisation alone. (The checker takes the wave's time
origin from the IC file — the seed run leaves `t_current = 1e-3` there, and
ignoring that offset puts a ~1.2e-3 phase-error floor under the fine grids.)

**(e) Pr sweep in a laminar channel** (`prsweep.ini`, Re 10, ny 64, t = 2,
against the analytic conduction series):

| Pr | D | steps | max\|s − series\| | centre s (analytic) | wall flux ratio |
|---|---|---|---|---|---|
| 0.1 | 1 | 5120 | 1.00e-05 | 0.990836 (0.990846) | 1.0011 |
| 1 | 0.1 | 512 | 1.41e-04 | 0.227799 (0.227909) | 1.0004 |
| 10 | 0.01 | 200 | 1.58e-03 | 0.000002 (0.000001) | 1.0041 |

The error grows with Pr exactly as expected: at Pr 10 the thermal layer at
`t = 2` is `√(Dt) ≈ 0.14`, i.e. ~4 cells. The step counts show the Peclet
limiter picking up the `1/Pr` scaling of the scalar diffusivity.

**(f) determinism** (`det.ini` = the Beltrami box with a 2:1 y slab,
manufactured scalar, 20 steps, **nofma builds**):

```
1 rank vs 4 ranks:  un vn wn pn s1   max_abs = 0.0   (EXACT)
CPU vs GPU:         un vn wn pn s1   max_abs = 0.0   (EXACT)
```

(With the default FMA-contracted builds the same comparison reads ~1e-14 on
the velocities and 7.8e-16 on the scalar — host/device contraction, not the
scalar path.)

Informative, NOT a gate: the same run's 2:1-interface conservation residual
is `1.05e-05` of `max|s|·V` over 20 steps (vs `6.8e-19` for the single-level
box) — the flux form telescopes exactly only across matching faces.

**(g) S0 restart round-trip** (`smoke.ini`):

```
restart WITH the dataset:    max|s - file| = 0.000e+00  (the ini's initial = 99 is ignored)
restart WITHOUT the dataset: value = 99.0 + "warning: restart file has no dataset for scalar 1"
```

**`[scalar] count = 0` bit-exactness** (`run_bitexact.sh`, nofma builds,
7-case suite, every dataset at tolerance 0):

```
CPU, 4 ranks (les_ibm cases pin [mpi] dims = 1 1 1, so those run on 1):
  min_channel  les_ibm  les_ibm_refine  beltrami_yslab  turb180  wf180_y30  lam30t
  -> ALL PASS, max_abs = 0 on un/vn/wn/pn (+ nut, k, omega, gamma, rethetat)
GPU, 1 rank: same 7 cases -> ALL PASS, max_abs = 0
```

Both sides are nofma builds (`compile_nofma.sh`); the reference binaries were
built from the pre-change tree. The gate covers S0 **and** S1 together:
`count = 0` never allocates a scalar plane, never calls a scalar kernel and
never issues the extra exchange, so the argument is by construction, not by
cancellation.

## Results — S2 (2026-08-03, branch `scalar`) — ALL PASS

Local runs (RTX 3060 GPU + nvhpc 25.9 CPU builds, 24 cores; the GPU was
shared with another job throughout, so the timings are not representative).

**(m) Kays–Crawford correlation** (`mpirun -n 1 build_cpu/scalar_test`):

```
scalar_test: ALL PASS
```

29 tabulated values at `Prt_inf` = 0.85 / 0.5 / 1.0 / 0.4 / 1.2 matched to
1e-13 relative, covering the `Pe_t = 0` guard (`Pr_t = 2 Prt_inf` exactly),
the direct expression, the small-`x` series, both sides of the `x = 1/2`
crossover (7.2310152534 / 7.2310152679 — the two branches join to 1.3e-10,
their own slope), `Pe_t = 1e300` (`Pr_t = 0.85` exactly, where the direct
expression overflows `a²` and its `1 - exp(-x)` rounds to 0), and a
600-point monotonicity sweep over `Pe_t` = 1e-4…1e8.

**(i) steady k-ω SST, resolved walls** (`turbsst.ini`, 186711 steps to
`t = 200`, CPU 4 ranks — a true fixed point, no time averaging):

| | `theta` (constant `Pr_t`) | `theta_kc` (Kays–Crawford) |
|---|---|---|
| `theta_tau` | 0.053413 | 0.049468 |
| total flux `J(y)` constancy | **5.2e-08** | **2.6e-08** |
| `theta` vs the `nut`-integral prediction | **0.096 %** of the wall difference | **0.094 %** |
| predicted `theta_tau` | 0.267 % | 0.262 % |
| `theta+/(Pr y+)`, `y+ ≤ 2` | 1.0000 | 1.0000 |
| `Pr_t` realised over the faces | 0.8500 | **0.8832 … 1.7000** |

The sharp gate here is the third row: because steady RANS has no resolved
fluctuations, the scalar equation reduces to one ODE whose only coefficient
is the `nut` the solver itself wrote, so integrating
`dθ/dy = −J/(1/(Re Pr) + nut/Pr_t(y))` from the snapshot's own `nut` and
comparing with the snapshot's own `theta` tests the S2 face diffusivity and
the `Pr_t` model **with the turbulence model divided out**. It reproduces the
solver to 0.1 % of the wall difference (the trapezoid-vs-face-flux
discretisation difference), for both `Pr_t` models.

`theta_kc` exercises `prt_kays` in the kernel over its whole range: the
correlation runs from 1.7000 = `2 Prt_inf` at the wall (`Pe_t → 0`) to 0.8832
in the core (`Pe_t = 15`), i.e. it DAMPS the near-wall eddy diffusivity, and
the wall flux duly drops 7.4 % below the constant-`Pr_t` scalar in the same
run. (The Python mirror of the correlation in `check_scalar_turb.py` is what
makes its flux constant to 2.6e-08; with the constant-`Pr_t` formula the same
field reads 0.19 — an incidental confirmation that kernel and mirror agree.)

Log-layer slope, reported for context, NOT gated against `Pr_t/kappa`: the
scalar carries a CONSTANT flux while the channel's momentum flux falls
linearly, so the reference is `Pr_t/(kappa (1 − y/h))`, which is 2.07 only as
`y/h → 0`. Measured over y+ ∈ [30,60] (y/h 0.17–0.33): 3.84 (`theta`) /
4.16 (`theta_kc`) against the window's 2.72.

**(h) developed turbulent channel under LES** (`turbles.ini`, WALE, Re_tau
180, Pr 0.71, 10 snapshots over t = 30.2 … 33.8):

```
theta_tau = 0.050860       (first half 0.051120, second half 0.050599 -> 1.0% drift)
total flux J(y): 0.049161 .. 0.052360,  max|J - J_wall|/theta_tau = 0.0334
theta+ vs Kader, y+ in [1,35]:   max rel dev 0.2067 (at y+ 6.3),  mean 0.0682
viscous sublayer theta+/(Pr y+), y+ <= 3:   0.9919 .. 1.0081
theta+/U+ : first cell 0.7233 (Pr = 0.71),  y+ 20-40  0.8558 (Pr_t = 0.85)
resolved <v'theta'>/J at y+ = 5:0.07  15:0.52  30:0.80  60:0.92  120:0.94
   peak +0.048494 = 0.953 theta_tau u_tau at y+ 155
```

Reading these:

- **`theta+/U+` is the sharpest physical statement.** At the first cell it is
  0.7233 against the molecular `Pr` = 0.71 (1.9 %), and across the log layer
  0.8558 against the configured `Pr_t` = 0.85 (0.7 %). The two halves of the
  face diffusivity are therefore each delivering what they were set to, in a
  genuinely turbulent field.
- **The turbulent heat flux** rises from 7 % of the total at y+ 5 through
  52 % at y+ 15 to 92–94 % beyond y+ 60, with its peak at 0.95 `theta_tau
  u_tau` — the textbook split, and the reason the closure only has to supply
  the remainder.
- **The constant-flux residual (3.3 %)** is simultaneously the convergence
  measure and a conservation check: the setup makes `J(y)` exactly constant
  in statistical steady state, and 10 snapshots over 3.6 t.u. leave that much
  statistical noise. `theta_tau` moves 1.0 % between the two halves.
- **Kader** is matched to 6.8 % in the mean over the wall layer. The 20.7 %
  maximum sits at y+ = 6.3, on **Kader's own blend kink**: the correlation
  jumps 2.49 → 3.63 → 5.61 across three points where the computed profile
  runs smoothly 2.64 → 4.29 → 5.86. The window stops at y+ 35 deliberately —
  Kader describes a constant-flux WALL layer, and beyond y/h ≈ 0.2 the
  channel's momentum flux has fallen enough that `U+` flattens into the wake
  while the constant-flux scalar keeps a log-like slope. That divergence is
  the setup, not the closure (the SST case above shows the same thing at a
  fixed point).

**(j) 2:1 wall-band-refined channel — NO spurious scalar band**
(`turbslab.ini`, 8 snapshots against the 10-snapshot uniform control):

```
2:1 interfaces at y+ = 89.2 and 270.8 (= 89.2 from the far wall)
theta'_rms band ratio (refined/control): adjacent 1.0014, core 1.0002, excess +0.0012
   per-row ratios from the interface outward:  1.011  0.992  1.000  1.001  1.004  1.001
<theta> footprint |refined - control|: adjacent 0.0320, core 0.0286 (wall difference 2.0)
```

The interface-adjacent coarse rows carry the same scalar fluctuation energy
as the unrefined control to **0.1 %**, and no row within six of the interface
departs from the control by more than 1.1 %. This is the metric that caught
the momentum-reflux u'/v' band in `../channel_interface/` — here there is
nothing to catch. (The comparison is made on the COARSE rows, where refined
run and control share cells exactly; the fine rows resolve more scales by
construction, so a difference there would be physics, not a band.)

**Cost of the S2 kernel change with NO turbulence model.** In a DNS run the
eddy branch is never taken (`turb%nut` does not even exist), but the kernel
still initialises the six face diffusivities and carries the six face-`nut`
registers, so the dead code was measured rather than assumed: the S1 and S2
**nofma** binaries were run alternately on the same DNS channel (`turbles.ini`
with `[les]` removed, 64x48x64 + one scalar, 100 steps, 5 repetitions each,
interleaved to cancel machine drift), together with the `count = 0` twin whose
code path is provably identical and therefore measures the noise floor:

| | S1 s/step | S2 s/step | S2/S1 |
|---|---|---|---|
| CPU 4 ranks, `count = 1` | 0.2102 | 0.2102 | **0.9998** |
| CPU 4 ranks, `count = 0` (identical code = noise floor) | 0.2006 | 0.1999 | 0.9962 |
| GPU, `count = 1` | 0.2120 | 0.2128 | **1.0038** |
| GPU, `count = 0` (identical code = noise floor) | 0.1950 | 0.1963 | 1.0068 |

The `count = 1` change is smaller than the `count = 0` noise on both paths, so
S2 costs nothing measurable in DNS. (The GPU was shared with an unrelated job
throughout, which is why its absolute numbers are high and its noise floor is
0.7 %; the CPU box was idle.) For reference, the scalar transport as a whole —
S1's cost, not S2's — is **5.2 %** of the step for one scalar on this grid.

**(k) wall functions + scalars rejected** (`wferr.ini`):

```
ERROR STOP [scalar] with [rans] wall_treatment = wall_function is not implemented
           (a thermal wall function is increment S5); use resolved
```

(SUPERSEDED by S5a, which implements the thermal wall function: the same ini
is now the positive control that the combination is accepted — see the S5a
results below.)

**(l) determinism on a scalar + turbulence case** (`detles.ini`, 20 steps,
nofma builds):

```
1 rank vs 4 ranks:  un vn wn pn nut theta   max_abs = 0.0   (EXACT)
CPU vs GPU:         un vn wn pn nut theta   max_abs = 0.0   (EXACT)
```

**`[scalar] count = 0` bit-exactness, S2 vs the S1 binaries**
(`run_bitexact.sh`, nofma, 7-case suite, every dataset at tolerance 0):

```
CPU, 4 ranks: min_channel les_ibm les_ibm_refine beltrami_yslab turb180 wf180_y30 lam30t
              -> ALL PASS, max_abs = 0 on un/vn/wn/pn (+ nut, k, omega, gamma, rethetat)
GPU, 1 rank:  same 7 cases -> ALL PASS, max_abs = 0
```

**Scalar run with turbulence OFF, bit-exact vs the S1 binaries**
(`run_bitexact_s1.sh`, nofma, tolerance 0) — the other half of the
by-construction argument: every face diffusivity keeps the molecular value
exactly, the eddy term sitting behind `if (useNut)`:

```
CPU, 1 rank: uniform3 (2 scalars, 3-level refined) det conduction prsweep conserve
             -> ALL PASS, max_abs = 0 on un/vn/wn/pn and every scalar
GPU, 1 rank: same 5 cases -> ALL PASS, max_abs = 0
```

**S1 gates re-run with the S2 binary** (`run_gates.sh`, all seven groups, CPU
1 rank) — the same numbers as the S0/S1 results above, to every digit:

```
uniform      max|theta - const| = max|phi - const| = 0.000e+00, uniform-flow gate exact
conserve     relative drift -6.838e-19
cond_lin     L2 = 0.000000e+00                         (CPU build -> exactly 0)
cond_16/32/64  L2 = 2.814767e-03 / 7.054808e-04 / 1.764823e-04   (order 2.00 / 2.00)
wave16/32/64   L2 = 1.554700e-02 / 3.909091e-03 / 9.786703e-04   (order 1.99 / 2.00)
pr 0.1/1/10  max|s - series| = 9.995e-06 / 1.415e-04 / 1.584e-03
             wall-flux ratio  1.001092 / 1.000366 / 1.004115
restart      dataset present -> max|s - file| = 0.000e+00; absent -> warn + 99.0
```

## Results — S3 (2026-08-04, branch `scalar`) — ALL PASS

Local runs (RTX 3060 GPU — shared with an unrelated job throughout — and
nvhpc 25.9 CPU builds).
`./run_gates_s3.sh [solid|conserve|balance|prep|missing|refine|det|cyl]`.

**(n) dirichlet body mode — the solid cell IS the body value** (`ibmwavy.ini`,
analytic wavy wall, 32×32×8 over a 1 × 0.25 × 0.25 box, `ibm_value = 1`):

```
cold start, 200 steps      704 solid cells,  max|theta - 1| = 0.000e+00
after the seeded restart   704 solid cells,  max|theta - 1| = 0.000e+00
refine_body, 2 levels     5984 solid cells,  max|theta - 1| = 0.000e+00
through the FILE path      704 solid cells,  max|theta - 1| = 0.000e+00
```

An EQUALITY, not a tolerance: inside the body `coef_p = SOLID/Re`, so
`mu_s = 1/(1 + dt_gamma coef_p/Pr) ≈ 4e-25` falls far below 2⁻⁵³,
`(1 - mu_s)` is exactly 1 and `ibm_value + mu_s·ŝ` rounds back to
`ibm_value`. (With `ibm_value = 0` the same cells would keep the residual
`mu_s·ŝ ~ 1e-25` instead of exactly 0 — the velocity penalization's `~1e-26`
behaviour, CLAUDE.md Phase 2. That residual is also what the heat-release
diagnostic needs and does not get; see the FINDING under (o2).) The last
line uses the PREPARED
`coef_p` tiles to select the solid cells instead of the analytic indicator;
the two select **identical** 704-cell sets.

**(o) adiabatic body mode — conservation and the sealed cells** (the `phi`
scalar of the same case, seeded with `2 + sin(k·x)`, 200 steps):

```
int phi dV: 1.2500000000000000e-01 -> 1.2500000000000000e-01
relative drift (per max|s| V, V = 0.0625) = 0.000e+00
sealed (all six staggered faces solid) cells: 624, max|delta phi| = 0.000e+00
```

Both halves matter. The drift is EXACTLY zero because the face mask is
symmetric across a face, so the flux form telescopes with the body in place
exactly as without it. The second line is the positive control — conservation
alone would also hold if the mask never fired, and 624 cells sealed on all
six sides are frozen to the last bit, which only happens if it does.

**(o2) the dirichlet body's heat release closes the energy budget** — and the
A2 penalization integral does NOT measure it (`ibmwavy.ini` again; the case
has no boundary flux at all, so the heat the body releases must be exactly
the rate at which the domain stores it):

```
body heat release  1/2 [Q(t1) + Q(t2)] = 8.1028873772e-02   (Q: 8.197512e-02 -> 8.008263e-02)
storage rate  d/dt int theta dV        = 8.0997088108e-02   (dt = 0.002, 10 steps)
relative difference = 3.924e-04
   [the penalization integral alone would give 5.104632e-02 = 63 % of the truth]
```

The closure to **4e-4** (the residual is the trapezoid's curvature error —
over a 200-step window, where `Q` falls 29 %, it grows to 3.1 % exactly as
`O(dt²)` predicts) is a full discrete energy-budget check on the Dirichlet
body: everything the penalization injects is accounted for, and nothing else
enters.

**FINDING (this is why `Q` above is not the penalization integral).** The A2
force diagnostic `F = ∫coef·u dV` does NOT transpose to a Dirichlet scalar.
A solid cell's stored scalar is the body value **to the last bit** — that is
gate (n) — so `coef_p (s_body − s)` evaluates to `1e28 × 0 = 0` there, while
the cell is in fact re-heated every substage by exactly the flux it loses to
its fluid neighbours. Measured split on this case: solid cells contribute
`0.000000e+00`, graded fluid cells `3.74e-02`, truth `6.95e-02`. The force
version survives the same arithmetic only because `u_body = 0` and the stored
velocity keeps a ~1e-26 residual whose product with `coef` is O(1) — the
scalar's residual is flushed to zero by the rounding. The heat release is
therefore measured as **staircase-interface flux + graded-cell
penalization** (`check_scalar_ibm.py surface`), which is cancellation-free,
and cross-checked against a control volume in the fluid (`cv`).

**(p) `moby_prepare`: `coef_p_blocks` and nothing else**

```
h5same.py ibmwavy_case_ns.h5 ibmwavy_case.h5 --ignore coef_p_blocks   -> IDENTICAL
(without --ignore: dataset lists differ by exactly coef_p_blocks)
prepare 1 rank == 4 ranks                                             -> IDENTICAL
solve from the case file vs the inline analytic solve, tolerance 0:
   un vn wn pn theta phi   max_abs = 0.0
coef_p_blocks vs the independent transcription of the graded formula:
   level 0: 3456 cells checked (144 graded, 252 solid), worst rel dev 1.399e-16
```

The case file prepared from the `[scalar]`-stripped twin is dataset-for-dataset
and bit-for-bit identical, so no existing case file changes. The transcription
check is independent of the solver: it re-derives `Σ((d0−d)/d)/d0²/Re` in
Python, bisection included, from the analytic wall.

**(p2) a case file without `coef_p_blocks` + `[scalar]`: hard error**

```
error: [scalar] is configured but the coefficient file
  carries no coef_p_blocks (cell-centred scalar coefficients);
  re-run moby_prepare with [scalar]: ibmwavy_case_ns.h5
ERROR STOP
```

This is not hypothetical: it caught the inherited S1 `uniform3` gate, whose
zero-force twin predates S3. `../multilevel_body/make_uniform_twin.py` now
writes a zeroed `coef_p_blocks` beside the zeroed `coef_blocks` (the scalar
analogue of "the body exerts no force" is "the body penalises no scalar").
The twin is a GENERATED artifact, not a tracked file — a stale local copy is
repaired by re-running `../multilevel_body/setup.sh`, NOT by the "re-run
moby_prepare with [scalar]" the error message suggests, which is the right
advice only for a real case file.

**(q) `refine_body`: the coefficients are evaluated at each leaf's own level**
(`ibmwavyr.ini`, 72 leaves = 8 level-0 + 64 level-1):

```
level 0: 1728 cells checked (0 graded, 0 solid), worst rel dev 0.000e+00
level 1: 13824 cells checked (432 graded, 2496 solid), worst rel dev 2.173e-16
solve from the multi-level case file vs inline analytic, tolerance 0:
   un vn wn pn theta phi   max_abs = 0.0
```

Level 0 carries no body cells at all under `refine_body` (everything touching
the wall is refined), which is why its deviation is trivially 0; level 1 is
where the graded coefficients live and they match the transcription to
round-off. The bit-exact solve is the second, stronger statement: the file's
per-level `coef_p` tiles and the solver's own inline `VAR_P` pass agree to
the last bit at every level, or `theta` would differ.

**(r) determinism on a scalar + IBM case** (`ibmwavy.ini`, 50 steps, nofma
builds, tolerance 0 on `un vn wn pn theta phi`):

```
1 rank vs 4 ranks   (analytic path)                        max_abs = 0.0
CPU vs GPU          (analytic path)                        max_abs = 0.0
file path vs analytic path (CPU)                           max_abs = 0.0
CPU vs GPU          (file path, coef_p_blocks both sides)  max_abs = 0.0
CPU vs GPU          (+ WALE LES, ibm_aware; incl. nut)     max_abs = 0.0
```

The last line is the three code paths at once — eddy diffusivity, body
coefficients and the scalar — on both devices.

**(s) heated cylinder, Re 40, Pr 0.71** (`cylheat.ini`, 512×512×8, restarted
from the A2 steady field at t = 100 with `pn` zeroed; `C_D` holds at
1.690 ± 0.003 and `|C_L| < 6e-4` throughout, i.e. the A2-validated wake is
undisturbed). Nusselt from the body-local method, snapshot by snapshot:

```
t = 105   Nu = 3.5199        (staircase flux 0.1825 + graded penalization 0.2069)
t = 110   Nu = 3.4149
t = 115   Nu = 3.3797
t = 120   Nu = 3.3655        drift 0.42 % over the last 5 t.u.
```

and the INDEPENDENT Gauss/CV border flux on the t = 120 field (storage term
included; it has fallen to 2.7 % of the outflow, which is the second
convergence measure):

| control volume | Nu (CV) | vs the body-local 3.3655 |
|---|---|---|
| `[4,8] × [6,10]` | 3.3958 | **0.90 %** |
| `[3,9] × [5,11]` | 3.4321 | 1.98 % |
| `[2,10] × [4,12]` | 3.4733 | 3.20 % |

The two measurements share nothing but the `theta` field: one integrates the
flux across the staircase solid/fluid faces plus the penalization delivered
into the graded band, the other the advective + diffusive flux through a
border several diameters away. Their spread with box size is the same
signature the A2 momentum cross-check showed (0.1 % / 2.4 % / 6.5 % there) —
longer faces accumulate more collocation and gradient error.

Against the literature, `Nu = 3.37` sits on **Churchill–Bernstein's 3.35**
(0.5 %) and inside the numerical band 3.2–3.5 for Re = 40, Pr = 0.71. The
setup's known biases push the same way as they do for `C_D`: ~6 % blockage
from the Dirichlet far field at 16 D, and a first-order staircase body whose
effective diameter is ~D + h.

**`[scalar] count = 0` bit-exactness, S3 vs the S2 binaries**
(`run_bitexact.sh`, nofma, 7-case suite, every dataset at tolerance 0):

```
CPU, 4 ranks: min_channel les_ibm les_ibm_refine beltrami_yslab turb180 wf180_y30 lam30t
              -> ALL PASS, max_abs = 0 on un/vn/wn/pn (+ nut, k, omega, gamma, rethetat)
GPU, 1 rank:  same 7 cases -> ALL PASS, max_abs = 0
```

**Scalar run WITHOUT AN IMMERSED BODY, bit-exact vs the S2 binaries**
(`run_bitexact_s3.sh`, nofma, tolerance 0) — the third by-construction
property: with `useIbm` off, `adiab` is `.false.` for every scalar, the six
face masks collapse to the FACE_CLOSED flags and the penalization statement
is never entered, so the S2 arithmetic survives byte for byte. It is also the
cheap, sharp form of "every S1/S2 gate still reads the same number": the S2
gate cases are campaigns costing hours to days (turbsst alone is 186711
steps), and if the S3 binary reproduces the S2 binary's fields at tolerance 0
on those very inis, no statistic derived from them can have moved.

```
CPU, 1 rank: uniform3 det conduction prsweep conserve      (the S1 cases)
             turbles turbslab turbsst detles               (the S2 cases)
             -> ALL PASS, max_abs = 0 on un/vn/wn/pn, nut, k, omega and every scalar
```

`uniform3` is the one case in that list that DOES enter the penalization
statement (its zero-force twin has `[ibm] enabled = true` with all
coefficients zero), so its exactness rests on the IEEE identity
`x*1.0 + 0.0*v = x` rather than on the branch not being taken — the same
class of argument as the IDDES `fd_force = 0` blend identity, and measured
rather than assumed.

**S1 gates re-run with the S3 binary** (`run_gates.sh`, CPU 1 rank) — the
same numbers as the S0/S1 results above, to every digit:

```
uniform      max|theta - const| = max|phi - const| = 0.000e+00, uniform-flow gate exact
conserve     relative drift -6.838e-19
cond_lin     L2 = 0.000000e+00                         (CPU build -> exactly 0)
cond_16/32/64  L2 = 2.814767e-03 / 7.054808e-04 / 1.764823e-04   (order 2.00 / 2.00)
wave16/32/64   L2 = 1.554700e-02 / 3.909091e-03 / 9.786703e-04   (order 1.99 / 2.00)
pr 0.1/1/10  max|s - series| = 9.995e-06 / 1.415e-04 / 1.584e-03
             centre s = 0.990836 / 0.227799 / 0.000002
restart      dataset present -> max|s - file| = 0.000e+00; absent -> warn + 99.0
```

All seven groups, every digit identical to the S0/S1 results table above.

**S2 gates re-run with the S3 binary**: the Kays–Crawford unit test
(`scalar_test: ALL PASS`, 29 tabulated values) and the wall-function
rejection reproduce identically. The three turbulent-channel CAMPAIGNS
(`sst`, `les`, `band`) are not re-run wholesale — `turbsst` alone is 186711
steps — because `run_bitexact_s3.sh` proves the stronger statement on the
very same inis: the S3 binary reproduces the S2 binary's fields at
TOLERANCE 0, so no statistic derived from them can have moved.

## Results — S4 (2026-08-04, branch `scalar`) — ALL PASS

Local runs (RTX 3060 GPU — shared with an unrelated job throughout — and
nvhpc 25.9 CPU builds). `./run_gates_s4.sh [group]`.

**(t) (t2) the in-solver rows ARE the snapshot's rows** (`s4stats.ini`, the S2
LES channel + two scalars, `check_scalar_stats.py profile`; the deviation is
the summation order alone — the solver reduces with atomics over 64×64 cells
per row, the checker with numpy):

```
ONE sample  (step 67610):    theta 2.881e-14, theta_kc 2.780e-14  (worst of 7 columns)
FOUR samples (67610..67640): theta 2.596e-14, theta_kc 2.683e-14
   columns: <s>  <s^2>  <u_c s>  <v s>|lo  J|lo  <v s>|hi  J|hi
```

`theta_kc` is the same scalar with `prt_model = kays`, so the SECOND column
block also exercises the Kays–Crawford branch of the statistics kernel's face
diffusivity in-kernel, and the pair gates the multi-scalar column layout.

**(t3) the plane (boundary-layer) layout**, same case, `stats_layout = plane`
(3072 rows of the global 64 × 48 plane): worst deviation **4.8e-16** /
4.5e-16 (`theta` / `theta_kc`) — smaller than the profile layout's because a row is one z line, not a
whole x-z plane. It is gated on the CHANNEL deliberately: the layout is chosen
by `[scalar] stats_layout`, not by the flow case, so a boundary-layer campaign
would exercise the same kernel far more slowly, and the plane statistics of a
channel are perfectly well defined (x-inhomogeneous rows, z averaged).

**(t3b) the 2:1 wall-band refined case** (`turbslab.ini` + statistics): the
row tables are PER LEVEL, so the run writes `s4slab.h5` (level 0) and
`s4slab_l1.h5`. Checked with `check_scalar_stats.py rows`, which accumulates
straight from the leaves (no global box) and therefore reaches what `profile`
cannot:

```
level 0: 16 of 48 rows covered (the unrefined core),   worst 2.563e-14
level 1: 64 of 96 rows covered (the two wall bands),   worst 2.122e-13
```

**(t4) restart continuation** — 2 samples, restart, 2 more samples, versus one
40-step run over the same four sample steps:

```
count / raw_sum / profile:  max|a - b|/max|a| = 0.000e+00   (EXACT)
(and the FIELDS of the restarted run are bit-identical to the continuous one)
```

**(t5) determinism of the statistics** (a tolerance, not an equality: the
sampling kernel reduces with atomics and the write reduces across ranks):

```
1 rank vs 4 ranks   2.346e-14
CPU vs GPU          9.383e-16
```

**(t6) the statistics do not touch the solution** — statistics ON vs OFF,
same binary, `un vn wn pn nut theta theta_kc`: **max_abs 0** on every dataset.

**(u) the body heat release: the solver reproduces the validated Python form**
(`ibmwavy.ini` + `heat_interval`, 20 samples over 200 steps, against
`check_scalar_ibm.py surface` on the very snapshots the same run wrote):

```
worst relative deviation over all samples and both terms:  2.8e-15
energy budget with the SOLVER's own Q (10-step window):
   1/2[Q(t1) + Q(t2)] = 8.5104733941e-02
   d/dt int theta dV  = 8.5066583783e-02      rel 4.5e-04
```

The second block is the physics: the wavy case has no boundary flux at all, so
the heat the body releases must be the rate at which the domain stores it, and
it is — to the trapezoid's `O(dt²)`, the same 4e-4 class as the S3 Python
measurement. (The runtime file is written at `ES24.16` for exactly this
comparison: `ES16.8` put a 1e-9 floor under it.)

**(u2) an adiabatic scalar exchanges nothing — the positive control**: the
seeded `phi` spans **1.004815 … 2.995185** (manifestly non-zero) and its three heat
columns read **0.000e+00** exactly. No penalization is applied and every body
face is masked, so any other number would be a flux the solver never applied.

**(u3) heated cylinder Re 40, Pr 0.71: the RUNTIME Nusselt number** (20 steps
from the S3 campaign's converged `t = 120` field):

```
staircase/Lz 1.739023e-01 + graded/Lz 1.983667e-01 = Q/Lz 3.722690e-01
Nu = 3.3653     (S3, post-processed at t = 120: 3.3655; Churchill-Bernstein 3.35)
solver vs check_scalar_ibm.py surface on the same snapshot:  1.3e-15
```

**(v) tooling.** `compare_fields.py` with no dataset arguments discovers
`un vn wn pn nut theta theta_kc` on a scalar+LES snapshot (canonical
variables first, the rest alphabetically; the `blocks` table and the node
lines drop out because they do not have `un`'s rank), and
`tools/scalar_stats.py profile` reads the statistics files this session
wrote:

```
wall flux: y_min +5.088135e-02   y_max +4.885805e-02  -> theta_tau = 0.049870
total flux J(y): max|J - J_wall|/theta_tau = 0.139        (FOUR samples over
   0.02 t.u.; the S2 10-snapshot average over 3.6 t.u. of the same case reads
   0.033 -- this is the statistical noise of a short window, not an error)
Nusselt (q_w H / (D dT), H = 2, dT = 2) = 6.50
```

**`[scalar] count = 0` bit-exactness, S4 vs the S3 binaries**
(`run_bitexact.sh`, nofma, 7-case suite, every dataset at tolerance 0):

```
CPU, 4 ranks: min_channel les_ibm les_ibm_refine beltrami_yslab turb180 wf180_y30 lam30t
              -> ALL PASS, max_abs = 0 on un/vn/wn/pn (+ nut, k, omega, gamma, rethetat)
GPU, 1 rank:  same 7 cases -> ALL PASS, max_abs = 0
```

**Scalar runs vs the S3 binaries** (`run_bitexact_s3.sh`, nofma, tolerance 0)
— the S1 + S2 + S3 scalar cases, i.e. the cheap sharp form of "every earlier
gate still reads the same number": with the statistics intervals off (the
default) nothing is allocated and no kernel is called, so the S3 arithmetic
survives byte for byte.

```
CPU, 1 rank: uniform3 det conduction prsweep conserve turbles turbslab turbsst detles
             -> ALL PASS (9/9), max_abs = 0 on un/vn/wn/pn, nut, k, omega and every scalar
```

The S4 statistics are therefore bit-exact-by-construction on TWO counts: the
`[scalar] count = 0` suite above (no scalar at all) and this one (scalars, no
statistics). The third — statistics ON leaving the fields untouched — is
gate (t6).

**Earlier increments re-run with the S4 binary** — the same numbers as the
tables above, to every digit:

- `run_gates.sh` (S1, all seven groups): `cond_lin` L2 = 0.000000e+00;
  `cond_16/32/64` 2.814767e-03 / 7.054808e-04 / 1.764823e-04;
  `wave16/32/64` 1.554700e-02 / 3.909091e-03 / 9.786703e-04; Pr 0.1/1/10
  max|s − series| 9.995e-06 / 1.415e-04 / 1.584e-03 with wall-flux ratios
  1.001092 / 1.000366 / 1.004115; restart present/absent as before.
  NOTE the `det` group must be given the **nofma** binaries
  (`BIN=../../build_cpu_nofma/moby_solve GBIN=../../build_gpu_nofma/moby_solve`):
  it compares CPU vs GPU at tolerance 0, and with the default FMA-contracted
  builds it reads ~1e-14 on the velocities and 6.7e-16 on the scalar — the
  documented contraction difference, not a regression. On the nofma pair:
  1 rank == 4 ranks and CPU == GPU both **max_abs 0**.
- `run_gates_s3.sh` (the body gates, all groups except the expensive `cyl`
  campaign): **19/19 PASS**.
- `run_gates_s2.sh kays` (`scalar_test: ALL PASS`, 29 tabulated values) and
  `wferr` (the wall-function rejection) reproduce identically. The three
  turbulent-channel CAMPAIGNS are not re-run wholesale — `run_bitexact_s3.sh`
  proves the stronger statement on the very same inis (tolerance 0).

## Results — S5a (2026-08-04, branch `scalar`) — ALL PASS

Local runs (nvhpc 25.9 CPU builds, 1 rank for the channels — they are
8 x ny x 8 fixed points — and the RTX 3060 for the GPU legs; the GPU was
shared with an unrelated job throughout). `./run_gates_s5.sh [group]`.

**(m2) the correlations** (`mpirun -n 1 build_cpu/scalar_test`, the S2
driver extended):

```
scalar_test: ALL PASS
```

Jayatilleke's `P` at `Pr/Pr_t` = 0.835 / 1 / 2 / 0.1 / 7 / 100 (`P(1) = 0`
exactly — the Reynolds analogy), the thermal sublayer thickness `y+_T` for
five (Pr, Pr_t) pairs against mpmath's own root (air `Pr = 0.71` gives
**12.178**, slightly THICKER than the momentum `y+_lam = 11.530`; oil
`Pr = 7` gives 6.815; a liquid metal `Pr = 0.025` is conductive out to
284), and `wall_diffusivity` at `y+` = 5 / 12 / 30 / 45 / 100 — including
that it is **exactly 0** below `y+_T` and continuous across the switch
(6e-15 just above it).

**(w0) the resolved reference** (`wfsst.ini` = `turbsst.ini` restarted from
its converged `t = 200` field, 200 steps, one sample):

```
wall flux: y_min +5.341323e-02   y_max +5.341323e-02  -> theta_tau = 0.053413
total flux J(y): max|J - J_wall|/theta_tau = 0.0000        (a true fixed point)
Nusselt (q_w H / (D dT), H = 2, dT = 2) = 6.8262           theta+_c = 18.7220
```

`theta_tau = 0.053413` is the S2 number to every digit (S2 results table,
gate (i), post-processed by `check_scalar_turb.py rans`) — so the S4
statistics and the S2 analysis agree exactly on the reference, before the
wall-function cases are compared against it.

**(w) the DELIVERED-FLUX identity — the sharp gate.** The thermal wall
function is installed as a wall-cell eddy diffusivity `nu(y+/theta+ - 1/Pr)`
whose ghost copy makes the wall face see `D = nu y+/theta+`, so the
statistics' discrete wall flux must BE `u_tau* (theta_w - theta_1)/theta+`,
with `theta+` from an independent Python transcription of
Kader/Jayatilleke. Measured on all four grids, cell by cell (the solver
averages the FLUX over the wall plane and the flux is nonlinear in the
cell's own k):

```
y+_1 = 5    rel 1.3e-15 / 1.5e-15        y+_1 = 15   rel 1.3e-16 / 6.4e-16
y+_1 = 30   rel 1.4e-15 / 1.6e-15        y+_1 = 45   rel 1.1e-15 / 7.3e-16
```

This is the thermal analogue of T3's "u_tau from the delivered wall stress
= 1.0000": what the solver applies at the wall IS the closure, not an
approximation to it. (`u_tau` from the delivered wall stress reads
**1.0000** on every grid here too.)

**(x) the y+_1 sweep: graceful degradation** (`theta_tau` and the implied
centreline `theta+_c = 1/theta_tau`, against the resolved reference; the
`theta+_1` column is the first cell against Kader at the VISCOUS `y+`):

| y+_1 | branch | theta_tau | dev vs resolved | theta+_1 | Kader | dev |
|---|---|---|---|---|---|---|
| 5  | conduction (`y+_k` 1.4) | 0.054155 | **+1.39 %** | 3.550 | 3.041 | 16.8 % |
| 15 | log (`y+_k` 14.6) | 0.054629 | **+2.28 %** | 9.279 | 7.924 | 17.1 % |
| 30 | log (`y+_k` 28.9) | 0.056424 | **+5.64 %** | 10.832 | 10.819 | **0.1 %** |
| 45 | log (`y+_k` 43.1) | 0.057272 | **+7.22 %** | 11.771 | 11.854 | **0.7 %** |

Reading these:

- **Monotone degradation, no dip.** The wall flux drifts up by 1.4 % to
  7.2 % as the first cell walks from `y+ 5` to `y+ 45` — the same shape and
  the same magnitude class as T3's momentum sweep (centreline `U+` off by
  −3.1 % / +2.8 % / +3.0 % / 1.2 % / 0.7 %). A double-counting bug would
  show as a DEFICIT, and there is none.
- **On its own ground the closure is exact to Kader.** At `y+_1` = 30 and
  45 — the range wall functions are for — the first cell's `theta+` lands
  on Kader's correlation to 0.1 % and 0.7 %. The two buffer-layer grids sit
  ~17 % above Kader at `y+` 5 and 15, which is Kader's own blend region
  (the S2 LES gate measured up to 20.7 % there against a resolved DNS-grade
  profile) and the documented T3 behaviour of the first cell below `y+ 30`.
- **The two branches both fire.** `y+_1 = 5` puts the k-based `y+_k` at
  1.4, below `y+_T = 12.18`, so that case runs the CONDUCTION branch, where
  the wall diffusivity is exactly zero and the treatment degenerates to the
  resolved one — and it is duly the closest to the resolved reference.
- **The fixed points are converged**: the total-flux constancy
  `max|J - J_wall|/theta_tau` reads 7.6e-12 / 1.0e-12 / 1.6e-14 / 1.2e-14.

**(x2) the Kays–Crawford scalar of the same runs** (`theta_kc`, index 2;
the last column is `theta_tau(theta_kc)/theta_tau(theta)` WITHIN one run):

| case | theta_tau | dev vs resolved | s2/s1 |
|---|---|---|---|
| resolved reference | 0.049468 | — | **0.9261** |
| y+_1 = 5  | 0.050369 | +1.82 % | **0.9301** |
| y+_1 = 15 | 0.053397 | +7.94 % | 0.9775 |
| y+_1 = 30 | 0.055728 | +12.65 % | 0.9877 |
| y+_1 = 45 | 0.056812 | +14.85 % | 0.9920 |

This is the DESIGN, measured, not a defect: the thermal wall function is
defined with a CONSTANT `Pr_t` (P and the log branch are), so a wall cell
treats `theta` and `theta_kc` identically and Kays–Crawford can only act in
the interior, where `Pe_t` is large and the correlation sits at its
`Prt_inf` asymptote. Hence `s2/s1 -> 1` as the first cell moves out into
the log layer. The `y+_1 = 5` case, whose wall cells take the conduction
branch (i.e. the resolved arithmetic), reproduces the resolved reference's
own ratio 0.9261 to **0.4 %** — which is the positive control that the
branch structure, not a coincidence, is what produces the trend.

**(k, flipped) wall functions + scalars are now ACCEPTED**:

```
thermal wall function: P = -1.4915  y+_T =  12.178
```

**(z) determinism** (`wfs180_y30.ini`, 20 steps, nofma builds, datasets
`un vn wn pn nut k omega theta theta_kc`):

```
1 rank vs 4 ranks ([mpi] dims = 1 1 4)   max_abs = 0.0   (EXACT, all nine)
CPU vs GPU                               max_abs <= 3.3e-13 (pn), 4.6e-14 (un),
                                         1.3e-15 on both scalars
```

The CPU-vs-GPU spread is the T3 wall-function class: the `log()` intrinsic
differs by an ulp between host and device libm, and both wall functions
call it (resolved mode stays exactly CPU == GPU — see the bit-exactness
block below).

**FINDING (pre-existing, NOT S5a) — RESOLVED 2026-08-05.** The rank
comparison pins `[mpi] dims = 1 1 4` because this channel was **not
rank-independent under an x split**, with or without scalars:
`../rans_sst/wf180_y30.ini` with no `[scalar]` section at all, run with the
**S4 reference binary**, reproduced the same deviation to the last bit
(`un` max_abs 8.800384e-03 after 20 steps — and already 7.6e-04 after ONE
step, so it is the initial state or the first substage on an x-split rank
box, not an accumulation; `k` 9.0e-02, `omega` 1.2e+00). A z split was
EXACT on the same case.

The cause was in the RANS layer, exactly as suspected — but neither of the
two suspects named: `init_rans_transport`'s `k = 1.5 (tu/100 |u|)^2`
initial condition interpolates the cell velocity from the two staggered
faces, so at a block's last interior cell it reads the halo `q(nb+1)`,
which `moby_solve.f90` did not fill until AFTER the whole init block. `k`
came out a factor 4 low on the last plane of EVERY block — one bad plane
per block, hence the decomposition dependence — and only in x here because
this channel's `v` and `w` are zero in the IC. The fix moves `apply_bc` +
`exchange_halos` ahead of the `[rans]` init; the full write-up, the
measured initial-state numbers and the re-gate are in
`../rans_sst/README.md`. The `dims = 1 1 4` pin below is kept but is no
longer load-bearing: every decomposition (`4 1 1`, `1 1 4`, `2 1 2`) is now
max_abs 0 at 20 steps. Cold-started RANS snapshots written before
2026-08-05 no longer reproduce bit-for-bit; their converged physics does
(`wf180_y30` / `wf180_y45` re-run to `t_final` reproduce the T3 gate to
every printed digit).

**`[scalar] count = 0` bit-exactness, S5a vs the S4 binaries**
(`run_bitexact.sh`, nofma, 7-case suite, every dataset at tolerance 0):

```
CPU, 4 ranks: min_channel les_ibm les_ibm_refine beltrami_yslab turb180 wf180_y30 lam30t
              -> ALL PASS, max_abs = 0 on un/vn/wn/pn (+ nut, k, omega, gamma, rethetat)
GPU, 1 rank:  same 7 cases -> ALL PASS, max_abs = 0
```

**Scalar runs vs the S4 binaries** (`run_bitexact_s3.sh`, nofma, tolerance
0) — the resolved-wall path must be untouched, which is by construction:
`wallfn` is false, the new branch is not entered and the S4 arithmetic is
reproduced byte for byte (the `wfYplus` array is then a 1-cell dummy).

```
CPU, 4 ranks: min_channel les_ibm les_ibm_refine beltrami_yslab turb180 wf180_y30 lam30t
              -> ALL PASS, max_abs = 0 on un/vn/wn/pn (+ nut, k, omega, gamma, rethetat)
GPU, 1 rank:  same 7 cases -> ALL PASS, max_abs = 0_S3
```

**Earlier increments re-run with the S5a binary** — the same numbers as the
tables above, to every digit:

```
run_gates.sh    (S1, all seven groups -> "no failures" x7)
  uniform      max|theta - const| = max|phi - const| = 0.000e+00; uniform-flow gate exact
  conserve     relative drift -6.838e-19
  conduction   cond_lin L2 = 0.000000e+00 (CPU build);
               cond_16/32/64 L2 = 2.814767e-03 / 7.054808e-04 / 1.764823e-04  (order 2.00)
  wave         wave16/32/64 L2 = 1.554700e-02 / 3.909091e-03 / 9.786703e-04   (order 1.99/2.00)
  pr           max|s - series| 9.995e-06 / 1.415e-04 / 1.584e-03;
               wall-flux ratios 1.001092 / 1.000366 / 1.004115
  restart      dataset present -> max|s - file| = 0.000e+00; absent -> 99.0
  det (nofma)  1 == 4 ranks and CPU == GPU both max_abs 0; 2:1 residual 1.047e-05

run_gates_s3.sh (the body gates, every group except the cyl campaign) -> 7/7 ALL PASS
  solid 704/704/5984/704 cells at max|theta - 1| = 0.000e+00; adiabatic drift
  0.000e+00 with 624 sealed cells frozen exactly; energy budget 3.924e-04;
  coef_p_blocks vs the transcription 1.399e-16; prepare/file gates identical

run_gates_s4.sh (the statistics, every group except cyl) -> 10/10, no failures
  rows vs snapshots 2.881e-14 / 2.780e-14 (one sample), 2.596e-14 / 2.683e-14
  (four), plane 4.808e-16 / 4.503e-16, levels 2.563e-14 / 2.122e-13, restart
  0.000e+00, 1 vs 4 ranks 2.346e-14, CPU vs GPU 8.873e-16, statistics-off twin
  max_abs 0, body heat vs the Python form 2.817e-15, an adiabatic scalar exactly
  zero, tools theta_tau = 0.049870

scalar_test    ALL PASS (the S2 kays values + the S5a correlations)
run_gates_s2.sh wferr -> the flipped gate above; kays is scalar_test
```

Every number is the one the S1 / S3 / S4 result tables above record, to every
digit. The three turbulent-channel CAMPAIGNS (`sst`, `les`, `band`) are not
re-run wholesale — `run_bitexact_s3.sh` proves the stronger statement on the
very same inis (tolerance 0), which is the argument the S3 and S4 sessions
used.

## Re-gate after the 2026-08-05 solver fixes — ALL PASS

Two source changes landed on 2026-08-05: the RANS cold-start IC fix
(`moby_solve.f90`, write-up in `../rans_sst/README.md`) and the body-heat
double-count fix (`scalar_stats.f90`, below). The first CHANGES every
cold-started RANS run's initial condition, so the S2/S5a gate numbers had to
be re-measured rather than assumed. They are unmoved:

| re-gate | result |
|---|---|
| `run_gates_s5.sh unit` (the two correlations vs mpmath) | `scalar_test: ALL PASS` |
| `run_gates_s5.sh ref` (the RESOLVED reference `theta_tau`) | **0.053413**, flux constancy **0.0000** — the recorded S5a number to every digit |
| `run_gates_s5.sh sweep`, `theta_tau` at `y+_1` = 5/15/30/45 | 0.054155 / 0.054629 / 0.056424 / 0.057272 = **+1.39 % / +2.28 % / +5.64 % / +7.22 %** vs the reference — every one the recorded S5a number |
| the same, first-cell `theta+` vs Kader | **17.1 % / 0.1 % / 0.7 %** at `y+_1` = 15 / 30 / 45 — as recorded |
| the same, flux-constancy convergence | 7.3e-12 / 1.1e-12 / 4.6e-14 / 4.4e-14 (recorded 7.6e-12 / 1.0e-12 / 1.6e-14 / 1.2e-14) |
| `run_gates_s2.sh kays` (Kays–Crawford unit test) | `scalar_test: ALL PASS` |
| `run_gates_s2.sh wferr` (the S5a-flipped acceptance gate) | **PASS** |
| `run_gates_s2.sh sst` (steady SST, resolved walls, 186711 steps) | `theta`: `theta_tau` **0.053413**, flux constancy **5.225e-08**, vs the nut-integral prediction **0.0962 %**; `theta_kc`: **0.049468**, **2.613e-08**, **0.0941 %**, `Pr_t` spans **0.8832 … 1.7000**, ratio to the constant-`Pr_t` twin **0.92615** — every one the recorded S2 number |
| S5a gate (z) determinism, `wfs180_y30` 20 steps, nofma | 1 rank == 4 ranks **max_abs 0** on `un vn wn pn nut k omega theta theta_kc` — **for the x split too**, which is the decomposition the fixed bug broke (the `[mpi] dims = 1 1 4` pin is now genuinely unnecessary) |
| the same, CPU vs GPU (GPU leg on istmcetus) | **0.0** on eight datasets, **5.6e-17** on `theta_kc` — tighter than the recorded ≤3.3e-13 |
| `run_bitexact_s3.sh` (scalar without a body, 9 cases, nofma CPU) | **max_abs 0** on 8 of 9; `turbsst` MOVED (un 1.6e-02, k 3.2e-01, theta 2.0e-03) because it is a COLD-STARTED RANS case — the intended signature, same class as turb180/wf180_y30/lam30t |
| `run_bitexact.sh` MODE=gpu (7-case suite, GPU legs on istmcetus, compared locally) | non-RANS 4/4 **max_abs 0**; the 3 cold-started RANS cases moved, `turb180`'s deviation IDENTICAL to the CPU's to every digit |

LANDMINE for re-runs: `turbles.ini`, `turbslab.ini` and `turbsst.ini` pin
`[mpi] dims = 1 1 1`, so `run_bitexact_s3.sh` must be given `RANKS=1` for
them — with `RANKS=4` they error-stop on "MPI Cartesian dimensions do not
match the number of ranks", which looks like a failure of the gate and is
not. Also: istmcetus has **no h5py**, so the GPU suite runs its solver legs
there (`env_cetus.sh`) and the comparisons run locally on the shared files.

## Results — the IMMERSED-WALL thermal wall function (2026-08-05) — ALL PASS

S5a's remaining open item. Case `ibmwf180.ini` (see the gate table above),
prepared with `moby_prepare` so the case file carries `coef_p_blocks`, run
to a steady state on 1 CPU rank; commands:

```
mpirun -n 1 ../../build_cpu/moby_prepare ibmwf180.ini ibmwf180_case.h5
mpirun -n 1 ../../build_cpu/moby_solve  .ibmwf180_solve.ini     # coeff_file form
python3 check_scalar_ibmwf.py wall   ibmwf180_ransgeom.h5 ibmwf180_<step>.h5
python3 check_scalar_ibmwf.py budget ibmwf180_<step>.h5 ibmwf180_case.h5 \
        --heat ibmwf180_heat.txt
```

| gate | result |
|---|---|
| (y1) the regime EXISTS: classified IBM wall cells | **128** cut cells (centre inside the solid, one fluid staggered face); `wallcell == 2` (fully solid): 0 |
| (y1) the penalization acts on those same cells | `u` = 2.1e-27, `theta` = 1.6e-29 there (`ibm_value` = 0), i.e. pinned, while the wall-function machinery reads them |
| (y1) the LOG branch DURING THE TRANSIENT (t ≈ 2.1, k still high) | `y+_k` **12.74 … 13.58** — past the thermal `y+_T` = **12.1776** (Jayatilleke `P(0.71/0.85)` = −1.491461) and the momentum `y+_lam` = 11.5301: fires on **128/128** wall cells, and the wall-cell `nu_t` reproduces an INDEPENDENT transcription of `nu(y+ kappa/ln(E y+) − 1)` to **2.7e-15** (`nu_t` 4.55e-04 … 7.70e-04) — so the log-branch arithmetic at an IBM wall cell is verified |
| (y1) the LOG branch AT CONVERGENCE | **NOT REACHED**: `y+_k` falls to **5.46 … 6.17** (mean 5.81), i.e. back on the conduction branch, where the wall diffusivity is exactly 0 (`nu_t` 0.000e+00 on all 128) — see the note below. This grid therefore does NOT yet gate the converged log-branch regime. |
| (y2) the steady budget, closed form: the body heat release must be `source * V_fluid` = 4.626377063010637e-01 | solver runtime heat file **−4.626377063010643e-01**, rel dev **1.3e-15** (stable to 15 digits from `t ~ 130` on). `check_scalar_ibmwf.py budget` also checks the underlying identity from the SNAPSHOT — `sum coef_p (s_body − s) dV/Pr` over all interior cells = `−source*V_total`, rel dev **3.6e-16** — and the solver's `graded` column against the snapshot's unpinned sum, **1.5e-15** |
| (y3) determinism: cold start, 20 steps, 1 rank vs 4 ranks (nofma) | **max_abs 0** on `un vn wn pn k omega nut theta` — which also re-exercises the RANS cold-start IC path fixed the same day (see `../rans_sst/README.md`) |
| (y4) the diagnostic fix does not touch the solution | the fixed binary's fields vs the pre-fix run: `un vn wn k omega nut` **max_abs 0** |
| (y2) `ibm_value` invariance: the same physics at `ibm_value` = 0 and 1 | a CONTROLLED pair, verified as such: `un vn wn k omega nut` are **max_abs 0** between the two runs and `theta_1 − theta_0 = 1` everywhere to **1.0e-14**, so the only difference is the body-value convention. The staircase column is bit-identical to **13 digits** (−4.3874549474718e-01); after the fix below the totals agree to **1e-15** (before it they differed by 128 %) |
| (y1) the converged `y+` is convention-independent | both runs read `y+_k` 5.4568 … 6.1706 (mean 5.8137), as they must — the velocity and `k` fields are bit-identical |

**DEFECT FOUND AND FIXED (`src/modules/scalar_stats.f90`, 2026-08-05): the S4
body-heat diagnostic double counted whenever `ibm_value = 0`.** The
staircase/penalization split rests on "a solid cell contributes exactly 0 to
the penalization sum, because it holds the body value to the last bit" (the
S3 FINDING). That cancellation is an **artefact of the value**: it is bitwise
only when `ibm_value` is large enough to swallow the O(1e-29) penalization
residual under its own ulp. At `ibm_value = 0` the residual survives,
`coef_p*(0 − 1.6e-29)` is O(0.1) per cell, and the cell's heat is counted
twice — once as its own penalization, once as the staircase flux into it.
Measured on ONE physical problem run with both conventions (`theta` differs
by a constant, so the staircase column is bit-identical):

```
ibm_value = 1.0  ->  total  -4.62637706301056e-01   = source*V_fluid  (rel 1.7e-14)
ibm_value = 0.0  ->  total  -1.05559576981527e+00   = 128 % HIGH
                     graded -6.16850275068085e-01   = -source*V_TOTAL exactly:
                     the whole source input, including what is deposited inside
                     the body and never crossed the interface
```

The fix excludes solid cells from the penalization accumulator EXPLICITLY
instead of relying on the cancellation. It is a **no-op at `ibm_value = 1`**,
which is measured, not asserted: re-running the S4 heat gates with the fixed
binary reproduces their recorded numbers — solver-vs-Python worst rel dev
**1.315e-15** (tol 1e-12), energy budget **4.485e-04** (tol 1e-2), and the
adiabatic control **exactly 0.000e+00** — and the fixed binary at
`ibm_value = 0` now reproduces the `ibm_value = 1` answer to 1e-15. The
diagnostic is now independent of `ibm_value`, which is the invariance it must
have. The change touches the diagnostic only; solver fields are unaffected.

**WHEN the old code was wrong — pinned down 2026-08-05 by an `ibm_value = 0`
sweep of the S3 body gates.** It takes BOTH conditions, which is why the
defect hid for three increments:

| `ibm_value` | solid-cell rhs | solid cell holds | old diagnostic |
|---|---|---|---|
| 1.0 | anything | `1.0` BITWISE (the residual is under 1.0's ulp) | correct |
| 0.0 | zero (no source; `ibmwavy.ini` as committed) | **exactly 0.0** (`mus ~ 5e-30` drives it to zero) | correct |
| 0.0 | nonzero (a volumetric `source`) | **1.311e-29** — measured | **double counts** |

So the S3 gate (n) "solid cell == `ibm_value` EXACTLY" is NOT resting on the
value in its committed form: re-run at `ibm_value = 0` it still reads
**0.000000e+00** on all 704 cells. Add a `source` and the same gate reads
1.311e-29 and FAILs its exact-equality test. The scheme's real invariant is
`|s − ibm_value| <~ 1e-28` ABSOLUTE, not bitwise equality; bitwise is what you
get when the residual falls under `ibm_value`'s ulp, or when there is nothing
to drive it. Worth knowing before anyone tightens that gate or reads it as a
guarantee.

The fixed diagnostic is `ibm_value`-independent on a CURVED body too, not
only on the flat slabs: the wavy case run as a controlled pair (`ibm_value` 1
with `initial` 0, vs `ibm_value` 0 with `initial` −1 — the same problem shifted
by exactly 1, both with `source = 0.05`) reports body heat
**8.3661018093697e-02** and **8.3661018093128e-02**, agreeing to **6.8e-12**
(not bitwise: shifting the field changes the rounding, and the staircase and
graded terms partly cancel).

LESSON (worth carrying): never use a floating-point cancellation as a
classification. The S3 FINDING was right that a solid cell's penalization
integral vanishes — but "vanishes" was true of the cases it was measured on,
not of the arithmetic.

**The converged log branch: CLOSED at Re_tau 1000 (`ibmwf1000.ini`,
2026-08-05, istmcetus A6000 GPU).** Coarsening the grid is NOT the lever at
an immersed wall. At a DOMAIN wall the first cell's `y+` grows with the grid
spacing, which is how the S5a sweep reached `y+_1 = 45`. At an IMMERSED wall
the classified wall cell is a CUT cell whose velocity is penalized, so its
`k` stays small however coarse the grid gets, and
`y+_k = C_mu^(1/4) sqrt(k) y_eff/nu` only picks up the (bounded) growth of
`y_eff ~ dwall <= dy/2` — `ibmwf180` already has `y_eff` = 0.103, about 6x
T3's `ibm180wf`, and still converges to `y+ ~ 5.8`. THE LEVER IS `nu`: at
fixed `u_tau = 1`, `y+ ∝ 1/nu ∝ Re`. `ibmwf1000.ini` is `ibmwf180.ini` with
`re = 1000` and NOTHING else changed (the case file is re-prepared because
the stored coefficients carry the `1/Re` scaling; `dwall`, `yeff` and
`wallcell` come out **bit-identical**, max_abs 0 — so `y+` moves only through
`nu` and `k`, which is the controlled comparison the argument needs).

| gate (Re_tau 1000, converged at `t = 993.8`) | result |
|---|---|
| wall-cell `y+_k` | **38.29 … 41.03** (mean **39.66**) vs 5.81 at Re_tau 180 — past `y+_T` = 12.1776 and `y+_lam` = 11.5301 |
| the LOG branch at CONVERGENCE | fires on **128/128** wall cells, for BOTH wall functions — the regime S5a never gated |
| wall-cell `nu_t` vs the independent transcription of `nu(y+ kappa/ln(E y+) − 1)` | **0.000e+00** (exact, on the GPU); `nu_t` 1.648e-03 … 1.805e-03 |
| the closed-form budget, `source * V_fluid` | total **−4.6263770630e-01**, rel dev **2.3e-15**; `graded` vs the snapshot's unpinned sum **8.0e-16** |
| the same case run with the PRE-FIX GPU binary (an accident that became a control) | total **−1.0751501972e+00** = **132 % high**, exactly the reconstruction `check_scalar_ibmwf.py budget` prints — the double count reproduced independently at a second Reynolds number and on a second device |

So the Dirichlet penalization and the thermal wall function DO coexist
correctly on the same cut cell in the log branch: the penalization holds
`u`/`theta` at the body value while the wall function sets the cell's eddy
diffusivity, and the resulting wall heat flux is the exact closed-form one.

## Notes and limitations

- **Order of the substage calls.** `scalar_transport` runs immediately BEFORE
  `momentum(...)`, a deliberate deviation from `docs/next_session_scalar.md`
  §2 (which says "right after"): `momentum` ends by copying `qs → q` and the
  velocity halos are only refreshed by the exchange that follows it, so
  "after" would advect the scalar with the non-solenoidal predicted velocity
  AND with stale halo values. See the doc's STATUS header.
- **Conservation across a 2:1 interface is not exact** and is reported as a
  measurement, not a gate: the coarse and fine sides of an interface face
  compute their own fluxes from their own face velocities and transferred
  halo values, so the flux form telescopes exactly only across matching
  faces (the same structure that made the momentum reflux an artifact —
  `validation/channel_interface`).
- **Central convection is unbounded**: sharp fronts over/undershoot until the
  shared TVD/van-Leer increment lands (deliberately no upwind fallback).
- **What `[flow] convection = skew` means for a scalar.** The plan's formula
  (`conv_div − s·(div u)`) is the ADVECTIVE form: it preserves a uniform
  scalar exactly for any advecting field, but gives up the divergence form's
  exact global conservation. Subtracting HALF would be the skew-symmetric
  form the momentum kernel uses (neutral in `∑s²`, but leaving `½ s div u` on
  a uniform field). The implemented behaviour is the plan's; the gates above
  all run the DEFAULT divergence form, so neither is exercised by them.
- The diffusive flux is masked at `FACE_CLOSED` faces (blocks removed inside
  an immersed body hold a zeroed halo). Immersed-body scalar coefficients
  themselves are increment S3.
- **The thermal wall function uses the CONSTANT `Pr_t` at wall cells**, even
  under `prt_model = kays`: Jayatilleke's P and the log branch are defined
  with a constant `Pr_t`, and Kays–Crawford is a correlation for the
  resolved interior, where it keeps running. Gate (x2) measures exactly what
  that costs and why the trend is the design.
- **Under wall functions the face eddy diffusivity is the average of the two
  CELL diffusivities**, not `eddy_diffusivity` of the averaged `nu_t`: a
  wall cell's diffusivity is not a function of its `nu_t` at all. Away from
  wall cells and at constant `Pr_t` the two are the same expression in a
  different order; the wall-function branch is separate code, so resolved
  runs keep their arithmetic byte for byte.
- **IBM wall cells take the thermal wall function too** (the same
  `wallcell` marker `rans_assemble_nut` uses), but the S5a gates are DOMAIN
  walls, and no IBM wall-function scalar case was run here. On the T3 IBM
  channel the k-based `y+` is 2-3 (`../rans_sst/README.md`, gate (c)), i.e.
  the conduction branch of this closure, where the wall diffusivity is
  exactly zero and nothing new happens — but that is T3's measurement of
  the momentum side, not a scalar gate. A body whose wall cells reach the
  log layer is UNGATED: the Dirichlet penalization and the wall function
  would then both act on the same cell, which needs its own increment (and
  its coefficient file re-prepared with `[scalar]`, for `coef_p_blocks`).
- `check_scalar_stats.py` (the S4 snapshot recomputation) does NOT know
  about the wall function: it rebuilds the face diffusivity as
  `1/(Re Pr) + nut/Pr_t`, which is not what the solver applies at a wall
  cell in wall-function mode. Use `check_scalar_wf.py wall`, whose identity
  is sharper anyway.
