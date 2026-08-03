# Passive scalars — S0/S1 gates

Gates for increments **S0** (layout, config, io, halos, BCs — the scalar is
carried but not advanced) and **S1** (the transport kernel, molecular
diffusivity only) of `docs/next_session_scalar.md`. S2 (turbulent closure,
`nut/Pr_t` + Kays–Crawford) and S3 (cell-centred IBM coefficients) are NOT
part of this directory yet.

```bash
./compile_nofma.sh cpu && ./compile_nofma.sh gpu   # bit-exactness builds
./run_gates.sh [uniform|conserve|conduction|wave|pr|det|restart]
REF=<pre-change nofma binary> MODE=cpu RANKS=4 ./run_bitexact.sh
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

`run_bitexact.sh` is the **`[scalar] count = 0` bit-exactness** gate: the
standard 7-case suite (min_channel, les_ibm ± refine_body, Beltrami y-slab,
turb180, wf180_y30, lam30t) run with a pre-change nofma binary and the new
one, compared at tolerance 0 on every field dataset (`un vn wn pn` plus
`nut`, `k`, `omega`, `gamma`, `rethetat` where the case has them).

## Results (2026-08-03, branch `scalar`) — ALL PASS

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
