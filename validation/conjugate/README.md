# Conjugate heat transfer at an immersed interface — C1 gates

Gates for increment **C1** of `docs/next_session_conjugate.md`: the baseline
conjugate scheme, in which the solid stops being a boundary condition and
becomes a real unknown, and a face whose two cell centres straddle the
interface takes the **distance-weighted harmonic mean** of the two materials'
diffusivities, built on the level-set fraction `w = φ_L/(φ_L − φ_R)`.

```bash
./run_gates_c1.sh [slab|converge|weight|capacity|contact|peclet|limits|conserve|guard|refine|stats|det|all]
```

Environment: `BIN` (default `../../build_cpu/moby_solve`), `PREP` (the
`moby_prepare` next to it), `NBIN`/`NPREP`/`GBIN` (the nofma pair, `det`
group), `RANKS`. Helpers:

| file | role |
|---|---|
| `make_slab_stl.py` | the flat solid slab as an **ASCII** STL — `moby_prepare` parses ASCII vertices to float64, so the interface plane is exact and the analytic reference needs no quantisation dance. |
| `seed_slab_ic.py` | writes the exact two-material steady profile into a snapshot, for the fixed-point form of gate 1. |
| `check_conjugate.py` | `slab` / `weight` / `conserve` / `limit`. |
| `flux_limit.py` | the `κ_s → ∞` / `κ_s → 0` limits as a **rate**, which is the sharp half of gate 2. |

`slab.ini` is a template (`@STL@ @CASE@ @KAPPA@ @CAP@ @PREFIX@ @NSTEPS@
@WRITE@`); `wavy.ini` is the analytic-geometry case.

## Cases

| ini | gate |
|---|---|
| `slab.ini` | **1D two-material slab.** A flat solid slab fills `y < y_wall` (STL → `moby_prepare`, so this is the FILE geometry path); the fluid fills the rest; no flow anywhere. `T(0) = 0` on a Dirichlet face *inside the solid*, `T(L) = 1` in the fluid. The steady solution is piecewise linear with the exact series resistance, and a piecewise-linear field is an exact fixed point of the discrete operator **only if the cut face's series resistance is exact** — i.e. only if `w` is the true cut fraction. Nothing else in the case can absorb an error in `w`. |
| `wavy.ini` | **oblique analytic geometry.** ibm.f90's wavy bottom wall (`0.010 … 0.035` over `ly = 0.25`, `ny = 32`) — nowhere grid-aligned, so every cut face carries a direction cosine `a < 1` and the level-set weight, the grazing guard and the masked convection are all exercised on a real surface. It is also the second dwall producer (`fill_body_distance_analytic`, not the case file's tiles). Used for conservation, the config guards, the 2:1 precondition and determinism. |

## Results (2026-08-27/28, CPU `build_cpu`, GPU `build_gpu_nofma`)

PROVENANCE: every number below was produced with all four binaries
(`build_{cpu,gpu}`, `build_{cpu,gpu}_nofma`) rebuilt from the FINAL source.
That is worth stating because the Peclet-rate fix below changes `dt`: a gate
run against a stale nofma pair silently gates a different code, and the
`det` group in particular compares three binaries against each other, so it
passes happily when all three are equally stale. The same lesson as the
`~/s5b_ref_binaries` landmine in CLAUDE.md, one step earlier in the chain.

### (1) the exact two-material profile is a fixed point — **PASS**

Cut position swept through a full cell, `δ_L/h ∈ {0.05, 0.20, 0.35, 0.50,
0.65, 0.80, 0.95}` (the arm here also happens to sit **on a block boundary**,
so the ghost-inclusive tiles are exercised), `κ_s ∈ {10⁻², 1, 10, 10³}` with
`C_s = κ_s`. 28 runs, 500 steps from the seeded exact profile:

```
max |theta - exact|  <= 9.66e-15   over all 28 (w, kappa_s) pairs
                        0.0        exactly, for every kappa_s = 1 case
```

The gate STARTS at the fixed point rather than converging to it: any error in
`w` changes `k_face` and the profile leaves immediately, at a rate set by the
error instead of by the slowest eigenmode. At `κ_s = 10⁻²` a cold start needs
`O(10⁵)` steps, which is what the `converge` group pays.

### (1a) cold start reaches the same profile — **PASS**

`δ_L/h = 0.5`, `κ_s ∈ {10⁻², 1, 10}`, 200 000 steps from a uniform field:

```
kappa_s = 0.01   max|theta - exact| = 2.387e-14   residual (last 10k steps) = 0.0
kappa_s = 1      max|theta - exact| = 7.438e-15   residual                  = 0.0
kappa_s = 10     max|theta - exact| = 4.774e-15   residual                  = 0.0
```

The residual is reported alongside the error, so "converged" is a
measurement, not an assertion — here it is exactly zero.

### (1b) the level-set weight itself — **PASS**

Rebuilt from the case file alone — `φ = ±dwall_blocks` signed by
`coef_p_blocks` — and compared with the analytic cut position of the STL
plane, i.e. without the solver in the loop:

```
delta_L/h = 0.05   72 cut y-arms   max |w - w_exact| = 2.226e-14
delta_L/h = 0.50   72 cut y-arms   max |w - w_exact| = 0.000e+00
delta_L/h = 0.95   72 cut y-arms   max |w - w_exact| = 2.220e-14
```

The residual is the floating-point floor of the STL point-triangle distance,
not a scheme error. **It scales with the STL box's size**: padding the slab to
±4 instead of ±0.5 loses ~64× more in the `d²` cancellation (measured
`|dwall − exact| = 3.2e-13` vs `~1e-14` at a cut 0.003 from the plane), which
is why `make_slab_stl.py` keeps the box tight.

### (1c) the steady state does not depend on the solid capacity — **PASS**

`κ_s = 10`, `C_s ∈ {0.5, 1, 8}`, 150 000 steps: the same exact profile to
`8.85e-14` / `4.37e-14` / `4.88e-15`. Capacity is irrelevant at steady state
(strategy doc §6), so this also pins the capacity division to the flux
divergence rather than into the conductivity.

### (1d) the Peclet limiter over materials — **PASS**

`α_s/α_f = κ_s/C_s = 200`. With the limiter on the run is bounded
(`max|θ| = 8.69e-01` after 4000 steps); the control with `pecletmax = 0` and
`dt = 10⁻²` goes to NaN, so the gate is testing something.

### (1e) contact resistance — **PASS**

`R_c` adds to the series resistance of the cut face, so the steady solution
stays piecewise linear with a **jump `q·R_c`** at the interface, and stays an
exact fixed point. `κ_s = 10`, `δ_L/h = 0.35`:

```
R_c = 0     q = 1.276426007180   max|theta - exact| = 6.94e-17
R_c = 0.25  q = 0.967644390686   max|theta - exact| = 0.00e+00
R_c = 4     q = 0.209054680865   max|theta - exact| = 4.34e-19
```

`R_c = 0` reproducing the plain harmonic mean is the arithmetic identity that
makes the key safe to ship.

### (2) the two limits — **PASS**

`κ_s → ∞` against the S3 `dirichlet` mode and `κ_s → 0` against `adiabatic`.

```
kappa_s ->  inf   max|theta_conj - theta_dirichlet| over the fluid
  1e3   3.193e-04     1e5   3.194e-06     1e7   3.191e-08
kappa_s ->  0     max|theta_conj - theta_adiabatic| over the fluid
  1e-3  2.866e-03     1e-5  2.875e-05     1e-7  2.875e-07

interface flux (flux_limit.py), q_Dirichlet = 1.333333333333:
  kappa_s = 1e3   q = 1.332889036988   |q - q_inf| = 4.443e-04
  kappa_s = 1e5   q = 1.333328888904   |q - q_inf| = 4.444e-06   order 1.000
  kappa_s = 1e7   q = 1.333333288889   |q - q_inf| = 4.444e-08   order 1.000
  kappa_s = 1e-3  |q|/kappa_s = 3.988036
  kappa_s = 1e-5  |q|/kappa_s = 3.999880                          order 0.999
  kappa_s = 1e-7  |q|/kappa_s = 3.999999                          order 1.000
```

`|q|/κ_s → 4` is the closed form (`q = 1/(y_wall/κ_s + L − y_wall) → κ_s/0.25`).

**Better than §10 asked for.** The gate is written as a limit, not an identity
— the S3 modes place their effective boundary on the staircase while the
conjugate interface sits at its true position, so an O(h) floor was expected.
There is none: the difference from the `dirichlet` twin falls like `1/κ_s`
straight down to `3e-8`, five decades below any `h²` here. That is §5's
algebraic claim confirmed — the S3 `dirichlet` mode, with its second-order
graded coefficient, **is** Luchini's λ, and `κ_s → ∞` of the cut-face
coefficient reproduces it rather than merely approaching it.

LOAD-BEARING in the low group: `C_s = κ_s` (so `α_s = 1`). With `C_s = 1` the
solid has `α_s = κ_s` and equilibrates over `t ~ L²/α_s = 6e3`, so a
30k-step run measures the solid's **charging** flux — measured `|q|` came out
7× the closed form and the rate read 0.58 instead of 1. The physics was
right; the run simply was not steady.

### (3) conservation — **PASS**

Insulated composite box on the oblique wavy wall, flow ON (so the masked
convective flux is part of the statement), step initial condition (fluid 1,
solid 0). The capacity map in the checker is the analytic wavy wall, so it is
independent of the solver's own marker:

```
sum(C theta dV):  5.7128906250000007e-02 -> 5.7128906250000000e-02
drift = -6.94e-18    relative = 1.22e-16      (100 steps)
```

### (3b) config guards — **PASS** (5/5 rejected)

`remove_solid = true`; `ibm_value` with `conjugate`; a `solid_*` key on a
non-conjugate scalar; no immersed body; `solid_k <= 0`.

### (3c) the 2:1 precondition — **PASS**

The cut-face coefficient is a same-level two-point arm, so a cut face may
never sit on a coarse/fine block face; `refine_body`'s one-block 26-neighbour
buffer is what guarantees that, and the solver **checks** it.

```
hand-placed refinement box across the wall : rejected (224 bad faces)
refine_body + keep_buried                  : runs (5984 solid cells, 1216 cut faces)
refine_body without keep_buried            : rejected
```

NOTE for anyone writing a similar negative test: the box has to put a 2:1 face
**through** the surface, and this surface is nearly horizontal — an x-plane
interface slicing the whole domain carries no cut arm at all (measured: 0 bad
faces). The working form refines only the bottom block row so a 2:1 *y*-face
lands inside the wall's height range.

### (3d) the statistics branch — **SMOKE ONLY at C1**

`scalar_stats.f90`'s y-face diffusivity and convective mask now take the same
conjugate branch, with the same helper, as the transport kernel — that is the
invariant S4 exists for: a row's `J` must be the flux the kernel *applied*.
C1 checks only that the branch executes on CPU **and** GPU and that every row
is finite (it does; 4 datasets, 0 non-finite on both). The **quantitative**
gate on the conjugate flux columns is C3's Nusselt increment, which is where
`docs/next_session_conjugate.md` §10 puts it. Do not read this row as more
than it says.

### (4) determinism — **PASS**

1 rank == 4 ranks and CPU == GPU at **tolerance 0** on both geometry paths
(`un vn wn pn theta`, all `max_abs = 0`): the analytic wavy case and a
prepared slab case.

### (5) the bit-exactness protocol — **PASS**

```
REF=~/s5c_ref_binaries/moby_solve_cpu_nofma MODE=cpu ../scalar/run_bitexact.sh     # ALL PASS
REF=~/s5c_ref_binaries/moby_solve_gpu_nofma MODE=gpu ../scalar/run_bitexact.sh     # ALL PASS
REF=~/s5c_ref_binaries/moby_solve_cpu_nofma MODE=cpu ../scalar/run_bitexact_s3.sh  # ALL PASS
REF=~/s5c_ref_binaries/moby_solve_gpu_nofma MODE=gpu ../scalar/run_bitexact_s3.sh  # ALL PASS
```

`max_abs 0` on every case and every dataset, CPU and GPU: the 7-case standard
suite (`[scalar] count = 0`) and the 9-case scalar suite (`ibm_wall` =
`dirichlet` / `adiabatic`). This is by construction — a conjugate-free run
takes neither the `conjug` branch nor the geometry build, `sc%phi` is a 1-cell
dummy, and `precompute_peclet_rate` is untouched.

## What gate 1 taught about the time step

The cut-position sweep produced NaNs twice before the explicit limit was
right, and neither cause is visible from a per-material argument. Both are
recorded in `scalar_conjugate_peclet_rate`'s header:

1. **The limit is not `max` over materials of `α = κ/C`.** A cut face carries
   `k_face` up to `max(κ_L, κ_R)` — that bound *is* strategy doc §7's argument
   that the interface is not stiff — but it feeds the cell on the **other**
   side, whose capacity belongs to the other material. A fluid cell against a
   `κ_s = 1000` solid therefore sees 1000× the fluid rate even when
   `α_s = α_f` exactly. Building the rate from the **actual** face
   coefficients is both correct and far less conservative than
   `max(κ)/min(C)`: at `w = ½` and `κ_s = 1000` the true penalty is 2×.

2. **A cut cell pays the Gershgorin factor the uniform interior never does.**
   The spectral radius is bounded by `2 A_ii/C_i`, not `A_ii/C_i`, and the
   shipped `pecletmax` convention is a factor ~1.9 short of that — uniform
   runs survive only because their extreme modes are never excited. An
   isolated cut cell's row is strongly asymmetric, so its worst mode is local
   and *is* attained. Measured: `(κ_s, C_s) = (0.01, 0.01)` at `w = 0.95`
   blows up at `pecletmax = 0.3` and is stable at 0.2; `(1000, 1000)` at
   `w = 0.80` blows up at 0.4 and is stable at 0.2. Doubling the rate at cut
   cells makes the nominal 0.4 behave as the measured-stable 0.2 in both,
   which is exactly the factor Gershgorin predicts. (This RK3's real-axis
   stability limit is 2.5 — computed from `rk_alpha`/`rk_beta` directly.)
