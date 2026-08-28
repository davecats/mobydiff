# Conjugate heat transfer at an immersed interface — C1 and C2 gates

Gates for increments **C1** and **C2** of `docs/next_session_conjugate.md`.

**C1** is the baseline conjugate scheme, in which the solid stops being a
boundary condition and becomes a real unknown, and a face whose two cell
centres straddle the interface takes the **distance-weighted harmonic mean**
of the two materials' diffusivities, built on the level-set fraction
`w = φ_L/(φ_L − φ_R)`.

**C2** is the measurement of what that baseline gives up — it drops the
tangential term of the exact cut-face flux — and of what putting the term
back costs. Its gates are in [C2 below](#c2--measuring-the-tangential-term);
`./run_gates_c2.sh`.

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

---

## C2 — measuring the tangential term

**VERDICT: the correction ships DISABLED; the indicator ships enabled.** The
term the C1 baseline drops costs ~3.6 % rms of the local interface flux at the
DNS-like tangential ratio and moderate contrast, and §3's closed form for it is
confirmed per face to every digit. Putting it back helps only above a measured
crossover ratio (`r* ≈ 0.027` at κ_s = 10 on a plane; a LOCAL ratio ≈ 0.1 on a
body), makes the SOLUTION worse even where it makes the face flux exact
(because a finite-difference cell balance wants the face AVERAGE, not the
midpoint value), and — the decisive one — is **40–65× worse than dropping the
term** at high contrast on a curved interface, where `|∇_tT| ∝ 2/(1 + κ_s)`
vanishes in the isothermal limit while the discrete estimate of it does not.
Details below, in gate order.

```bash
./run_gates_c2.sh [flux|indicator|bvp|cylinder|dt|c1|all]
```

Same environment as the C1 driver. New helpers:

| file | role |
|---|---|
| `make_geometry_stl.py` | ASCII STLs for the two C2 geometries: a **tilted half-space** (a rotated box whose top face is the interface plane) and a **faceted cylinder**. Both are padded out of range, and both stay TIGHT for the same reason `make_slab_stl.py` does — the BVH distance loses precision quadratically in the vertex magnitude. |
| `check_oblique.py` | `flux` (the measurement) and `field` (the BVP). Rebuilds `φ` from the case file's own `dwall_blocks`/`coef_p_blocks` and evaluates the scheme's face quantities on the ANALYTIC field, so it measures the scheme rather than a transient — and needs no solver at all. |
| `check_cylinder.py` | gate 2, on the same face machinery. |
| `seed_manufactured.py` | writes either manufactured solution into a snapshot (the `seed_slab_ic.py` idiom). |

and one new solver-side entry point, `check_oblique.py residual`, which
evaluates the LOCAL TRUNCATION ERROR of both schemes on the exact field — the
measurement that explains gate 1b.

`oblique.ini` is the shared template (`@THETA@`-free: the geometry lives in the
STL; the ini carries `@NX@ @NY@ @NZ@ @LX@ @LY@ @LZ@ @NB@ @KAPPA@ @CAP@ @GX@
@GY@ @DT@ @IND@ @TANG@ …`).

### The manufactured solutions

**Oblique plane** (gate 1). With `n` the plane's unit normal, `t` the in-plane
tangential direction and `ξ = n·(x − x₀)`,

```
T = q_n ξ/k(side) + A (t·x),     k = 1 (fluid) | κ_s (solid)
```

is an exact conjugate solution for any `κ_s` — linear in each material,
continuous, and with `k ∂_nT = q_n` continuous across the interface. The
**controlled ratio** of §3 is `|∇_tT|/|∂_nT| = A/q_n = r` on the fluid side.

**Cylindrical shell** (gate 2). A solid cylinder of radius `a` with a uniform
volumetric source `S` gives the exact radial (log) solution and a constant
`q_n = −S a/2`. Because the field is radial, **`s_t = 0` identically**: the
correction must be inert, and what is left is the level-set weight's curvature
error.

### (1) oblique plane: the interface-flux error against the ratio, and against h — **PASS**

Per cut face, from the case file's `φ` and the analytic field (no solver in the
loop), the relative error of the **interface-NORMAL flux** each scheme implies:

```
q_n^base = (T_R − T_L)/(a R)              q_n^corr = (T_R − T_L − h_d s_t)/(a R)
```

θ = 30°, `h ∈ {1/64, 1/128, 1/256}`, 404–1616 cut faces. **rms over cut faces**:

| κ_s | r | baseline | corrected | κ_s | r | baseline | corrected |
|---|---|---|---|---|---|---|---|
| 10 | 0 | **7e-11** | 9.6e-2 | 10³ | 0 | **4e-10** | 1.4e-1 |
| 10 | 0.01 | **3.6e-2** | 9.6e-2 | 10³ | 0.01 | 1.3e-1 … 5.6e-1 | **1.4e-1** |
| 10 | 0.1 | 3.6e-1 | **9.6e-2** | 10³ | 0.1 | 1.3 … 5.6 | **1.4e-1** |
| 10 | 1 | 3.6 | **9.6e-2** | 10³ | 1 | 13 … 56 | **1.4e-1** |
| 10 | ∞ | 3.6e-1 (flux) | **1e-11** | 10³ | ∞ | 3.6e-1 (flux) | **1e-11** |

Four things this says, all of them load-bearing for the verdict:

1. **§3's closed form is exact, per face.** The measured baseline error and
   `e/(1 − e)` agree to every digit printed (`1.5080e+00` against
   `1.5080e+00` at κ_s = 10, r = 0.1; `9.9029e-01` against `9.9029e-01` at
   κ_s = 10³, r = 0.01). The prediction is not an estimate.
2. **The baseline error is strictly linear in r and INDEPENDENT of h** — §3's
   "set by the flow, not by `h`", now measured over a 4× refinement:

   ```
   κ_s = 10, r = 0.01   q_n rel err rms   3.559e-02  3.774e-02  3.586e-02   order -0.08 +0.07
                        face flux rms     3.284e-02  3.331e-02  3.318e-02   order -0.02 +0.01
                        e_face rms        3.857e-02  4.138e-02  3.893e-02   order -0.10 +0.09
   ```

   It cannot converge: `e_face = h_d s_t/(T_R − T_L)` is a RATIO OF TWO
   GRADIENTS and both scale like `h`. So §10's expectation that the baseline
   "stalls near first order" is optimistic by a full order — in the usual
   accounting the cut-face flux is **zeroth order**, i.e. inconsistent
   (`F_num − F_exact = s_t(k_face − k_loc)`, O(1) absolute). What makes the
   scheme usable is not consistency but CONSERVATION plus sign alternation:
   `k_loc` is a staircase in `w` that flips as the interface sweeps past the
   face midpoint, so neighbouring faces carry opposite-signed errors that
   cancel in the cell balances. The one genuine first order in this suite is
   a DIFFERENT error — gate 2's curvature error in `w`.
3. **The conductivity contrast amplifies it.** §3's `≈ tanθ |∇_tT|/|∂_nT|` is
   the equal-material estimate; the normal term it competes against is
   `a q_n R`, and `R` is reduced by the high-conductivity leg. At r = 0.01 the
   rms error is 3.6 % at κ_s = 10 but 13–56 % at κ_s = 10³.
4. **The corrected error is flat in r** (9.6e-2 / 1.4e-1 rms) — it is no longer
   the tangential term but the residual bias of the discrete `s_t`. Hence a
   CROSSOVER ratio, which is the number the ship decision turns on:

   ```
   r* = 0.027            (κ_s = 10,  all three h)
   r* = 0.003 … 0.011    (κ_s = 10³)
   ```

   Below `r*` the correction makes the interface flux WORSE than dropping the
   term; above it, better. §3's DNS-like ratio `r ~ 10⁻²` sits **at or below**
   `r*`.

**The note's own construction is not what ships, and this is why.** The raw
projection of §7.3's construction 1 — ordinary central differences, no
straddling correction — is reported alongside (`RAW s_t`, and the third
column of `c2_flux.dat` set):

| κ_s | r | baseline rms | raw `s_t` rms | shipped rms |
|---|---|---|---|---|
| 10 | 0.01 | **3.6e-2** | 3.6e-1 | 9.6e-2 |
| 10³ | 0.01 | 1.3e-1 | 1.7 | **1.4e-1** |

It is 10× worse than dropping the term. The bias is `μ n_d (½ − w)(1 − n_d²)`
with `μ` the jump in the normal derivative, so it GROWS with the contrast the
correction is supposed to help with — the note's argument that the projection
removes the jump holds for the two tangential differences (which weight the
sides ½, ½) but not for the arm difference (which weights them `w, 1−w`).
What ships removes it in closed form; see `conjugate_tangential`.

The `r = ∞` row is the sharp positive statement: for a purely tangential field
the corrected FACE flux is exact to `1e-11` (the STL distance floor divided by
`h`) where the baseline is 36 % rms.

θ = 45° is also run (`c2_flux.dat`) but is **lattice-degenerate**: a 45° plane
on a uniform grid gives every cut face of a direction the same `(w, a)`, so
`max = rms` and the numbers track the offset `mod h` rather than a statistic.
It is kept because §10 asks for it; the 30° numbers are the ones to read.
LANDMINE recorded on the way: with `y₀` ON the lattice, a whole diagonal of
cell centres lands EXACTLY on a 45° plane, `w` collapses to 0 or 1, and the
measured flux is off by the full contrast (`κ_s − 1`). That is a degenerate
geometry, not a scheme error — the driver offsets `y₀` by 0.5117.

### (1a) the solver's own `e_face` against the checker's — **PASS**

The indicator is a deliverable, so it is shown to compute what it claims. On
the `r = ∞` field (the one whose Neumann data is exactly consistent on every
face, so the solver's ghosts ARE the analytic field), `e_face = 1` exactly by
hand — `s_t = ∇_dT` when `∇T` is tangential:

```
solver : conjugate indicator step 1 'theta': cut faces 404
         max|e_face| = 1.0000E+00  rms = 1.0000E+00  (0 below the denominator floor)
checker: e_face                    max = 1.0000e+00  rms = 1.0000e+00
```

### (1b) the one global statement: the `r = ∞` BVP — **PASS**, and the number that decides

Every FINITE ratio has a jump in `∇T`, so its exact boundary data differs
BETWEEN THE TWO MATERIALS on the two domain faces the plane cuts, and the
solver's per-face rows are constants — a finite-ratio BVP with an analytic
solution is therefore not available, and adding per-point boundary data for a
gate was rejected as solver surface no production case would use. `r = ∞` is
the one ratio with no jump: `∇T` is then the same in both materials, all six
faces carry the exact constant Neumann value, and the BVP is well posed up to
the pure-Neumann constant (which the checker removes).

Seeded at the exact solution, relaxed to the discrete steady state
(`κ_s = 10`, `n = 32`, `t = 0.52`, the error at half the time given as the
convergence measurement, not a claim):

```
tangential_correction = false   L2 = 3.327519e-03   Linf = 1.124148e-02   (half time 3.263e-03)
tangential_correction = true    L2 = 8.744530e-03   Linf = 3.949333e-02   (half time 8.536e-03)
```

**The correction is 2.6× WORSE — at the ratio most favourable to it, and with
its face flux exact to 1e-11.** The error is not localized either (rms in the
band `|ξ| < 2h` vs away: `2.5e-3 / 3.3e-3` baseline, `9.5e-3 / 8.4e-3`
corrected): the cut faces act as a source and the response fills the domain.

Why, measured rather than argued (`check_oblique.py residual` — the discrete
flux divergence of the EXACT field, whose true value is zero):

```
C1 baseline   k_face   cut-cell |div| max = 2.5319e+02   rms = 1.2987e+02   (uncut cells max 6.82e-12)
C2 shipped    k_loc    cut-cell |div| max = 2.8800e+02   rms = 1.7406e+02   (uncut cells max 6.82e-12)
PROPOSED      k_area   cut-cell |div| max = 1.5883e-07   rms = 1.3324e-08   (uncut cells max 6.82e-12)
```

**Making the pointwise face flux exact makes the cell balance 34 % worse** —
and the third row says what would fix it; see "the way out" at the end of
this section.
§4 derives `F` as the flux AT THE FACE MIDPOINT, and that is what the
escalation makes exact; but a finite-difference divergence sums face fluxes as
if they were face AVERAGES, and at a cut face the exact average is
`f k_f + (1−f) k_s` with `f` the face's fluid AREA fraction. `k_face`, being
smooth in `w`, approximates that better than the exact midpoint value
`k_loc`, which is a staircase in `w`. An escalation that would actually pay
has to correct the face average — which needs the area fraction, i.e. exactly
the per-face geometric data the C1 baseline exists to avoid.

(One grid. The statement is a comparison, not an order — gate 1 carries the
h dependence — and finer grids are priced out: `dt = 0.4 h²/κ_s` against a
response that spreads over `L²/α_f = 1` costs 1e5 steps per leg at `n = 64`.)

### (2) cylindrical shell: the curvature error in `w` — **PASS**

`a = 0.25`, `S = 1`, `h ∈ {1/64, 1/128, 1/256}` (`κ_curv h = 0.0625 … 0.0156`),
512–2048 cut faces. `s_t = 0` exactly here, so this isolates `w`:

| h | baseline rms | baseline max | corrected rms | raw `s_t` rms |
|---|---|---|---|---|
| 1/64 | 1.61e-2 | 2.59e-2 | 4.11e-2 | 8.35e-2 |
| 1/128 | 8.38e-3 | 1.47e-2 | 8.82e-2 | 2.97e-1 |
| 1/256 | 4.09e-3 | 7.21e-3 | 5.80e-2 | 1.02e-1 |

The baseline converges at **order 1.0** (1.92, 2.05 in the rms ratios), which
is exactly the `O(κ_curv h/a)` of §8 — the curvature caveat, measured. The
numbers are the same to three digits at κ_s = 10 and 10³, so the error is
geometric, not material, as predicted.

**The correction does not belong here.** It has no tangential term to supply,
and its discrete `s_t` is not zero on a curved interface (the local model
behind it is a PLANE), so it injects 4–9 % where the baseline had 0.4–1.6 %,
and it destroys the first-order convergence.

This gate says nothing about a curved interface that DOES carry a tangential
gradient, because a radial field has none. That case — the cylinder in a
uniform far-field gradient — is measured at the end of this section, and it is
where the verdict's strongest argument comes from.

### (3) every C1 gate with the correction ON — **PASS**

```bash
TANG=true ./run_gates_c1.sh all       # or ./run_gates_c2.sh c1
```

`TANG=true` injects `tangential_correction = true` into the two C1 TEMPLATES
(`slab.ini`, `wavy.ini`), so the whole suite re-runs with the correction
active. The grid-aligned cases have `s_t = 0` BY CONSTRUCTION — a plane wall
normal to `e_y` gives `∇φ = e_y`, so the projection removes everything — and
they must therefore come back to round-off, not merely close:


```
slab, all 28 (w, κ_s) pairs   max|θ − exact| ≤ 9.658940e-15   (C1 recorded 9.66e-15)
cold start, κ_s = 1           max|θ − exact| =  7.438494e-15   (C1 recorded 7.438e-15)
                              residual over the last write interval = 0.000e+00
```

i.e. **C1's own numbers, digit for digit** — as do the `κ_s → ∞ / → 0` limits
(`2.866401e-03` at κ_s = 10⁻³ against C1's `2.866e-03`; `q = 1.332889036988`,
`|q − q_∞| = 4.443e-04`, rate **order 1.000**, all identical) and the level-set
weight (`2.226e-14 / 0.000e+00 / 2.220e-14` over the three cut positions).

The oblique `wavy.ini` half of the suite is a genuine re-measurement rather
than a reproduction, and it also comes back unchanged with the correction
active:

```
sum(C theta dV)  drift = 6.939e-18   relative = 1.215e-16   (C1: -6.94e-18 / 1.22e-16)
config guards                        5/5 rejected
2:1 precondition                     hand-placed box rejected; refine_body + keep_buried
                                     runs (5984 solid cells, 1216 cut faces); without
                                     keep_buried rejected
scalar_stats branch                  4 datasets, 0 non-finite, CPU and GPU
```

**`ALL C1 GATES PASS`, 61 checks, zero failures** — the whole suite, with the
correction on.

### determinism with the correction ON — **PASS**

On the nofma pair, `wavy.ini`, 20 steps, `un vn wn pn theta`:

```
1 rank vs 4 ranks   max_abs 0     (every dataset)
CPU vs GPU          max_abs 0     (every dataset)
```

`max_abs 0` across ranks is the real check on the correction's stencil: it
reaches the EDGE ghosts of the 3³ neighbourhood, and those are filled by the
26-neighbour exchange, not by `apply_scalar_bc`. The **indicator** is
bit-identical too — `max|e_face| = 9.1116E-01`, `rms = 3.5031E-01` at step 0
from all three runs — so its reduction is deterministic across ranks and
devices as well.

(That `rms ≈ 0.3` on `wavy.ini` is worth reading as what it is: the wavy wall
there is resolved by 1.3–4.5 cells, so the indicator is reporting a warning
about the CASE, not about the scheme.)

### (4) the time-step penalty of the correction — **PASS**, and a C1 finding

The correction is an explicit spatial operator at the cut faces whose size
grows with the contrast, so `scalar_conjugate_peclet_rate` adds its Gershgorin
row sum. Oblique plane, θ = 30°, `κ_s = 10³`, `h = 1/64`, 2000 steps from the
seeded field:

```
pecletmax = 0.4  tangential_correction = false  dt = 5.750066e-08  max|theta| = NaN
pecletmax = 0.4  tangential_correction = true   dt = 4.271588e-08  max|theta| = 0.676
pecletmax = 0.2  tangential_correction = false  dt = 2.875033e-08  max|theta| = 0.679
pecletmax = 0.2  tangential_correction = true   dt = 2.135794e-08  max|theta| = 0.677
```

**The correction costs 26 % of the time step** at κ_s = 10³ (`dt` ratio 0.743),
far less than its worst-case row sum suggests — the geometric factors
`(1 − n_d²)` and the tangential metric keep it well below the principal term.
That is a small price; it is not what decides the verdict.

**FOUND HERE, AND IT IS A C1 PROPERTY, NOT A C2 ONE: the shipped
`pecletmax = 0.4` is marginal at an OBLIQUE, high-contrast interface.** The
baseline run above goes to NaN at 0.4 and is stable at 0.2 (measured, and at
0.1). The arithmetic says why: `scalar_conjugate_peclet_rate` reports
`diag/(3 C)` at a cut cell, so the step it grants is `dt = pecletmax·3C/diag`,
while Gershgorin plus this RK3's real-axis limit of 2.5 allows
`dt ≤ 2.5 C/(2 diag) = 1.25 C/diag`. At `pecletmax = 0.4` the step is
`1.2 C/diag` — **96 % of the bound**. C1's own gates never saw it because
their interface is grid-aligned, where the bound is not attained; obliquity
excites it.

Nothing was changed here: touching the rate moves `dt` for every conjugate
run and would invalidate C1's recorded numbers, which is a C3-sized decision.
**Until then: run a conjugate case with a strongly oblique interface and
`κ_s ≳ 10²` at `pecletmax ≤ 0.2`.** The correction happens to mask this (its
own rate contribution shrinks `dt` by 1.35×, which is why the `true` run
survives at 0.4) — do not read that as the correction being more stable.

### (5) the bit-exactness protocol — **PASS**

```
REF=~/s5c_ref_binaries/moby_solve_cpu_nofma MODE=cpu ../scalar/run_bitexact.sh     # 7/7
REF=~/s5c_ref_binaries/moby_solve_gpu_nofma MODE=gpu ../scalar/run_bitexact.sh     # 7/7
REF=~/s5c_ref_binaries/moby_solve_cpu_nofma MODE=cpu ../scalar/run_bitexact_s3.sh  # 9/9
REF=~/s5c_ref_binaries/moby_solve_gpu_nofma MODE=gpu ../scalar/run_bitexact_s3.sh  # 9/9
```

**32 case-runs, every one `max_abs 0` on every dataset, zero failures**, CPU
and GPU: the 7-case standard suite (`[scalar] count = 0`) and the 9-case
scalar suite (`ibm_wall` = `dirichlet` / `adiabatic`). By construction, and
the construction is worth naming: the correction is added as its OWN flux
divergence after the C1 term (`if (tang) diff = diff + …`) rather than folded
into the fused expression, so a run without it executes the same arithmetic it
did before — no `+ 0.0`, no FMA-reassociation argument. The indicator is a
separate routine that is never called unless asked for.

…and the S-suite physics groups, which the same protocol asks for, all with
`rc = 0` and no failures:

```
BIN=build_cpu_nofma/moby_solve GBIN=build_gpu_nofma/moby_solve ./run_gates.sh   #  8 checks
./run_gates_s2.sh   #  6   (scalar_test: ALL PASS)
./run_gates_s3.sh   # 36   (S3 gates (all): ALL PASS)
./run_gates_s4.sh   # 33
./run_gates_s5.sh   #  8   (S5a gates: ALL PASS)
```

LANDMINE honoured: `run_gates.sh`'s `det` group compares CPU against GPU at
TOLERANCE 0, so it gets the nofma pair.

**PROVENANCE.** Every number in this C2 section was produced with binaries
built from the FINAL source, and that was checked rather than assumed: after
the last edit all four (`build_{cpu,gpu}`, `build_{cpu,gpu}_nofma`) were
rebuilt and the protocol re-run on a representative case —
`min_channel` back to `max_abs 0` against `~/s5c_ref_binaries`, and the
correction-ON determinism pair (1 vs 4 ranks, CPU vs GPU) back to `max_abs 0`
with the indicator reporting `9.1116E-01 / 3.5031E-01` identically from all
three. The C1 lesson one step earlier in the chain: a gate run against a stale
binary silently gates different code.

### the way out, measured — what a working escalation would have to be

The three schemes differ in **one number**: the coefficient multiplying `s_t`
in the face flux. For the manufactured field, `k e_d·∇T = a q_n + k s_t` on
each side and `a q_n` is CONTINUOUS, so

```
F = a q_n + K s_t          K = k_face   C1 baseline (harmonic in w)
                           K = k_loc    C2 as shipped (the midpoint, a staircase in w)
                           K = k_area   the exact face AVERAGE:  f k_f + (1-f) k_s
```

with `f` the face's fluid AREA fraction. `k_area` is the one the cell balance
wants, and the argument is two lines: with exact face averages the discrete
divergence IS the exact surface integral `∮ k∇T·n dS`, which equals
`∫ ∇·(k∇T) dV = 0`. (Equivalently: apply the divergence theorem to the FLUID
part of the cell — the interface patch contributes `(G·n)|Γ∩cell| = 0`.)

**It costs no new data.** For a plane, `f` follows in closed form from `φ` at
the face centre — the mean of the two cell values, exact for a plane — and the
in-plane components of the face-centred `∇φ`, both of which `s_t` already
forms (`face_area_fraction`, the 2D plane-in-rectangle form; validated against
brute force to the quadrature error, 9e-4 at 400² samples).

**One thing has to change besides the multiplier**: it must be applied at every
face the interface CLIPS, not only where the two cell markers disagree. Those
are different sets — the interface can clip a face whose two centres agree —
and the difference is not small:

```
cut-cell |div| rms, exact field, h = 1/64, κ_s = 10, r = ∞
  k_area at marker-cut faces only    3.92e+01
  k_area at every face               1.33e-08      <- exact
```

Measured across ratios and contrasts (`h = 1/64`, θ = 30°, rms):

| κ_s | r | C1 `k_face` | C2 `k_loc` | `k_area` |
|---|---|---|---|---|
| 10 | ∞ | 1.30e+02 | 1.74e+02 | **1.33e-08** |
| 10 | 1 | 1.30e+02 | 1.72e+02 | **8.51e+00** |
| 10 | 0.1 | 1.30e+01 | 1.72e+01 | **8.51e+00** |
| 10 | 0.01 | **1.30e+00** | 7.27e+00 | 8.51e+00 |
| 10³ | ∞ | 1.92e+04 | 1.93e+04 | **1.55e-06** |
| 10³ | 1 | 1.92e+04 | 1.91e+04 | **1.21e+03** |
| 10³ | 0.01 | **1.92e+02** | 1.08e+03 | 1.21e+03 |

Read it as three statements:

1. `k_area` **removes the paradox**: it is exact wherever `s_t` is exact
   (`r = ∞`, both contrasts, ~10¹⁰ down), so "making the flux exact makes the
   divergence worse" was an artefact of the WRONG multiplier, not of the
   escalation idea.
2. At finite `r` it hits a **floor that does not depend on `r`** (8.5 at
   κ_s = 10, 1.2e3 at κ_s = 10³) — that floor is the discrete `s_t`'s own
   residual error, not the multiplier.
3. So the decision moves, but does not flip: the floor still crosses the
   baseline at `r ≈ 0.065` (both contrasts — `8.51/1.30e2` and
   `1.21e3/1.92e4`), which is **still above the DNS-like 10⁻²**. `k_area`
   alone would not earn the correction its default-on.

What the floor is, and the obvious next step: the de-bias in
`conjugate_tangential` models the straddle of the ARM difference only. The two
TANGENTIAL differences straddle as well wherever the interface crosses their
stencils, and the same closed-form argument applies to them, with their own
weights. That is the third correction, and it is the one that would decide
whether the escalation can ever pay at DNS ratios.

### …and the same fix on a CURVED interface, which is where it stops

The plane is the friendly case: `s_t` is exact there by construction at
`r = ∞`. The curved counterpart is the **cylinder in a uniform far-field
gradient** (`check_cylinder.py dipole`) —

```
T_out = -G cos t (r + β a²/r),  β = (1-κ)/(1+κ)        T_in = -G cos t γ r,  γ = 2/(1+κ)
```

— harmonic on both sides with no source, so the exact divergence is zero and
the truncation test applies unchanged; `[T] = 0` and `[k ∂_nT] = 0` hold at
`r = a`; and unlike the radial log solution it **carries a tangential
gradient** whose local ratio `|∇_tT|/|∂_nT| = |tan t|` sweeps `0 → ∞` around
the body, so one run samples every ratio the decision turns on. The nearest
interface point of a face centre is its radial projection, so the face
centre's polar angle is the exact angle to compare at — no O(h) error enters
the reference.

**(a) Where on a body the correction pays** (κ_s = 10, h = 1/64, 512 cut
faces; interface-flux error rms, binned by the LOCAL ratio):

| ratio band | faces | baseline | corrected | winner |
|---|---|---|---|---|
| 0.0 – 0.1 | 32 | **1.26e-2** | 1.37e-2 | baseline |
| 0.1 – 0.3 | 56 | 9.58e-2 | **1.95e-2** | corrected |
| 0.3 – 1.0 | 168 | 1.31e-1 | **4.92e-2** | corrected |
| 1.0 – 3.0 | 160 | 2.88e-1 | **4.68e-2** | corrected |
| 3.0 – ∞ | 96 | 1.08e+0 | **2.32e-2** | corrected |

Position-resolved, the picture is *friendlier* to the correction than the
plane's global-ratio sweep: the crossover sits at a local ratio ≈ 0.1, and
most of a body's surface is above it.

**(b) But `s_t` itself is badly wrong on a curved interface**, and this is the
number that governs everything else (rms over cut faces, against the exact
`s_t`, which is known here):

| h | exact rms | raw projection | shipped de-bias |
|---|---|---|---|
| 1/64 | 7.39e-2 | 58 % | **43 %** |
| 1/128 | 7.41e-2 | 77 % | **54 %** |
| 1/256 | 7.42e-2 | 57 % | **46 %** |

It does **not converge**. The de-bias is derived for a PLANE; curvature breaks
its model, and no refinement repairs that.

**(c) So `k_area` helps but no longer fixes** (cut-cell truncation rms; the
uncut background is the ordinary O(h²) Laplacian truncation):

| κ_s | h | `k_face` | `k_loc` | `k_area` | background |
|---|---|---|---|---|---|
| 10 | 1/64 | 1.66e+1 | 2.29e+1 | **1.01e+1** | 3.6e-2 |
| 10 | 1/128 | 3.07e+1 | 4.12e+1 | **2.28e+1** | 1.1e-2 |
| 10 | 1/256 | 6.70e+1 | 8.51e+1 | **4.14e+1** | 3.0e-3 |
| 10³ | 1/64 | **2.24e+1** | 8.55e+2 | 1.45e+3 | 4.4e-2 |
| 10³ | 1/256 | **9.40e+1** | 4.38e+3 | 6.25e+3 | 3.6e-3 |

At moderate contrast `k_area` is consistently ~1.6× better than the baseline
and ~2.2× better than the shipped `k_loc`. **At κ_s = 10³ both corrections are
40–65× WORSE than the baseline** — and the reason is physics, not arithmetic:

```
|∇_tT| at the interface = G |sin t| · 2/(1 + κ_s)     ->  0  as κ_s -> ∞
```

The isothermal limit **removes the tangential term the correction exists to
add**, while the discrete estimate of it stays O(5e-2) — and that noise is
then multiplied by `(k_loc − k_face) ≈ κ_s`. Measured: the exact `s_t` rms
falls 90× from κ_s = 10 to 10³ (7.4e-2 → 8.1e-4, matching `2/(1+κ)` to two
digits), while the estimate's error does not fall at all.

**So the fix works for the defect it targets and cannot rescue the feature.**
`k_area` is right, and provably so on a plane; the binding constraint is the
`s_t` estimate on curved interfaces, and at high contrast there is no signal
left to estimate. **None of it is implemented** — it is what the measurement
says a fix would have to be, recorded so the next increment starts from the
number rather than from the idea. `check_oblique.py residual` and
`check_cylinder.py dipole` report all three multipliers.
