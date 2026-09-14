# Conjugate heat transfer at an immersed interface — C1, C2 and C3 gates

Gates for increments **C1**, **C2** and **C3** of
`docs/next_session_conjugate.md`.

**C1** is the baseline conjugate scheme, in which the solid stops being a
boundary condition and becomes a real unknown, and a face whose two cell
centres straddle the interface takes the **distance-weighted harmonic mean**
of the two materials' diffusivities, built on the level-set fraction
`w = φ_L/(φ_L − φ_R)`.

**C2** is the measurement of what that baseline gives up — it drops the
tangential term of the exact cut-face flux — and of what putting the term
back costs. Its gates are in [C2 below](#c2--measuring-the-tangential-term);
`./run_gates_c2.sh`.

**C3** is the TRANSIENT (the fluid-fraction-weighted capacity), the interface
heat (the Nusselt diagnostic) and the time step. Its gates are in
[C3 below](#c3--the-fraction-weighted-capacity-the-nusselt-diagnostic-and-the-time-step);
`./run_gates_c3.sh`. **Read the time-step subsection first**: the cut-cell
share changed, so every conjugate number recorded before C3 was measured
against a different `dt`, and the C1/C2 tables above have been re-measured
rather than left standing.

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

---

## Stage 0 of the route-B escalation — a one-sided `s_t` CONVERGES

**KILL GATE 1 PASSES.** The plan is `docs/next_session_tangential.md`; this is
its first stage, and it overturns the finding that closed C2 — that the `s_t`
estimate on a curved interface cannot be repaired.

```bash
./run_gates_c2.sh stzero          # needs no solver; the case files above suffice
```

### What was measured

`s_t = e_d·∇_tT` and **T is continuous across Γ**, so its surface gradient is
single-valued: only `∂_nT` jumps, and either side estimates the same number.
The shipped estimator instead differences `T` ACROSS the interface and removes
the straddle in closed form — a de-bias derived for a PLANE, which is exactly
why it does not converge on a curved one. A same-side stencil never straddles,
so there is nothing to de-bias.

**It fits the existing one-deep halo**, which is the reason C2 did not try it:
`conjugate_tangential`'s DEVIATION comment rules out same-side least squares
because a full 3D fit needs the far cell's far neighbour. It is not needed. One
material layer gives the two coordinate-tangential derivatives within it, and
the missing third component closes against the side's own normal derivative
`∂_nT = q_n/κ`:

```
s_t = A q_n - B,    A = (1 - n_d^2)/(n_d kappa),   B = (n_1 t1 + n_2 t2)/n_d
```

with `q_n` from either the face balance (`flux`, linear in `s_t`, denominator
provably `>= 1`) or from demanding the two sides agree (`agree`, which
eliminates `q_n` outright: `q_n = (P_f - P_s)/[(1 - n_d^2)(1/k_f - 1/k_s)]`).
Both are in `check_oblique.onesided`.

### The numbers (rms error as % of the exact `s_t`; `check_st.py`)

The plane is the unit test, not the gate: the field is piecewise LINEAR and φ
is exact, so a one-sided estimator must return the geometry floor — **it does,
7.4e-10 absolute against a signal of 0.658**, where the shipped de-bias leaves
4.3 % and the raw projection 7.1 %, neither converging.

The gate is the curved interface. `multipole --mode 2` exists because **the
dipole's interior is exactly LINEAR** (`gamma r cos t = gamma x`), which would
flatter any solid-side score; at m = 2 neither side is linear.

| quadrupole, κ_s = 10 | h=1/64 | 1/128 | 1/256 | order |
|---|---|---|---|---|
| raw projection | 37.3 % | 51.2 % | 34.5 % | 0.06 |
| **shipped de-bias** | **36.0 %** | **25.6 %** | **21.8 %** | **0.37** |
| one-sided, fluid side | 42.3 % | 22.1 % | 11.1 % | 0.97 |
| **one-sided, solid side** | **1.96 %** | **1.04 %** | **0.44 %** | **1.08** |
| one-sided, `agree` | 11.6 % | 5.92 % | 2.97 % | 0.98 |

and at κ_s = 10³, where C2 concluded there is "no signal left to estimate":
shipped 4014/2800/2461 %, **solid side 2.99/1.45/0.95 %, order 0.83**.

**The side rule, swept over six decades** (m = 2, h = 1/256):

| κ_s | 0.001 | 0.01 | 0.1 | 0.5 | 2 | 10 | 1000 |
|---|---|---|---|---|---|---|---|
| shipped | 1.65 % | 1.65 % | 1.60 % | 1.05 % | 2.27 % | 21.8 % | 2461 % |
| fluid | 5.19 % | 3.23 % | 1.23 % | 0.41 % | 1.75 % | 11.1 % | 1164 % |
| **solid** | **0.65 %** | **0.64 %** | **0.59 %** | **0.50 %** | **0.38 %** | **0.44 %** | **0.95 %** |

The solid (body-INTERIOR) side is uniformly best and essentially
contrast-independent, and the fluid side's error is `~κ_s ×` the solid's at
high contrast — which is the whole of C2's defect 3. Neither conductivity nor
linearity explains it (it holds at κ_s < 1 and at m = 2); the working
hypothesis is that the interior harmonic `r^m` is the smoother of the two near
Γ, while the exterior carries the decaying `r^-m` multipole whose normal
derivatives grow. **That hypothesis is not yet tested on a CONCAVE body, where
"interior" and "solid" come apart** — do that before hardwiring the side.

### Two defects in the C2 gate found on the way

1. **`Dipole.ratio` returned `|tan t|`, the κ_s = 1 form.** The tangential and
   fluid-normal gradients at `r = a` carry `(1 + β)` and `(1 - β)`, and
   `(1+β)/(1-β) = 1/κ_s`, so the fluid-side ratio — the quantity `Plane`
   controls directly, and the one the whole decision is defined on — is
   `|tan t|/κ_s`. FIXED. **Consequence: the position-resolved crossover in
   table (a) above was binned a DECADE too high.** Its "local ratio ≈ 0.1"
   is really ≈ 0.01, i.e. at the DNS-like ratio rather than a decade above it.
   The C2 verdict was reached partly on that number.
2. **The dipole's linear interior** would have made any solid-side estimator
   look exact for free. `Multipole` (verified: `[T] = 0`, `[k ∂_nT] = 0`,
   `grad_outer` against finite differences) is the de-flattered replacement.


### Stage 1 — the multiplier, and where the escalation now stops

`k_area` + the stage-0 `s_t`, scored on the CUT-CELL TRUNCATION RESIDUAL (the
exact divergence is zero, so anything left is the scheme's; this is what the
cell balance sees, and it is the measurement that decided C2).

| case (h = 1/256 unless noted) | C1 baseline | `k_area` + shipped `s_t` | `k_area` + 1-sided | gain |
|---|---|---|---|---|
| plane κ_s=10, r=0.01 | 5.18 | 34.4 | **1.8e-6** | 2.8e6 |
| plane κ_s=10, r=1 | 518 | 34.4 | **1.8e-6** | 2.9e8 |
| plane κ_s=10, r=∞ | 518 | 3.9e-7 | **3.7e-7** | 1.4e9 |
| plane κ_s=10³, r=0.01 | 764 | 4875 | **1.4e-4** | 5.4e6 |
| quadrupole κ_s=10, n=64/128/256 | 7.6 / 13.5 / 31.0 | 5.5 / 10.8 / 19.2 | **1.07 / 1.40 / 0.92** | 7 / 10 / **34** |
| dipole κ_s=10, n=64/128/256 | 16.6 / 30.7 / 67.0 | 10.1 / 22.8 / 41.4 | **1.13 / 1.30 / 1.06** | 15 / 24 / **63** |
| quadrupole κ_s=10³, n=64/128/256 | **10.2 / 23.2 / 45.4** | 751 / 1500 / 2703 | 114 / 177 / 87 | 0.09 / 0.13 / 0.52 |
| dipole κ_s=10³, n=64/128/256 | **22.4 / 47.4 / 94.0** | 1448 / 3412 / 6247 | 146 / 185 / 130 | 0.15 / 0.26 / 0.72 |

Three statements, and the third is the one that matters:

1. **On a plane the scheme becomes EXACT** — at every ratio and every contrast,
   six to nine orders below the baseline. The `r`-independent floor the README
   attributed to "the discrete `s_t`'s own residual error" was exactly that,
   and removing it removes the floor completely. `k_area` was the right
   multiplier all along; it was waiting on an `s_t` that converged.
2. **On a curved interface at moderate contrast it is 7-63× better AND
   BOUNDED**, where the C1 baseline's cut-cell residual grows like `1/h`
   (7.6 → 13.5 → 31.0). The gain IMPROVES with refinement, which is the
   signature of removing a term the baseline does not converge on.
3. **At high contrast on a curved interface it is 1.4-11× WORSE**, and
   **`s_t` is no longer the cause**: with the SAME one-sided `s_t`, `k_loc`
   sits at baseline parity (13.7/22.4/48.7 against 10.2/23.2/45.4) while
   `k_area` is 2× above it. So the residual defect is in the MULTIPLIER on a
   curved high-contrast interface — a different defect from the one C2
   identified, and the single remaining blocker.

**What it is not**: the plane assumption inside `face_area_fraction` is real
(|f_plane − f_exact| rms 2.2e-2, max 0.25 against exact quadrature on the
n = 256 cylinder) but too small to explain the excess — it predicts ~2.7 of
the ~38 by which `k_area` exceeds `k_loc` at κ_s = 10³. **Unresolved; start
here.** The likely candidate is the derivation's own premise: `F_exact =
k_area s_t + n_d q_n` assumes `∂_dT` is piecewise CONSTANT over the face, and
its error carries the same `(κ_s − 1)` amplification.

**Two limiters were tried and neither is never-worse.** The free confidence
indicator does discriminate — spread/|s_t| is 0.11 where the correction wins
34× and 11.6 where it loses, a 100× separation — but converting that into a
rule does not close the gap: a hard gate needs a FITTED threshold (c = 0.3
recovers high contrast but discards good corrections at coarse h, and still
loses on one row), and parameter-free soft shrinkage
`s_t → sign(s_t) max(0, |s_t| − spread)` keeps the moderate-contrast win
(4-19×) but not the high-contrast loss. Both ship in the checker
(`area1sg`, `area1soft`, `MOBY_SPREAD_GATE`) as measurements, not proposals.

**GATE 2 VERDICT: tier 2 of the plan** — a large, real improvement with a
documented regime boundary in (contrast × curvature), not the never-worse
tier 1. A plane interface of any contrast, and a curved one up to
κ_s ~ 10, are solved. Do NOT implement in Fortran until statement 3 is
resolved: the Fortran cost is the same either way, and a multiplier that is
wrong at high contrast on a curved body would have to be gated on exactly the
quantity that is not yet understood.


### Stage 1b — which premise breaks, and it is none of the ones we suspected

Stage 1 left the high-contrast curved case unexplained. `check_cylinder.py
correction` decomposes it: per clipped face it forms what the correction
SHOULD be, `C_exact = <k d_dT>_quadrature - k_face gtd`, with no model in it,
and compares the scheme's `C = s_t (K - k_face)` against the same formula fed
EXACT inputs.

```bash
./check_cylinder.py correction cyl_256.h5 --radius 0.25 --kappa 1000 --mode 2
```

**The formula is right.** At κ_s = 10³, h = 1/256: `|C_ideal − C_exact|` =
4.1e-3 against `|C_exact|` = 0.223 — **1.8 %**. So the
piecewise-constant-`∂_dT` premise holds, and the `Cov(k, s_t)` hypothesis (the
cell balance wanting the face average of the PRODUCT) is REFUTED, despite the
excess scaling linearly in κ_s, which is that hypothesis's signature.

**The area fraction is not it either.** Exact `s_t` + the plane `f` gives
3.9e-3; estimated `s_t` + the exact `<k>` gives 3.1e-1. The whole error is the
`s_t` estimate, even though stage 0 measured it at 0.95 %.

**Because stage 0 measured it on the wrong set.** Split the faces `k_area`
acts on:

| face set (κ_s = 10³, h = 1/256) | count | `s_t` error | correction error |
|---|---|---|---|
| marker-cut | 608 | **5.8e-6** (1.1 %) | 5.8e-3 |
| clipped-only | 424 | 6.2e-3 (1134 %) | **4.7e-1** |

and split the clipped-only faces again, by which material the ARM lies in:

| κ_s | arm in SOLID (296) | arm in FLUID (128) |
|---|---|---|
| 10 | **1.0 %** | 18.2 % |
| 10³ | **1.0 %** | **1905 %** |

**128 faces out of 1032 carry the entire high-contrast defect**, and the
mechanism is stage 0's rule a third time: `|∇_tT| ∝ 2/(1 + κ_s)` at the
interface, so `s_t` scales with the SOLID-side amplitude, and extracting it
from an O(1) fluid field is catastrophic cancellation.

**THE FIX: extend `s_t` from the faces that have a solid-side view.** `s_t` is
continuous along the interface, so the value exists nearby. Sweeping it in
over the 6 face neighbours:

| | κ_s = 10 | κ_s = 10³ |
|---|---|---|
| clipped-set `s_t` error, n = 64 / 128 / 256 | 13.9 / 14.2 / 6.0 % | 1421 / 1492 / 628 % |
| **after 2-4 extension sweeps** | **4.78 / 2.65 / 2.73 %** | **4.73 / 2.49 / 3.45 %** |

Contrast-independent, as it must be. **Two caveats, both for the next round:**
it PLATEAUS near 3 % rather than converging, and it needs 2-4 sweeps — which
in the solver is a halo-DEPTH question, not a free local gather. And the
residual table has NOT been re-measured with it; that is the immediate next
step, not a result.

**Four things that did not work, recorded so they are not retried:**

1. *Extend tangentially.* ZERO clipped-only faces have a marker-cut tangential
   neighbour — they are a different family. A marker-cut face has the
   interface roughly PARALLEL to it (|n_d| median 0.87); a clipped-only face
   has it roughly perpendicular (0.23), i.e. grazing. The good neighbour lies
   along `d`, not in-plane.
2. *Apply the one-sided estimate only at marker-cut faces*, on the argument
   that a same-material arm straddles nothing: **24× worse**. The shipped
   de-bias is not a fallback there — with `k_L = k_R` its coefficient vanishes
   identically, so it leaves the STRADDLING central tangential differences
   completely uncorrected (10760 % error, bit-identical to the raw projection).
3. *A direct projection at clipped-only faces* (`s_t = gtd(1−n_d²) − n_d P`,
   no `1/n_d`, which the closure needs and which is ill-conditioned at exactly
   those grazing faces): right in principle, 10× better than the shipped
   estimator, but still 1091 % — the cancellation, not the conditioning, is
   what dominates.
4. *A limiter.* Covered under Stage 1; with this diagnosis it is clear why
   neither could work — the defect is confined to an identifiable 12 % of the
   faces, so it wants a different ESTIMATE there, not a smaller correction
   everywhere.


### Stage 1b, THE TEST — every row beats the baseline, planes exactly

`./run_gates_c2.sh residual`. Cut-cell |div| rms on the exact field, whose
exact divergence is zero. This is what the cell balance sees, and it is the
measurement that decided C2.

| plane, h = 1/256 (theta = 20, 30, 35 deg; r = 0.01, 1, inf) | C1 baseline | `k_area` + shipped | **STAGE 1b** | gain |
|---|---|---|---|---|
| kappa_s = 10, 18 rows | 4.1 - 598 | 0.16 - 38 | **1.1e-7 - 1.6e-6** | 3e6 - 3e9 |
| kappa_s = 10^3 | 672 - 8.1e4 | 1.9e-5 - 5.8e3 | **1.6e-6 - 1.4e-4** | 6e6 - 2e9 |

| curved (dipole + quadrupole) | n = 64 | 128 | 256 | gain at n = 256 |
|---|---|---|---|---|
| kappa_s = 0.01: C1 | 9.7 / 20.0 | 21.7 / 42.5 | 41.3 / 85.7 | |
| **STAGE 1b** | **1.08 / 1.03** | **1.07 / 0.97** | **1.11 / 1.02** | **37 / 84** |
| kappa_s = 10: C1 | 7.6 / 16.6 | 13.5 / 30.7 | 31.0 / 67.0 | |
| **STAGE 1b** | **0.71 / 0.84** | **0.72 / 0.82** | **0.71 / 0.84** | **44 / 80** |
| kappa_s = 10^3: C1 | 10.2 / 22.4 | 23.2 / 47.4 | 45.4 / 94.0 | |
| **STAGE 1b** | **0.93 / 1.63** | **1.13 / 1.65** | **1.08 / 1.89** | **42 / 50** |

**NEVER WORSE, on all 36 rows** — three plane angles, two curved geometries,
five decades of contrast, three grids. A plane interface is solved EXACTLY at
any angle, contrast and tangential ratio. On a curved interface the residual is
**BOUNDED** (0.7 - 1.9, flat in both h and kappa_s) where the C1 baseline grows
like 1/h, so the gain IMPROVES with refinement -- 11 -> 21 -> 42 at kappa_s =
10^3. That is a complete reversal of C2's verdict, whose best variant was
40-65x WORSE than dropping the term.

**What it is made of** (`face_flux_field` mode `area1sext`), in the order the
measurements forced:

1. `k_area`, the face AREA-weighted multiplier, applied at every CLIPPED face
   (C2's own proposal -- it was never the problem);
2. the stage-0 one-sided `s_t` at MARKER-CUT faces, closed against the face's
   own normal flux -- inside the existing one-deep halo;
3. a DIRECT projection at clipped-only faces, which need no closure (the arm
   lies in one material) and must not have one (they are the grazing faces,
   where 1/n_d is ill conditioned);
4. the stage-1b extension sweep, carrying `s_t` into faces whose arm lies in
   the FLUID and which therefore cannot see a quantity scaling as
   `2/(1 + kappa_s)`.

**Two caveats, unchanged and honest.** The curved residual is BOUNDED, not
converging: ~1 against an uncut background of 3.6e-3 at n = 256, so cut cells
still carry a far larger truncation error than the interior -- the 1/h growth
is removed, pointwise consistency is not achieved. And the extension needs up
to 4 sweeps, i.e. a stencil reaching ~4 cells for ~12 % of clipped faces; most
need one. Since the closure is LINEAR in T with purely geometric coefficients
and phi is static, the whole thing collapses to a per-face (cell, weight) list
precomputed at init -- one sparse gather at runtime, no iteration -- but the
reach is a halo-DEPTH question for the solver.

**FOUR BUGS OF MINE, all of which read as "the new mode equals the baseline"
or as a silent regression** -- recorded because the checker now has eleven
modes and the next one will hit them too: a `detail` block dedented out of
scope so the correction was never applied; `"area1sg".endswith("1s")` false;
a gate applied AFTER the value it gates; and `usable` requiring BOTH sides'
closures when only the selected side is used, which disqualified faces whose
solid-side estimate was perfect and cost the curved high-contrast case a
factor 30. The last two interact: promoting every solid-side-viewable face
into `usable` is load-bearing (1.08 vs 32.5) but is only SAFE once `usable`
tests the selected side alone -- promoting with the two-sided test turned an
exact plane result at 20 degrees into 2.19. Isolate one edit at a time; the
bisect that found this took four runs and no guessing.


### THE CONVERGENCE TEST — the baseline is FIRST ORDER on a curved interface

`./run_gates_c2.sh converge`. Everything above measures a truncation residual
or an `s_t` error. Neither is a statement about the SOLUTION: the residual
lives on a codimension-one set, and a conservative scheme can damp it. This
solves the discrete conjugate BVP and refines.

**The problem.** A conducting cylinder in a harmonic far field (`Multipole`):
exact solution known in both materials, `div(k grad T) = 0`, Dirichlet data
from the exact solution on the two outer cell layers, interface far from the
boundary. **Solved EXACTLY, not iteratively** -- every scheme here is linear in
T with purely geometric coefficients (phi is static), so the matrix is
recovered by graph colouring and solved directly; no Krylov tolerance enters
the measured error. The recovery is checked against the operator itself
(2e-16), because a stencil wider than the colouring radius would otherwise
attribute entries to the wrong column *silently*.

**Observed L2 / Linf order, grids 32-256** (quadrupole):

| κ_s | C1 baseline | | STAGE 1b | | L2 error ratio at n=256 |
|---|---|---|---|---|---|
| | L2 | Linf | L2 | Linf | |
| 0.01 | 1.11 | 0.94 | **2.00** | 1.60 | 16.6x |
| 0.1 | 1.10 | 0.92 | **2.04** | 1.62 | 41.0x |
| 2 | 1.07 | 0.86 | **1.87** | 1.60 | 23.2x |
| 10 | **1.00** | 0.87 | **1.87** | 1.56 | 22.8x |
| 100 | 1.36 | 1.34 | **1.97** | 1.47 | 4.6x |
| 1000 | 1.89 | 1.65 | **1.99** | 1.49 | 1.4x |

At κ_s = 10 over `n = 64...512` -- an EIGHTFOLD refinement range -- the
baseline reads **0.999** and stage 1b **1.92**, with the gap widening as it
must: 23x at n = 256, **42x at n = 512**. The DIPOLE reproduces all of it
(baseline 1.03-1.09, stage 1b 2.09-2.11).

**Three things follow.**

1. **The C1 baseline as shipped is FIRST ORDER at a curved conjugate
   interface**, across the whole physically interesting contrast range. The
   as-built note said the solution order there had not been measured; it has
   now, and this is a correction to it (`docs/conjugate/conjugate_ibm_asbuilt`
   §6).
2. **The conservative cell balance damps the non-converging face error by
   exactly one power of h.** The face flux plateaus at O(r) with observed
   order -0.01; the solution is first order, not stalled. That is worth
   knowing on its own -- it is the quantitative content of "the flux is
   conservative, so partial cancellation is plausible".
3. **The baseline's WORST case is moderate contrast, not high**, which is the
   opposite of what the face-flux and residual measurements suggest. It
   recovers second order as κ_s → ∞ for a physical reason:
   `|grad_t T| ~ 2/(1 + κ_s)` at the interface, so in the isothermal limit the
   dropped term vanishes on its own. Stage 1b is second order at every
   contrast, so its gain is largest exactly where the baseline is weakest.

**A degeneracy check that passes**: at κ_s = 1 there is no contrast, the exact
field is `x^2 - y^2`, the discrete Laplacian reproduces a quadratic exactly,
and both schemes sit at round-off (2.5e-16) and are BIT-IDENTICAL. An
"observed order" fitted to round-off is noise, and the sweep reports one
(-1.30) -- which is why the row is a control and not a data point.

---

## C3 — the fraction-weighted capacity, the Nusselt diagnostic, and the time step

Three items that share one cost: the capacity and the time-step convention
both move `dt` at cut cells, so they land together and the re-gate is paid
once.

```bash
./run_gates_c3.sh [fraction|transient|nusselt|budget|channel|all]
```

Same environment as the C1/C2 drivers, plus `CBIN`/`CPREP` — the archived C2
binaries (`~/c2_ref_binaries`, commit `c243e56`), which are the POINTWISE-
capacity control of the transient gate. New helpers:

| file | role |
|---|---|
| `check_transient.py` | the two-material slab's lowest DECAYING EIGENMODE: the eigenvalue by bisection on the exact transcendental condition, the field and the interface flux in closed form, and the seeder. |
| `check_nusselt.py` | the interface-heat diagnostic three ways: `sum` (an independent Python transcription of the same discrete sum), `slab` (a closed form), `budget` (the control-volume cross-check). |
| `transient.ini`, `chan_conj.ini` | the two new cases. |

### THE TIME-STEP DECISION, and what it cost

**The cut-cell share is now 2, not 3** (`scalar_conjugate_peclet_rate`). C2
measured that `dt = pecletmax·3C/diag` at a cut cell sits at 96 % of the
`1.25 C/diag` that Gershgorin plus this RK3's real-axis limit allow, and found
the case that attains it. `share = 2` grants `0.8 C/diag` at the default
`pecletmax = 0.4` — **64 % of the bound, a 1.56× margin** — and the case that
went to NaN now runs:

```
oblique plane, theta = 30, kappa_s = 1e3, h = 1/64, 2000 steps
  pecletmax = 0.4  tangential_correction = false  dt = 3.833377e-08  max|theta| = 0.6794849   (C2: NaN)
  pecletmax = 0.4  tangential_correction = true   dt = 2.847725e-08  max|theta| = 0.6768589
  pecletmax = 0.2  tangential_correction = false  dt = 1.916689e-08  max|theta| = 0.6789875
  pecletmax = 0.2  tangential_correction = true   dt = 1.423863e-08  max|theta| = 0.6773903
```

The alternative — leave the rate and make the combination a documented config
error — was rejected: it asks the user to know a bound the solver can compute,
and the bound is not one anybody would guess (it is not a per-material `α`,
and it is attained only at cut cells and only when the interface is oblique).

**What it cost: `dt × 2/3` at cut cells, and nothing anywhere else.** The
baseline `dt` above is exactly 2/3 of C2's `5.750066e-08`. The correction's
own penalty is unchanged at `0.743` (`2.847725/3.833377`), which is C2's
recorded number — the two effects are independent, as they should be.

Every conjugate number recorded before this decision was measured against a
different `dt`, which is why C1's and C2's suites are re-run below.

### (0) the fluid volume fraction itself — **PASS**

`plane_box_fraction` is the increment's one new piece of arithmetic, so it is
checked where a mistake could not be shared with its own derivation
(`build_cpu/scalar_test`, the S2/S5a unit-test driver):

```
vfrac vs brute force: max deviation 2.9457E-06   (300^3 midpoints; the quadrature error is O(1/N))
scalar_test: ALL PASS
```

plus hand-derivable exact values in each of the three degeneracy regimes — the
plane parallel to two axes (`f = clip(½ + φ/h)`, which is the grid-aligned wall
and therefore the COMMON case, not a special one), to one axis (the triangular
corner, `f = s²`), and to none (the tetrahedral corner, `f = s³√3/2`) — the
symmetry `f(φ) + f(−φ) = 1` over 40 random planes, the clipping to exactly 0/1
outside the cut band, and the medial-axis fallback.

TWO THINGS THE CLOSED FORM NEEDED THAT THE PAPER FORM DOES NOT SAY:

1. **Degenerate directions are the common case, not an edge case.** A
   grid-aligned wall has two of the `a_i = |n_i| h_i` at zero (up to the float
   noise of a distance field), and the 3D inclusion-exclusion then divides a
   numerator that has cancelled to nothing by a denominator that is nothing.
   The limits `a₃ → 0` and `a₂ → 0` ARE the 2D and 1D forms, so the routine
   drops negligible directions and evaluates the reduced form.
2. **Clip first, and exactly.** Outside the cut band the inclusion-exclusion
   is an identity that holds analytically and cancels numerically — its terms
   grow like `s³` while the answer stays `6a₁a₂a₃` — so a cell a few `h` from
   the interface came back as `1 − 6e-13` instead of `1`. Nothing downstream
   would have broken, but "which cells are cut" is a CLASSIFICATION, and the
   2026-08-05 body-heat lesson is that a classification must not be decided by
   round-off. Measured on `wavy.ini` before and after: 4080 "partial" cells
   became **304**, which is the right number (a 32-column wall crossing one to
   two cell rows over 8 spanwise planes).

### (0b) the one new config guard — **PASS**

`heat_interval` together with `tangential_correction` is rejected:

```
error: [scalar] heat_interval with [scalar.N] tangential_correction = true:
 the interface-heat diagnostic reports the baseline cut-face flux, which is
 not the one the kernel applies with the correction on
```

The diagnostic reports the BASELINE cut-face flux. With the correction on the
kernel applies a different flux at those faces, and a diagnostic that does not
report the flux the kernel applied is worse than no diagnostic — that is the
invariant the S4 accumulators exist for. The correction ships disabled by
measurement (C2's verdict), so rather than carry a second, ungated copy of its
six-face stencil in the statistics, the combination is a hard error that says
which two keys are in conflict and why.

### (1) the transient two-material slab — **PASS**, and the order the old capacity cost

The gate the capacity change exists for. C1 gated the STEADY state, where the
capacity is irrelevant by construction (it divides an rhs that vanishes); this
is the transient, and the reference is the one two-material transient
available in closed form — the lowest DECAYING EIGENMODE,

```
T(y, t) = exp(-mu t) X(y),   X = A sin(k_s y) | B sin(k_f (L - y)),
kappa_s k_s cot(k_s y_w) + k_f cot(k_f (L - y_w)) = 0
```

(`check_transient.py`; the eigenvalue by bisection between 0 and the first
pole, to `|F| ~ 5e-14`, and `κ_s = C_s = 1` returns `mu = pi^2` and `B = 1`,
which is the sanity check on the solver of the transcendental equation).

`y_w = 0.31`, `κ_s = 10`, `C_s = 4`, `n_y ∈ {16, 32, 64}`, `dt = h²/100` FIXED
(`pecletmax = 0`) so the RK3 temporal error — order 3, i.e. `h⁶` here — cannot
be mistaken for the spatial one, and so that the control leg runs the same
`dt`. `t_end = 0.05`, at which the mode has decayed to 0.397.

| | `n_y` | field `L₂` | field `L∞` | decay rate `mu` |
|---|---|---|---|---|
| **C3, fraction-weighted** | 16 | 6.918e-3 | 1.119e-2 | 6.557e-3 |
| | 32 | 1.739e-3 | 2.828e-3 | 1.654e-3 |
| | 64 | 4.350e-4 | 7.080e-4 | 4.139e-4 |
| | *order* | **1.99, 2.00** | **1.98, 2.00** | **1.99, 2.00** |
| C2, pointwise (the control) | 16 | 7.313e-3 | 1.166e-2 | 7.377e-3 |
| | 32 | 2.230e-3 | 3.346e-3 | 3.035e-3 |
| | 64 | 1.049e-3 | 1.587e-3 | 2.082e-3 |
| | *order* | 1.71, **1.09** | 1.80, **1.08** | 1.28, **0.54** |

**The pointwise capacity costs a full order, and the control measures it
rather than the argument asserting it.** The mechanism is the one the plan
predicted: `mu` is a global functional of the capacity distribution, so a cut
cell carrying one material's whole heat capacity is an O(1) error on an O(h)
band. The control leg is the archived C2 binary (`~/c2_ref_binaries`, commit
`c243e56`) on the SAME generated inis — same geometry, same `dt`, same
everything but the line this increment changed.

**DEVIATION from Section 10's wording.** It asks for "second order in the
time-resolved interface flux". The flux's TIME DEPENDENCE — its decay rate,
which is the time-resolved content — is second order, as the table shows. Its
INSTANTANEOUS AMPLITUDE is not, and cannot be: the cut-face flux
`(T_R − T_L)/(δ_L/k_L + δ_R/k_R)` is the flux at the resistance-weighted point
of the arm, which is O(h) away from the interface, so it carries
`q'·(δ_R²/k_R − δ_L²/k_L)/2R` — first order, with a coefficient whose SIGN
depends on where in the cell the interface falls. Measured, and mixed with the
O(h²) eigenvalue error times `t` (which enters the amplitude at fixed time):
`1.38e-3 / 1.13e-3 / 1.75e-3` relative — bounded and small, but not an order.
That is the same class of statement as C2's deviation 2: what converges at the
scheme's order is the SOLUTION, not the pointwise cut-face flux.

**FOUND HERE, and it is why the diagnostic needed its own branch:** the
PRE-C3 interface-heat diagnostic on a conjugate case reported a meaningless
number. The control leg's heat file reads

```
C2 binary   theta_staircase -3.108e-01   theta_graded -7.658e-01   total -1.077e+00
C3 binary   theta_staircase -4.480e-01   theta_graded  0.000e+00   total -4.480e-01
```

— the `graded` column is the S3 penalization term `coef_p(s_body − s)`, which
in a conjugate run is a large spurious contribution from graded FLUID cells
(nothing pins their value), and the staircase column carried the molecular
diffusivity instead of the cut-face coefficient. C1's deviation 6 recorded the
columns as smoke-gated; they were also wrong. (The `mu` column of the control
leg survives that, and is still meaningful, because both terms are
proportional to `theta` and therefore decay at the same rate.)

### (2) the interface-heat (Nusselt) diagnostic — **PASS**

C1 left the conjugate flux columns of `scalar_stats.f90` smoke-gated. Three
independent statements, none of which needs a reference run.

**(2a) a closed form.** A slab whose outer face is INSULATED and whose solid
carries a volumetric source must, at steady state, deliver every watt it
generates across the interface: `H = C_s S y_w A`. It is exact rather than
approximate because the fluid fraction of a cell cut by a PLANE is exact, so
the fraction-weighted source integrates to the true solid volume — i.e. the
gate tests the fraction and the diagnostic at once. `y_w = 0.3125`,
`κ_s = 4`, `C_s = 1`, `S = 2`, `A = 0.0625`, `n_y = 16`:

```
solver 3.9062499999982139e-02   closed form 3.9062500000000000e-02   relative 4.572e-13
```

at `t = 15.6` (400 000 steps); the residual is the run's remaining transient,
not the diagnostic — the previous sample differs by `1.164e-09` relative, so
the approach to the closed form is exponential (a factor ~2e4 per write
interval: at half the steps it stood at `1.165e-09`), and the gate reports
that alongside so "converged" is a measurement rather than an assertion.

**(2b) the solver against an independent Python transcription** of the same
discrete sum, rebuilt from the case file's `φ` and the snapshot's `θ` (the S4
idiom — it catches a transcription slip on either side):

```
slab     16 cut faces    python 3.9062499999982153e-02  solver ...139   relative 3.553e-16
channel  8192 cut faces  python -4.1139642840667688e-05 solver ...776  relative 2.141e-15
```

The channel line is the LES/IBM case of gate 3, so the transcription check
also covers the branch where `ν_t` exists and must NOT enter a cut face.

**(2c) the control-volume cross-check, the way A2 validated `C_L`/`C_D`** —
but stated so it holds at any TIME rather than only at a steady state, which
is what makes it affordable. Over the cells the solver classifies fluid the
discrete flux form telescopes and leaves exactly the interface faces, so

```
d/dt sum_fluid C theta dV  =  H  -  (flux out through a plane in the fluid)
```

with the plane flux (convective AND diffusive) evaluated in Python from the
snapshot. C1's insulated wavy box with the flow ON — oblique geometry, the
analytic `dwall` path, `κ_s = 5`, `C_s = 2`:

```
                                   dE/dt          H            plane flux     residual   relative
niter =  40, dt = 2e-4, whole    -3.303044e-2  -3.303022e-2       --          -2.131e-7  6.453e-06
niter =  40, dt = 2e-4, plane    -3.263943e-2  -3.303022e-2   -3.909909e-4    -2.001e-7  6.058e-06
niter = 200, dt = 2e-4, whole    -3.303044e-2  -3.303022e-2       --          -2.131e-7  6.453e-06
niter = 200, dt = 2e-4, plane    -3.263943e-2  -3.303022e-2   -3.909909e-4    -2.001e-7  6.058e-06
niter =  40, dt = 1e-4, whole    -3.310271e-2  -3.310266e-2       --          -5.355e-8  1.618e-06
```

Three readings:

1. **The residual is the centred difference's own truncation, measured, not
   assumed**: halving `dt` divides it by `3.98`. It is not the diagnostic's
   error.
2. **It does not move with `niter`** — identical to four digits at 40 and 200.
   The A2 landmine (a border flux is only as good as the projection's
   divergence) is real for a FORCE, which reads the pressure; the scalar
   border flux reads the velocity and the residual here is dominated by
   something else entirely. Recorded because the next person will ask.
3. **The plane term is doing work**: it is `3.9e-4`, three orders above the
   residual, so a gate that dropped it would fail by 1.2 %.

### (3) the conducting channel wall — **PASS**

The scheme on the full production stack — WALE LES + file-based IBM + (the
`refine` leg) 2:1 block refinement — on `validation/channel_interface/les_ibm`'s
geometry. That geometry is chosen for what its walls are: flat planes sitting
MID-CELL (`y = 8.3 dy` and `72.3 dy` on a uniform `dy = 0.03125`), so every cut
cell has a non-trivial fluid fraction — the thing this increment changed —
while the interface stays grid-aligned, where the C1 cut-face coefficient is
EXACT. Anything the gate sees is therefore the capacity, the diagnostic or the
stack, not the interface scheme's own truncation.

The coefficient file is prepared here (`moby_prepare`, the same two wall STLs)
rather than reused: the committed `ibm_coeff.h5` predates S3 and carries
neither `coef_p_blocks` nor `dwall_blocks`. The refined leg COLD-STARTS —
`refine_body` + `keep_buried` keeps 3328 leaves where the committed
`IC_refine.h5` (prepared without `keep_buried`, which a conjugate run may not
do) has 2560, so their block tables do not match. The budget identity does not
care what the flow is.

```
flat    640 leaves    65536 solid cells    8192 cut faces
  fluid-cell budget   dE/dt -4.094375e-05   H -4.094306e-05   residual 1.681e-05 relative
  vs Python           8192 cut faces        relative 2.141e-15
refine  3328 leaves  524288 solid cells   32768 cut faces   (2:1 interfaces present)
  fluid-cell budget   dE/dt -4.821949e+01   H -4.821926e+01   residual 4.817e-06 relative

determinism, BOTH legs, tolerance 0, on un vn wn pn nut theta vfrac:
  1 rank vs 4 ranks   max_abs 0     CPU vs GPU   max_abs 0
```

Three things this covers that nothing else does. **`ν_t` is active**, so the
budget also states that the eddy diffusivity enters the diagnostic and the
kernel identically — and that it enters NEITHER at a cut face (if it did, the
two would disagree, since only one of them would have added it). **The refined
leg carries 2:1 interfaces**, so the conjugate scheme, the capacity and the
diagnostic are exercised on the multi-level path; its 32768 cut faces are 4×
the flat leg's. And **`vfrac` is in the determinism comparison**, which is the
check that the fraction is a deterministic function of the geometry rather
than of the decomposition — it reads `φ`'s ghost layer, which on the refined
leg is filled per leaf at that leaf's own level.

The budget legs run on the GPU (they are a physics statement, not a
bit-exactness one, and the refined case is 50 minutes of one CPU rank against
under a minute of device time); the determinism trio is the one that needs the
nofma pair.

**A driver bug worth recording, because it produced a plausible wrong number
rather than an error.** The two-stage run first chained its stage-B restart by
INSERTING a `[restart]` section before `[output]` — while the template already
carried one. Two sections, the last one wins, so stage B silently restarted
from the ORIGINAL initial condition and reset the step counter; the snapshots
the budget wanted did not exist and the checker died on a missing file, which
is the lucky outcome. Had the names collided instead, the gate would have
measured a budget across a discontinuity. The line is now a REPLACEMENT.

### (4) every C1 gate, re-run under the new `dt` and the new capacity — **PASS**

`./run_gates_c1.sh all`, 61 checks, `ALL C1 GATES PASS`. The recorded numbers
move only where they must:

| C1 gate | C1 recorded | C3 re-measured |
|---|---|---|
| slab fixed point, 28 `(w, κ_s)` pairs | ≤ 9.66e-15 | **≤ 8.47e-16** |
| cold start, `κ_s = 10⁻²/1/10` | 2.387e-14 / 7.438e-15 / 4.774e-15 | 3.975e-14 / 1.044e-14 / 7.772e-15, residual **exactly 0.0** |
| level-set weight `w` | 2.226e-14 / 0.0 / 2.220e-14 | **identical** (geometry) |
| capacity irrelevant at steady state | 8.85e-14 / 4.37e-14 / 4.88e-15 | 1.154e-13 / 6.550e-14 / 9.104e-15 |
| contact resistance, `q` | 1.276426007180 / 0.967644390686 / 0.209054680865 | **identical**; errors 5.55e-17 / **0.0** / **0.0** |
| `Σ C θ dV` drift, insulated wavy box | −6.94e-18, rel 1.22e-16 | 6.939e-18, rel **1.206e-16** |
| `κ_s → ∞`: `q`, `|q − q_∞|`, rate | 1.332889036988, 4.443e-04, order 1.000 | **identical** |
| `κ_s → 0`: `|q|/κ_s`, rate | 3.988036 / 3.999880 / 3.999999, order 1.000 | **identical** |
| Peclet limiter on / control | bounded 0.869 / NaN | bounded **0.839** / NaN |
| guards, 2:1 precondition, stats smoke | 5/5, both ways, 4 datasets 0 non-finite | **identical** |

Two of these are the C3 items showing through, and both are as designed:

- **the `converge` group still converges.** Its 200 000 steps buy 2/3 of the
  physical time they used to, and the cold start still reaches the exact
  profile with a residual between the last two writes of exactly 0.0. This was
  checked, not assumed — it is the one C1 gate the time-step decision could
  have broken.
- **the conservation total MOVED, and the invariant did not.**
  `Σ C θ dV` is now `5.7536963634433440e-02` where C1 recorded
  `5.7128906250000007e-02`, because `C` is a different (correct) map; the
  DRIFT is unchanged at round-off. The checker had to be changed with it: an
  analytic pointwise capacity map would leave `Σ (C_ref − C) Δθ dV`, which is
  ~1e-3 relative here, so the gate would have measured the map mismatch
  instead of conservation. It now reads the solver's own `vfrac` from the
  snapshot — which is what that dataset is for — and reports the deviation
  from the pointwise marker (0.472 at a cut cell) as information.

### (5) the bit-exactness protocol — **PASS**

```
REF=~/s5c_ref_binaries/moby_solve_cpu_nofma MODE=cpu ../scalar/run_bitexact.sh     # ALL PASS
REF=~/s5c_ref_binaries/moby_solve_gpu_nofma MODE=gpu ../scalar/run_bitexact.sh     # ALL PASS
REF=~/s5c_ref_binaries/moby_solve_cpu_nofma MODE=cpu ../scalar/run_bitexact_s3.sh  # ALL PASS
REF=~/s5c_ref_binaries/moby_solve_gpu_nofma MODE=gpu ../scalar/run_bitexact_s3.sh  # ALL PASS
```

**32 case-runs, every one `max_abs 0`** — the 7-case standard suite
(`[scalar] count = 0`) and the 9-case scalar suite (`ibm_wall` =
`dirichlet` / `adiabatic`), CPU and GPU, zero non-zero deviations anywhere in
the four drivers' output. By construction, and the construction is three
things: `sc%vfrac` is a 1-cell dummy without a conjugate scalar (so the
capacity expression is never evaluated), the `vfrac` dataset is gated on
`size > 1` exactly as `iddes_fd` is (so no other file gains a byte), and the
diagnostic's conjugate branch is mode-gated per scalar.

### (6) every C2 gate, re-run under the new `dt` — **PASS**

`./run_gates_c2.sh all` → **`ALL C2 GATES PASS`**, and its `c1` group →
**`ALL C1 GATES PASS`** (the whole C1 suite with `tangential_correction =
true`). What moved and what did not is the useful part:

- **The analytic-field measurements did not move at all**, as they cannot:
  `check_oblique.py flux`, `check_oblique.py residual` and
  `check_cylinder.py` evaluate the schemes on a MANUFACTURED field with `φ`
  from the case file, so no solver and no `dt` enters them. The solver's own
  indicator still reads `max|e_face| = 1.0000E+00`, `rms = 1.0000E+00` on 404
  cut faces against the checker's `1.0000e+00 / 1.0000e+00`.
- **The `r = ∞` BVP moved, and its conclusion did not.** It is a relaxation to
  a discrete steady state, so a different `dt` reaches a slightly different
  point at the same step count: `L2 = 3.310e-03` (correction off, C2 recorded
  3.328e-03) against `8.682e-03` (on, C2 recorded 8.745e-03) — **still 2.6×
  worse with the correction**, at the ratio most favourable to it.
- **The C1 suite with the correction on reproduces C1's own numbers**: slab
  fixed point ≤ `8.465451e-16` (identical to the correction-off re-run, as it
  must be — `s_t = 0` by construction on a grid-aligned interface), the
  level-set weight `2.226e-14 / 0.0 / 2.220e-14`, the `κ_s → ∞ / → 0` rates
  1.000, conservation drift `1.388e-17` (relative `2.412e-16`), and 24
  `max_abs = 0` dataset comparisons with zero non-zero ones in the
  determinism group.

The `dt` group is the time-step decision's own measurement, at the top of this
C3 section.

…and the S-suite physics groups, all with `rc = 0` and no failures:

```
BIN=build_cpu_nofma/moby_solve GBIN=build_gpu_nofma/moby_solve ./run_gates.sh
                                          scalar S1 gates: no failures
./run_gates_s2.sh                         scalar S2 gates: no failures
./run_gates_s3.sh                         S3 gates (all): ALL PASS
./run_gates_s4.sh                         10/10 PASS, no failures
./run_gates_s5.sh                         S5a gates: ALL PASS
```

LANDMINE honoured: `run_gates.sh`'s `det` group compares CPU against GPU at
TOLERANCE 0, so it gets the nofma pair.

**PROVENANCE.** Every number in this C3 section was produced with all four
binaries (`build_{cpu,gpu}`, `build_{cpu,gpu}_nofma`) rebuilt from the FINAL
source before the gate runs — which matters more here than in C1 or C2,
because this increment MOVES `dt`: a run against a stale nofma pair would gate
different code, and the `det` groups compare three binaries against each other
and pass happily when all three are equally stale. The C2 control binaries
live in `~/c2_ref_binaries` (commit `c243e56`, cut from a COMMIT and recorded
in its `PROVENANCE.txt`, the `~/s5b_ref_binaries` lesson).
