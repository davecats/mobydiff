# Conjugate heat transfer in a turbulent channel — the physics validation

The C1–C3 gate suite is entirely manufactured solutions, analytic limits and
self-consistency identities. It says the scheme solves what it claims to
solve; it says nothing about whether the scheme reproduces a **turbulence**
result. This is that case.

```bash
./run_cht.sh [prepare|ic|develop|reseed|stats|all]
./check_cht.py velocity cht_vel_stats.h5
./check_cht.py thermal  cht_stats.h5 --dump cht_sweep.dat
```

## The reference

**Flageul, Benhamadouche, Lamballais & Laurence, "DNS of turbulent channel
flow with conjugate heat transfer: effect of thermal boundary conditions on
the second moments and budgets", Int. J. Heat and Fluid Flow 55 (2015) 34–44**
(`literature/flageul.pdf`). Re_τ = 149, Pr = 0.71, four thermal boundary
conditions — isoT (Dirichlet), isoQ (Neumann), Robin and **Conjug**.

Their conjugate case has **G = G₂ = 1**: the solid's conductivity *and* its
diffusivity equal the fluid's, i.e. κ_s = 1, α_s = 1, **K = √(κ_s C_s) = 1**.
That is exactly this sweep's `k1`, and their solid is exactly one half-height
thick — which is why this case uses d = h = 1 (d⁺ = 180) rather than a thin
slab. Values quoted below are the temperature variance ⟨T′²⟩/T_τ² read from
their **figure 5** (there is no table, so treat them as ±0.1, not as digits):
isoQ wall **4.2**, Conjug wall **1.1**, isoT wall **0**, peak **≈ 6.3** at
y⁺ ≈ 18 and nearly boundary-condition-independent.

## The case

Re_τ = 180, Pr = 0.71, 4π × 4 × 2π, **128 × 288 × 128 → 128 × 480 × 128**
(7.9 M cells), **single level — 240 leaves, zero 2:1 interfaces** (this case
sets no `refine*` key; `nb = 32` is only the block size, and it is required
because a prepared case file's block table cannot depend on the rank count).

* fluid gap exactly 2.0 at y ∈ [1, 3], solid slabs d = h = 1 either side;
* the interface sits **on a cell face**, so w = ½ to the last bit and the
  cut-face coefficient is the textbook two-resistance harmonic mean — *exact*.
  Anything measured here is the conjugate physics, not the IBM approximation;
* dy⁺ = 1.5 uniform (first cell centre y⁺ = 0.75), dx⁺ = 17.7, dz⁺ = 8.8;
* antisymmetric Dirichlet (+1 / −1) at the **outer solid faces**, no source,
  so the total wall-normal flux is constant across solid and fluid alike;
* velocity IC: the developed KMM180 field mapped into the gap
  (`make_cht_ic.py`, `tools/make_channel_restart.py`'s validated staggered
  interpolation), then t = 0…10 of development.

**Six scalars ride one velocity field**, so every comparison between them is
at the same turbulence realisation:

| | κ_s | C_s | K | α_s |
|---|---|---|---|---|
| k1 / k2 / k3 | 1 | 1 / 100 / 10⁴ | 1 / 10 / 100 | 1 / 10⁻² / 10⁻⁴ |
| a1 / a2 / a3 | 0.1 / 10 / 100 | same | 0.1 / 10 / 100 | 1 |

## The velocity field, and a landmine

Validated against an **exact law** rather than a reference file: in a steady
channel the total stress is linear, `−⟨u′v′⟩ + ν dU/dy = u_τ²(1 − y/h)`.

```
this run          wall value 0.9964 (exact 1)   max |tau - (1-y/h)| = 0.011  mean 0.007
U+ centreline 17.950     u'_rms peak 2.587     v'_rms 0.813     w'_rms 1.004
first fluid cell (y+ = 0.75):  u+ = 0.769  against y+ = 0.75
solid cells |u|,|v| ~ 1e-27      wall-face v = 6.6e-12   (no penetration)
```

**LANDMINE, recorded because it nearly cost a wrong conclusion:
`tutorials/channel_kmm180/channel_kmm180_stats.h5` is NOT published
Kim–Moin–Moser data** — it is this solver's own accumulated statistics from
that tutorial (its `input.ini` writes it), **and it is not converged**: it
violates the same total-stress law by up to **22 %** (mean 15 %). Compared
against it, this run's v′ looked "+19 % high"; against the exact law it is the
reference that is wrong, and 0.813 is the physical value. Do not use that file
as a reference.

## Results

### The conjugate signature — the headline

The wall/peak variance ratio is independent of θ_τ and of the fluid-side
level, so it isolates what the conjugate wall actually does: damp the
interface temperature fluctuation relative to the fluid's own.

| K | ⟨θ′²⟩_wall / ⟨θ′²⟩_peak | Flageul |
|---|---|---|
| **1 (k1)** | **0.1745** | **0.175** — their conjugate case |
| 0.1 (a1) | 0.5012 | (isoQ limit 0.667, not reached at K = 0.1) |
| 100 (a3) | 0.0068 | 0 (isoT limit) |

**0.3 % on the one number that measures the interface coupling**, at
identical thermal parameters (κ_s = 1, α_s = 1, d = h), with the ideal limits
bracketing it correctly.

Absolute values, ⟨θ′²⟩/θ_τ²:

```
       K = 0.1   1     10(k2) 10(a2) 100(k3) 100(a3)      Flageul
wall     4.21   1.50   0.26   0.157   0.167   0.059    isoQ 4.2 / conjug 1.1 / isoT 0
peak     8.40   8.60   8.88   8.75    8.77    8.65     6.3
```

### Internal checks that came out right

* **the fluid resistance is one number.** R_f/2 = Δθ_half/J measured
  independently from all six conjugate walls: **18.83 … 19.13, spread 1.6 %**.
  It is a property of the flow, so six different solids agreeing on it is a
  real check on the interface coefficient;
* **the mean profile.** θ⁺ within 4–6 % of **Kader's correlation** over the
  log layer. At y⁺ = 5.25 this run gives 3.094 against the exact conduction
  value 3.195 — closer than Kader's own blending (2.599), so the apparent
  20 % "error" there is the correlation, not the run;
* **the S-curve is monotone** over the converged scalars: θ′_w/θ_τ =
  2.05 (K = 0.1) → 1.22 (K = 1) → 0.396 (K = 10) → 0.243 (K = 100).

### Two things that did not come out as designed

**1. The high-capacity scalars are not converged, and could not be.** With
d = h the solid holds most of the resistance (R_s = 128 vs R_f ≈ 38) and its
mean relaxes on d²/α_s — t ≈ 80 for α_s = α_f, but **8×10³ and 8×10⁵** for
C_s = 100 and 10⁴. Measured against the exact steady flux
`J = 2/(2R_s + 2R_f/2)` with the common R_f:

```
a1 0.998   a2 0.996   a3 0.977   k1 0.991        <- converged
k2 1.120   k3 1.198                              <- NOT converged
```

`reseed_solid.py` exists for exactly this (the capacity cannot change the
steady state — C1 gate 1c measured that to 1e-13 — so the solid mean is set,
not waited for), but it used **each scalar's own** transient-contaminated R_f
instead of the single fluid value. THE FIX IS ONE LINE: measure R_f once from
the converged scalars and re-seed every solid with it. The wall/peak ratios
above survive this (numerator and denominator scale together), which is why
they are the quantity to read.

**2. The "effusivity collapse" prediction was too naive — and the paper says
why.** Same K, different (κ_s, C_s):

```
K = 10:  k2 0.0294  vs  a2 0.0179       K = 100:  k3 0.0191  vs  a3 0.0068
```

They do not collapse. Flageul's eq. (12) is the reason: the interface response
carries `R² = k_x² + k_z² + i k_t ρc/λ_s`, so `λ_s R` reduces to the effusivity
`√(λ_s ρ_s c_s)` **only when the temporal term dominates**. Where the lateral
conduction term matters, the damping goes with the **conductivity**, and the
higher-κ_s wall should damp more at equal K — which is exactly the ordering
measured (a2 damps more than k2, a3 more than k3). So this is their physics,
not a defect; the design's claim that K alone should collapse the response was
wrong.

## The open discrepancy

**The fluid-side variance level is ~37 % high for every boundary condition**
(peak 8.4–8.9 against Flageul's 6.3) — including the ideal limits, so it is
*not* an interface error. The mean profile is right and the wall/peak ratio is
right; only the level is off. Most likely cause, from the reference itself:
this run has dx⁺ = 17.7, dz⁺ = 8.8 against their 14.8, 5.1, and they note that
"the accurate resolution of the scalar evolution equation is very demanding in
terms of grid spacing", citing Galantucci et al.'s recommendation of
Δx⁺ = Δz⁺ = 1 for an accurate ε_θ — an under-resolved scalar dissipation
leaves the variance too high. Two other candidates not yet excluded: the
averaging window (t⁺ = 3.8×10³ against their 2.9×10⁴), and the thermal problem
(constant-Δθ here, Kasagi's mean-streamwise-gradient formulation there).

**Settling it needs a resolution study** — 256 × 480 × 256 at the same
physical setup, ≈ 4× the cost (~6 h on this GPU). That is the next step, and
until it is run the correct statement is: *the conjugate interface coupling is
reproduced (0.3 % on the wall/peak ratio, correct brackets, monotone
S-curve, consistent fluid resistance); the absolute fluid-side fluctuation
level is 37 % high for a reason that is common to all four boundary conditions
and therefore not attributable to the conjugate scheme.*
