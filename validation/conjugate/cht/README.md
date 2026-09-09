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

Re_τ = 180, Pr = 0.71, 4π × 3.2 × 2π, **160 × 384 × 224** (13.8 M cells),
**single level — 420 leaves, zero 2:1 interfaces** (this case sets no
`refine*` key; `nb = 32` is only the block size, and it is required because a
prepared case file's block table cannot depend on the rank count).

* fluid gap exactly 2.0 at y ∈ [0.6, 2.6], solid slabs d⁺ = 108 either side —
  **sized from the measurement**: campaign 1 (d⁺ = 180) showed the interface
  variance down to 0.3–1 % of its interface value by d⁺ = 108 and to 1.6e-7 at
  d⁺ = 180, so the deepest 70 wall units bought nothing and their cells went
  into x–z resolution instead;
* the interface sits **on a cell face**, so w = ½ to the last bit and the
  cut-face coefficient is the textbook two-resistance harmonic mean — *exact*.
  Anything measured here is the conjugate physics, not the IBM approximation;
* dy⁺ = 1.5 uniform (first cell centre y⁺ = 0.75), **dx⁺ = 14.1, dz⁺ = 5.05**
  — comparable to Flageul's 14.9 and 5.0. Uniform y is forced: the solver has
  no piecewise or file-based node line, and every analytic family it does have
  refines at the DOMAIN edges, which here are the solid outer faces. It is
  also expensive in a way stretching would not fix — see the last section;
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
campaign 2        wall value 0.9915 (exact 1)   max |tau - (1-y/h)| = 0.0050  mean 0.0024
U+ centreline 18.199 (the textbook Re_tau 180 value)
u'_rms peak 2.569     v'_rms 0.827     w'_rms 1.052
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

Two campaigns and a wall-normal grid-convergence study, all at
`interfaces on a cell face`, single level:

| | dx+ | dz+ | dy+ | window (wall units) |
|---|---|---|---|---|
| campaign 1 | 17.7 | 8.8 | 1.5 | 3.8e3 |
| campaign 2 | 14.1 | 5.1 | 1.5 | **2.9e4** (= Flageul's) |
| y study | 14.1 | 5.1 | 1.5 / 1.0 / **0.75** | 4.5e3 each |
| Flageul | 14.9 | 5.0 | 0.49 - 4.8 | 2.9e4 |

### A COMPARISON ERROR, FOUND AND CORRECTED — read this first

Campaign 1 reported the temperature variance "37 % high" and a wall/peak ratio
matching Flageul to 0.3 %. **Both numbers were wrong, for the same reason.**
The checker took the peak as `max` over the fluid. In THIS case's thermal
problem the total flux is CONSTANT across the channel (antisymmetric walls, no
source: measured `J(y)/J_wall = 1.0000` everywhere), so production never
switches off and the variance keeps rising to the **centreline** — that is
where the global maximum is. Flageul's wall-flux problem has the flux falling
linearly to zero at the centre, so their variance decays and their global
maximum IS the near-wall peak. The two `max`es are different quantities.

```
   y+     <theta'^2>     J(y)/J_wall
    5.3      3.21          1.0001
   20.3      7.07          1.0002     <- the near-wall peak: THIS is Flageul's 6.3
  144.8      7.96          0.9998
  177.8      8.79          0.9992     <- the centreline: the old "peak" of 8.8
```

Corrected, the near-wall peak is **+11 %**, not +40 %; and the wall/peak ratio
agreement is 1 %, not 0.3 %, but is now grid-converged rather than accidental.
`check_cht.py` takes the peak in `5 < y+ < 40` and reports the centreline
value separately, with a comment saying why.

### The figure

![conjugate channel validation](cht_validation.png)

`./plot_cht.py` — (a) the mean against **Kader's correlation** (a published
analytic curve, so it is drawn as a curve; dashed where it no longer applies);
(b) the LOCAL SLOPE, which is the panel that says what the outer profile is
actually doing (see below); (c) the variance for the K ladder with **Flageul's
values as markers with an uncertainty bar** — they publish no table, those four numbers were read off
their figure 5 to about ±0.1, and drawing a reference *line* through eyeballed
points would dress a reading up as data; (d) the interface coupling against
wall-normal resolution.

Panel (c) also shows the comparison trap directly: our curves keep rising
toward the centreline while theirs decay, because our flux is constant across
the channel and theirs falls to zero. **Only the shaded band is comparable.**
The y axis is logarithmic because the wall values span two decades (3.83 down
to 0.051) — on a linear axis the two high-K cases sit on the axis and cannot
be told apart, which is precisely the quantity the figure exists to show.

### Why our flux is constant and theirs is not

The steady, x-z homogeneous mean scalar balance is

    dJ/dy = S ,      J = <v'theta'> - D dTheta/dy

so the flux profile is set entirely by the source. **Ours has S = 0** —
antisymmetric wall temperatures, no volumetric source — so J is CONSTANT:
heat passes straight through from the hot wall to the cold one and nothing in
the interior absorbs it. Measured `J(y)/J_wall = 1.0000` everywhere, solid
slabs included, which is also why the series-resistance model of `reseed_solid.py`
is exact.

**Theirs has S != 0.** Their eq. (2) carries Kasagi's `f_T u_x`: the mean
temperature grows downstream, so convection acts as a distributed sink
proportional to the LOCAL velocity, `dJ/dy = -f_T u(y)`. Both their walls are
heated identically, so symmetry forces J = 0 at the centreline.

That is a different physical problem, not a defect in either. It also means
that saying their flux "falls **linearly**" — as an earlier draft here and the
figure caption both did — is wrong as stated: that is the MOMENTUM stress,
which is linear because its source (the pressure gradient) is uniform. Their
SCALAR flux follows the cumulative flow rate, because its source is `u(y)`; it
lags a straight line near the wall, where u is small, and overtakes it in the
core. Computed from this run's own
mean velocity profile:

| y+ | J/J_wall, Kasagi form | 1 − y/h | ours |
|---|---|---|---|
| 10 | 0.984 | 0.946 | 1.000 |
| 60 | 0.744 | 0.671 | 1.000 |
| 100 | 0.506 | 0.446 | 1.000 |
| 179 | 0.000 | 0.004 | 1.000 |

### The outer profile is NOT logarithmic, and that is physics

An earlier draft of this README claimed the mean profile "keeps a log-like
slope where Kader flattens". **It does not**, and panel (b) exists because the
eye cannot tell on a semilog plot. Measured:

| y+ | dθ⁺/dy⁺ | 1/(Pr_t κ y⁺) | ⟨vθ⟩/J |
|---|---|---|---|
| 0.8 | **0.708** | 3.83 | 0.000 |
| 29 | 0.128 | 0.098 | 0.813 |
| 100 | **0.0456** ← minimum | 0.0288 | 0.936 |
| 179 (centre) | **0.0560** | 0.0160 | 0.921 |

The slope equals Pr = 0.71 in the conduction sublayer (exact), follows the log
law only loosely through the log layer, reaches a **minimum at y⁺ ≈ 97 and
then turns back up**. A logarithmic profile would keep falling as 1/y⁺.

The mechanism is the last column. This thermal problem holds the TOTAL flux
constant across the channel, so

    dθ⁺/dy⁺ = Pr / (1 + Pr·D_t/ν),

and the turbulent fraction peaks at 0.936 near y⁺ = 100 and *falls* to 0.921
at the centreline — the eddy diffusivity drops off there while the flux may
not, so the molecular part has to grow and the gradient steepens. It is not a
wake, and it is not a log region.

**This is the same structural difference that makes panel (c)'s outer region
incomparable**: Flageul's flux falls to zero at the centre, so their
profile has no such constraint. The two setups can only be compared inside the
wall layer, which is what the shaded band marks and what
`check_scalar_turb.py` has always said about the Kader window.

### The conjugate signature, grid-converged

| dy+ | near-wall peak (y+ ~ 20) | wall, K = 1 | **wall/peak** |
|---|---|---|---|
| 1.5 | 7.07 | 1.367 | 0.1933 |
| 1.0 | 7.02 | 1.290 | 0.1839 |
| **0.75** | 7.14 | 1.262 | **0.1769** |
| Flageul | 6.3 | 1.10 | **0.1746** |

**The interface coupling converges monotonically under y-refinement to 1 % of
the reference.** The ratio is the right quantity: it is free of theta_tau and
of the fluid-side level, so it isolates what the conjugate wall does — damp
the interface fluctuation relative to the fluid's own. The brackets are
ordered correctly around it (K = 0.1 -> 0.523; K = 100 -> 0.0075 against their
isoT 0).

Note that dy+ = 1.5 IS too coarse at the wall, exactly as expected: it costs
10 % on the ratio (0.193 vs 0.177). It does not show up in the near-wall peak,
which is grid-converged to ~2 % across the three levels.

### What is left, and it is physics not numerics

The near-wall peak sits at **7.0-7.1 against 6.3, +11 %**, and it is converged
in all three directions and in time (campaign 2 has Flageul's full averaging
window; the y study adds dy+ 0.75). Two structural differences remain, neither
of which is the conjugate scheme:

* **the thermal problem is not theirs.** Ours is constant-flux-across-the-
  channel (antisymmetric walls); theirs is Kasagi's wall-flux formulation with
  a mean streamwise gradient, whose flux falls to zero at the centre.
  That changes the whole outer profile — visible above at y+ = 145, where we
  have 7.96 and they have ~0.8 — and feeds back on the near-wall budget
  through turbulent transport;
* **Re_tau 180 against their 149.**

Reproducing their formulation needs the Kasagi source term `f_T u_x` (a mean
streamwise temperature gradient), which this solver does not have. That, not
resolution, is what a like-for-like comparison of the near-wall peak requires.

### Internal checks that came out right

* **the fluid resistance is one number.** R_f/2 = dtheta_half/J measured
  independently from all six conjugate walls: **19.06 .. 19.46, spread 4.2 %**
  (campaign 2 reseed). It is a property of the flow, so six different solids
  agreeing on it checks the interface coefficient;
* **capacity does not touch the steady state.** theta_tau across the three
  kappa_s = 1 scalars: spread **1.7 %** (campaign 2), against 20.9 % in
  campaign 1 — see the reseed note below;
* **the velocity is a Re_tau 180 DNS**, validated against an EXACT law rather
  than a reference file: the steady total stress is linear, and campaign 2
  holds it to **0.24 % mean, 0.50 % max**, with U+_c = 18.20;
* **theta+** within 4-6 % of Kader over the log layer.

### The effusivity collapse fails, and the paper says why

Same K, different (kappa_s, C_s), campaign 2 wall/peak: K = 10 gives 0.0313
(k2) against 0.0201 (a2); K = 100 gives 0.0212 (k3) against 0.0075 (a3). The
design predicted these should collapse. They do not, and Flageul's eq. (12) is
the reason: the interface response carries
`R^2 = k_x^2 + k_z^2 + i k_t rho c/lambda_s`, so `lambda_s R` reduces to the
effusivity `sqrt(lambda_s rho_s c_s)` ONLY when the temporal term dominates.
Where lateral conduction matters the damping follows the CONDUCTIVITY, so the
higher-kappa_s wall must damp more at equal K — which is exactly the ordering
measured. Their physics; the design's claim was wrong.

### The reseed, and why it needs ONE fluid resistance

With a thick wall the solid holds most of the resistance and its mean relaxes
on d^2/alpha_s — t ~ 8e3 and 8e5 for C_s = 100 and 1e4. It is therefore SET,
not waited for (the capacity cannot change the steady state at all: C1 gate 1c
measured that to 1e-13). Campaign 1 used each scalar's OWN fluid-resistance
estimate, which for a high-capacity wall is contaminated by its frozen
transient, and left C_s = 100 and 1e4 running 12 % and 20 % off their steady
flux for the whole window. Campaign 2 measures R_f once and shares it: the
three kappa_s = 1 scalars then get identical J_steady and theta_i, and the
theta_tau spread falls to 1.7 %.


## Why the grid is uniform, and what it costs

The wall-normal spacing is uniform at dy⁺ = 1.5 — simultaneously too coarse at
the wall (it costs 10 % on the wall/peak ratio, see above) and finer than
needed at the centre. Two reasons, and the second is the one that matters for
the solver's future:

1. **It is not expressible otherwise.** `[grid.y] distribution` offers only
   analytic families over the whole domain (`uniform`, `cosine`, `tanh`,
   `natural`, `blayer`, `geometric`); there is no piecewise or file-based node
   line. In an extended conjugate domain the domain edges are the SOLID OUTER
   FACES, so `natural`/`tanh` would refine there — exactly backwards. A
   `channel_slab` distribution (uniform through each slab, stretched inside
   the gap, a node pinned to each interface) is the fix.
2. **Stretching would save cells but not time steps.** The explicit diffusive
   limit sets `dt` from the smallest cell ANYWHERE, so `dt ~ dy_wall^2`
   whether the grid is stretched or not:

   | dy+ at the wall | dt | steps for one 2.9e4-wall-unit window |
   |---|---|---|
   | 1.5 | 1.4e-3 | 112 k |
   | 0.75 | 3.6e-4 | 447 k |
   | 0.49 (Flageul) | 1.6e-4 | ~1.0 M |

   Flageul reach dy⁺ = 0.49 only because their wall-normal diffusion is
   **Crank–Nicolson**; ours is fully explicit. This is why the study above
   refines y at a SHORT window rather than running one converged case at their
   spacing, and it is the strongest argument in this whole campaign for
   implicit (or direction-split implicit) diffusion in the solver.

## Reproducing

```bash
./run_cht.sh prepare ic develop reseed stats   # campaign 2 (3.5 h on an RTX 5090)
./run_yconv.sh 1.5 1.0 0.75                    # the y study (4.8 h)
./check_cht.py velocity cht_vel_stats.h5 --interfaces 0.6 2.6
./check_cht.py thermal  cht_stats.h5      --interfaces 0.6 2.6
```

`run_corax.sh` drives the same phases on istmcorax (which needs its own cc120
build and an explicit PATH — it has no modulefile, and `/tmp` is NOT shared
between the hosts even though the home filesystem is).

---

# Campaign 3 — BULK HEATING (`cht180_bulk.ini`)

The constant-flux problem above differs from Flageul's *structurally*: their
source makes the wall-normal flux fall to zero at the centreline, ours holds it
constant across the channel. That is what made their global variance maximum a
near-wall peak and ours a centreline value, and it is the root of the
comparison error recorded further up.

This campaign removes that difference with **no solver code**: `[scalar.N]
source` already exists and the conjugate kernel weights it by the fluid
fraction, so a uniform `S = 1` heats only the fluid. Both outer solid faces are
held at zero instead of ±1. The steady mean balance is then `dJ/dy = S`.

The remaining difference from Flageul is the source's SHAPE — theirs follows
`u(y)` (Kasagi), ours is uniform. Measured against this run's own mean velocity
that is ≤ 7 % anywhere and 2.4–3.8 % over the near-wall comparison band,
against the 100 % difference at the centreline that constant flux gives.

## What the configuration buys

* **`J_wall = S h` EXACTLY, independent of κ_s.** All the heat generated in a
  half-channel leaves through that wall. So θ_τ is identical for all six
  scalars *by construction* rather than being an outcome of the wall/fluid
  series resistance — and the normalised variance cannot be contaminated by a
  θ_τ that is still converging.
* **The solid mean is closed-form**, `θ_i = J R_s`. The seed is already the
  answer, so the reseed phase is GONE — and with it the hardest practical
  problem of campaigns 1 and 2, where the `C_s = 10⁴` wall could never relax on
  its `d²/α_s ≈ 8e5` timescale. `run_bulk.sh` has two phases, not four, and no
  velocity development at all (it seeds from a campaign-2 snapshot on the
  identical grid).

## Result

Grid, case file and conjugate walls are campaign 2's; only the thermal problem
differs. 14 000 develop + 64 000 statistics steps (t = 20 → 115, ≈ 18 000 wall
units), local RTX 3060, 0.68 s/step.

| | ours | Flageul | |
|---|---|---|---|
| θ_τ vs the exact `S h` | 1.0019–1.0135 | — | max dev 1.35 % |
| `dJ/dy = S` | 1.38 % | — | J(centreline) = −0.00002 `S h` |
| ⟨θ'²⟩_wall, conjugate K = 1 | 1.20 | 1.10 | +9 % |
| ⟨θ'²⟩_wall, isoQ bracket | 3.97 | 4.20 | −5 % |
| ⟨θ'²⟩_wall, isoT bracket | 0.043 | 0.0 | |
| near-wall peak ⟨θ'²⟩, K = 1 | 4.90 | 6.30 | −22 % |
| wall/peak, K = 1 | 0.244 | 0.175 | +40 % |

**The near-wall peak is now the GLOBAL maximum** for 5 of 6 scalars, so it is
finally the same quantity Flageul report. And the reference is BRACKETED by our
two thermal problems: constant flux gives a peak of 7.07 (+12 %), bulk heating
4.90 (−22 %), Flageul 6.30 — exactly as panel (b) of `bulk_vs_constflux.png`
predicts from the three flux profiles, since the Kasagi flux lies between a
constant and our linear one. The wall value moves the same way and lands closer
(+9 % against campaign 2's +24 %).

So the residual disagreement is NOT a defect in the conjugate scheme; it is the
thermal problem. Closing it needs the Kasagi `f_T u_x` source (a mean
streamwise temperature gradient), which this solver has no term for.

## Two findings

**a1 (κ_s = C_s = 0.1) is a different regime and must not be pooled.** Its
variance never peaks near the wall — it rises to 13.4 at y⁺ ≈ 60 and stays high
— because a nearly-insulating wall leaves the bulk-heated core with no thermal
sink for large-scale structure. Verified as physics, not a bug: the flux split
`J_total = J_conv + J_diff` closes to **0.000 %** against the mean gradient
measured independently from the profile, and a1 carries LESS of the same total
flux by convection (0.577 vs k1's 0.611 at y⁺ = 60) and correspondingly more by
molecular diffusion — its large-scale core structures are poorly correlated
with the wall-normal velocity, so its eddy diffusivity is genuinely lower and
its mean gradient steepens to compensate. Its instantaneous *within-plane*
variance is 6.6 against k1's 0.30, stable across snapshots, so it is spatial
structure and not a slow drift.

**A REPORTING TRAP OF THE SAME FAMILY AS THE `max()` ONE.** Averaging the peak
over the sweep gave 6.15 against their 6.30 — a 2 % "agreement" that was pure
coincidence, produced by a1's 13.17 dragging up five scalars that all sit at
4.6–4.9. A scalar with no near-wall peak returns its value at the SEARCH BAND
EDGE, which is not the same quantity. `check_cht.py` now excludes such scalars
by name and reports the range plus k1 alone (k1 *is* their case). The general
lesson, twice learned: **never aggregate across cases before checking they are
the same quantity.**

## Reproducing

```bash
./run_bulk.sh ic develop stats                       # ~15 h on an RTX 3060
./check_cht.py thermal bulk_stats.h5 --source 1.0    # switches on the two bulk gates
./plot_bulk.py --vel cht_vel_stats.h5                # bulk_vs_constflux.png
```

`--source` also switches gate (1) to the interface-referenced (offset-free)
form. The absolute solid level carries a slowly-decaying charging transient —
observed to sort by the EFFUSIVITY K, not by the solid diffusion time (k2/a2
agree to 1 % and k3/a3 to 0.2 % with their timescales a factor 100 and 10⁴
apart), i.e. Flageul's eq. (12) showing up in the mean. θ⁺ and the variance are
both offset-free, so nothing that is compared with the literature depends on
it. `check_cht.py thermal` still exits FAIL on gate (2), the effusivity
collapse (15.8 %/39.7 %) — the known and documented physics, unchanged from
campaign 2, not a bulk-heating result.

## The main figure

`cht_validation.png` (regenerated by `./plot_cht.py`) now carries BOTH thermal
problems, since the comparison only means something with the two together:

* **(a)** the mean. A result that came out of drawing them together: the
  bulk-heated profile follows **Kader to within 4.3 %** of the centreline value
  over the whole profile (2.2 % beyond y⁺ 30), while the constant-flux one
  departs by **34 %** in the outer region. That is the reverse of what Kader's
  derivation suggests — the correlation is FITTED to channel and pipe data,
  whose flux decays toward the centre, so it describes the bulk-heated channel
  rather than the literally constant-flux one.
* **(b)** why: the constant-flux slope has a minimum at y⁺ = 97 and turns back
  up, whereas the bulk-heating slope keeps falling because `J → 0` at the
  centre.
* **(c)** the variance, with the K ladder under constant flux and the same
  K = 1 wall under bulk heating. Flageul's peak marker falls **between** the
  two curves — the bracketing, drawn.
* **(d)** unchanged, and now labelled as what it is: the wall-resolution
  convergence of the CONSTANT-FLUX campaign.

---

# The reference, digitised (2026-09-01)

Everything above rested on FOUR VALUES read off Flageul's figure 5 by eye,
because the paper tabulates nothing. `digitize_flageul.py` now extracts the
black solid "Conjug" curve itself into `flageul_fig5a_conjug.dat`.

**How, and how it is checked.** The figure is vector art, but `pdftocairo`
drops the figure-5 paths, so the page is rendered at 600 dpi and the curves are
separated BY COLOUR: the other three series are saturated red/green/blue, and
the Conjug line is the only achromatic dark ink inside the axes. Frames, ticks
and the legend rectangle are masked; axes are calibrated on the TICK MARKS
(found on an axis the curves do not touch), whose uniform spacing is asserted.
**The validation is free and is the point:** figure 5 plots the same curve
twice, panel (a) on a log x-axis and panel (b) on a linear one. Both are
digitised independently — different axis type, different calibration, different
legend mask — and agree to **0.015 max / 0.004 rms**, with the peak landing at
6.208 in both. A calibration slip would not survive that; the script asserts it.

## TWO OF THE EYEBALLED NUMBERS WERE WRONG, AND ONE MATTERED

| | eyeballed | digitised |
|---|---|---|
| peak ⟨T'²⟩ | 6.3 | **6.208** at y⁺ 17.6 |
| "wall" ⟨T'²⟩ | 1.1 | **1.270** at y⁺ 0.75 |

The peak was close. The wall value was not, and the reason is a height
mismatch, not a reading error: 1.1 is roughly where their curve meets the axis,
but **their first data point is at y⁺ = 0.49 and the curve is still climbing**.
Our first cell is at y⁺ = 0.75, where their curve reads 1.270. The old
comparison put our y⁺ = 0.75 cell against their y⁺ ≈ 0.49 value.

**This retracts the campaign-2 headline.** The wall/peak ratio was reported as
converging monotonically under y-refinement to 0.1769 against "their 0.1746",
i.e. 1 %. Both sides of that were wrong. Their ratio *at their own first point*
is 0.1933, and at a MATCHED y⁺ = 0.75 it is **0.2045**. Ours at the matched
height are 0.1953 / 0.1947 / 0.1911 for Δy⁺ = 1.5 / 1.0 / 0.75 — i.e. **4.5 to
6.6 % below, and not converging toward them.** The apparent convergence was an
artefact: our first cell moved (0.75 → 0.5 → 0.375) down the rising near-wall
curve, so the number fell for a reason that had nothing to do with accuracy.
Panel (d) now reads both sides at y⁺ = 0.75 and plots the old first-cell
measure in grey, labelled as not comparable.

## The full-profile comparison

Against the digitised curve, over 0.75 < y⁺ < 137:

| | at y⁺ 0.75 | near-wall peak | profile rms | bias |
|---|---|---|---|---|
| Flageul (digitised) | 1.270 | 6.208 @ y⁺ 17.6 | — | — |
| constant flux | 1.37 (+8 %) | 7.07 (+14 %) | 2.02 (32.5 % of their peak) | +1.24 |
| bulk heating | 1.20 (−6 %) | 4.89 (−21 %) | 0.83 (13.3 %) | −0.63 |

The bracketing survives the correction and is now visible as a curve in panel
(c) of `cht_validation.png` rather than inferred from three points: their curve
lies between ours over essentially the whole profile, the constant-flux run
biased high and the bulk-heating one low, with bulk heating **2.4× closer
overall**. What the digitisation changes is the WALL value — bulk heating is
−6 % there, not the +9 % reported against the eyeballed 1.1.

## The two ideal brackets, digitised too

`digitize_flageul.py` now extracts all three series — the conjugate case from
its solid line, and the **isoQ** (green ×, the K → 0 limit) and **isoT** (blue
+, the K → ∞ limit) brackets from their symbol series. Each is digitised off
BOTH panels of their figure 5 and its own 5a-vs-5b agreement is measured and
written into the `.dat` header:

| series | peak | at y⁺ 0.75 | 5a vs 5b |
|---|---|---|---|
| Conjug (line) | 6.208 @ y⁺ 17.6 | 1.270 | 0.015 max / 0.004 rms |
| isoQ (× markers) | 6.365 @ y⁺ 16.3 | 4.294 | 0.209 / 0.048 |
| isoT (+ markers) | 5.777 @ y⁺ 18.5 | (see below) | 0.359 / 0.083 |

A LINE is located to a fraction of a pixel; a SYMBOL CENTROID is much coarser,
and its bias differs between the panels because the same ~20 px marker spans
very different y⁺ widths on a log and on a linear axis. So the brackets are
good to ~0.2–0.4 in ⟨T'²⟩ where the conjugate line is good to 0.015, and the
script enforces a per-series tolerance rather than one number.

**isoT is UNRELIABLE below ⟨T'²⟩ ≈ 0.3, and that is not fixable.** A `+` spans
about ±0.08 in ⟨T'²⟩, so where the curve is that small the symbol straddles the
axis and its lower half is **clipped in their own figure**. That ink is not in
the image and no estimator recovers it: the digitised values there go
non-monotone (0.094, 0.079, 0.122 at y⁺ 0.5, 0.75, 1.0) instead of following
the y² an ideal Dirichlet wall must. The exact limit is known analytically
anyway — isoT → 0 at the wall, by definition of the boundary condition — so the
checker compares isoT at its PEAK and states the wall limit, and the figures
draw the band's lower edge only where it is data.

(Found on the way: an earlier version applied the 16 px frame margin to the
colour masks too, which ate the lower half of every near-axis symbol and biased
isoT's wall value from 0.079 up to 0.163. The margin exists to drop the black
frame and ticks, which a colour mask already excludes — so coloured series get
a 2 px clip instead. The symptom that exposed it was the checker reporting our
K = 100 wall as damping *more* than an ideal isothermal one, which is
impossible.)

With the brackets in, the sweep lands where it should: our K = 0.1 wall sits on
the isoQ edge (3.97 vs 4.29 at y⁺ 0.75, −8 %) and our K = 100 wall on the isoT
edge (peak 4.82 vs 5.78, −17 %) — the same ≈20 % peak deficit the conjugate
case shows, uniform across the sweep, which is the bulk-heating thermal problem
rather than anything K-dependent.

---

# Campaign 4 — LIKE-FOR-LIKE: Kasagi's source at Re_tau = 149

Campaigns 1–3 left two structural differences from Flageul et al., and the
bracketing in campaign 3 said those differences were the whole residual. This
campaign removes **both at once** and tests that claim.

* **The source.** `[scalar.N] source_type = velocity` implements Kasagi, Tomita
  & Kuroda's fully-developed channel source. For a uniform wall flux the mean
  temperature rises linearly in x, so writing `T = beta x + theta` with theta
  homogeneous in x turns `u.grad T` into `beta u + u.grad theta`: the extra
  term is proportional to the **instantaneous** streamwise velocity, not the
  mean. The mean flux then follows the CUMULATIVE FLOW RATE — their profile
  exactly, where campaign 3's uniform source only approximated it (measured:
  the Kasagi flux sits 3–19 % above the linear one).
* **The Reynolds number.** `re = 149`, matching theirs. The case file is
  re-prepared at that Re; verified dataset by dataset that the geometry
  (`blocks`, `dwall_blocks`, node lines) is IDENTICAL and only the coefficients
  move, by exactly 180/149.

## Result

60 000 statistics steps (t = 36 → 127, ~16 000 wall units), local RTX 3060.

| | ours | Flageul | |
|---|---|---|---|
| theta_tau, kappa_s-independent | 0.994–1.006 | — | spread 1.16 % |
| J(centreline) | +0.0002 | 0 | |
| ⟨θ'²⟩ at y⁺ 0.62 (matched height) | 1.38 | 1.24 | +12 % |
| near-wall peak, K = 1 | 5.85 | 6.208 | **−6 %** |
| **full profile, 0.62 < y⁺ < 137** | | | **rms 0.32 = 5.1 % of their peak** |
| isoQ bracket at y⁺ 0.62 | 4.42 | 4.30 | +3 % |
| isoT bracket, at its peak | 5.69 | 5.78 | −2 % |
| wall/peak, K = 1 | 0.236 | 0.199 | +19 % |

**The full-profile agreement across the three campaigns is 2.02 → 0.83 →
0.32**, i.e. 32.5 % → 13.4 % → **5.1 %** of their peak. Removing the two
structural differences removed most of the residual, which is what campaign 3's
bracketing predicted and is the point of running this case. Both ideal brackets
now agree to 3 % and 2 %, so the agreement is not specific to K = 1.

## What is left

The residual is a slightly FLATTER near-wall profile: +12 % at the wall against
−6 % at the peak, i.e. wall/peak +19 %. The most likely remaining cause is
wall-normal RESOLUTION — Δy⁺ = 1.24 here against their 0.49. They can afford
0.49 because their wall-normal diffusion is Crank–Nicolson; ours is fully
explicit, so Δy⁺ is capped by the diffusive stability limit. That is the same
constraint recorded at the end of the y-refinement study, and it is now the
single identified obstacle to closing the comparison further — the strongest
argument in this whole campaign for implicit (or direction-split implicit)
diffusion.

Secondary: our solid is d⁺ = 89 against their 149, so our wall clips more of
the low-frequency interface response. That damps θ'_wall, i.e. it has the
WRONG SIGN to explain a +12 % wall value, so it is not the leading candidate.

## Reproducing

```bash
./run_kasagi.sh ic ; ./run_kasagi.sh develop ; ./run_kasagi.sh stats   # ~16 h
./check_cht.py thermal kasagi_stats.h5 --source 1.0 --re 149 --kasagi
./plot_cht.py             # cht_validation.png  -- the 4-panel comparison
./plot_kasagi.py          # kasagi_vs_flageul.png -- flux precondition + variance
```

**BOTH FIGURES NOW SHOW ONLY THIS ONE-TO-ONE COMPARISON.** The constant-flux
and bulk-heating campaigns bracketed the reference and are what motivated this
run, but they are a different thermal problem, so drawing them beside it
invites exactly the apples-to-oranges reading this README has had to correct
twice. They stay documented above and are no longer plotted;
`bulk_vs_constflux.png` is deleted and `plot_bulk.py` is now `plot_kasagi.py`.

`cht_validation.png`: (a) mean vs Kader, (b) the kappa_s sweep against their
digitised conjugate curve and isoT-isoQ band, (c) the same comparison as an
ABSOLUTE deviation -- the only way to read a 5 % agreement off a log plot --
and (d) the interface response vs effusivity, which is what shows the agreement
is not a K = 1 coincidence. `kasagi_vs_flageul.png`: (a) the flux profile,
which is the PRECONDITION -- our measured flux against A*int(<u> dy) from the
run's own mean velocity, max departure 0.005 -- and (b) the variance on a
linear axis, where 5 % is actually visible.

`--kasagi` relabels gate (0) and drops the linear `dJ/dy = S` test, which does
not apply when the flux follows the flow rate; `J(centreline) = 0` and a
kappa_s-independent theta_tau still do.

**TWO TRAPS, both found here.**

1. **A restart file's `re` WINS over the ini.** Seeding a Re = 149 run from a
   Re = 180 snapshot silently runs at 180. It surfaced only because the IBM
   coefficient file is Re-specific and was then rejected as mismatched —
   without that cross-check the run would have looked healthy at the wrong
   Reynolds number. `make_bulk_ic.py --re` now stamps it (and fixes `R_s` in
   the seed with it).
2. **The velocity `stats_file` was a HARD-CODED name** in all three cht inis,
   so this campaign APPENDED its samples to campaign 3's rather than starting
   fresh, mixing two Reynolds numbers in one file. Now `@PREFIX@_vel.h5`, per
   run. Campaign 3's velocity statistics were lost this way; they were a
   diagnostic, not a committed result (their conclusions are recorded above),
   and the Re = 149 flow is certified from snapshots instead: the exact steady
   total-stress law holds to 3.0 % over 3 instantaneous fields (campaign 3's
   0.24 % came from thousands of time samples) with U⁺_c = 18.3.

---

# Campaign 5 — FLAGEUL-MATCHED (`cht149_flageul.ini`)  [DONE]

Campaign 4 left one identified obstacle: wall-normal resolution. This campaign
removes it, and the two remaining non-resolution differences with it.

| | Flageul | campaign 4 | **campaign 5** |
|---|---|---|---|
| Δy⁺ wall → centre | 0.49 → 4.8 | 1.24 uniform | **0.490 → 4.87** |
| Δx⁺ | 14.8 | 11.7 | **8.36** |
| Δz⁺ | 5.1 | 4.18 | **4.18** |
| d⁺ (solid) | 149 | 89 | **149** |
| outer solid face | imposed flux | Dirichlet | **imposed flux** |

224³ cells, 343 leaves, nb = 32.

**The y line comes from a file** (`make_ynodes.py` → `ynodes_f149.dat`, read by
the new `[grid.y] nodes_file`): every built-in distribution clusters at the
DOMAIN ENDS, and this case needs the two interior fluid/solid interfaces
resolved AND landing exactly on cell faces. Solid 53 + fluid 118 + solid 53;
the fluid is a symmetric tanh solved on the REALISED first cell rather than the
analytic derivative — they differ by 3 % at this stretching, and the realised
value is what sets both the resolution and the diffusive time-step limit.

**The imposed-flux BC needed no solver code.** `apply_scalar_bc_q`'s
non-Dirichlet branch is already a general nonzero Neumann, `ghost = interior +
dn*value`, so `bcValue` is dθ/dy and an imposed flux q is `q/(kappa_s D_f)` —
per scalar, since kappa_s spans four decades. Measured: J at the outer face is
**−1.0000 for every scalar**, exactly the prescribed value.

That makes the scalar problem pure Neumann, so the temperature level is free.
The seed therefore puts **θ = 0 at the interface**, which keeps the
kappa_s = 0.1 scalar O(10) in the fluid instead of the O(600) the Dirichlet
version carried — the same physics, far better conditioned for a variance
computed as ⟨s²⟩ − ⟨s⟩².

## Cost, and why

dt = 3.07e-4 against campaign 4's 1.63e-3: the 0.49 wall spacing costs 6.4× in
dt, because the explicit diffusive limit is dy². Develop took 80 000 steps
(12.5 h at 0.561 s/step); the statistics window is 311 000 steps ≈ 48 h for the
SAME 95.5 time units campaign 4 averaged over. Local RTX 3060.

## Reproducing

```bash
./make_ynodes.py                                   # the y line
mpirun -n 8 ../../../build_cpu/moby_prepare .fprep_in.ini cht149_flageul.h5
./run_flageul.sh ic ; ./run_flageul.sh develop ; ./run_flageul.sh stats
./check_cht.py thermal flageul_stats.h5 --source 1.0 --re 149 --kasagi \
    --interfaces 1.0 3.0
./plot_cht.py ; ./plot_kasagi.py                   # both default to this case
```

The plot scripts now take `--interfaces`, `--re` and (plot_kasagi) `--source`
and `--snapshots`, so the earlier campaigns are still reproducible:

```bash
./plot_cht.py --stats kasagi_stats.h5 --re 149 --interfaces 0.6 2.6
./plot_kasagi.py --stats kasagi_stats.h5 --interfaces 0.6 2.6 --re 149 \
    --source 0.064134 --snapshots 'kstat_*0000.h5'
```

## TWO OPERATIONAL TRAPS, both hit here

1. **A stale binary silently ran a DIFFERENT GRID.** `build_cpu` was rebuilt
   after `nodes_file` was added and `build_gpu` was not, so the solve ignored
   the key and fell back to the channel case's natural distribution while
   `moby_prepare` used the file line. **Nothing caught it**: the
   coefficient-file check compares lx/ly/lz/re and the block table, but NOT the
   node lines. It surfaced only as an absurd time step (3.6e-6 instead of
   3.1e-4). The case file already stores x/y/z_nodes, so the solver could
   cross-check them at read — **that guard is not implemented**.
2. **mpirun's exit code is not a success test.** The develop leg finished
   cleanly, wrote its snapshot, and then failed at MPI teardown with PMIX
   `NO-PERMISSIONS` errors (a /tmp permissions problem). The driver's
   `|| exit 1` guard fired on that and the chained statistics leg never
   started — 12.5 h of correct work sat idle. `run_flageul.sh` now tests for
   `main loop ended` IN THE LOG instead.

## Result (2026-09-07)

311 000 statistics steps, 48.6 h at 0.563 s/step, t = 24.5 → 119.9 (the same
95.5 time units campaign 4 averaged over).

| | ours | Flageul | |
|---|---|---|---|
| theta_tau, kappa_s-independent | 1.0008–1.0034 | — | spread **0.25 %** |
| J(centreline) | −0.00036 | 0 | |
| **full profile, 0.49 < y⁺ < 137** | | | **rms 0.197 = 3.2 % of their peak** |
| near-wall peak, K = 1 | 5.99 | 6.208 | **−4 %** |
| ⟨θ'²⟩ at y⁺ 0.49 | 1.25 | 1.20 | +5 % |
| isoT bracket, at its peak | 5.69 | 5.78 | −2 % |
| isoQ bracket at y⁺ 0.45 | 3.58 | 4.28 | −16 % |
| wall/peak, K = 1 | 0.2096 | 0.1933 | +8 % |

**The four campaigns, full-profile rms against their digitised curve:**

| | rms | % of their peak | k1 peak |
|---|---|---|---|
| constant flux, Re_tau 180 | 2.02 | 32.5 % | +14 % |
| bulk heating, Re_tau 180 | 0.83 | 13.4 % | −21 % |
| Kasagi, Re_tau 149 | 0.31 | 5.1 % | −6 % |
| **Flageul-matched** | **0.197** | **3.2 %** | **−4 %** |

## What the resolution match bought, beyond the 5.1 % → 3.2 %

Two QUALITATIVE artefacts of the coarser grids disappeared, and both were
visible in the earlier figures:

* **Every K now peaks near the wall.** On the coarser grids K = 0.1 rose
  monotonically to the centreline and had to be excluded from the peak
  statistics by name. Here all 6 scalars peak at y⁺ 16–19 and the peak is the
  global maximum for 6/6.
* **K = 0.1's wall value is now BELOW the ideal-isoQ limit**, where campaign 4
  put it 3 % ABOVE. A finite-effusivity conjugate wall cannot exceed the
  K → 0 plateau, so the earlier +3 % was unphysical — an artefact of reading a
  under-resolved wall cell. The −16 % now reported is the finite-K deficit,
  which is the right side of the limit; whether −16 % is the RIGHT size is not
  settled, since K = 0.1 is an approximation to isoQ and not isoQ.

## Two open items

* **The residual is a slightly FLATTER profile**: +5 % at the wall, −4 % at the
  peak, wall/peak +8 %. It has survived every refinement, so it is no longer
  attributable to resolution — Δy⁺ is now theirs, Δx⁺ and Δz⁺ are FINER, and
  the outer BC and solid thickness match. What remains unmatched is the
  DOMAIN: ours is 1872 × 936 wall units against their 3814 × 1270, i.e. half
  the streamwise extent. Their own §5 attributes the near-wall conjugate
  behaviour to very large-scale thermal structures, and reports that the
  temperature variance is sensitive to the box length — so the domain is the
  first thing to try, not the grid.
* **4.8 % of the interface variance survives at our outer solid face**, where
  their figure 12 has ~1e-3 mid-solid. Same thickness, same BC now. The likely
  cause is statistical: they averaged 29 000 wall units and say explicitly it
  was to converge the DEEP SOLID, where the correlation times are longest;
  ours is 14 200. Not resolved.

`check_cht.py` still exits FAIL on gates (1) and (2): the absolute solid level
(the charging transient, offset-free form 6.0e-2 against a 5e-2 tolerance) and
the effusivity collapse (22 %/38 %). Both are documented physics, unchanged.

## The Kader comparison, and why the profile sits above it

Panel (a) shows θ⁺ about 9 % above Kader (1981) at the centreline, and the log
slope steeper (2.80 against Kader's 2.27 over y⁺ 30–100). That is NOT a scalar
error, and the check is the momentum profile of the same run:

| y⁺ | θ⁺ − Kader | U⁺ − (2.44 ln y⁺ + 5) |
|---|---|---|
| 29 | +0.67 | +0.78 |
| 51 | +1.07 | +1.11 |
| 100 | +1.33 | +1.25 |
| 147 | +1.04 | +0.96 |

**The two departures are the same curve to within ~0.1.** Kader is a HIGH-Re,
constant-flux correlation; at Re_τ = 149 the FLOW itself sits about one wall
unit above the high-Re log law (measured wake +0.92, dU⁺/d ln y⁺ = 2.74 against
2.44), and the thermal profile simply inherits it. So we agree with Kader to
the extent Kader can be expected to apply, and the scalar-specific residual is
**≤ 0.1 θ_τ**. The inset in panel (a) plots the two departures together; they
overlay.

Two smaller contributions, both real but subordinate:
* the conduction sublayer is EXACT — at y⁺ 0.7 and 1.9 the two agree to 0.007;
* the largest LOCAL departure is at y⁺ ≈ 5 (+0.51), the buffer layer, where
  Kader's two-branch blend is weakest;
* our flux DECAYS (Kasagi source) where Kader assumes it constant. That works
  in the opposite direction — a decaying flux lowers θ⁺ — so it partly cancels
  the low-Re excess rather than causing it.

## Figure improvements (2026-09-07)

* `kasagi_vs_flageul.png` panel (b) gains a RESIDUAL STRIP sharing its x axis.
  Two curves agreeing to 3 % are indistinguishable on the main axes, which is
  the point but also hides both the SIZE and the SHAPE of the disagreement;
  the strip shows both, and is where the sign structure (high at the wall, low
  at the peak) is read.
* `cht_validation.png` panel (a) gains the inset described above, so the
  Kader question is answered in the figure rather than only in this file.

---

# The RAW data (2026-09-07) — the comparison stops being digitised

Flageul et al. publish their statistics: **repo.ijs.si/CFLAG/incompact3d**,
folder `gxay`. `fetch_flageul.py` downloads and converts them (Scilab
`write_csv` writes COMMA decimals and TAB separators, which numpy reads as
nonsense if you let it). Column meanings are taken from their `conjug.sce`
lines 494–502, not guessed.

Their naming maps onto ours as: `gXaY` has X = G = fluid-to-solid DIFFUSIVITY
ratio → α_s = 1/G, and Y = 1/G₂ with G₂ the solid-to-fluid CONDUCTIVITY ratio →
κ_s = 1/Y. So C_s = G/Y and K = √G/Y. **Their `g1a1` is κ_s = α_s = C_s = K = 1
— exactly our k1.**

## FIRST: the digitisation was biased, and every earlier "we are low" was overstated

| | digitised (mine) | RAW (theirs) |
|---|---|---|
| peak ⟨T'²⟩ | 6.208 | **5.941** |
| at y⁺ 0.5 | 1.200 | **1.142** |

Systematically HIGH: rms 0.15, **bias +0.13**, 2.5 % of the peak. The
line-tracing centroid sat above the curve. Consequence: campaign 5's headline
"−4 % on the peak" was measured against a curve that was itself 4.5 % high —
the true figure is **+0.7 %**. The `.dat` files stay in the tree for
provenance; **nothing should be compared against them now.**

## The detailed comparison (`compare_flageul.py`, `flageul_detail.png`)

Ours vs their g1a1, over 0.5 < y⁺ < 137:

| quantity | rms | max | bias | % of that quantity's peak |
|---|---|---|---|---|
| U⁺ | 0.095 | 0.179 | +0.061 | **0.5 %** |
| θ⁺ | 0.165 | 0.287 | −0.158 | **1.1 %** |
| u'² | 0.066 | 0.114 | −0.054 | **1.0 %** |
| v'² | 0.014 | 0.031 | −0.009 | **2.2 %** |
| w'² | 0.046 | 0.080 | −0.037 | **4.3 %** |
| −u'v' | 0.004 | 0.009 | −0.003 | **0.6 %** |
| u'T' | 0.067 | 0.146 | −0.039 | **1.1 %** |
| **T'²** | 0.130 | 0.186 | +0.034 | **2.2 %** |

This is the first time the VELOCITY has been checked against a trustworthy
reference — the in-repo `channel_kmm180_stats.h5` is this solver's own
unconverged output (see the warning above). u'², v'², w'² and −u'v' agree to
0.6–4.3 %, so the flow is right, and the thermal agreement is not resting on a
compensating error.

## THEIR OWN DATA DOES NOT COLLAPSE ON THE EFFUSIVITY

Their nine cases include two pairs at equal K:

* K = 1.414: `g05a05` 0.781 vs `g2a1` 0.888 — **14 % apart**
* K = 0.707: `g05a1` 1.466 vs `g2a2` 1.561 — **6.5 % apart**

So "same effusivity ⇒ same interface response" is not exact in the REFERENCE
either. `check_cht.py` gate (2) has been failing a 10 % tolerance with 22 %/38 %
and calling it documented physics; this is the documentation. Our spread is
larger because our κ_s spans four decades where theirs spans a factor 4.

Across their nine cases our interface variance is **mean +1.0 %, rms 8.0 %,
max 13 %** — scatter comparable to their own non-collapse.

## THE ONE REAL DISCREPANCY: the variance INSIDE the solid

| depth −y⁺ | ours | theirs | ratio |
|---|---|---|---|
| 0 (interface) | 1.147 | 1.042 | 1.10 |
| 20 | 0.177 | 0.110 | 1.6 |
| 77 | 0.038 | 0.011 | 3.4 |
| 145 (outer face) | 0.059 | 0.007 | **8.0** |

Ours decays too slowly and then **turns back up** at the insulated outer face.
Two hypotheses are already eliminated:

* **not a charging transient** — the solid MEAN is at its exact steady linear
  profile, slope 0.70986 against their 0.71019 (**0.05 %**);
* **not the boundary condition or the thickness** — both are matched now.

What is left: our solid is 2nd-order FD on a grid coarsening to Δy⁺ 8.4, where
theirs is SPECTRAL wall-normal on a Chebyshev grid clustered at both ends; and
our averaging is 14 200 wall units against their 29 000, which their §2 says
was chosen specifically to converge the deep solid, where correlation times are
longest. The upturn at the Neumann face looks like a reflected low-frequency
mode that has not been averaged out. **Open.**

## The deep-solid discrepancy is CONVERGENCE, not the grid (2026-09-07)

Before refining the solid grid, the obvious question is whether the excess is
spatial structure a finer grid could resolve. It is not. Compare the
time-ACCUMULATED variance with the INSTANTANEOUS within-plane variance from the
snapshots at the same depth:

| depth y⁺ | accumulated | instantaneous (per snapshot) | Flageul | accum/inst |
|---|---|---|---|---|
| −5 | 0.606 | ~0.62 | 0.474 | 1.0 |
| −40 | 0.064 | ~0.054 | 0.034 | 1.2 |
| −102 | 0.046 | ~0.014 | 0.0084 | 3.3 |
| −145 | 0.059 | ~0.012 | 0.0073 | **4.9** |

**The spatially-resolved variance is fine at every depth** — at the outer face
it is 0.012 against their 0.0073, the same order. What grows with depth is the
GAP between it and the time-accumulated value, and that gap is the temporal
wander of the plane mean: the very-low-frequency modes whose correlation time
grows with depth and which a 14 200-wall-unit window cannot average out. A
finer solid grid cannot fix a temporal sampling problem.

This is also exactly what Flageul et al. warn about: their §2 says the
statistics were accumulated "over an interval significantly longer than in
previous studies (about five times longer)" specifically "to ensure a
satisfactory statistical convergence deep inside the solid domain".

**So the fix is a longer window, and `run_flageul.sh extend` is it**: it
restarts from the last snapshot WITHOUT deleting the statistics file, and
`read_stats_restart` reads that file, so the sums carry over. +311 000 steps
doubles the window to ~28 400 wall units — matching their 29 000 — for 48.6 h
against the 88 h a solid-grid refinement would cost.

**The grid refinement remains available and is scoped**: ny 224 → 320 (101
solid cells each side instead of 53). dt is UNCHANGED, because the minimum
spacing is at the INTERFACE and that does not move; only the cell count grows,
16.1 M against 11.2 M, so develop + statistics ≈ 88 h. Worth doing only if the
longer window fails to close the gap.

## IT WAS THE TRANSIENT, and the fix is a three-leg run (2026-09-09)

The "convergence" diagnosis above was only half right. Decomposing the
deep-solid variance into a SPATIAL part (within-plane) and a TEMPORAL part
(variance of the plane mean across snapshots) separates a sampling problem from
a transient:

| depth y⁺ | window | spatial | temporal | total | Flageul |
|---|---|---|---|---|---|
| −145 | all, t > 30 | 0.0109 | **0.0227** | 0.0336 | 0.0073 |
| −145 | t > 120 only | 0.0083 | **0.0005** | **0.0088** | 0.0073 |
| −100 | all | 0.0129 | 0.0145 | 0.0273 | 0.0085 |
| −100 | t > 120 only | 0.0091 | 0.0006 | **0.0097** | 0.0085 |

**The temporal term falls by a factor 45** when the transient is dropped, and
the outer-face value goes from 4.6× the reference to **1.2×**. Not 1/T sampling
noise, and not the box.

**Why it was there.** The solid diffusion time is d²/α_s = 1/D_f ≈ **106 time
units**, and the statistics window was t = 24.5 → 120 — LESS THAN ONE — started
from a develop leg of only 24.5. The solid was still relaxing throughout. The
snapshot means show it directly: at y⁺ = −145 the plane mean moves −102.26 →
−101.65 over t = 31 → 92, then creeps.

**A second, smaller cause, now fixed.** The residual creep after t = 120 is a
generation/removal imbalance: `source` had been set from the PREVIOUS run's
bulk velocity, so generation exceeded the imposed outflux by **0.30 %**. The
predicted level drift, −7.5e-4 per time unit, matches the measured −9.3e-4.
`source` is now 0.064932, from this run's own ∫⟨u⟩dy.

### The three legs

1. **`reseed_solid_exact.py`** — put every solid on its EXACT steady mean, in
   closed form. At steady state the solid profile is linear with slope
   J/(κ_s D_f) and its LEVEL is a pure-Neumann null-space mode carrying no
   physics, so the transient is removed by replacing the mean and leaving the
   fluctuations alone (subtracting a y-dependent mean leaves θ' untouched, and
   the operator is linear). *This is leg 1's goal reached without running it:*
   a capacity-accelerated leg would cost ~76 h and would still leave the
   κ_s = 1, C_s = 10⁴ scalar unequilibrated, since its α_s = 10⁻⁴ gives
   d²/α_s ~ 10⁶ time units. Measured shifts: 0.10–0.49 θ_τ.
2. **`run_flageul.sh settle`** — 65 000 steps (t → 20) at the CORRECT
   capacities, statistics off, with the corrected source.
3. **`run_flageul.sh clean`** — statistics FROM ZERO, 650 000 steps, a
   **200-time-unit window = 29 800 wall units**, matching their 29 000.

715 000 steps ≈ **112 h (4.7 days)**. The contaminated statistics are kept as
`flageul_stats.h5` / `flageul_stats_311k.h5` for the before/after.
