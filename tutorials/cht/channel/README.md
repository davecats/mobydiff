# Conjugate heat transfer in a turbulent channel

Validation of mobydiff's **conjugate heat-transfer scheme** on a **flat,
grid-aligned** solid–fluid interface, against

> S. Flageul, S. Benhamadouche, É. Lamballais, D. Laurence, *DNS of turbulent
> channel flow with conjugate heat transfer: Effect of thermal boundary
> conditions on the second moments and budgets*, Int. J. Heat Fluid Flow 55
> (2015) 34–44 — their published raw profiles, in [`asset/flageul_data/`](asset/flageul_data/).

This is the companion of [`../pipe/`](../pipe/), and the two are deliberately
complementary. **Here the interface sits on a cell face**, so the two
straddling cell centres are symmetric about it and the level-set weight is
w = ½ to the last bit: the cut-face coefficient is then the textbook
two-resistance harmonic mean and is *exact*. Anything this case measures is
therefore the conjugate **physics**, not the immersed-boundary approximation.
The pipe adds the curved interface on top of that.

---

## The case

| | |
|---|---|
| flow | turbulent channel, Re_τ = 149, Pr = 0.71, fixed pressure gradient |
| domain | 4π × 4.0 × 2π, fluid gap exactly 2.0 at 1.0 < y < 3.0 |
| solid | conjugate slabs d⁺ = 149 on both sides (= Flageul's), **imposed heat flux** on the outer faces |
| grid | 224³ = 11.2 M cells, 343 leaves, `nb = 32` |
| spacing | Δy⁺ 0.490 → 4.87 wall to centre (theirs 0.49 → 4.8), Δx⁺ 8.36 (14.8), Δz⁺ 4.18 (5.1) |
| heating | Kasagi streamwise reduction: source ∝ u_x, i.e. **their** source, not an approximation to it |
| cost | 0.562 s/step on an RTX 3060; the campaign is 715 000 steps ≈ 112 h |

**Six conjugate scalars ride one velocity field.** Three share κ_s = 1 with
C_s = 1, 100, 10⁴ (a capacity sweep, whose steady state must be identical);
three sweep κ_s = 0.1, 10, 100 with C_s = κ_s (a conductivity sweep). K is the
fluid-to-solid effusivity ratio.

| scalar | k1 | k2 | k3 | a1 | a2 | a3 |
|---|---|---|---|---|---|---|
| κ_s | 1 | 1 | 1 | 0.1 | 10 | 100 |
| C_s | 1 | 100 | 10⁴ | 0.1 | 10 | 100 |
| K | 1 | 10 | 100 | 0.1 | 10 | 100 |

`k1` is exactly their published case g1a1 (κ_s = α_s = C_s = 1).

### Three things the case needs, and how each is met

**The y line comes from a file.** Every built-in distribution clusters at the
*domain ends*; this case needs the two **interior** fluid/solid interfaces
resolved *and* landing exactly on cell faces. `make_ynodes.py` writes
`ynodes_f149.dat` (solid 53 + fluid 118 + solid 53) and `[grid.y] nodes_file`
reads it. The fluid part is a symmetric tanh solved on the **realised** first
cell rather than the analytic derivative — they differ by 3 % at this
stretching, and the realised value is what sets both the resolution and the
diffusive time-step limit.

**The imposed-flux outer face needed no solver code.** `apply_scalar_bc_q`'s
non-Dirichlet branch is already a general nonzero Neumann,
`ghost = interior + dn·value`, so `bcValue` is dθ/dy and an imposed flux q is
`q/(κ_s D_f)` — per scalar, since κ_s spans four decades. Measured: J at the
outer face is **−1.0000 for every scalar**, exactly the prescribed value.

**That makes the scalar problem pure Neumann**, so the temperature level is
free. The seed therefore puts θ = 0 *at the interface*, which keeps the
κ_s = 0.1 scalar O(10) in the fluid instead of the O(600) a Dirichlet version
carried — the same physics, far better conditioned for a variance computed as
⟨s²⟩ − ⟨s⟩².

---

## Results against Flageul et al.

Window t = 20 → 219 = **29 684 wall units**, matching their 29 000; 650 000
steps from a settled state (see "the three legs" below).

### The full profiles, against their raw data

`asset/compare_flageul.py`, over 0.5 < y⁺ < 137, as a percentage of each
quantity's own peak:

| quantity | rms deviation |
|---|---|
| U⁺ | **0.4 %** |
| θ⁺ | 1.2 % |
| u′² | **0.9 %** |
| v′² | 1.8 % |
| w′² | 3.8 % |
| −u′v′ | **0.4 %** |
| u′T′ | 1.3 % |
| T′² | 2.5 % |

**The near-wall variance peak, which is the conjugate quantity:** ours 5.981
at y⁺ 17.0 against their raw **5.941** at y⁺ 17.85 — **+0.7 %**.

> Use their **raw** data, not the digitised figure. The digitisation of their
> Fig. 5a was biased by rms 0.15 / bias +0.13, and the peak it gives (6.208)
> is 4.5 % high — every earlier "we are low by 4 %" in the campaign notes was
> that bias, not the solver. `asset/flageul_data/` is their published data;
> `asset/flageul_fig5a_*.dat` is the digitisation, kept only for the two
> **ideal** brackets (isothermal / isoflux), which their repository does not
> publish.

### Through the solid

θ′² as a fraction of the interface value, against depth, for k1 — and what a
too-short window does to it:

| depth y⁺ | contaminated window | **clean** | Flageul | clean/theirs |
|---|---|---|---|---|
| −5 | 0.618 | 0.593 | 0.485 | 1.22 |
| −40 | 0.064 | 0.044 | 0.034 | 1.28 |
| −77 | 0.038 | **0.0113** | 0.0113 | **1.00** |
| −102 | 0.047 | 0.0075 | 0.0083 | 0.90 |
| −145 (outer face) | 0.0586 | **0.0063** | 0.0073 | **0.86** |

At the outer face the ratio goes from **8.0 to 0.86**. What remains is the
solid within ≈ y⁺ 40 of the interface reading 22–28 % high, which is the
fluid side's near-wall excess conducted inward, not a separate problem.

### Internal consistency checks

`check_cht.py`, which needs no reference data:

* **θ_τ is κ_s-independent** — it must be, under the Kasagi source with h = 1 —
  to **0.019 %**, spread 0.023 % across a sweep in which κ_s spans four decades;
* **J(centreline) = +0.00023 S·h**, i.e. vanishes;
* the wall-normal flux lies 8.0 % above the linear profile, which is the
  Kasagi-vs-uniform-source gap and not an error;
* θ′_wall/θ_τ is **monotone decreasing in K** across the whole sweep
  (1.90 → 1.10 → 0.33 → 0.13), the conjugate S-curve.

### What the gate reports as FAIL, and why that is correct

`check_cht.py` prints **FAIL**, on two checks, both understood:

1. **The effusivity collapse.** Two walls with the same K but different
   (κ_s, C_s) do not collapse: 20 % apart at K = 10, 38 % at K = 100. **The
   design's premise was wrong and their paper says why** — Flageul's eq. (12)
   carries `R² = k_x² + k_z² + i k_t ρc/λ_s`, so `λ_s R` reduces to the
   effusivity √(λ_s ρ_s c_s) *only* when the temporal term dominates. Where
   lateral conduction matters the damping follows the **conductivity**, so the
   higher-κ_s wall must damp more at equal K — which is exactly the ordering
   measured.
2. **The capacity sweep's solid mean** differs by 7.3 % of θ_τ across
   κ_s = 1, C_s = 1/100/10⁴. The steady mean is capacity-independent (gated to
   1e-13 in the C1 suite), so this is a frozen transient: d²/α_s is ~8 × 10³
   and ~8 × 10⁵ time units for C_s = 100 and 10⁴. `reseed_solid_exact.py`
   sets those means in closed form rather than waiting for them; 7.3 % is what
   survives.

Neither touches the comparison above, which is k1 (C_s = 1).

### What is left

⟨T′²⟩ is **+11 % in the first fluid row** (y⁺ ≈ 0.45: ours 1.25 against their
1.13) and wall/peak is +10 %, i.e. a slightly flatter
near-wall profile, while the peak itself is +0.7 % and the whole profile 2.5 %.
It survived matching Δy⁺, the outer boundary condition, the solid thickness
and the averaging window, so the remaining candidate is the **domain**:
1872 × 936 wall units against their 3814 × 1270.

---

## Reproducing

```bash
module load toolkits/nvhpc/25.9
export PY=$HOME/ibmc/bin/python
export BIN=$PWD/../../../build_gpu/moby_solve
export PREP=$PWD/../../../build_cpu/moby_prepare

./make_ynodes.py                  # -> ynodes_f149.dat (committed; this rewrites it)
./run_flageul.sh prepare          # -> cht149_flageul.h5 from the two wall slabs
SEED=<a developed Re_tau 149 channel snapshot> ./run_flageul.sh ic
./run_flageul.sh develop          #  80 000 steps, statistics off
# put every solid on its EXACT steady mean (source: the last statistics
# snapshot of the contaminated window; the leg below restarts from IC_clean.h5)
$PY ./reseed_solid_exact.py Fstat_640000.h5 IC_clean.h5
./run_flageul.sh settle           #  65 000 steps at the correct capacities
./run_flageul.sh clean            # 650 000 steps, statistics from zero
```

**The `ic` leg needs a developed velocity field**, and the campaign took one
from the earlier Re_τ 149 Kasagi run on a different grid (160 × 384 × 224
uniform) — `make_flageul_ic.py` interpolates it onto this grid (the fluid gap
is the same size, so the map is a shift plus trilinear interpolation) and
**seeds the scalars analytically** rather than interpolating them. Any
converged Re_τ 149 channel field will do; there is no small one to ship, so a
user starting from nothing must develop one first.

### The three legs, and why there are three

The first statistics campaign was **contaminated**. The solid diffusion time
is d²/α_s ≈ 106 time units and the window was t = 24 → 120 — *less than one
such time*, started from a develop leg of only 24. Decomposing the deep-solid
variance into a spatial (within-plane) and a temporal (plane-mean-across-
snapshots) part separates sampling from transient, and **the temporal term
falls by a factor 45** when the transient is dropped.

1. `reseed_solid_exact.py` — put every solid on its exact steady mean, in
   closed form. At steady state the solid profile is linear with slope
   J/(κ_s D_f) and its *level* is a pure-Neumann null-space mode carrying no
   physics, so replacing the mean and leaving the fluctuations alone removes
   the transient exactly (subtracting a y-dependent mean leaves θ′ untouched,
   and the operator is linear). A capacity-accelerated leg would have cost
   ~76 h and *still* left the C_s = 10⁴ scalar unequilibrated.
2. `settle` — 65 000 steps at the **correct** capacities, statistics off.
   `source` was also corrected here: it had been set from the previous run's
   bulk velocity, so generation exceeded the imposed outflux by 0.30 % and the
   level drifted at −9.3e-4 per time unit (predicted −7.5e-4).
3. `clean` — statistics **from zero** over one uncontaminated window.

## Reproducing the comparison

```bash
cd asset
./compare_flageul.py --stats flageul_clean.h5 --vel Cstat_vel.h5   # the table above
./plot_cht.py        --stats flageul_clean.h5 --vel Cstat_vel.h5
cd .. && ./check_cht.py thermal asset/flageul_clean.h5 \
    --source 1.0 --re 149 --kasagi --interfaces 1.0 3.0
```

`asset/plot_wall_fluctuation.py` draws θ′ in the first fluid cell as a plan
view for three walls that differ only in effusivity. It is the one script that
needs a **snapshot** (989 MB), so it cannot run from the shipped assets;
`asset/figures/wall_fluctuation.png` is its output. The amplitude collapses —
rms 1.98 → 1.11 → 0.074 θ_τ from K = 0.1 to 100, a factor 27 — while the
streak *shape* does not change: the same turbulence writes the same footprint
on every wall, and only the wall's effusivity decides how much survives.

## Files

**Here** (to run the case): `run_flageul.sh` (one leg per subcommand),
`cht149_flageul.ini` (the prepare ini is generated from it, so the two cannot
drift), `make_ynodes.py` + `ynodes_f149.dat`, `make_flageul_ic.py`,
`reseed_solid_exact.py`, `wall_lo_f.stl` / `wall_hi_f.stl`, `check_cht.py`.

**In [`asset/`](asset/)**: `flageul_data/` (their published raw profiles,
fetched by `fetch_flageul.py`), `flageul_fig5a_*.dat` (the digitised ideal
brackets, `digitize_flageul.py`), our statistics (`flageul_clean.h5`,
`Cstat_vel.h5`, and `flageul_stats.h5` — the *contaminated* window, kept so
the before/after above stays reproducible), the analysis scripts, `figures/`,
and [`CAMPAIGN_NOTES.md`](asset/CAMPAIGN_NOTES.md) — the full lab notebook of
all five campaigns that led here, including the ones this directory no longer
carries.

## Landmines hit on the way

* **A stale binary silently ran a DIFFERENT GRID.** `build_cpu` was rebuilt
  after `nodes_file` was added and `build_gpu` was not, so the solve ignored
  the key and fell back to a built-in distribution while `moby_prepare` used
  the file. **Nothing caught it** — the coefficient-file check compares
  lx/ly/lz/re and the block table, but *not* the node lines. It surfaced only
  as an absurd time step (3.6e-6 instead of 3.1e-4). The case file already
  stores x/y/z_nodes, so the solver could cross-check them at read; **that
  guard is still not implemented.**
* **mpirun's exit code is not a success test.** A develop leg finished
  cleanly, wrote its snapshot, then failed at MPI teardown with PMIX
  `NO-PERMISSIONS` errors; the driver's `|| exit 1` fired and the chained
  statistics leg never started — 12.5 h of correct work sat idle.
  `run_flageul.sh` now tests for `main loop ended` **in the log**.
* **Block datasets are stored `(nbz, nby, nbx)`**, so an x–z slice must be
  transposed before placement. Without that each tile is internally transposed
  while the tiles land correctly — a plausible-looking field that is wrong. It
  does not fail loudly, so it needs a positive check: *the mean |∂/∂x| across a
  block boundary must match the mean inside a block*. It was 5–8× larger (15–38×
  in z); transposed, 0.9–1.2. `make_flageul_ic.py` had the same error, which is
  why the campaign's IC had x and z swapped within each block — harmless here,
  since 80 000 develop steps plus 715 000 more followed it.
