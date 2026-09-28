# Conjugate heat transfer in turbulent pipe flow

Validation of mobydiff's **immersed-boundary conjugate heat-transfer scheme**
on a **curved** solid–fluid interface, against a body-fitted spectral-element
DNS:

> J. Neuhauser, *Conjugate Heat Transfer in Turbulent Pipe Flows with
> Non-Uniform Heating Effects at Two Prandtl Numbers*, NekRS v26.0,
> DOI [10.35097/26za20q32xsz43yk](https://doi.org/10.35097/26za20q32xsz43yk).

Everything needed to **re-run the case** is in this directory. Everything
needed to **reproduce the comparison** — the reference profiles, the analysis
scripts, the figures and the full report — is in [`asset/`](asset/).

---

## Why this case

Every earlier validation of this scheme used a **flat, grid-aligned**
interface (Flageul's channel in [`../channel/`](../channel/), the C1–C3
manufactured slabs in `validation/conjugate/`), where the cut-face coefficient
is the textbook two-resistance harmonic mean and is *exact*. A pipe is the
first **curved** conjugate interface: the cut cells carry arbitrary level-set
weights, the wall is a staircase in the Cartesian mesh, and the solver's own
manufactured convergence test measures the scheme as **first order** in that
geometry. This campaign asks what that costs in a turbulence result, and
whether the error behaves as discretisation error should.

The answer, in one line: **on the production grid every conjugate statistic
matches the body-fitted DNS to 1–2 %**, the residual sits where the theory
says it must (the cut cells and the under-resolved core), and it moves with
the grid.

## The case

| | |
|---|---|
| geometry | pipe, D = 1, R = 0.5, axis z, streamwise periodic, L_z = 12.5 D |
| solid | conjugate shell 0.5 < r < 0.6 (d⁺ = 36), exactly insulated beyond |
| domain | Cartesian [0, 1.3]² × [0, 12.5], pipe centred at (0.65, 0.65) |
| driving | fixed pressure gradient at Neuhauser's **measured** u_τ: `forcing_z = 2u_τ²/R = 1.8702750564e-2` |
| flow | Re_τ = 181.2, δ_v = 2.7593e-3, ν = 1/5300; **Re_b is an outcome, not an input** |
| thermal | Pr = 0.71, Kasagi streamwise reduction ⇒ volumetric source ∝ u_z |
| grid | 256 × 256 × 896 = 58.7 M cells, uniform; Δ⁺ = 1.84, Δz⁺ = 5.06 |
| blocks | `nb = 32`, 1792 leaves, **single level** (see "no 2:1 refinement" below) |
| cut faces | 702 464 conjugate faces |
| cost | 0.52–0.56 s/step on one RTX 5090; the statistics leg is 168 000 steps ≈ 24–26 h |

**Seven scalars ride one velocity field**, so every comparison between them is
at the same turbulence realisation. The dataset's (K, λ_sf) map to the solver
as κ_s/κ_f = λ_sf and ρc_s/ρc_f = 1/(K²λ_sf), K being the **fluid-to-solid**
effusivity ratio (settled four independent ways — see the report):

| | c0 | c1 | c2 | c3 | c4 | mbc | isof |
|---|---|---|---|---|---|---|---|
| K | 1 | 1 | 1 | 4 | 0.25 | → 0 | → ∞ |
| λ_sf | 1 | 2 | 0.5 | 1 | 1 | — | — |
| `solid_k` | 1 | 2 | 0.5 | 1 | 1 | (Dirichlet body) | 0.01 |
| `solid_rhocp` | 1 | 0.5 | 2 | 1/16 | 16 | — | 0.01 |
| α_s/α_f | 1 | 4 | 0.25 | 16 | 1/16 | — | 1 |

`c1` and `c2` are the only cases with a **conductivity** contrast and are
therefore the real test of the conjugate face coefficient, which is a
conductivity harmonic mean and is the identity when κ_s = κ_f. `c3` and `c4`
test the fluid-fraction cut-cell capacity instead. `mbc` (isothermal
fluctuations) and `isof` (isoflux) bracket the sweep.

### Two solver features this case leans on

**The shell is a level set, not a mask.** `solid_thickness = 0.1` carves the
conjugate band out as −0.1 < φ < 0, where φ is the signed distance to the pipe
wall — so the shell *is* 0.5 < r < 0.6 exactly, at any orientation, with no
new dataset. Everything deeper is the inert jacket at `solid_outer_k = 0`, an
exact insulator, which makes the shell's outer surface exactly adiabatic at
depth 0.1. That matters: the reference solid is thermally **thin**, the
fluctuation reaches its outer surface at 38 % and reflects, so a jacket run
out to the box faces would give an azimuthally varying thickness and destroy
the premise.

**A nonzero Neumann condition at an immersed surface is not expressible**, so
the heat is removed by a volumetric **sink in the shell** instead. The
substitute is *exact for every fluctuation statistic* — a steady, azimuthally
uniform source lives entirely in the mean, and the solid's fluctuation
equation does not contain it — and changes only ⟨θ⟩ inside the solid, whose
profile is then known in closed form.

### No 2:1 refinement, deliberately

The conjugate cut-face coefficient is a **same-level arm**. Both the pipe wall
and the band boundary 0.1 deeper would have to lie inside the finest level,
which is nearly the whole annulus — so refinement would buy almost nothing
here and the solver checks the precondition at init rather than assuming it.
`asset/figures/fig8_blocks.png` shows the single-level lattice.

---

## Results against Neuhauser

Statistics window t = 55…245 D/u_b, 168 000 steps. Every number below is
printed by `asset/pipe_numbers.py` from the shipped accumulators; none is
typed by hand.

**Provenance of the shipped accumulators (2026-09-28).** They come from a
RE-RUN of the statistics leg with the solver at numerics-review step 2
(`docs/numerics_review_2026-09-26.md`, finding F2): a conjugate scalar now
runs the ADVECTIVE convective form, where the original 2026-09-20 campaign
ran the divergence form whose cut-face convective mask is inconsistent with
continuity on a curved wall. Same settled state (`p_pr_settle2_30000.h5`),
same ini (`niter = 6` Chebyshev, so the binary is the ONLY change), same
window; the run diverges from the original realisation chaotically, so the
two are different samples of the same statistics, not a bit-for-bit pair.
What moved, old → new (the old numbers are the 2026-09-20 campaign's, the new
are the tables below): interface heat balance 9.760 → 9.761 (c0; every
scalar within 0.04 %), q_w 0.24854 → 0.24856, ⟨θ⟩⁺ at the axis +0.1 → +0.7 %,
θ′⁺ at the axis −6.2 → −6.0 %, θ′⁺ for y⁺ < 50 +0.7 → −0.1 %, the interface
signature at y⁺ 7.25 by ≤ 0.004 (c4 0.903 → 0.907, mbc 0.862 → 0.866, isof
1.209 → 1.207, all toward the reference), and — the one systematic move —
the fluctuation through the solid, normalised at the interface: 0.798 /
0.638 / 0.528 / 0.420 / 0.385 / 0.374 → 0.802 / 0.633 / 0.521 / 0.404 /
0.365 / 0.353 against Neuhauser's 0.757 / 0.591 / 0.490 / 0.388 / 0.353 /
0.343, i.e. the deep-shell excess fell from +9 % to +3 %, and renormalised at
depth 0.1 d the decay now reads −2.0 … +0.4 % where it read +0.1 … +7.7 %.
That is the direction F2 predicts (the divergence-form drift at the fluid
cut faces was a spurious low-frequency source at the interface, and the
solid integrates exactly that), but one realisation cannot separate it from
sampling: the VELOCITY statistics, which F2 does not touch at all, moved by
comparable amounts between the two realisations (u_z′ in the core +3.1 →
−3.4 %, ⟨u_r′u_z′⟩ −1.2 → +1.0 %, bulk velocity 1.0075 → 1.0092), which
bounds what any single thermal difference can claim. The interface-heat and
budget tables are unchanged to their printed digits. The `asset/` figures
were regenerated from the new accumulators; `asset/pipe_report.html` keeps
the 2026-09-20 campaign's figures and tables (they are a four-grid study of
which this is one column) and carries the re-run as an addendum.

### Hydrodynamics — validated first, with no scalars present

| quantity | deviation | where |
|---|---|---|
| ⟨u_z⟩ | **+1.1 %** | core, r < 0.1 |
| u_z′ rms | −3.4 % | core |
| u_r′ rms | −0.6 % | 0.35 < r < 0.47 |
| u_θ′ rms | −0.1 % | " |
| ⟨u_r′u_z′⟩ | +1.0 % | " |
| bulk velocity | 1.0092 vs 1.0042 (**+0.5 %**) | Re_b 5349 vs 5322 |

**An exact law, independent of the reference data.** Integrating the axial
momentum balance of a steady, fully developed pipe once gives

    −ν d⟨u_z⟩/dr + ⟨u_r′u_z′⟩ = u_τ² r/R

and the measured deviation is **max 3.1 %, mean 0.49 %** of u_τ² — confirming
both statistical convergence and that the immersed wall delivers the stress it
is driven with. *Sign note:* with r measured from the axis the mean shear is
negative, so the down-gradient Reynolds stress is **positive** and adds to the
viscous term — the opposite of the channel's ⟨u′v′⟩ < 0.

### Thermal field, case c0

q_w = 0.24856 from the interface-heat balance |Q|/(2πRL); θ_τ = q_w/u_τ =
3.6351.

| quantity | deviation |
|---|---|
| ⟨θ⟩⁺ at the axis | **+0.7 %** |
| θ′⁺ for y⁺ < 50 | **−0.1 %** |
| θ′⁺ at the axis | −6.0 % |

The core θ′ deficit is the one residual, and §8 of the report traces it: it is
**turbulent transport of variance into the core**, carried by z-elongated
structures, and it converges in the **axial** spacing (−16.6 % at Δz⁺ 14.2,
−14.5 % at 10.1, −6.0 % at 5.1) while being nearly independent of the radial
spacing. It is not the source/sink treatment — see below.

### The conjugate signature — the actual test of the scheme

θ′ at the interface relative to case c0, quoted at r = 0.48 (y⁺ 7.25), just
outside the cut cells:

| | c0 | c1 | c2 | c3 | c4 | mbc | isof |
|---|---|---|---|---|---|---|---|
| present | 1.000 | 0.997 | 1.001 | 1.104 | 0.907 | 0.866 | 1.207 |
| Neuhauser | 1.000 | 0.995 | 1.002 | 1.092 | 0.921 | 0.880 | 1.190 |
| deviation | — | +0.2 % | −0.1 % | +1.1 % | −1.5 % | −1.6 % | +1.4 % |

The whole effusivity bracket — from isothermal to isoflux, a factor 1.4 in
θ′ — is reproduced to ≤ 1.6 %. **The conductivity pair c1/c2 agrees to
0.2 %**, which is the face coefficient itself passing.

**The quoting radius matters, and the caveat is measured rather than
asserted.** The outermost fluid bin is a **cut cell**, where what is stored is
a penalization blend over a cell straddling the wall; comparing that against a
body-fitted DNS is not like-for-like, and it is where the profiles are
steepest. At r = 0.4977 (y⁺ 0.83) the same table reads c3 1.492 vs 1.444 and
mbc 0.118 vs 0.202 — that is what the caveat is worth, and `pipe_numbers.py`
prints the cut-cell row every time so the reader can see it.

### Through the solid

θ′ normalised by its interface value, against depth (r − R)/d:

| (r−R)/d | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|---|
| present | 0.802 | 0.633 | 0.521 | 0.404 | 0.365 | 0.353 |
| Neuhauser | 0.757 | 0.591 | 0.490 | 0.388 | 0.353 | 0.343 |

The profile does not decay to zero — it **flattens over the last 30 %**, the
signature of the adiabatic outer surface reflecting the fluctuation, and the
solver reproduces that shape. The +6 % (shallow) to +3 % (deep) offset is the
cut-cell denominator: normalising one cell deeper (at depth 0.1, as Fig. 4a
does) leaves −2.0 … +0.4 %. (The 2026-09-20 campaign read a uniform +9 %
here; see the provenance note above.)

### The ⟨θ′²⟩ budget closes on both sides

Wall units, each side by its own θ_τ (`asset/variance_budget.py`):

| y⁺ | | P | −ε | T_turb | T_mol | source | sum |
|---|---|---|---|---|---|---|---|
| 110.5 | present | 0.026 | −0.038 | 0.012 | 0.001 | 0.002 | 0.002 |
| 110.5 | Neuhauser | 0.025 | −0.038 | 0.011 | 0.001 | 0.002 | 0.001 |
| 38.1 | present | 0.148 | −0.154 | 0.007 | 0.002 | 0.006 | 0.009 |
| 38.1 | Neuhauser | 0.154 | −0.168 | 0.000 | 0.003 | 0.007 | −0.005 |
| 19.9 | present | 0.334 | −0.221 | −0.070 | −0.027 | 0.010 | 0.025 |
| 19.9 | Neuhauser | 0.338 | −0.237 | −0.075 | −0.038 | 0.010 | −0.001 |

(Three snapshots; the residual at y⁺ 19.9 is the sampling floor of that
count.)

The `source` column is the **fluctuating** part of the Kasagi heating: it is
proportional to u_z, so it feeds the variance directly at 2S⟨u_z′θ′⟩. It
agrees to 1–3 % at every radius and is ~4 % of the budget — which is why the
core deficit cannot be blamed on the heating treatment. Two further grounds:
the sink difference lives entirely in the solid's **mean** field, and the
treatment was identical on all four grids while the deficit moved from
−16.6 % to −6.0 % with Δz⁺ alone. A modelling difference cannot respond to
grid spacing.

### Interface heat balance

Per scalar, time mean over the statistics window (should be equal across
scalars — each carries the same wall flux by construction):

| c0 | c1 | c2 | c3 | c4 | mbc | isof |
|---|---|---|---|---|---|---|
| 9.761 | 9.768 | 9.741 | 9.772 | 9.723 | 9.822 | 9.756 |

Spread 1.0 %, which is the sampling scatter of the slowest scalars.

### Two methodological points that changed the answers

**(a) The shell sink must balance the DISCRETE generation.** Sizing it from
the continuum `u_b = 1` leaves a 2.79 % imbalance, because the discrete
`∫u_z dV` gives u_b = 0.9755 on the coarse grid. The field then cools
uniformly. The mean profile is immune — a uniform ramp cancels when referenced
to the interface temperature — **but the time-averaged variance is not**: the
drift adds ≈ (λT)²/12, spatially uniform, hence worst where the true
fluctuation is smallest. Measured: θ′ inflated by **+32 % at the axis**, and,
because λ ∝ 1/C_total, a *different* bias per scalar that inverted the
conductivity sweep. `measure_generation.py` measures each grid's own factor
before its statistics leg (production: 1.029210).

**(b) θ_τ must come from the heat balance, not from the near-wall slope.** A
slope differenced across the cut cells reads 6.6 % low on the coarse grid and
0.0 % on the fine one, so normalising by it *flatters* the coarse grid by
cancelling a too-small scale against a too-low profile. Every wall unit above
uses q_w from the interface-heat diagnostic.

### What is not claimed

* The Cartesian mesh leaves a **4-fold azimuthal signature** at the wall at
  ~1 % of ΔT (Fig. 6b,c). Its **convergence rate was withdrawn**: the
  measurement band 0.49 ≤ r < 0.4999 holds a different set of cell layers on
  each grid, so the amplitudes are not a convergence sequence.
* The **isoflux** case shows high-frequency noise inside the solid
  (Fig. 7, bottom row). It is the near-insulating shell (κ_s = 0.01) resolving
  the staircase of the immersed wall as a real temperature pattern; it does
  not contaminate the fluid-side statistics, which agree to 1.6 %.
* Statistical convergence, not the scheme, sets the floor on several numbers
  above: the window is 190 D/u_b.

---

## Reproducing the run

Requires: an MPI build of `moby_solve` (GPU recommended — 22.7 GB of device
memory) and `moby_prepare` (**CPU build**: the GPU build computes the
coefficients on the device and differs by libm ulps), plus a Python with
numpy, h5py, matplotlib and trimesh.

```bash
module load toolkits/nvhpc/25.9
export PY=$HOME/ibmc/bin/python
export BIN=$PWD/../../../build_gpu/moby_solve
export PREP=$PWD/../../../build_cpu/moby_prepare
export GRID=prod TAG=p_pr

./run_pipe.sh geom        # the annulus STL (grid independent)
./run_pipe.sh prepare     # -> pipe_prod.h5, the case file
./run_pipe.sh check       # gate: the level set -0.1 < phi < 0 IS 0.5 < r < 0.6

# a one-step cold start, only to get a snapshot with the right layout and
# metadata for make_pipe_ic.py to fill
DEV_STEPS=1 DEV_FIELD=1 ./run_pipe.sh dev
$PY ./make_pipe_ic.py p_pr_dev_1.h5 IC_prod.h5   # see the warning below

IC=IC_prod.h5 ./run_pipe.sh dev                  # leg A: hydrodynamics
./run_pipe.sh vel                                # leg A statistics window
./run_pipe.sh therm                              # leg B: 7 scalars, c4 accelerated
./run_pipe.sh settle                             # leg C: c4 capacity restored
./run_pipe.sh settle2                            # leg C': until the interface heat is flat

# size each grid's sink from ITS OWN discrete generation (see (a) above)
SET=$(ls -1 p_pr_settle2_[0-9]*.h5 | sort -t_ -k4 -n | tail -1)
SCALE=$($PY ./measure_generation.py "$SET" pipe_prod.h5 --current 2.21110542)
FROM=$SET SRCSCALE=$SCALE STAT_STEPS=168000 STAT_FIELD=12000 \
    ./run_pipe.sh stats                          # leg D: the measurement
```

**Do not cold-start the pipe.** A plug profile plus white noise
**relaminarises** — a plug carries no shear to feed the turbulence.
`make_pipe_ic.py` writes a 1/7 power-law mean plus a divergence-free
perturbation `u′ = ∇ × (g(r)A)` with `g = (1 − (r/R)²)²`, which does sustain.

**The actual campaign took a shortcut** and the README states it rather than
hiding it: the production run was restarted from the converged **fine** grid
(256² × 448, same radial spacing) through `make_pipe_refine.py` +
`merge_scalars.py`, so only z was refined and the slow solid-side state — which
relaxes on d²/α_s and is the reason the z-fine run's shell profile was the
worst of the three — mapped across without being interpolated in r at all.
Running the chain above from scratch reproduces the statistics to the quoted
sampling scatter, not bit for bit.

Then turn the snapshots into radial profiles and refresh the caches:

```bash
$PY ./pipe_stats.py p_pr_stat_*.h5 --case pipe_prod.h5 --out p_pr_snaps.npz
$PY ./pipe_stats.py --stats p_pr_stats.h5 --case pipe_prod.h5 --out p_pr_statsD.npz
cp p_pr_statsD.npz asset/pipe_prod_statsD.npz
cp p_pr_snaps.npz  asset/pipe_prod_snaps.npz
cp p_pr_stat.heat.txt asset/pipe_prod_heat.txt
$PY asset/make_caches.py .            # the figure caches
```

## Reproducing the comparison

No download, no large files — everything below runs from `asset/`:

```bash
cd asset
./pipe_numbers.py                     # every table above
./plot_pipe.py                        # figures 1-7
./plot_pipe.py --blocks               # figure 8
./compare_neuhauser.py velocity pipe_prod_snaps.npz
./compare_neuhauser.py thermal  pipe_prod_statsD.npz --scalar c0 --wall-flux 0.24856
```

The reference profiles come from `asset/neuhauser_profiles.npz`, a 250 kB
reduction of the 11.8 GB published archive. **It is checkable, not merely
convenient**: `extract_neuhauser.py` rebuilds it from the archive, and figures
1–5 come out **byte-identical** either way.

## Files

**Here** (to run the case):

| | |
|---|---|
| `run_pipe.sh` | the campaign driver: one leg per subcommand, documented inline |
| `pipe_prep.ini` | preprocessing; `[scalar] count = 1` is what makes `moby_prepare` write `coef_p_blocks` |
| `pipe_hydro.ini` | leg A, no scalars |
| `pipe_thermal.ini` | legs B–D, the seven scalars |
| `../../../tools/make_geometry_stl.py` | the annulus STL, shared with the C2 gate suite (`--axial-segments`: see the landmine below) |
| `make_pipe_ic.py` | the turbulent initial condition |
| `make_pipe_refine.py`, `merge_scalars.py` | the grid-to-grid restart chain |
| `../../../tools/check_annulus.py` | the geometry gate on the prepared case file, shared with the F5 gate |
| `measure_generation.py` | the discrete-generation sink factor |
| `pipe_stats.py` | snapshots / solver statistics → radial profiles |
| `combine_stats.py` | merge consecutive statistics legs |

**In [`asset/`](asset/)** (to reproduce the comparison): the reference cache,
the analysis scripts, `figures/`, and [`pipe_report.html`](asset/pipe_report.html)
— the full publication-grade report, which carries the derivations, the
figure discussion and the findings summarised above.

## Landmines hit on the way

* **Sliver triangles.** The annulus skins spanning the whole pipe length made
  the solver's point-triangle distance lose precision: φ error 2.09e-6.
  `--axial-segments 8` fixes it (verified 1.88e-11). Keep it.
* **Radial bins narrower than a cell** empty near the axis and alias onto
  m = 4, 8, 12 — exactly where a Cartesian artifact would land.
  `pipe_stats.azimuthal_mean` therefore uses an **adaptive** ring width
  w(r) = max(w_min, NΔ²/2πr).
* **Rings straddling r = 0.6** mix shell cells with inert jacket cells held at
  a constant value: residual 667 against 4–6 inside the shell. Excluded
  everywhere it is measured or drawn.
* **rms definition.** Accumulating the variance over a radial bin adds the
  azimuthal variance of the steady staircase pattern (inflation 22–36 % for
  c0, 72–113 % for isoflux). `pipe_stats.py` accumulates the per-cell
  **temporal** variance instead and keeps the bin-moment version under
  `*_rms_binmoment` for comparison.
* **cos φ under transposition.** Built on an (x, y) mesh and broadcast against
  a (z, y, x) field, cos φ silently becomes sin φ and the radial flux comes
  out 5× too small. The field layout is (k, j, i); build the mesh as (y, x).
