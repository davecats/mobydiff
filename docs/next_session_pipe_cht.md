# Conjugate pipe flow vs Neuhauser (NekRS, body-fitted) — session handout

STATUS: PLANNED. Written 2026-09-15 on branch `scalar`. Nothing implemented.
Everything below marked *measured* was read out of the dataset itself in the
session that wrote this; everything marked *estimate* is not.

Data: `validation/conjugate/neuhauser_data/10.35097-26za20q32xsz43yk/data/dataset`
(23 GB, BagIt). Jonathan Neuhauser, *Conjugate Heat Transfer in Turbulent Pipe
Flows with Non-Uniform Heating Effects at Two Prandtl Numbers*, thermal DNS in
NekRS v26.0, DOI 10.35097/26za20q32xsz43yk.

**Read this first, because it sets expectations**: a pipe wall is a CURVED
conjugate interface, and `validation/conjugate/README.md` ("THE CONVERGENCE
TEST") measures the shipped scheme as **FIRST ORDER** there — 0.999 over an
eightfold refinement range — where every previous conjugate validation
(Flageul channel, slab gates) was a flat wall and therefore *exact*. This
campaign is the first test of the scheme on the geometry where it is weakest.
That is the point of running it, and it is also why **two grids are
mandatory** (§5): at first order a single grid cannot separate scheme error
from physics.

---

## 1. The case to run — and why it is the cheapest

### What the dataset contains (measured)

| | |
|---|---|
| geometry | pipe, D = 1, R = 0.5; conjugate solid shell 0.5 < r < 0.6 |
| bulk | Re_b = u_b D/ν = 5300, u_b = 1, **ν = 1/5300 = 1.8868e-4** |
| measured wall | **u_τ = 0.068379, Re_τ = 181.2, δ_v = 2.7593e-3** |
| solid thickness | d = 0.1 = 0.2 R, **d⁺ = 36.2** |
| Prandtl | 0.025 and 0.71 |
| meshes | `cht_long` (50 D), `cht_short` (12.5 D), `cht_fine` (12.5 D, finer) |
| statistics | z-averaged; `.interp.zarr` is polar (r, phi), r = 1e-4 … 0.59994 |

`z` is a scalar coordinate in every case (6.13 / 6.18 / 24.88 = the domain
mid-points), i.e. **the data are streamwise-averaged and the problem is
fully-developed and streamwise-periodic** — structurally the same reduction as
the Flageul channel, not a thermal entrance problem.

### The thermal problem (measured, not assumed)

The mean fluid temperature is **higher on the axis** than at the wall
(+2.51e2), and the mean solid profile is **pure logarithmic**,
`t − t_w = −1881.8 ln(r/R) + 0.18` — the 0.18 offset against 1882 says there
is no volumetric generation in the solid and *all* the mean heat crosses it
radially. So: heat is generated in the fluid by the streamwise-periodic
reduction `T = βz + θ ⟹ source = −βw` (Kasagi's form, exactly what
`source_type = velocity` implements) and removed through the outer solid
surface at a constant flux.

**The outer surface is load-bearing.** θ'_rms through the solid, as a fraction
of its interface value: 1.00, 0.77, 0.61, 0.52, 0.43, 0.393, **0.382** at
(r−R)/d = 0, 0.1, 0.2, 0.31, 0.51, 0.70, 1.00. It does not decay to zero, and
it *flattens* over the last 30 % — the signature of a constant-flux (i.e.
adiabatic-to-fluctuations) boundary reflecting the fluctuation. **The solid is
thermally THIN at d⁺ = 36, so its thickness must be 0.1 everywhere.** This one
measurement rules out the otherwise-obvious Cartesian shortcut of letting the
solid run out to the box faces (§3, F2).

### The selection

**`cht_short`** (12.5 D) — `cht_fine` is the same domain refined, `cht_long` is
4× the cells. And **`bccode = 0`** (azimuthally uniform): codes 2–5 are
sinusoidal azimuthal heating and need an azimuthally varying immersed flux,
which is out of scope.

**Pr = 0.71, not 0.025.** α = ν/Pr, so Pr = 0.025 makes the scalar diffusivity
**40×** the momentum one and the explicit scalar-diffusion limit dominates the
time step. Pr = 0.71 gives α = 2.6575e-4, close enough to ν that convection
sets dt.

**Run all five CHT variants as five SCALARS in ONE simulation.** `solid_k`,
`solid_rhocp` and `ibm_wall` are per-`[scalar.N]`, and all of
`bccode = 0` shares one velocity field — exactly the structure the Flageul
channel campaign used with six scalars. The five are `isc` 0–4:

| isc | K | λ_sf | α_s/α_f (measured) | note |
|---|---|---|---|---|
| **0** | 1 | 1 | 1 | **no contrast at all** — start here |
| 1 | 1 | 2 | 2 | |
| 2 | 1 | 0.5 | 0.5 | |
| 3 | 4 | 1 | 1 | the conductivity-contrast case |
| 4 | 0.25 | 1 | 1 | |

plus `isc` 5 (`MBC`, K→0) and 6 (`IF`, isoflux) as non-conjugate controls,
which map to `ibm_wall = dirichlet` and `adiabatic`.

**Note what isc = 0 is and is not.** K = λ = 1 means κ_s = κ_f and α_s = α_f:
the materials are identical and the conjugate interface has *no contrast*. It
is the cheapest and the right first target — it validates geometry, forcing,
statistics and the whole pipeline — but it is the *weakest* test of the
conjugate scheme (our own gates show κ_s = 1 is a degeneracy where every
scheme is exact). **isc = 3 (K = 4) is the real conjugate test.**

**PIN DOWN FIRST (unresolved here):** the mapping from (K, λ_sf) to
(κ_s/κ_f, ρc_s/ρc_f). `diffusionCoeffSolid/diffusionCoeff` is unambiguous and
gives α_s/α_f = λ_sf (measured, table above). With K the effusivity ratio
`K² = (κ_s/κ_f)(ρc_s/ρc_f)` this forces **κ_s/κ_f = K √λ_sf** and
**ρc_s/ρc_f = K/√λ_sf**, which is consistent with every row — but
`transportCoeffSolid` in the file equals `1/λ_sf`, not that, so the convention
is not self-evident. **Check it against the thesis before setting `solid_k`.**
Getting it wrong silently changes the physics, not the numerics.

---

## 2. Reading the data

`zarr` v3, and **no venv here had it** — `zarr 2.17` cannot open these files.
Done in the writing session: `~/ibmc/bin/pip install "zarr>=3" xarray dask`
(now zarr 3.3.0 / xarray 2026.7.0). Untar into scratch, then

```python
ds = xarray.open_zarr("joined_datasets.interp.zarr", consolidated=False)
sel = (ds.type=='CHT') & (ds.K==1.0) & (ds.lambda_sf==1.0) & (ds.bccode==0)
d = ds.isel(isc=int(np.flatnonzero(sel.values)[0])).sel(Pr=0.71).mean("time")
```

`dataset/datawrapper.py` is the authors' own accessor and carries their
conventions (`u_tau` from `d(w)/dr` at r = 0.5, `y+ = (0.5−r)/δ_v`,
`Re_τ = 0.5/δ_v`). **Use their definitions, not ours**, for every comparison.
`w` is the streamwise velocity; `time` is a batch axis for uncertainty
estimation, so quote error bars from it rather than presenting single numbers.

---

## 3. Feature gaps

**ONE real feature (F2), plus two small ones.** F3 turned out to fold into F2
rather than being separate, and its substitute is exact for every fluctuation
statistic; F4 is not a gap at all.

**F1 — the Kasagi source is hardwired to `u_x`.** `scalar.f90`:
`srcVal = source*0.5*(uw + ue)`, the x-faces. `stats_layout = plane` averages
over **z** and returns the (x,y) plane — which is exactly a pipe cross-section
*if the axis is z*. So the two constrain the orientation in opposite
directions, and one must be generalised. **Generalise the source direction**
(a `source_dir` key; a few lines, and a natural completion of the feature) and
put the pipe axis along **z**. `forcing_z` already exists.

**F2 — the solid must be an ANNULUS, and `solid_k` is a per-scalar CONSTANT.
This is the ONE feature the campaign needs; F3 folds into it.**
A Cartesian box cannot exclude the corner region r > 0.6. Letting the solid run
to the box faces makes its thickness vary from 0.1 (face mid-points) to 0.35
(corners), and §1 measured that the fluctuation reaches the outer surface at
38 % and reflects off it — so an azimuthally varying thickness gives
azimuthally varying reflection and destroys the `bccode = 0` premise. The
region beyond r = 0.6 therefore needs to be *insulating*, and there is
currently no way to say so: `solid_k` / `solid_rhocp` / `solid_source` are
scalars, not fields.
*Minimal fix*: let them be modulated by a per-leaf mask written by
`moby_prepare` beside `coef_p_blocks` (the `dwall_blocks` precedent), so the
face coefficient reads κ_s(x) instead of a constant. Contained, and it also
buys F3.
*Also check*: φ = ±dwall in an annulus is the distance to the NEAREST of the
two surfaces. For cut faces at r = 0.5 both adjacent cells are within Δ of the
inner surface, so `w` is unaffected — but verify rather than assume, and note
the outer surface at r = 0.6 becomes a *second* conjugate interface with
whatever lies beyond it.

**F3 — no nonzero Neumann on an immersed surface. NOT A SEPARATE FEATURE: it
folds into F2, and the substitute is exact for everything this campaign
measures.** `ibm_wall = neumann` is only an alias for `adiabatic` (zero flux,
six-face masking); the faithful outer BC is a constant nonzero flux at r = 0.6.
Replace it by a `solid_source` sink in the annulus with an adiabatic outer
surface, sized so the total heat matches.

*Why that is not an approximation where it matters.* The scalar is PASSIVE and
its equation is LINEAR, so split θ = ⟨θ⟩ + θ'. A steady, azimuthally uniform
source lives ENTIRELY in the mean: the solid's fluctuation equation is
`C_s ∂θ'/∂t = κ_s ∇²θ'`, in which `S` does not appear at all, and θ'_s is
driven only by the interface conditions and the outer boundary — adiabatic to
fluctuations in Neuhauser's constant-flux case and in ours alike. In the fluid,
θ' is driven by `u'·∇⟨θ⟩_f`, and ⟨θ⟩_f is fixed by the mean interface flux,
which is unchanged by construction. **So θ' everywhere — θ'_rms through the
solid, the turbulent heat fluxes, the whole conjugate signature — is EXACTLY
unaffected.** The only quantity that differs is ⟨θ⟩ *inside the solid*, and
since the source distribution is known, its profile is known in closed form:
compare our solid mean against its own analytic prediction, and compare
everything else against Neuhauser directly.

*And if the solid mean profile is wanted too*, concentrate the source in the
outermost solid cells rather than spreading it: with ~18 cells across the
solid, a one-cell layer approaches a surface flux to O(Δ/d) ≈ 5 %. That needs a
field-valued `solid_source` — **which is the same mechanism F2 already
requires**, so it costs nothing extra once F2 lands. It is in fact *mandatory*
once F2 lands: a per-scalar constant source would also fire inside the
insulating jacket, which cannot conduct the heat away and would run away.

*Why the "faithful" route is the worse one.* An annulus presents the solver
with TWO cut surfaces, at r = 0.5 and r = 0.6, and the IBM marker is binary —
locally they are indistinguishable, so a surface-flux BC would need a
per-surface tag *on top of* cut-face area machinery that exists only in the
Python checkers (`face_area_fraction`), not in the solver. That is strictly
more work than F2, and it needs F2-like fields anyway.

*Implementation note*: give the jacket κ small but NONZERO (e.g. 1e-6, not 0) —
the face coefficient is a harmonic mean `dm/(w/κ_L + (1−w)/κ_R)` and κ = 0
divides by zero — and C_jacket = 1 so the time-step limiter stays finite.

**F4 — radial statistics: NOT missing.** `stats_layout = plane` gives the
z-averaged (x,y) cross-section per scalar; radial binning is post-processing.
Velocity statistics are the gap — the channel stats module is y-profiles — so
take those from snapshots.

**F5 — no annulus STL generator.** Extend `validation/conjugate/
make_geometry_stl.py` (it already writes a faceted cylinder). Use enough facets
that the chord error is well below h: at 16384 facets on R = 0.6 it is ~1e-8.

---

## 4. Setting the case up

```
axis            z (periodic), L_z = 12.5
cross-section   [-0.65, 0.65]^2  (r = 0.6 solid, thin insulating jacket beyond)
nu              1.8868e-4                 [flow] re = 5300
forcing_z       1.8703e-2   = 2 u_tau^2/R   (gives Re_tau = 181, u_b ~ 1)
scalars         7   (5 conjugate + dirichlet + adiabatic controls)
Pr              0.71
source_type     velocity (Kasagi), source = beta, direction z   [F1]
ibm_wall        conjugate;  solid_k / solid_rhocp per scalar    [see PIN DOWN]
required        keep_buried = true, remove_solid = false, no wall functions
```

The driving force is *derived from the measured u_τ*, so Re_b is an outcome;
check it lands near 1 and say so rather than tuning. Grid (estimate):
Δ⁺ ≈ 2 ⟹ Δ = 5.5e-3 ⟹ 224² in cross-section; Δz⁺ ≈ 10 ⟹ 448 in z;
**≈ 22 M cells**. `refine_body` with a fine wall band could cut that ~3× and is
a pass-2 optimisation — it needs `keep_buried` and the 2:1 precondition (no cut
face on a coarse/fine block face), so do not take it on in pass 1.

---

## 5. Running it

**istmcorax (RTX 5090, 32 GB) is the target** — measured 0.111 s/step on a
13.8 M-cell conjugate channel with 6 scalars, ~6× the local 3060 per cell, and
usually idle. It needs its own `build_gpu_corax/` (`-gpu=cc120`) and has no
modulefile; `export PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.9/{compilers/bin,
comm_libs/12.9/hpcx/latest/ompi/bin}`. istmcetus (2× A6000) is the fallback —
**check `nvidia-smi` first, it is often someone else's**. Home is shared across
the hosts; **`/tmp` is not**, so driver scripts and logs live in the repo.
See the `remote-hosts` memory.

Cost (estimate): convection sets dt ≈ 1.7e-3 (CFL 0.4, u_max ~ 1.3); the
scalar-diffusion limit at Pr = 0.71 is ~10× looser. At ~0.2 s/step for 22 M
cells × 7 scalars, **≈ 400 k steps ≈ 22 h** buys ~700 D/u_b. Neuhauser's own
`cht_short` window is t = 517 → 4942, i.e. **≈ 4400 D/u_b** — six times that.
So plan the **three-leg structure that worked for the channel**: a cheap
transient (optionally with an inflated solid heat capacity to accelerate the
solid's approach to equilibrium), a settle leg at the correct properties, and
only then a statistics leg, extended until the batch-means error bars from
Neuhauser's own `time` axis are matched.

---

## 6. Comparison plan

Their definitions throughout (`datawrapper.py`). Bin our z-averaged (x,y) plane
into r, and compare against `.interp.zarr` azimuthally averaged.

Their radial grid is clustered **at the interfaces**, not in z: of 226 points,
48 sit in 0.45 < r < 0.5 and 41 in 0.5 < r < 0.55, with dr down to 1e-4 =
**0.04 wall units**. (That also settles the technical remarks' phrase "near the
start of the heated region" — it is the near-wall RADIAL region where the mean
gradient is steep, not a streamwise entrance; z is a scalar coordinate equal to
each domain's mid-point, so the data are streamwise-averaged.) Their radial
resolution is far finer than ours will be, so interpolate THEIR profile onto
OUR bins, never the reverse.

1. **Hydrodynamics first, before any thermal claim** — U⁺(y⁺), u'/v'/w'_rms,
   −⟨u'v'⟩, u_τ, Re_τ, u_b. If the velocity does not match, nothing thermal
   means anything. This also isolates the IBM pipe from the conjugate scheme.
2. **Mean temperature** ⟨θ⟩⁺(y⁺) in the fluid; θ_τ; Nusselt from the
   interface-heat diagnostic (conjugate branch, already gated three ways).
3. **The conjugate signature** — θ'_rms at the interface and through the solid,
   against the profile in §1 (1.00 → 0.382). This is the measurement the whole
   campaign exists for, and it is where the K and λ_sf variants separate.
4. **Turbulent heat fluxes** ⟨θ'u'⟩, ⟨θ'w'⟩; and ⟨θ'²⟩ budgets if the
   derivative products (`d(t)dx*d(t)dx` etc. are in the dataset) are worth it.
5. **Azimuthal uniformity** — our own check, not theirs: at `bccode = 0`
   everything must be φ-independent, so the azimuthal spread of ⟨θ⟩ and θ'_rms
   at the interface is a direct measure of the Cartesian-grid/IBM artefact on a
   curved wall. **This is the number that quantifies the first-order interface
   error**, and no body-fitted code can produce it.
6. **Two grids.** Repeat 2–5 at Δ⁺ ≈ 2 and ≈ 3 (or 4). At first order the
   interface error halves rather than quarters; showing that it moves *at the
   expected rate* converts a discrepancy from "unexplained" into "resolution",
   which is the whole reason for running two.

---

## 7. Risks, in the order they are likely to bite

1. **F2 is a genuine feature, not a config change**, and it is now the only
   one on the critical path. If it slips, the fallback is a thick-solid run
   that is *not* Neuhauser's d⁺ = 36 — report it as a different case, do not
   present it as a match.
2. **The (K, λ_sf) → (κ_s, ρc_s) mapping** (§1). Wrong here = wrong physics,
   silently.
3. **First order at a curved interface.** Expect a visible interface error;
   the two-grid study is what makes it interpretable rather than embarrassing.
4. **Averaging window.** 700 D/u_b against their 4400. Quote batch-means error
   bars from their `time` axis and do not over-claim agreement inside them.
5. **The `t_final` landmine** — a `t_final`-terminated run takes one extra step
   whose dt is round-off and whose final snapshot is a bad restart. Fixed in
   the solver, but prefer `nsteps` for campaign legs anyway.
