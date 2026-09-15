# Conjugate pipe flow vs Neuhauser (NekRS, body-fitted) — session handout

STATUS: written 2026-09-15 on branch `scalar` as PLANNED; **every feature gap
is now closed and the campaign can run.**
UPDATE 2026-09-15: (a) the §1 **PIN DOWN is resolved** — and the reading this
handout was written with was WRONG; its consequences for dt and for the
transient legs are folded into §5. (b) **F1 (`source_dir`) and F5 (the annulus
STL) are implemented and gated** — §3 and §3a. (c) **F2 is DONE** — §3 and the
new §3b — but NOT as this handout proposed: the per-leaf mask field turned out
to be unnecessary, because φ already carries the shell. Read §3's F2 entry
before writing any case ini. **Nothing on the critical path is left; §4 and §5
are the next session's work.**
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

| isc | K | λ_sf | `solid_k` = κ_s/κ_f | `solid_rhocp` = ρc_s/ρc_f | α_s/α_f | note |
|---|---|---|---|---|---|---|
| **0** | 1 | 1 | 1 | 1 | 1 | **no contrast at all** — start here |
| 1 | 1 | 2 | **2** | 0.5 | 4 | conductivity contrast |
| 2 | 1 | 0.5 | **0.5** | 2 | 0.25 | conductivity contrast |
| 3 | 4 | 1 | 1 | **1/16** | **16** | pure CAPACITY contrast; **sets dt** |
| 4 | 0.25 | 1 | 1 | **16** | 1/16 | pure CAPACITY contrast; **sets settling time** |

plus `isc` 5 (`MBC`, K→0) and 6 (`IF`, isoflux) as non-conjugate controls,
which map to `ibm_wall = dirichlet` and `adiabatic`.

**Note what isc = 0 is and is not.** K = λ = 1 means κ_s = κ_f and α_s = α_f:
the materials are identical and the conjugate interface has *no contrast*. It
is the cheapest and the right first target — it validates geometry, forcing,
statistics and the whole pipeline — but it is the *weakest* test of the
conjugate scheme (our own gates show κ_s = 1 is a degeneracy where every
scheme is exact). **isc = 1 and 2 (λ_sf = 2 and 0.5) are the real test of the
conjugate face coefficient**, because they are the only two with κ_s ≠ κ_f and
the coefficient is a conductivity harmonic mean — it is IDENTITY at κ_s = κ_f
whatever the capacity does. isc 3 and 4 (the K sweep) carry no conductivity
contrast at all and instead test the C3 capacity `C_cell = f + (1−f)C_s`.

**PIN DOWN — RESOLVED 2026-09-15, and the earlier reading in this handout was
WRONG.** The mapping is

> **κ_s/κ_f = λ_sf** and **ρc_s/ρc_f = 1/(K² λ_sf)**, hence α_s/α_f = K² λ_sf².

so, directly, **`solid_k` = λ_sf** and **`solid_rhocp` = 1/(K² λ_sf)** (the
solver's κ and C are both ≡ 1 in the fluid, `scalar.f90` "kappa = 1 in the
fluid"). λ_sf is simply the CONDUCTIVITY ratio — λ as in the German λ for
thermal conductivity, "sf" = solid/fluid — and **K is the FLUID-to-SOLID
effusivity ratio**, K² = (κρc)_f/(κρc)_s, i.e. the INVERSE of the usual Tiselj
activity ratio. That inversion is why the dataset's own labels (§technical
remarks) read backwards at first sight and are in fact right: K → ∞ is the
inert solid = isoflux, K → 0 the infinitely effusive solid = isothermal
fluctuations.

The superseded reading (λ_sf = α_s/α_f, κ_s/κ_f = K√λ_sf) predicts κ_s/κ_f = 4
for isc 3; the true value is 1. **It would have put a conductivity contrast on
the two cases that have none and understated α_s by 16×.**

*Four independent confirmations, all on `cht_short`:*
1. **The file's own coefficients.** `diffusionCoeff`/`transportCoeff` are Nek's
   `vdiff`/`vtrans`, i.e. CONDUCTIVITY and ρc_p — not diffusivity. The fluid
   pair is exactly (ν/Pr, 1) for both Pr, which fixes the identification;
   `diffusionCoeffSolid/diffusionCoeff` is then κ_s/κ_f, and it equals λ_sf on
   all five rows, while `transportCoeffSolid` = 1/(K² λ_sf) on all five.
2. **Flux continuity in the data.** The φ- and time-averaged ∂T/∂r just inside
   and just outside r = 0.5 jumps by exactly λ_sf (measured 1.0002, 2.0005,
   0.5001, 1.0002, 1.0003 for isc 0-4 — the residual is the polar
   interpolation), never by K√λ_sf. isc 3 and 4 show NO jump.
3. **The authors' own script.** `mean_dtdr_comparison.py:39` builds the flux as
   `where(r < 0.5, 1, lambda_sf) * (−∂T/∂r)/(Pr/ν)`, i.e. it multiplies the
   solid gradient by λ_sf to make the flux continuous. Only consistent with
   λ_sf = κ_s/κ_f.
4. **The controls bracket the K sweep.** θ'_rms at r = 0.5⁻ (Pr = 0.71),
   normalised by isc 0: MBC 0.004, K = 0.25 → 0.446, λ sweep 0.95-1.03,
   K = 4 → 1.511, IF 2.078. Monotone in K between the two limits, exactly as
   K = e_f/e_s requires.

*This also CONFIRMS the control mapping below*: MBC carries θ'_rms ≈ 0 at the
interface (0.4 % of the conjugate case) ⇒ `ibm_wall = dirichlet` — it is the
Kasagi "mixed" idealisation, constant mean flux with isothermal fluctuations —
and IF carries the largest θ'_rms with zero fluctuating flux ⇒ `adiabatic`.

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

**ALL CLOSED (2026-09-15).** F1 and F5 were the two small ones; F2 was
expected to be the one real feature and turned out not to be, for the reason
its entry gives. F3 folded into F2 as predicted, and its substitute is exact
for every fluctuation statistic; F4 was never a gap.

**F1 — the Kasagi source was hardwired to `u_x`. DONE 2026-09-15.**
`[scalar.N] source_dir = x | y | z` picks the component the VELOCITY source
rides; default x, so the arithmetic of an ini that does not name one is the
original expression operand for operand. `stats_layout = plane` averages over
**z** and returns the (x,y) plane — exactly a pipe cross-section *if the axis
is z* — which is why the source had to move rather than the statistics: the
pipe axis is **z**, and `forcing_z` already existed. A `source_dir` on a
`source_type = uniform` scalar is a hard config error, not a silent no-op.
Gates in `validation/conjugate/run_gates_pipe.sh` (`source_dir`, `guard`),
CPU and GPU, plus the standard bit-exactness suites — see §3a.

**F2 — the solid must be an ANNULUS, and `solid_k` is a per-scalar CONSTANT.
DONE 2026-09-15 — but by a DIFFERENT and much smaller mechanism than this
handout proposed. F3 folds into it as predicted.**

*The problem, unchanged.* A Cartesian box cannot exclude the corner region
r > 0.6. Letting the solid run to the box faces makes its thickness vary from
0.1 (face mid-points) to 0.35 (corners), and §1 measured that the fluctuation
reaches the outer surface at 38 % and reflects off it — so an azimuthally
varying thickness gives azimuthally varying reflection and destroys the
`bccode = 0` premise. The region beyond r = 0.6 must be *insulating*.

*The mechanism that shipped: `[scalar.N] solid_thickness`, a LEVEL-SET BAND.*
The handout's plan was a per-leaf κ_s(x) mask written by `moby_prepare`. That
is not needed, and the reason is the same one that gave the C1 baseline its
geometry for free: **φ is already the signed distance to the body surface**,
so a shell of uniform thickness d wrapped around that surface is *exactly*
the level set −d < φ < 0, and everything deeper is φ ≤ −d. The shifted level
set ψ = φ + d is itself a distance function, so the shell/outer interface is
handled by the SAME distance-weighted harmonic mean on the SAME obliquity
lemma, with ψ in place of φ — the grazing guard included, since it tests
|φ_L − φ_R| = |ψ_L − ψ_R|. No mask, no `moby_prepare` stage, no case-file
dataset, no re-preparing any existing case file.

    [scalar.N] solid_thickness    = 0.1   ! shell depth; 0 (default) = one material
               solid_outer_k      = 0.0   ! default 0 = an exact insulator
               solid_outer_rhocp  = 1.0

`solid_source` fires in the SHELL only, by construction — which the handout
correctly identified as mandatory (a source in the insulating jacket has
nowhere to send its heat and would run away). `solid_outer_k = 0` is handled
by returning an exact zero rather than by the small-κ surrogate the handout
suggested: the harmonic mean would divide by it. An insulated outer band is
then wholly inert — no flux, no convection, no source, no contribution to the
explicit time-step limit — so it simply holds `solid_init`, and **the shell's
outer surface is an exactly adiabatic one at depth d.**

*Where it is weaker than a mask field, stated plainly.*
1. The band is a uniform-thickness OFFSET of the body surface. Exact for a
   pipe, and fine for any skin whose offset does not self-intersect; it is
   NOT a general region mask, and on a medial axis it follows φ rather than
   the geometry a user had in mind.
2. The outer surface is a STAIRCASE: band membership is decided at the cell
   centre, so the effective thickness carries the usual ±h/2 — at Δ⁺ ≈ 2 that
   is ±1.0 wall unit out of d⁺ = 36, i.e. ±2.8 % of d = 0.1. There is no point
   refining that with a volume fraction while the conduction geometry stays
   staircase — and a fraction-blended capacity would be *worse*, since it
   would dilute a straddling cell with a fictitious material. The fluid-side
   cut cells keep the C3 fluid-fraction capacity exactly as before: `vfrac`
   is 0 throughout the band region, so the two never interact.
3. `tangential_correction` (C2) with a band is a hard config error. C2's
   stencil, its Gershgorin rate and its indicator are all written against the
   single κ_s of the C1 baseline, and C2 ships disabled by measurement
   anyway.
4. The band boundary is a same-level arm like the body surface, so it is
   subject to the SAME 2:1 precondition (C1 gate 3c), which now tests both.
   That matters for the pass-2 `refine_body` optimisation of §4: its
   one-block 26-neighbour buffer keeps the SURFACE inside the finest level,
   and nothing yet guarantees the same for a band boundary 0.1 deeper.
   Expect to have to widen the refinement region, and note that the solver
   will say so rather than run.

*The one precondition, CHECKED at init rather than assumed.* φ is
1-Lipschitz, so an arm of length h changes it by at most h and a face can only
ever join ADJACENT bands — unless the band is thinner than the arm, in which
case the shell has holes and a face joins the fluid straight to the outer
material. `check_scalar_bands` counts exactly that, on the real φ field, and
hard-errors. It is stated as the thing itself, not as a proxy on the spacing,
so it holds on stretched and refined grids with no margin argument.

*And the premise itself was measured, not asserted*: on the prepared 64²×8
pipe case, the level set −0.1 < φ < 0 is the annulus 0.5 < r < 0.6 with
**0 flips** out of 64 000 ghost-inclusive cells, against the ANALYTIC polygon
distance (§3b).

**F3 — no nonzero Neumann on an immersed surface. NOT A SEPARATE FEATURE: it
folded into F2 exactly as predicted, and the substitute is exact for
everything this campaign measures.** `ibm_wall = neumann` is only an alias for `adiabatic` (zero flux,
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
solid, a one-cell layer approaches a surface flux to O(Δ/d) ≈ 5 %. **NOT
IMPLEMENTED** — the band confines `solid_source` to the shell, which is what
made the runaway impossible, but it does not concentrate it. It would be a
second band (`−d < φ < −d + t`) on the machinery that now exists, and it buys
only the solid MEAN profile, which is known in closed form for a uniform
source anyway. Do it only if that profile turns out to be wanted.

*Why the "faithful" route is the worse one.* An annulus presents the solver
with TWO cut surfaces, at r = 0.5 and r = 0.6, and the IBM marker is binary —
locally they are indistinguishable, so a surface-flux BC would need a
per-surface tag *on top of* cut-face area machinery that exists only in the
Python checkers (`face_area_fraction`), not in the solver. That is strictly
more work than F2, and it needs F2-like fields anyway.

*Implementation note, SUPERSEDED*: the draft said to give the jacket κ small
but nonzero (1e-6) because the harmonic mean `dm/(w/κ_L + (1−w)/κ_R)` divides
by κ = 0. The shipped code branches on it instead and returns an exact zero,
so `solid_outer_k = 0` is a true insulator rather than a nearly-one, and the
default `solid_outer_rhocp = 1` keeps the limiter finite as intended (it never
multiplies a non-zero rhs there anyway).

**F4 — radial statistics: NOT missing.** `stats_layout = plane` gives the
z-averaged (x,y) cross-section per scalar; radial binning is post-processing.
Velocity statistics are the gap — the channel stats module is y-profiles — so
take those from snapshots.

**F5 — no annulus STL generator. DONE 2026-09-15.**
`make_geometry_stl.py annulus --axis {x,y,z} --r-inner R --box-half H |
--r-outer R2 --facets M --a0/--a1` writes the closed shell. The pipe's body is
everything OUTSIDE r_inner — the fluid is the hole — so the outer surface is
pure padding; a SQUARE one (`--box-half`) clears the domain corner at a smaller
vertex magnitude than a cylindrical one, which the module's own precision
argument (padding costs precision quadratically) makes the right default for a
square cross-section. `--domain-half` asserts the property F2 needs and this
handout said to verify rather than assume: that inside the domain the pipe
wall, not the padding, is the nearest surface. 16384 facets on R = 0.5 gives a
chord error 9.19e-09, as predicted. Gates in §3a.

---

## 3a. What F1 and F5 were gated on (2026-09-15)

`validation/conjugate/run_gates_pipe.sh [source_dir|guard|annulus|ranks|all]`
— ALL PASS, and the `source_dir` group reproduces its numbers exactly on GPU.

| gate | measured |
|---|---|
| F1 closed form, all three branches | uniform (u,v,w) = (1,2,3) in a periodic box is an exact steady solution and a uniform scalar has no gradient, so θ = source·u_dir·t exactly: 0.04 / 0.08 / 0.12, max err **5.6e-17**. Three DIFFERENT answers, so a branch reading the wrong component lands on another one's number. |
| F1 same arithmetic in x and z | channel driven along the x–z diagonal, u(y) = w(y) (checked: max\|u−w\| **0.0**) ⇒ `source_dir` x and z **BIT-IDENTICAL**, every dataset 0.0 |
| F1 config guards | `source_dir` on a uniform source, and a bad direction name, both hard-error |
| F1 inert by default | `validation/scalar/run_bitexact{,_s3}.sh` vs `~/f1_ref_binaries` (cut from commit `0e0e225` with the edit stashed): 7-case and 9-case suites **max_abs 0, CPU AND GPU** |
| F5 classification | 64000 ghost-inclusive cells, solid ⇔ outside the 16384-gon, **0 flips** |
| F5 φ | max\|φ − φ_exact\| **5.6e-13** (z axis) / **5.9e-13** (x axis) against the analytic distance to the FACETED polygon |
| F5 nearest surface | pipe wall nearest everywhere inside the domain, margin **0.21** (box outer) / **0.29** (cylindrical outer) |
| F5 rank independence | prepared case file 1 == 4 ranks, dataset-identical |

TWO THINGS FOUND while building those gates, both recorded in the gate inis
because both cost real time:

1. **`dns%ibm_enabled` defaults to `.true.`** (`init.f90`), so an ini with NO
   `[ibm]` section silently runs the default analytic WAVY WALL. A triply
   periodic box that should preserve uniform flow exactly instead showed an
   8 % spread with a y-edge deficit modulated sinusoidally in x — which is
   precisely a wavy wall, and for an hour looked like a solver defect. With
   `[ibm] enabled = false` the spread is **0.0**. Any gate ini that does not
   want a body must say so.
2. **`tgv3d` is a pure GRADIENT**, u = ∇[−cos(kx)cos(ky)cos(kz)/k], so the
   projection annihilates it (measured: velocity 4e-8 instead of O(1)). It is
   the one x↔z-symmetric initial condition on offer and is therefore useless
   for a rotation test — hence the diagonal-channel bit-identity above
   instead. The module comment already says it is "NOT an NS solution"; this
   is what that means in practice.

Also corrected: `check_conjugate.py`'s docstring claimed case-file tiles are
`(k, j, i)`. They are `(i, j, k)` — `check_oblique.py` had already found this
and says so; the field datasets are the ones that are `(k, j, i)`. A pipe can
tell them apart (invariant in z, not in x or y): with the axes swapped the φ
error above reads 0.18 instead of 5.6e-13.

---

## 3b. What F2 was gated on (2026-09-15)

`validation/conjugate/run_gates_pipe.sh [band|insulate|bandguard|bandannulus|banddet]`
— ALL PASS, on binaries rebuilt from the final source (the README's own
PROVENANCE lesson). Full write-up and the "where the band is weaker than a
mask field" list in `validation/conjugate/README.md`.

| gate | measured |
|---|---|
| `band` — three-material slab | the C1 slab one layer deeper (jacket κ_o / shell κ_s = 2 / fluid 1, `T(0)=0`, `T(L)=1`). The steady profile is piecewise LINEAR in three layers and is an exact fixed point ONLY if both internal faces carry the true series resistance — the body-surface one on φ and the band one on ψ. Band boundary swept through a full cell × κ_o ∈ {0.01, 1, 100}: `max\|θ − exact\|` **0.0**, residual **0.0**, on all 12. |
| …and it is LIVE | mutation control: give the SOLVER a band depth wrong by a tenth of a cell and seed the same exact profile — it leaves at once, residual **1.9e-2**. |
| `insulate` — κ_o = 0 | the campaign's actual configuration. Shell + fluid reach the `y = L` Dirichlet 1; the jacket holds `solid_init = 0.5`, a value NEITHER domain face carries, so a leak either way shows. **4.4e-14**, residual 0.0. |
| `bandannulus` — the premise | on the prepared 64²×8 pipe, `−0.1 < φ < 0` vs the ANALYTIC polygon distance: 13520 shell / 20960 jacket cells, **0 flips** of 64000. The level set IS the annulus. |
| `bandguard` | 5/5 rejected: band thinner than the grid (caught at INIT on the real φ, not by a config rule); negative thickness; a `solid_*` key on a non-conjugate scalar; `tangential_correction` with a band; **a band boundary on a 2:1 block face** — C1 gate 3c's precondition, which tested the sign of φ only and so could not see the band. Paired with a CONTROL (same geometry, same refine box, band off) that must run, so the probe cannot pass by catching the body surface instead. |
| `banddet` | banded conjugate run: 1 == 4 ranks == **GPU**, `max_abs 0`. |
| inert by default | `run_bitexact{,_s3}.sh` vs `~/f1_ref_binaries`: 7-case and 9-case suites **max_abs 0, CPU AND GPU**; the whole C1 gate suite re-run unchanged; `scalar_test` gained 17 band assertions. |

**One time-step note**, for whoever re-derives §5's `dt`: a material interface
inside the solid excites the same extreme Gershgorin mode a cut cell does, so
a band cell takes the cut-cell `share = 2` — but ONLY when
`solid_outer_k > 0`. An insulated jacket can only remove face coefficients,
never amplify one, so charging it the cut-cell share would cost the whole run
a factor ~2.5 in `dt` for nothing. The campaign's jacket is insulated, so it
costs nothing.

---

## 4. Setting the case up

```
axis            z (periodic), L_z = 12.5
cross-section   [-0.65, 0.65]^2  (r = 0.6 solid, insulating jacket beyond)
geometry        make_geometry_stl.py annulus --axis z --r-inner 0.5
                    --box-half <clears the corner> --facets 16384   [F5, done]
nu              1.8868e-4                 [flow] re = 5300
forcing_z       1.8703e-2   = 2 u_tau^2/R   (gives Re_tau = 181, u_b ~ 1)
scalars         7   (5 conjugate + dirichlet + adiabatic controls)
Pr              0.71
required        keep_buried = true, remove_solid = false, no wall functions
```

and each of the five conjugate scalars, with (K, λ_sf) from §1's table:

```
[scalar.N]
pr                = 0.71
source_type       = velocity        ; Kasagi
source            = <beta>
source_dir        = z               ; F1 -- the pipe axis
ibm_wall          = conjugate
solid_k           = <lambda_sf>     ; = kappa_s/kappa_f       [see PIN DOWN]
solid_rhocp       = <1/(K^2 lambda_sf)>
solid_thickness   = 0.1             ; F2 -- the shell; beyond it, the jacket
solid_outer_k     = 0.0             ; (the default) an exact insulator
solid_source      = <-beta u_b/(2 d) ...>   ; the F3 substitute sink; SHELL only
```

`solid_source` must remove exactly the heat the fluid source puts in, or the
case has no steady state. Size it from the closed balance, then CHECK it
against the interface-heat diagnostic (`[scalar] heat_interval`, the C3
Nusselt branch) rather than trusting the arithmetic.

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
scalar-diffusion limit at Pr = 0.71 is ~10× looser *in the fluid*. **But the
pin-down changed this: isc 3 has α_s = 16 α_f** (κ_s = κ_f with C_s = 1/16),
and our limiter builds the rate from the actual face coefficient over the LOCAL
capacity — so those solid cells run 16× the fluid scalar-diffusion rate and the
scalar limit lands at ≈ 1.1e-3, i.e. **BELOW convection**. Expect dt ≈ 1.1e-3,
a factor ~1.6 under the estimate below, on top of the C3 cut-cell share. At
~0.2 s/step for 22 M cells × 7 scalars, ≈ 400 k steps ≈ 22 h buys ~700 D/u_b at
the old dt — so budget **≈ 35 h**, and re-derive rather than trust this.
Neuhauser's own `cht_short` window is t = 517 → 4942, i.e. **≈ 4400 D/u_b** —
six times that.

**And isc 4 sets the SETTLING time, symmetrically**: C_s = 16 gives
α_s = α_f/16 and a solid diffusion time d²/α_s ≈ **600 D/u_b** — comparable to
the entire planned statistics window, so its solid will still be equilibrating
while everything else is converged. So plan the **three-leg structure that
worked for the channel**: a cheap transient, a settle leg at the correct
properties, and only then a statistics leg, extended until the batch-means
error bars from Neuhauser's own `time` axis are matched.

*On accelerating the transient — the direction matters and the earlier draft
had it backwards.* To reach the solid's equilibrium faster you **REDUCE** C_s
(raising α_s), never inflate it; inflating the capacity slows the approach.
This is legitimate and already gated: C1's "capacity-independent steady state"
gate says the fixed point does not depend on C_s at all, so a transient leg may
run isc 4 at, say, C_s = 1 and restore C_s = 16 for the settle leg. isc 3 needs
no such help (its solid equilibrates in ~2 D/u_b) — it is the one paying in dt.

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

1. ~~**F2 is a genuine feature, not a config change**~~ — **DONE 2026-09-15**
   (§3, §3b), and much smaller than this handout expected, because φ already
   carries the shell. Nothing is left on the critical path. What survives of
   this risk is its second half: the shell's outer surface is a STAIRCASE at
   ±h/2, i.e. d = 0.1 ± 0.0028 azimuthally at Δ⁺ ≈ 2 (±1.0 wall unit out of
   d⁺ = 36, i.e. ±2.8 %). That is a real, if small, departure from Neuhauser's exact
   d⁺ = 36 — quote it, and note that the two-grid study of §6.6 moves it by
   the same factor as everything else.
2. ~~**The (K, λ_sf) → (κ_s, ρc_s) mapping**~~ — **RESOLVED 2026-09-15**, four
   ways, §1. `solid_k` = λ_sf, `solid_rhocp` = 1/(K² λ_sf). The reading this
   handout shipped with was wrong; if anything downstream still assumes
   κ_s/κ_f = K√λ_sf, fix it there too.
3. **First order at a curved interface.** Expect a visible interface error;
   the two-grid study is what makes it interpretable rather than embarrassing.
4. **Averaging window.** 700 D/u_b against their 4400. Quote batch-means error
   bars from their `time` axis and do not over-claim agreement inside them.
5. **The `t_final` landmine** — a `t_final`-terminated run takes one extra step
   whose dt is round-off and whose final snapshot is a bad restart. Fixed in
   the solver, but prefer `nsteps` for campaign legs anyway.
