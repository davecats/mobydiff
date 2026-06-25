# 2:1 interface spurious band — handout

Branch `claude/jacobi-interface`. Latest relevant commit `fdfd477`.

## The problem

A turbulent Re_tau=180 channel with a 2:1 wall-band refinement
(`validation/channel_interface/refined_y110`) **blows up at the interface**
(~step 180; div_l2 and global mass explode). Root symptom found earlier:

> The **clean** field (long before any blow-up) already carries a **spurious
> band on the COARSE-INTERIOR cell of every 2:1 interface, in all variables**
> (u,v,w,p). The fluctuation rms spikes there (u' rms ~1.14 -> **1.55** -> ~1.23
> across y, with the same bump at both wall bands), and the channel advection
> grows it into the blow-up. Cross-section: `/tmp/chanval/xsection_reflux.png`,
> zoom `/tmp/chanval/clean_interface_zoom.png`.

Davide's requirement: **a 2:1 interface method that does NOT limit the stability
of the underlying scheme and is not too sensitive to the interface position.**

## What is RULED OUT (do not re-investigate)

- **Not the projection accelerator.** Jacobi and Chebyshev blow up identically at
  the same step (300-step comparison `/tmp/chanval/divergence_comparison_refined.png`).
  In the stable regime Chebyshev holds ~3.4x lower divergence (so it is the better
  default, already set), but it does not change stability.
- **Not the reflux.** reflux ON vs OFF: same location/mode; reflux only HALVES the
  magnitude and DELAYS onset (163->181), it does not cure it. (Off is worse.)
- **Not (mainly) the velocity prolong halo.** A real artifact WAS found and FIXED
  there (below), but it is NOT this band (see "Key distinction").
- **Not the niter / sor.** sor=0.8 niter=6 Jacobi is fine on the UNIFORM channel
  (div decays, mass round-off, 6000 steps). The earlier "uniform blow-up" was the
  restart overriding sor to 1.5 -> fixed (commit cc31dc2). NB nb is the BLOCK SIZE
  (cells/block), not blocks-per-dim.

## What WAS fixed (commit fdfd477) — and why it is NOT the band

Gate 1 (linear-field exactness = `MOBY_HALO_AUDIT`) localized a genuine defect:
the coarse->fine **velocity prolong** injected a **tangential checkerboard** for a
smooth field. The face blend `(2C+F)/3` is normal-2nd-order but injected the coarse
value tangentially (pair-constant) while the fine-interior term is per-cell -> a
linear field gained `+-slope*dx/3` alternating. Fixed by a **two-pass exchange**
(same-level copies fill coarse halos, THEN cross-level prolong/restrict) plus a
**2-point coord-weighted tangential interpolation** of the coarse value
(`copy_local_entries`, node-line weights, handles staggering). Validated: audit
checkerboard GONE, bit-exact with no interface, channel stable.

**KEY DISTINCTION (do not re-conflate):** the velocity prolong fills fine-block
HALOS; the channel band sits on the COARSE-block INTERIOR cell (and the reassembled
global field shows interiors, not halos). So this fix could not, and did not, move
the band. The phi scalar prolong has NO blend (pure injection -> a smooth step, not
a checkerboard) so it is not it either. **The band is a separate, coarse-side
defect** -- in how the COARSE interface cell's value is set by the predictor
(differencing the restrict-filled halo across the 2:1 face) and/or the projection
(2:1 divergence/pressure consistency).

## The validation pipeline (gate suite)

Spec: `validation/interface_suite/README.md`. The shared metric the old gates
MISSED is **interface roughness**: RMS of the per-cell discrete Laplacian
(checkerboard/high-k proxy), split BAND vs INTERIOR. A good interface keeps
`roughness(band) ~ roughness(interior)`. Old order gates used smooth fields (mode
tiny); old conservation gates were global (mode is locally div-free) -> both blind.

Ordered gates (simplest/most-diagnostic first; position-sensitivity dropped):

| # | gate | how | status |
|---|------|-----|--------|
| 0 | roughness metric | `interface_diagnostics.py` (error-based, Beltrami); reference-free FIELD variant TO BUILD | partial |
| 1 | linear-field exactness | `MOBY_HALO_AUDIT` (manufactured linear, checks every exchanged halo vs the transfer design; per-var/op/un-written breakdown) | RUN: velocity prolong checkerboard fixed; residual = periodic-wrap artifact of the non-periodic test field |
| 2 | bug localization | toggle pieces: `MOBY_NORECON`, reflux on/off, projonly vs predonly, copy vs prolong vs restrict | partial (prolong localized; coarse-cell side OPEN) |
| 3 | projection-only | Beltrami base + `MOBY_MANUF` (globally conserving, not div-free), `MOBY_PROJONLY MOBY_RESLOG MOBY_DIVDUMP`, niter sweep {6,12,50,200,1000} jac+cheb, band patch no edges, 2-3 res. Gates: div->0 with niter; projected==base; roughness(band)~int at every niter | TO RUN |
| 4 | predictor-only, per term | Beltrami, `MOBY_PREDONLY MOBY_TERMDUMP MOBY_RHSDUMP`. Gates: each term O(h^2) bands (O(h) edges) via `rhsband.py`/`rhsterms.py`; roughness(band)~int for EVERY term. bands first | TO RUN |
| 5 | interface-mode growth | a BAND interface on a CONTROLLED advected base (uniform/laminar at U), seed noise, growth-vs-U -> the stability limit, isolated from turbulence (generalizes `interface_decay`) | TO BUILD |
| 6 | channel 200 steps | `refined_y110`, {niter 6,12} x {jac,cheb}, `MOBY_STEPDIV` + field dumps. Gates: stability (bounded div/mass, onset step); reference-free roughness over time | RUN: both blow up ~181, Chebyshev ~3.4x lower until then |
| 7 | regression | bit-exact no-interface + CPU==GPU | velocity-prolong fix passes both |

### Tools / env flags
- `tools/interface_diagnostics.py FIELD.h5 [FINE.h5] [--band B]` -- band/interior
  error L2/Linf, ORDER (2 files), momentum drift, **roughness (high-k) proxy**.
- `tools/rhsband.py`, `rhsterms.py` -- per-term order, band vs interior.
- `tools/divsum.py`, `compare_fields.py` (`--export-global` to reassemble a
  block-table field to a global grid), `check_beltrami.py`.
- env (off by default): `MOBY_PROJONLY MOBY_PREDONLY MOBY_MANUF=<amp> MOBY_DIVDUMP`
  `MOBY_RHSDUMP MOBY_TERMDUMP MOBY_RESLOG MOBY_STEPDIV MOBY_NORECON MOBY_HALO_AUDIT`
  `MOBY_PHASETIME`.
- `MOBY_HALO_AUDIT` now prints: bad per var, per op (copy/prolong/restrict), max
  FINITE error per op, UN-WRITTEN (sentinel) count + face/edge/corner split, and
  the first 40 bad / 16 un-written cells. Run with the refined IC (block layout).

### Reproducing the band
```
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
./compile.sh gpu
# IC (256x128x256 finest; band layout) once:
python3 tools/make_channel_restart.py --mode refined --band-cells 24 \
   --source tutorials/channel_kmm180/channel_kmm180_restart.h5 --out IC.h5
# run ~150 steps, dump field, export to global, measure interface rms by y-row
# (the y-idx-78 spike in u/v/w/p is the band). See the scratch python used in the
# session (per-y rms over x,z of the x-mean-removed field).
```

## Next session focus: localize + fix the COARSE-cell band

Goal: find which interface operation elevates the COARSE interface cell's rms, then
fix it so `roughness(band) ~ interior` for a smooth field AND the channel is stable.

Concrete plan:
1. **Coarse-cell gate.** Build a coarse-side roughness/error gate: with a smooth
   manufactured field on a band patch, measure the COARSE interface cell's error and
   discrete-Laplacian roughness vs a uniform-COARSE reference (and vs uniform-fine).
   The halo-side audit (gate 1) does not see this; this is the missing piece of
   gate 0/2.
2. **Split predictor vs projection.** `MOBY_PREDONLY` (coarse cell after the
   predictor only) vs `MOBY_PROJONLY` (coarse cell after the projection only) on the
   same manufactured field -> attribute the coarse-cell band to one.
3. **Predictor suspect:** the coarse predictor differences its interface halo (the
   RESTRICT of the fine flux/velocity) against its coarse interior. The restrict
   VALUES are exact for a linear field (audit), but the GRADIENT across the 2:1 face
   (coarse interior vs fine-averaged halo at mismatched locations) may be
   inconsistent -> spurious high-k at the coarse cell. Check d/dn and d2/dn2 of the
   coarse cell across the interface vs a uniform-coarse stencil.
4. **Projection suspect:** the 2:1 divergence/pressure consistency at the coarse
   cell (the phi prolong is injection -> a tangential step in the pressure ghost;
   the divergence operator across the 2:1 face). Check whether the projection alone,
   on a div-free-except-grad manufactured field, leaves a coarse-cell roughness.
5. Whatever it is, the fix must keep: bit-exact no-interface, CPU==GPU, conservation
   round-off, and `roughness(band) ~ interior`; validate with the channel (band gone,
   stable) and gates 3/4/6.

## Session 2026-06-25b: localized (not yet fixed). What is now KNOWN.

Reproduced the band end-to-end and ran the plan above. Tools built this session
(committed): `tools/channel_band_profile.py` (per-y-row fluctuation rms of a
block-table channel field, reassembled per level; FLAGS the band) and
`tools/interface_coarse_gate.py` (per-coarse-row tangential-Laplacian roughness on
the Beltrami slab). Fast inputs: `validation/beltrami/slab_y_diag.ini` (32^3 1-step
y-band, ~5 s CPU), `uniform32.ini` (no-refine reference), `slab_y_diag500.ini`.

**GROUND TRUTH (GPU, refined_y110 IC from make_channel_restart, 150 steps, stable):**
the band sits on the COARSE interior cell touching each 2:1 interface, BOTH
orientations, **strongest in streamwise u and in p**:
- bottom interface (coarse-ABOVE-fine, `physLow(2)==FACE_FINE`), y=0.643, gj=24:
  u_rms 1.06 (fine) -> **1.38** (coarse) -> 0.93. p likewise elevated.
- top interface (coarse-BELOW-fine, `physHigh(2)==FACE_FINE`), y=1.357, gj=39:
  u_rms 1.14 -> **1.55** -> 1.23 (this is the handout's "1.14->1.55->1.23").
v,w bands are minor. The RAW interpolated IC is SMOOTH across the interface
(u_rms 0.955->0.925->0.891, monotone) -> **the band is SOLVER-CREATED, not an IC /
make_channel_restart interpolation artifact.** It develops over the run (u-band
~+16% by step 80, ~+40% by step 150).

**ATTRIBUTION (clean, on the real signal):** restart from the developed step-150
field and run ONE step. u_rms at gj=24 is IDENTICAL after PREDONLY (1.3875) and
after FULL (1.3875) -- **the projection does not touch the streamwise band**; the
predictor sets it (it even nudges it up, 1.3778->1.3875/step). The projection DOES
set p (full lowers p(gj=24) 1.111->1.021 in one step, but it re-equilibrates
banded). So: **u-band = predictor; p-band = projection (p is only ever written by
the projection)**, and the standing u-band is consistent with being driven by the
banded pressure gradient that the predictor reads each step.

**RULED OUT this session (do NOT re-try):**
- *Projection under-convergence.* niter 6 vs 50 give a BIT-IDENTICAL band at step 80
  (u 1.1178 vs 1.1176; channel) and the Beltrami coarse-cell roughness is
  niter-independent 50 vs 500. The band is STRUCTURAL (the consistent discrete
  interface operators), not an unconverged residual. (Matches the handout's
  Jacobi==Chebyshev stability finding.)
- *Projection of a smooth divergent field.* `MOBY_PROJONLY MOBY_MANUF=0.3` niter=500
  on the slab leaves NO coarse-cell roughness (flat profile). The projection's 2:1
  divergence/pressure operators are consistent for smooth input; they do not
  manufacture the band. (The channel p-band comes from projecting the TURBULENT
  predictor divergence, which k=1 Beltrami/MANUF does not excite.)
- *Coarse-side TANGENTIAL deep-halo reconstruction (the obvious predictor fix).*
  Extending inc-5 to reconstruct the coarse cell's tangential u,w deep halos (cubic
  `q(0)=3q(1)-3q(2)+q(3)`) **BLOWS THE CHANNEL UP** -- unlimited at ~step 150,
  slope-LIMITED at ~step 150 too (stable to ~100 then diverges). Seed is the FINE
  cell across the interface (gj=47): perturbing the coarse cell feeds back through
  the prolong->fine-halo->restrict loop. This confirms the repo's prior choice
  (inc 5 = normal-only). **The band is NOT an order defect** (coarse tangential is
  already ~2nd order on smooth fields) and raising the halo order is the wrong,
  destabilizing lever. Reverted; baseline restored & rebuilt.
- *reconstruct_interface_halos overall (`MOBY_NORECON`).* Does not change the
  coarse-above-fine band (it only adds the normal-v ghost there).

**WORKING HYPOTHESIS for the fix (next session).** The band is the classic AMR
coarse-fine ENERGY-PILEUP / reflection: the fine band carries more resolved high-k
tangential energy; the coarse cell adjacent to it cannot dissipate the high-k
content fed across the interface (advectively + via the injected/restricted phi),
so it piles up as a localized streamwise+pressure band that is 2nd-order-consistent
(vanishes as h->0) but unphysical at finite h. Since RAISING order destabilizes,
the productive direction is a BOUNDED, targeted high-k TANGENTIAL filter/dissipation
applied ONLY at the coarse interface band (stable, unlike extrapolation; vanishes
for smooth flow; tune to bring roughness(band)~interior). Secondary lead: the
pressure side -- the phi/pressure PROLONG is still pure injection (a tangential step
into the fine pressure halo; the VELOCITY prolong got the two-pass tangential interp
in fdfd477 but the SCALAR phi prolong did not) -- check whether a tangential-interp
phi prolong reduces the p-band (and thus the predictor-read pressure-gradient that
sustains u). Validate any fix with `channel_band_profile.py` (band gone, 200+ steps
stable), bit-exact no-interface, CPU==GPU.

## Session 2026-06-25c: the band is LOW-k; two fixes tried separately (gated, off by default)

**Spectral character (decisive).** FFT of the per-coarse-row fluctuation at the
interface (channel_band_profile reassembly): the band is **low-wavenumber** in
ALL variables. High-k energy fraction (top quarter of tangential wavenumbers) at
the coarse interface cell is ~0.000-0.002 (= interior), 2dx mode ~0; only p has a
small high-k-in-x component (0.014-0.025). The dominant band (u +58%) is a
LARGE-SCALE amplitude amplification, NOT a tangential checkerboard. (The
"checkerboard" in this doc's earlier sections was the velocity-PROLONG fine-halo
artifact fixed in fdfd477, a different thing.) This explains why both fixes below,
which target high-k / the pressure coupling, do not clean the band.

Both fixes are committed but GATED OFF by default (baseline bit-exact; verified
the channel baseline still gives u(gj=24)=1.3778, div=2.46147e-3 exactly):

**(i) Targeted tangential filter at the coarse interface band** (`MOBY_IFFILT=<alpha>`,
`step.f90 filter_interface_band`): a mean-zero tangential discrete-Laplacian
smoothing of the TANGENTIAL velocity components on the coarse FACE_FINE cell-rows
(normal component left alone for conservation; race-free via a frozen buffer;
re-exchange after). STABLE (div ~3.3e-3, mass round-off) -- unlike the cubic
reconstruction. But it OVER-DAMPS: u(gj=24) 1.378 -> 0.514 (alpha=0.05) / 0.341
(alpha=0.15), BELOW the neighbours (~0.85-1.1). Because the band is low-k (same
scales as the physical streaks), a filter cannot separate the band from the real
turbulent energy -- any alpha that reduces the band also smears the streaks. Not a
clean fix (and tuning-/position-sensitive, against Davide's requirement).

**(ii) Tangential-interp phi/pressure prolong** (`MOBY_PHIINTERP`, comm.f90
`copy_local_scalar_entries` now two-pass + weighted gather on the VAR_P node line,
mirroring the velocity prolong fix fdfd477; threads `blk` through
`exchange_scalar_halos`). STABLE (div 2.459e-3 == baseline). A genuine correctness
improvement (interp > injection for the fine pressure halos) but **negligible on
the band**: u(gj=24) 1.3778 -> 1.3770 (-0.06%), p(gj=24) 1.1108 -> 1.1075 (-0.3%).
The low-k coarse u-band barely couples to the pressure-prolong defect. (NOTE: the
MPI pack/unpack scalar path was NOT updated for interp -- multi-rank with an
interface coinciding with a rank boundary would inject there; single-rank tested.)

**Where this leaves the fix.** The band is a low-k, predictor-momentum-driven
energy amplification at the coarse cell next to the better-resolved fine band
(an AMR resolution-jump / energy-reflection effect at the RESOLVED scales). It is
neither high-k (so filtering/hyperviscosity over-damps) nor pressure-prolong-driven
(so (ii) is inert on it) nor an order defect (so cubic reconstruction is wrong and
unstable). Remaining ideas for next session: (a) an ENERGY-CONSERVATIVE (Galerkin)
restrict/prolong pair so the coarse-fine transfer neither creates nor reflects
resolved-scale energy at the interface (the current restrict=simple-average /
prolong=inject pair is not energy-consistent); (b) a skew-symmetric / energy-
conserving advection discretization for the coarse interface cell so its advective
coupling to the fine streaks does not pump energy; (c) accept the band as intrinsic
and instead gate STABILITY (Davide's actual requirement) -- find what tips the
~step-180 blow-up and bound only that. (a)/(b) are the principled directions.

## Session 2026-06-25d: Galerkin (adjoint) transfer pair -- tested, not the lever

The velocity prolong is linear INTERPOLATION (fdfd477 default) while the restrict
is simple AVERAGE -> a NON-adjoint pair. The two adjoint alternatives are
injection-prolong+average-restrict (0th order) and interp-prolong+full-weighting-
restrict (1st order). Tested the prolong half cheaply (`MOBY_VELINJECT`, comm.f90:
injection instead of tangential interp -> the injection+average ADJOINT pair):
channel 150 steps, STABLE (div 2.457e-3), band u(gj=24) 1.378 -> 1.322 (~4%),
gj=39 unchanged. Baseline (velInject off) bit-exact (1.3778) -- the toggle is
inert by default.

So making the transfer pair adjoint barely moves the band. The remaining half
(full-weighting restrict adjoint to the interp prolong) is provably no better:
the band is LOW-k, and EVERY consistent linear restrict (average, full-weighting)
is the identity on low-k modes -- they differ only at O((k h)^2), i.e. on the
high-k content the band does not have. A low-k band is immune to the choice of
(consistent) grid-transfer operator. Hence the full Galerkin pair was NOT built
(intricate gather-map + MPI-pack widening for a provably-null change); the
`MOBY_VELINJECT` diagnostic is kept.

CONCLUSION across all four attempts (cubic reconstruction, high-k filter,
phi-prolong interp, Galerkin transfer pair): the band is a LOW-k, predictor-
ADVECTION-driven energy pileup at the coarse cell next to the fine band, immune to
grid-transfer and filtering. The one untried principled lever is the
DISCRETIZATION: a skew-symmetric / kinetic-energy-conserving advection at the
coarse interface cell so its advective coupling to the resolved fine-side streaks
neither produces nor reflects energy (Morinishi/Verstappen-style energy-conserving
convection, restricted to the interface-adjacent coarse cells). Alternative if that
also fails: accept the band as intrinsic to the 2:1 resolution jump and gate only
STABILITY (Davide's real requirement) -- bisect what tips the ~step-180 blow-up
(CFL at the interface? the band amplitude crossing a threshold?) and bound that.

## Session 2026-06-25e: bug hunt (task i) -- NO bug feeds the coarse-cell band

Done before any energy-conserving implementation, to rule out a bug (Davide's
request). All checks point to the band being a real nonlinear effect, not a defect:

- **Periodic borders: NOT a bug.** The band is UNIFORM in x,z (u' vs x at gj=24:
  mean 1.35, no spike at the x=0/Lx or z=0/Lz seam). MOBY_HALO_AUDIT bad-count is
  IDENTICAL with x,z periodic vs x,z walls (292848) -> the bad cells are NOT at the
  seams (else walls would skip them). So the periodic halo wrap is not implicated.
- **Conservative transfers feeding the coarse cell: EXACT for a linear field.**
  MOBY_HALO_AUDIT (linear manufactured field, every exchanged halo vs the transfer
  design): COPY = 0 bad, RESTRICT = 0 bad. The coarse cell's tangential deep-halo
  inputs and the divergence interface face are restrict/copy-filled -> no
  first-order transfer bug at the coarse cell. (The PROLONG "bad" cells are
  fine-block halos that reconstruct_interface_halos OVERWRITES -- the exchange-only
  audit cannot validate them, and they feed the FINE block, not the coarse cell.)
- **Reconstruction (incl. its uniform-spacing cubic on the stretched grid): ~3.5%.**
  MOBY_NORECON on the channel: u(gj=24) 1.378 -> 1.330, gj=39 1.548 -> 1.536; band
  persists. So the predictor deep-halo reconstruction is not the cause. (There IS a
  real minor stretched-grid inexactness -- the cubic 3q1-3q2+q3 assumes uniform
  spacing; worth fixing for accuracy but it is not the band.)
- **Not the stretched grid.** The band SEED appears on a UNIFORM grid too:
  refined-vs-uniform-coarse Beltrami at IDENTICAL dt (5 steps, dtmax=1e-3) deviates
  most at the interface-adjacent coarse cell (j=24: dv 3.2e-3, du 1.6e-3), decaying
  inward -- on a uniform grid.
- **Not projection convergence** (niter 6==50, earlier) ; **not high-k** (FFT: low-k).

WHICH TERM / WHY-CHANNEL-NOT-MANUFACTURED: the band is the PREDICTOR's advection at
the coarse interface cell (projection doesn't touch u; reconstruction ~3.5%; the
remainder is the coarse cell's nonlinear advective coupling to the resolved fine-
side streaks). It DOES appear in a manufactured field -- in the FULL Beltrami step
(predictor+projection), interface-localized at matched dt -- just small, because
k=1 Beltrami has little interface energy/shear; the channel has a lot, so the band
is large. Since the coarse cell's LINEAR inputs are exact, the band is a NONLINEAR
(kinetic-energy) effect, not a linear transfer bug.

CONCLUSION (task i): no bug. The hypotheses hold -- low-k, predictor-advection-
driven, kinetic-energy pileup at the coarse cell adjacent to the better-resolved
fine band. => proceed to (ii): skew-symmetric / KE-conserving advection at the
interface-adjacent coarse cells (Morinishi/Verstappen), so the coarse cell's
advective coupling to the fine-side resolved streaks neither produces nor reflects
kinetic energy. (Diagnostic kept: MOBY_VELINJECT in comm.f90 = injection vs interp
velocity prolong, the adjoint-pair toggle.)

## Session 2026-06-25f: literature check -- the approach is sound AND sharpened

Confirmed against the literature before implementing (Davide's request). The band
is a KNOWN phenomenon and the energy-conserving fix is canonical:

- **Cevheri & Stoesser, Int. J. Numer. Methods Fluids 82 (2016)** (local mesh
  refinement for LES): explicitly reports "significant energy accumulation when the
  grid is suddenly coarsened" in turbulent flow -- exactly our coarse-side band.
- **Verstappen & Veldman, J. Comput. Phys. 187 (2003) 343-368** (symmetry-
  preserving discretization): the convective operator is a SKEW-SYMMETRIC matrix,
  the diffusive a symmetric pos-def one => the scheme conserves mass, momentum and
  (inviscid) kinetic energy and is STABLE ON ANY GRID. It explicitly covers LOCAL
  GRID REFINEMENT. On coarse grids it "behaves nicely" where truncation-minimizing
  schemes make the resolved fluctuations "too high" (their §6) -- our band.
- **Morinishi et al., JCP 143 (1998)**: fully-conservative finite-difference
  convection (divergence/advective/skew-symmetric forms), uniform grids.

PRECISE CONDITION (V&V §2.1.2, decisive): the discrete convective operator is
skew-symmetric (energy-conserving) IFF (a) the velocity-to-face interpolation uses
CONSTANT 1/2 weights -- NOT mesh-size/metric-dependent weights -- and (b) the
convective flux through a shared face is SINGLE-VALUED ("computed independent of
the control volume in which it is considered"). Minimizing local truncation error
with metric weights gives nonzero diagonal entries on nonuniform grids => breaks
skew-symmetry => energy not conserved. V&V deliberately choose constant weights
(energy conservation + stability) over truncation order on nonuniform grids.

This is GLOBAL (Davide's point): the WHOLE convective operator is skew-symmetric;
the refinement interface must simply not violate it. The cross-grid transfer fits
the same frame -- the discrete gradient = transpose of the divergence (their
§2.1.3), i.e. restriction = volume-weighted transpose of prolongation (the adjoint/
Galerkin condition), AND constant-1/2 interpolation.

TENSION THIS EXPOSES IN OUR CODE (key for next session): our accuracy-motivated
interface fixes use GRID-DEPENDENT weights that VIOLATE V&V's constant-1/2
condition -- the fdfd477 velocity prolong (node-line/metric tangential interp) and
the cubic deep-halo reconstruction (3q1-3q2+q3). They were chosen for truncation
order; V&V say that is exactly what destroys interface energy conservation. The
energy path likely wants constant-1/2, single-valued interface convective fluxes
(lower formal order at the interface, but stable + no energy pileup -- V&V's
explicit trade). The MOBY_VELINJECT proxy moved the band only ~4% because adjoint
transfer ALONE is insufficient: V&V need BOTH constant-1/2 interpolation AND a
skew-symmetric convective FORM at the interface, not just an adjoint transfer pair.

VERDICT: the skew-symmetric / symmetry-preserving direction is well-founded and is
the literature-standard cure for this exact band. Implement it as a GLOBAL
symmetry-preserving convective discretization (the interior central scheme is
already 1/2-weighted = skew-symmetric for div-free; the work is the interface:
single-valued constant-1/2 convective flux + adjoint transfer), per V&V 2003.
Next-session prompt: docs/next_session_skew_symmetric.md.

## Session 2026-06-25g: MOBY_KEBAL gate BUILT + V&V attribution MEASURED (commit 47961fe)

Built the kinetic-energy-balance gate (next_session_skew_symmetric.md step 1) and
used it to measure the V&V prescription. The gate is `print_step_ke_balance` in
main.f90 (mirror MOBY_STEPDIV), off by default (`MOBY_KEBAL`), bit-exact unset.

**What it reports.** Per RK step, the convective KE production `Sum vol*u*C(u)` (no
pressure/diffusion/forcing, stencil identical to `momentum`), split BAND (cells in a
FACE_FINE *or* FACE_COARSE interface row -- both sides) vs INTERIOR, in two forms:
- **DIV** = `vol*u*C_div`: the divergence-form production the scheme actually uses.
  It telescopes to `-Sum KE*(div u)`, so it is **contaminated by the projection's
  residual divergence** (large at niter=6) -- not a clean interface signal.
- **SKEW** = `vol*( u*C_div + 1/2 u^2 * Div_cv )`, `Div_cv` = the velocity divergence
  on the component's control volume (same face neighbours as the convective fluxes).
  This is the production of the skew-symmetric form `C_skew = C_div + 1/2 u(div u)`;
  it is **identically zero in the interior for ANY field** (constant-1/2 telescoping,
  divergence-state independent), so the band SKEW value is the **pure interface energy
  defect** -- the pass/fail. (Sign is PLUS: a minus gives exactly 2*DIV; verified by
  the uniform interior going to round-off only with the plus sign.)

**Validation.** Uniform (no interface): interior SKEW = **2.4e-14** (round-off), band
empty. Beltrami y-slab (`slab_y_diag.ini`, 1 step, ~5 s CPU): band SKEW = **+1.27e-1
PRODUCTION** -- the defect, on the clean uniform-grid testbed.

**Slab toggle sweep (band SKEW, 1 step) -- the key result:**

| config | band SKEW | vs base |
|---|---|---|
| baseline (no reflux, recon, metric) | +0.127 | -- |
| reflux | +0.124 | -2% (inert on energy) |
| reflux + NORECON (no cubic) | +0.056 | **-56%** |
| reflux + NORECON + VELINJECT (const-1/2 prolong) | +0.045 | **-65%** |
| NORECON only | +0.054 | -58% |

=> the **cubic deep-halo reconstruction** (the increment-5 `3q1-3q2+q3`, a non-
constant-1/2 metric extrapolation) is the **DOMINANT interface energy producer**; the
**metric prolong** is secondary; the single-valued momentum **reflux is ~inert on
ENERGY** (it conserves momentum, not energy). This is *exactly* V&V 2003: constant-1/2
weights conserve energy, metric/cubic weights destroy it -- and it explains why all
four prior fixes (none on the constant-1/2 energy path) failed. The residual **+0.045**
(~35%) is the genuine div-form 2:1 skew-defect that needs the skew FORM, not just
constant-1/2 inputs.

**Channel (GPU, refined_y110, IC `make_channel_restart --mode refined --band-cells 24`,
nsteps=150).** Ran STABLE to step 150 (wrote channel_field_150.h5). KEBAL band SKEW
~ -95 (steady, |.|~1.75e4) vs interior SKEW ~ -59 (|.|~1.33e5): band ~25x higher
RELATIVE skew-imbalance -- a real localized interface signal -- BUT the channel's
near-wall **y-stretching contaminates SKEW** (the div form is not energy-conserving on
a stretched grid either, concentrated exactly where the wall band is) and flips its
sign vs the uniform slab. So **KEBAL is cleanest on the UNIFORM slab**; the channel's
decisive gate remains `channel_band_profile.py` (the u_rms gj=24/39 band).

**NOT done (next session):** (1) on the channel, run NORECON+VELINJECT(+reflux) and
measure the PHYSICAL band with `channel_band_profile.py` (gj=24/39) + stability past
step 200 -- the slab energy sweep predicts -56..-65% energy; test whether the u_rms
band actually shrinks AND stays stable (V&V's order-for-energy trade; the prompt says
gate the cubic/metric prolong OFF on the energy path). (2) For the residual +0.045,
implement the skew FORM at the interface band (add `1/2 u(div u)` at band cells to turn
div->skew there; caveat: it breaks LOCAL momentum conservation but -> 0 as the
projection converges -- weigh against the reflux). Scratchpad files (regenerable):
`IC.h5`, `chan_ke.ini` (refined_y110 + that IC + nsteps=150), `s1_reflux.ini`.
RUN GPU CASES ONE AT A TIME and pkill stray `build_*/main` between runs (they
accumulate and starve everything).

## Session 2026-06-26: skew-symmetric interface convection + STABILITY achieved

Implemented the V&V energy-conserving interface correction and ran the channel
stability test. **Key result: the 2:1 interface no longer limits stability** --
Davide's actual requirement is met -- and it is the **constant-1/2 condition**
(dropping the cubic deep-halo reconstruction + the metric velocity prolong), not
the skew-form term, that does it.

**`MOBY_KESKEW` (commit 4c43741)** -- `skew_interface_correction` in step.f90 adds
`1/2 u (div u)` at the interface-band cells (both sides) into `refluxCorr` so the
existing `reflux_apply` lands it on the predictor RHS, making the band cells'
convection skew-symmetric (`C_skew = C_div + 1/2 u (div u)`) while the interior
keeps the bit-exact divergence form. Gated (requires `momentum_reflux`); inert
without a 2:1 interface => bit-exact.

**Slab band-SKEW sweep (MOBY_KEBAL, the energy each lever removes):**
`baseline 0.127 | reflux+NORECON 0.054 | +VELINJECT 0.045 | +KESKEW 0.036`. Each
V&V lever chips the band energy; the residual **0.036 is the cross-interface
area-mismatch** term the LOCAL skew correction cannot remove (it restores the
within-block telescoping; the 2:1 face needs the volume-weighted ADJOINT transfer,
restrict = transpose of prolong -- the remaining unbuilt piece).

**Channel stability (refined_y110, GPU, 250 steps, MOBY_STEPDIV):**
| run | result |
|---|---|
| BASELINE (cubic recon + metric prolong) | **blows up ~step 200** (div_l2 -> 2.3e8; div_max climbs 0.09->0.24 from step ~180) |
| NORECON+VELINJECT (constant-1/2, no skew) | **STABLE to 500** (2.5x past blow-up); div_max DECREASES 0.09->0.056, div_l2 <= 6.4e-3 -- the channel SETTLES, not just bounded |
| + KESKEW (local skew) | **STABLE to 250**, identical to NORECON+VELINJECT (div_l2 <= 0.062) |

=> the truncation-optimal **cubic/metric weights break interface energy
conservation and DESTABILIZE**; **constant-1/2** (V&V's explicit order-for-energy
trade) cures it. KESKEW adds energy accuracy (-20% more band energy) but no extra
stability here. The coarse-side **amplitude band PERSISTS** (u_rms ~1.5, intrinsic
to the resolution jump -- the coarse side genuinely cannot resolve the fine side's
content) but is now **bounded/stable instead of growing to blow-up**: "band gone"
was the wrong target; "band does not limit stability" is the achievable + correct
one (matches V&V -- energy conservation bounds the resolved fluctuations, it does
not erase the resolution jump).

**Production path for Davide (decision needed).** The stabilizing levers are
currently env-var DIAGNOSTICS (`MOBY_NORECON` = drop the cubic increment-5
reconstruction; `MOBY_VELINJECT` = inject vs metric-interp the velocity prolong).
To ship: make constant-1/2 the DEFAULT on refined runs -- either disable the cubic
reconstruction or, better, REPLACE it with a constant-1/2 (linear-average) ghost so
some interface accuracy survives while satisfying V&V; and make the velocity prolong
constant-1/2. This changes the no-flag refined behaviour (NOT bit-exact vs current),
so it is a deliberate choice -- single-level runs stay bit-exact (no interface).
`MOBY_KESKEW` can stay an optional energy-accuracy refinement. NOT yet done:
(a) the adjoint (volume-weighted transpose) transfer for the residual 0.036 band
energy; (b) turbulence-statistics validation (t=5..25 stats leg) that the bounded
band does not corrupt the mean profile. (Longer-run stability is settled:
NORECON+VELINJECT runs 500 steps with div_max DECREASING to 0.056.)

## Memory
`interface-validation-suite` (the gates + this state), `corner-reconstruction-todo`,
`restart-overrides-config-sor`, `refinement-perf-profile`, `chebyshev-jacobi-plan`.
