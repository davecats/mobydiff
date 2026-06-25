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

## Memory
`interface-validation-suite` (the gates + this state), `corner-reconstruction-todo`,
`restart-overrides-config-sor`, `refinement-perf-profile`, `chebyshev-jacobi-plan`.
