# 2:1 interface — edge/corner validation (embedded core patch)

Validates the 2:1 refinement at **edges** (two interface faces meeting) and
**corners** (three) in turbulence — the flat-face (y-band) case validated on
2026-06-29 has only y-face interfaces. See `docs/next_session_edges_les.md`
open item 1. Production interface settings: `interface_constant_half = true`
(default), `momentum_reflux = false`.

## The case

A single level-1 box floating in the **core** of an Re_tau=180 channel (base
128×64×128, nb=8): x-blocks 6–9, z-blocks 6–9, y-rows 2–5 (a 4×4×4 base-block
patch → 2496 leaves, 512 fine). It exercises **all** interface orientations: 6
faces, 12 edges, 8 corners, all genuine fine↔coarse (no degenerate wall edges).

Why the core, not a wall patch (the decision):
- **Clean attribution.** The patch sits where the base grid is already adequate
  (the y-interfaces land at y+≈54.5, log layer; the patch straddles the
  centreline), so a clean run gives refined ≈ uniform-128/256 and *any*
  localised excess at a face/edge/corner is unambiguously the interface
  treatment — not a coarse-resolution deficit (the wall-patch confound).
- **Full orientation coverage + symmetry.** y-centred ⇒ the two y-interfaces are
  mirror images about the centreline (orientation-symmetry check), with real
  mean shear across them.
- Caveat: the core is a *gentle* shear test. A buffer-region off-wall patch is
  the follow-up stress case before trusting a true wall deployment.

## Files

- `core_patch.ini` — the case (restart from `IC.h5`, 250 steps, writes
  `channel_field_250.h5`).
- `IC.h5` — interpolated from the KMM180 restart, leaf layout matching the
  solver's enumeration (verified bit-identical `blocks` table):
  `python3 tools/make_channel_restart.py --mode patch --refine-box <the 6
  numbers from [blocks] refine> --out IC.h5`.
- `uniform_patch.ini` — gate 2 (uniform-flow exactness), a periodic cube with
  the same central-patch geometry and the production interface settings.

## Gates run (2026-06-29, GPU)

1. **Bit-exact no-interface** — N/A this round: the change is Python-only
   (`make_channel_restart.py` gained a box leaf table); no Fortran diff, so the
   single-level path is byte-identical by construction.
2. **Uniform-flow exactness through the patch** — PASS. `uniform_patch.ini`
   (const (1.0, 0.5, 0.25), periodic, no forcing = exact steady state): after 10
   steps `max|u−1|=max|v−0.5|=max|w−0.25|=0.0`, p constant. The transfer
   (RESTRICT/PROLONG/copy) is exact for constants at every face/edge/corner under
   const-½ + reflux-off.
3. **250-step turbulent stability** (`MOBY_STEPDIV=1`) — PASS. Step 1
   `div_max=9.24` (interpolated-IC transient, absorbed by the first projection),
   then post-transient `div_max` peaks 0.078 and **decays to ~0.012**; `div_l2`
   monotonically decreasing; `max|mass|=1e-15`; no NaN. Well under the benchmark
   threshold (`div_max < 0.5`); the old `const-½=false` path blows up to ~1e11
   around step 200.

4. **Edge/corner banding metric** — `tools/patch_interface_diff.py` (built
   2026-06-29). Reference-based (a localised patch has no homogeneous direction to
   average, and the band is low-k so a reference-free roughness proxy misses it):
   the patch run's coarse cells outside the patch start bit-identical to the
   base-128 control (`base.ini` from `--mode base`), so the cell-by-cell
   difference after the same steps, classified by each coarse leaf's adjacency to
   the patch (interior/face/edge/corner), isolates the interface footprint. Gate:
   edge/corner not worse than face. FIRST READING (250-step snapshot, GPU):

   ```
   class    nleaf    u_rmsΔ     v_rmsΔ     w_rmsΔ     p_rmsΔ
   interior  1832   7.00e-03   4.63e-03   5.65e-03   2.18e-01
   face        96   1.28e-01   1.26e-01   1.29e-01   3.44e-01
   edge        48   4.74e-02   3.81e-02   4.28e-02   2.21e-01
   corner       8   1.60e-02   7.44e-03   1.55e-02   1.19e-01
   ```
   (shell 96+48+8 = 152 = 6·16 + 12·4 + 8·1 ✓). Edge/corner are 0.30–0.37× /
   0.06–0.12× the face level — progressively SMALLER (contact-area-natural: face >
   edge > corner). A band pathology would make edge/corner anomalously HIGH; the
   opposite is seen. PASS (first reading). CAVEAT: t=0.0025 (IC-dominated) — the
   magnitudes are not converged turbulence; the relative ordering is the signal.
   The quantitative statement needs the developed run (gate 5).
   ```bash
   python3 tools/patch_interface_diff.py \
       validation/channel_interface/core_patch/channel_field_250.h5 \
       validation/channel_interface/core_patch/base_field_250.h5
   ```

## Gate 5 — developed time-averaged banding (run on the cluster)

The converged version of gate 4. For the CORE patch the per-y-row
`channel_stats` are useless (the patch breaks x,z homogeneity, so a row average
smears its localised effect), so this uses **3D time-statistics**: dump field
snapshots over the averaging window for both the patch run and a base-128
control, then compute per-cell time-mean + fluctuation-rms and aggregate by
interface-adjacency class.

Driver (two legs per case, mirrors `../developed/run_developed.py`):

```bash
cd validation/channel_interface/core_patch
python3 run_gate5.py --arch gpu --ranks 2          # runs patch + base
# (or --case patch / --case base; --t-transient 5 --t-average 20 --snap-interval 320)
```

It writes `runs/gate5/{patch,base}/stats/*.h5` (snapshots; ~200 per case at
`dtmax`, ≈7–8 GB/case — raise `--snap-interval` to use less). Then the metric:

```bash
python3 ../../../tools/patch_interface_stats.py \
    --patch 'runs/gate5/patch/stats/channel_field_*.h5' \
    --base  'runs/gate5/base/stats/base_field_*.h5'
```

Two diagnostics, both classified interior/face/edge/corner and matched
cell-by-cell (patch and base share the base lattice outside the patch):
- **MEAN footprint** `rms(<q>_patch − <q>_base)` — interior ≈ 0; face/edge/corner
  = the mean interface transport footprint. Gate: edge/corner not worse than
  face; small.
- **BAND ratio** `<fluct_rms>_patch / <fluct_rms>_base` — interior ≈ 1.0; a
  band shows as a ratio > 1 at face/edge/corner. Gate: ratios ≈ 1, edge/corner
  not worse than face, consistent with the flat-face ~5% small-scale v'/w' loss
  (ratio ≲ 1 if anything).

The uniform-256 reference (`../reference.ini`, `../run_reference.sh`) and the
uniform-128 isolation control (this `base` case) complete the picture, as for the
flat band.

### Result (2026-06-29, 201 snapshots/case over t = 5..25) — PASS

`gate5_metric.log`. No edge/corner band in developed turbulence:

- **BAND ratio** (`<fluct_rms>_patch / <fluct_rms>_base`, matched cells):
  interior 1.00 (method floor), face/edge/corner all within **1–2%** for u,v,w
  (max corner v' 1.022, on 8 cells). Contrast the old flat-face u' spike ~1.5.
- **MEAN footprint** (single-run deviation from the x,z-homogeneous row mean):
  the apparent face/edge/corner excess over interior is a few-cell SAMPLING
  artifact — the patch-free `base` run shows the SAME pattern, stronger (base
  corner w 1.42× interior vs patch 1.20×), and the patch run's class deviations
  are ≤ the base run's in every class. No mean distortion above the noise floor.
- p is uninformative (floating pressure datum; the base run's global p drift
  inflates the band-ratio denominator → p shows as ~0).

Conclusion: const-1/2 + `momentum_reflux=false` is clean at edges and corners,
matching the flat-face result. Open item 1 of `docs/next_session_edges_les.md`
(edge/corner validation) is DONE; open item 2 (LES across the interface) remains.

### Cross-section figures

`tools/plot_patch_slice.py PATCH.h5 BASE.h5 --out slice.png [--axis z|x|y]`
reassembles both runs onto the finest lattice (base upsampled, so it renders
uniformly coarse) and plots u,v,w,p without vs with refinement on the true
node coordinates, patch outline dashed. `slice_xy.png` (z-normal) and
`slice_zy.png` (x-normal) through the patch centre, from the developed
(t≈25) snapshots.

### Side observation — base-128 pressure null mode (NOT an interface issue)

In the developed base-128 control the PRESSURE grew a large-amplitude mode
(pn std ≈ 2e5, ±3e6) while the velocities stayed healthy (un/vn/wn std
identical to the patch run). It is a velocity-decoupled pressure null mode of
the under-converged niter=6 Chebyshev-Jacobi projection on the uniform run; the
refined-patch run does NOT show it (pn ±13). It does not affect the velocity-
based gate-5 result, and it is why the gate-5 p column is uninformative (the band
ratio denominator is huge → p ≈ 0). Worth a separate look at pressure pinning /
niter on the uniform channel, independent of the 2:1 interface.
