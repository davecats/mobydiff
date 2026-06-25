# Next session — localize and fix the 2:1 interface COARSE-cell band

Read `docs/interface_band_handout.md` first (full state, what is ruled out, the
validation pipeline). Memory: `interface-validation-suite`. Branch
`claude/jacobi-interface`, HEAD around `fdfd477`.

## The one problem to solve

A 2:1 wall-band refined turbulent channel (`validation/channel_interface/
refined_y110`) blows up at the interface (~step 180). The cause is a **spurious
high-k band on the COARSE-INTERIOR cell of every 2:1 interface** that is present
already in the clean field (all variables; u' rms ~1.14 -> 1.55 -> 1.23 across y at
the interface) and is advected into the blow-up. The velocity-prolong checkerboard
(fine halos) was a separate artifact and is already FIXED (commit fdfd477) -- do
NOT re-chase it; the band is on the COARSE cell, set by the predictor/projection,
not the halo exchange. Reflux, the projection accelerator, niter/sor are all ruled
out (handout).

Davide's requirement: a 2:1 interface method that does NOT limit the stability of
the underlying scheme and is not too sensitive to the interface position.

## Do this, in order

1. **Build the coarse-side gate** (the missing metric). On a band-refined patch
   with a SMOOTH manufactured field (Beltrami, or a smooth periodic field; NOT a
   non-periodic linear field -- it breaks at the periodic wrap), measure the COARSE
   interface cell's error and discrete-Laplacian roughness, BAND vs INTERIOR, vs a
   uniform-coarse reference. Reuse `tools/interface_diagnostics.py` (it already does
   error + roughness band/interior for Beltrami) and add the reference-free FIELD
   roughness variant for the channel. The gate must FLAG the current band.

2. **Attribute it: predictor vs projection.** Same manufactured field, one RK
   sub-step: `MOBY_PREDONLY` (coarse cell after predictor only) vs `MOBY_PROJONLY`
   (after projection only). Which one creates the coarse-cell roughness?

3. **Predictor branch:** the coarse predictor differences its 2:1 interface halo
   (the RESTRICT of the fine field) against its coarse interior. Restrict VALUES are
   exact for linear (audit), but the wall-normal gradient/Laplacian across the 2:1
   face mixes the coarse interior with the fine-averaged halo at MISMATCHED metric
   locations -> check d/dn and d2/dn2 of the coarse cell across the interface vs a
   uniform-coarse stencil; the fix is the coarse-side analogue of the fine-side
   reconstruction (a consistent interface stencil / metric on the coarse cell).

4. **Projection branch:** the 2:1 divergence/pressure consistency at the coarse
   cell (the phi prolong is pure injection -> a tangential step in the pressure
   ghost). On a manufactured globally-conserving-but-not-div-free field
   (`MOBY_MANUF`), `MOBY_PROJONLY` with large niter must drive div->0 AND leave NO
   coarse-cell roughness; if it does not, the 2:1 divergence/gradient operators at
   the coarse cell are inconsistent.

5. **Fix it** so a smooth field gives `roughness(band) ~ interior`, then validate:
   bit-exact with NO interface, CPU==GPU, conservation round-off, the channel band
   GONE (re-measure the y-row rms), and the channel STABLE past step ~200 (the old
   onset). Run gates 3, 4, 6 from the suite (`validation/interface_suite/README.md`).

## Don't

- Don't re-chase the velocity-prolong checkerboard (fixed) or the reflux/accelerator
  (ruled out). Don't use a non-periodic linear field for the channel-side roughness
  (periodic-wrap artifact). Don't declare the band fixed from the audit alone -- the
  audit is halo-side; re-measure the channel coarse-cell rms.

## Validate / tools

`MOBY_PROJONLY MOBY_PREDONLY MOBY_MANUF MOBY_RHSDUMP MOBY_TERMDUMP MOBY_DIVDUMP
MOBY_RESLOG MOBY_STEPDIV MOBY_NORECON MOBY_HALO_AUDIT`. Tools:
`interface_diagnostics.py` (roughness band/int + order), `rhsband.py`/`rhsterms.py`
(per-term order), `divsum.py`, `compare_fields.py --export-global`,
`check_beltrami.py`. Build BOTH `-Mnofma` (CPU) / `-Mnofma -gpu=nofma` (GPU) for
bit-exact comparisons; ALWAYS run through `mpirun` even single rank.
