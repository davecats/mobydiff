# Next session — energy-conserving (skew-symmetric) 2:1 interface convection

Read `docs/interface_band_handout.md` FIRST (full state; the "Session 2026-06-25b..f"
sections are this investigation). Memory: `interface-validation-suite`. Branch
`claude/jacobi-interface`, HEAD around `87a97bd` (checkpoints 1-3 of this work).

## The settled diagnosis (do NOT re-litigate)

The 2:1 wall-band refined channel (`validation/channel_interface/refined_y110`)
grows a **spurious band on the COARSE interior cell of every interface**, strongest
in streamwise u and pressure (u_rms ~1.14->1.55->1.23 across y; gj=24 and gj=39 on
the reproduction). It is:
- SOLVER-created (the raw interpolated IC is smooth across the interface), grows
  over the run, blows up ~step 180;
- LOW-wavenumber (FFT: high-k fraction ~interior) -- a large-scale amplitude
  amplification, NOT a checkerboard;
- PREDICTOR-driven (1-step-from-developed-field: projection leaves u unchanged;
  the predictor sets it);
- NOT a bug: periodic borders clean; COPY/RESTRICT halos exact for a linear field
  (MOBY_HALO_AUDIT); reconstruction ~3.5% (MOBY_NORECON); not stretch (appears on a
  uniform grid); not projection convergence (niter 6==50);
- a **kinetic-energy pileup at the coarse-fine interface** -- a documented LES/AMR
  phenomenon (Cevheri & Stoesser 2016; Verstappen & Veldman 2003).

Four fix classes are RULED OUT (gated toggles kept as diagnostics, all off by
default, baseline bit-exact): cubic coarse-tangential reconstruction (unstable);
high-k tangential filter `MOBY_IFFILT` (over-damps a low-k band); phi-prolong
tangential interp `MOBY_PHIINTERP` (inert); adjoint transfer pair `MOBY_VELINJECT`
(~4% -- adjoint transfer ALONE is insufficient).

## The approach (literature-confirmed): Verstappen-Veldman symmetry-preserving

Make the discrete CONVECTIVE operator skew-symmetric across the 2:1 interface so it
conserves kinetic energy (V&V 2003, JCP 187:343, which explicitly covers local grid
refinement; Morinishi 1998 for the forms). This is GLOBAL: the whole convective
operator is skew-symmetric. The interior 2nd-order central scheme is ALREADY
1/2-weighted = skew-symmetric for a div-free field, so the work is the interface.

PRECISE conditions (V&V §2.1.2-2.1.3), the spec for the fix:
1. **Single-valued interface convective flux**: the advective momentum flux through
   the 2:1 shared face must be computed IDENTICALLY from the coarse and fine sides
   ("independent of the control volume in which it is considered"). The momentum
   reflux matches the momentum flux but the ENERGY-conserving requirement is that
   the flux FORM be single-valued at the face.
2. **Constant 1/2 interpolation weights** for velocity-to-face interpolation at the
   interface -- NOT mesh-size/metric weights. This is the crux: metric weights
   minimize truncation error but make the convective matrix's diagonal nonzero ->
   break skew-symmetry -> energy not conserved (V&V show this explicitly).
3. **Adjoint transfer**: restriction = volume-weighted TRANSPOSE of prolongation
   (the discrete gradient = -divergence^T structure), so cross-grid transfer
   preserves skew-symmetry.

## Do this, in order

1. **Build a discrete KE-conservation gate.** With viscosity off (or measured
   separately) and the interface present, the convective operator must change total
   discrete kinetic energy Sum(vol * 1/2 |u|^2) by ~0 per step (round-off). Build a
   `MOBY_KEBAL` monitor (mirror MOBY_STEPDIV): the convective contribution to dKE/dt,
   split interface-band vs interior. TODAY it should show net KE PRODUCTION localized
   at the coarse interface band (that IS the band). This is the pass/fail gate.
   Cheap testbed: the Beltrami slab (`validation/beltrami/slab_y_diag.ini`, CPU,
   ~5 s) full step -- the band reproduces there (small) without turbulence.
2. **Read `momentum` in `step.f90`** (the divergence-form convection; see
   `compute_momentum_terms` in main.f90 for the exact flux stencils -- it already
   uses (u_i+u_{i+1})/2, i.e. 1/2 weights, in the interior). Identify the
   interface-adjacent flux terms (the coarse cell's d(u v)/dy etc. reading the
   restricted/reconstructed halo).
3. **Make the interface convective flux single-valued + constant-1/2.** Replace the
   coarse cell's interface-normal advective flux so it equals the (constant-1/2,
   energy-conserving) flux the fine side computes, summed conservatively. Likely
   reuses the reflux machinery (restrict the fine flux to the coarse) but in the
   ENERGY-conserving form, and must NOT use the metric-weighted prolong / cubic for
   the velocities entering this flux (those break the constant-1/2 condition --
   gate them off on the energy path; `MOBY_VELINJECT` already toggles the prolong).
4. **Verify**: KEBAL band production -> round-off; then channel band re-measured
   (`tools/channel_band_profile.py` -- gj=24/gj=39 spike gone), stable past step
   ~200 (the old onset), conservation round-off (MOBY_STEPDIV), bit-exact with NO
   interface, CPU==GPU.

## Watch out
- V&V trade local truncation ORDER for energy conservation on nonuniform grids;
  the interface will be lower formal order but stable + pileup-free. Don't chase
  order at the interface -- chase the KE balance (that is the whole point).
- The existing accuracy fixes (fdfd477 metric prolong, the reconstruction cubic)
  fight this; expect to gate them off (or replace) on the energy path and re-verify
  the momentum-interface gates (`validation/momentum_interface`) still hold for the
  conservative quantities.
- Keep the conservative RESTRICT of the interface FACE (mass/divergence) -- it is
  exact and must stay; the change is the CONVECTIVE FLUX form, not the divergence.
- CPU `-Mnofma` reference, GPU `-Mnofma -gpu=nofma` bit-identical; run ONE GPU case
  at a time; the band reproduction recipe + tools are in the handout.

## Tools / env (all off by default)
`tools/channel_band_profile.py` (the band gate, reassembled per level),
`tools/interface_coarse_gate.py` (Beltrami coarse-row roughness),
`tools/make_channel_restart.py --mode refined --band-cells 24` (the IC),
`MOBY_STEPDIV MOBY_NORECON MOBY_HALO_AUDIT MOBY_IFFILT MOBY_PHIINTERP MOBY_VELINJECT`
+ build the new `MOBY_KEBAL`.
