# NACA 0012 validation case — alpha = 5 deg, Re_c = 4e5, vs OpenFOAM

The finalized validation case of the Re 4e5 campaign
(docs/next_session_naca_re4e5.md carries the full investigation trail):
a NACA 0012 at 5 degrees incidence, chord Re 4e5, k-omega SST with the
OpenFOAM-matched configuration, compared against a body-fitted
simpleFoam kOmegaSST reference (48400-cell O-mesh, wall y+ 1.3).

## Configuration (c11_aoa5.ini)

- Grid: 128c x 96c chord-lift domain, nose at (50, 48), span y
  (0.1875c, ny 8), [blocks] refine_dims = xz, 11 levels: base c/3,
  finest c/6144 (first-cell y+ 1.4-2.1), concentric wake-skewed
  per-level boxes + refine_body wall bands. 16042 leaves, 8.2M cells.
- Momentum convection: skew-symmetric (the production form, hardwired;
  docs/next_session_skew_convection.md).
- Turbulence: SST, resolved wall. Ambient MATCHES OpenFOAM: NO sustain
  (their decayControl false), inlet k/omega = their inlet pair
  (tu 21.35 %, nut_ratio 1.0936 -> nut_in = 1.09 nu, free decay).
- Forced transition MATCHES OpenFOAM's fvOptions: [rans] kpin_box
  (k = 0 for x < LE + 0.09c) + two ktrip_box strips at x/c 0.10-0.12
  (rate 0.385 = their 2e-5 absolute over their set volume).
- dt = dtmax = 5e-5 (CFL at the finest band), niter 18 chebyshev.

## Running it

1. Geometry/case file (4-7 min, CPU build):
   mpirun -n 4 --oversubscribe ../../build_cpu/moby_prepare .prep_c11.ini ibm_coeff_c11.h5
2. Either restart from a converged state (c11_aoa5ab_370013.h5 if you
   have it; any converged snapshot works via [restart] file = ...) and
   run a few t.u., or from scratch with the STAGED protocol (3-4x
   cheaper than all-fine): run the same physics on the L10 twin
   (refine_levels 10, dt 1e-4) to t ~ 12, then
   tools/../tutorials/naca/interp_restart.py onto the L11 layout and
   finish here (t ~ 5 more). See the cost notes in
   docs/next_session_naca_re4e5.md.
3. mpirun -n 1 ../../build_gpu/moby_solve c11_aoa5.ini

## Validation results (CONVERGED, plain run, t = 30; stationary t 29-30)

|                | mobydiff C11 | mobydiff C10 (cheap) | OpenFOAM |
|----------------|--------------|----------------------|----------|
| C_L (CV)       | 0.514 +- 0.009 | 0.519 +- 0.002     | 0.5142 |
| C_D (CV)       | 0.0130 +- 0.0007 | 0.0158 +- 0.0006 | 0.0134 |
| Cp_min         | -1.796       | -1.820               | -1.780 |
| trip (both sides) | yes       | yes                  | (fvOptions) |

Converged state: c11_aoa5_450013.h5 (plain marching; the BoostConv-
assisted leg was DISCARDED — the recombination suppressed the
pressure-side trip; docs/next_session_boostconv.md CAVEAT).
C10 VARIANT (c10_aoa5.ini, refine_levels 10, dt 1e-4, ~5x cheaper):
sufficient for LIFT and Cp (1-2 %); NOT for drag (+18 %, the y+ 3-4
resolved-wall penalty — L11 or wall functions needed when C_D is a
target). POST-PROCESSING NOTE: when t_final is not an exact step
multiple, the final write's dt-clipped micro-step inflates the stored
pn scale — extract pressure-based quantities from regular-cadence
snapshots.

- Cp curves overlay within extraction accuracy on BOTH sides
  (cpcf_c11_final_vs_openfoam.png); at full convergence the peak and
  both coefficients MATCH OpenFOAM within the extraction/CV accuracy.
- Cf: laminar dip + trip-jump structure aligned with OF; the laminar
  zone carries a 1.5-1.7x Thwaites excess = staircase roughness (step
  height ~15 % of the laminar delta) — the known first-order-IBM
  signature; turbulent zone matches.
- Extraction: surface_cp_cf.py (k-gated hybrid Cf estimator, Cp
  extrapolation depth 12h — the converged depth), cv_forces.py
  (flux-exact finite-volume borders, wind axes via --aoa 5),
  compare_openfoam.py for the overlay.

## References on disk

- OpenFOAM case: ~/auswertung/20251027_MA_JannikWeber/run (their
  printed kOmegaSST coefficient dict is IDENTICAL to ours, c1 limiter
  and nut clip included).
- XFOIL polars: xfoil_re4e5_n{1,9}.dat (Debian xfoil SIGFPEs on any
  second viscous point: one ALFA per process, parse stdout).
- BoostConv steady-state accelerator study (module committed, gated):
  docs/next_session_boostconv.md.
