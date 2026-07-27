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

## Validation results — UNDER REVALIDATION (2026-07-28)

CAUTION: the t = 29 numbers below came from the BoostConv-assisted leg,
which was found to KILL the pressure-side trip (strip k 7.6e-3 ->
7.9e-6; lower Cf halved — the recombination overwrites the localized
k source; user diagnosis from the Cf overlay). A PLAIN rerun from the
pre-boost t = 26 state is producing the replacement numbers; the table
and figures will be regenerated from it.

## Results table (t = 29, boosted leg — transition corrupted, see above)

|                | mobydiff (IBM) | OpenFOAM (body-fitted) | delta |
|----------------|----------------|------------------------|-------|
| C_L (CV)       | 0.516 +- 0.005 | 0.5142                 | +0.4 % |
| C_D (CV)       | 0.0134 +- 0.0004 | 0.0134               | exact |
| Cp_min (depth-converged) | -1.787 | -1.780               | +0.4 % |
| TE separation  | x/c ~0.99      | 0.996                  | -- |

(The earlier t = 26 numbers -- C_L 0.506, Cp_min -1.763 -- were the
still-converging circulation; the peak followed the circulation-scaling
prediction exactly. Converged state: c11_v2_430013.h5. The final leg
ran WITH the BoostConv accelerator: [rans] boostconv alpha 0.02, N 20,
p 25; overhead +2.4 % s/step; plateau reached within ~1 t.u. of
activation.)

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
