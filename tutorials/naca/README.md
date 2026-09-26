# NACA 0012 tutorials

**The live case is [`rans/`](rans/README.md)** — the converged alpha = 5,
Re_c = 4e5 k-omega SST validation against a body-fitted OpenFOAM reference
(C_L 0.5199 vs 0.5142, Cp_min matching to four digits), with its drivers
(`run_case.sh`, `postprocess.sh`), prepare inis and post-processing scripts
(`rans/postProcess/`: control-volume forces, surface Cp/Cf, the OpenFOAM
overlay). Start there.

> **History (2026-09-25/26).** This directory used to hold an earlier
> polar-sweep generation (aoa -2..5 at Re 4e5 against XFOIL/OpenFOAM, on the
> R2D xz-quadtree stack with `span = y`) and its post-processing scripts. They
> were REMOVED: `rans/` replaces them, and their C_L/C_D came from the
> penalization integral, which the solver no longer computes — the runtime
> force is the control-volume momentum budget over `[case.airfoil] cv_box`
> (`docs/next_session_cv_forces.md`). The live versions of the Cp/Cf and
> force scripts are in `rans/postProcess/`. Recover the sweep from git
> history if it is wanted again; it will need a `cv_box` and a re-measured
> reference.

## Modelling notes that still apply

- **Resolution.** A resolved-wall SST needs y+_1 of order 1-3 over the
  turbulent chord (cf ~ 4e-3 at Re 4e5 puts that at a surface spacing of a
  few 1e-4 c). The xz quadtree (`[blocks] refine_dims = xz`, span along y
  via `[case.airfoil] span = y`: chord x, LIFT z) makes such grids cheap —
  the same resolution in 3D octree mode is two orders of magnitude more
  cells. `rans/` reaches c/12288 at the nose with `refine_body_levels` +
  `refine_body_box`.
- **Fully turbulent vs transitional.** A fully-turbulent run (`[rans]
  transition` off) is what XFOIL with forced transition at the LE
  (`vpar -> xtr 0.01 0.01`) and a plain OpenFOAM SST compare against; it
  also avoids the gamma-Re_theta_t fine-grid unreliability found with the
  first-order-upwind transition scalars (the front smears over ~100 cells
  and separation-induced transition stops firing on finer grids — see
  CLAUDE.md, R2D-3 and A3 increment 3). `rans/` instead reproduces
  OpenFOAM's fvOptions forced transition with `[rans] kpin_box` /
  `ktrip_box`.
- **Forces.** Use the runtime control-volume budget (`cv_box`, a TIGHT box
  around the profile) and a clean stored pressure (a higher projection
  `niter`; long niter = 6 IBM runs carry a velocity-neutral pn mode that
  pollutes the border integrals). Buried interior blocks may be removed —
  the budget does not see them.
- **Blockage.** A Dirichlet far field at ~12c costs ~10-15 % of the 2 pi
  lift slope; compare slopes with that in mind or use a larger box (`rans/`
  uses 128c x 96c).
- **Surface Cp/Cf** from an immersed boundary: see the method notes at the
  top of `rans/postProcess/surface_cp_cf.py` (Cf from the penalization band,
  Cp by wall extrapolation).
