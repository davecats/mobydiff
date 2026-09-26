# A0 freestream gates — inlet/outlet patch types + Dirichlet-pressure outlet

Validation for phase A0 of the airfoil plan (`docs/next_session_airfoil.md`):
the `[boundary] <dir>_<side>_patch = wall | patch | inlet | outlet` face
concept (`resolve_face_bcs` derives the per-variable BC rows) and the
Dirichlet-pressure outlet in the projection (outlet face counted `2*d1f*mu`
in the Jacobi denominator, corrected with `d1f` against the MIRRORED phi
ghost, `apply_bc` writing the zero-gradient outflow copy at the predictor
stage only).

Run everything (one solver job at a time) with a built `build_cpu`/`build_gpu`:

    ./run_gates.sh              # or one group: oblique | pois | vortex | ranks | config

Groups map to the phase gates:

- **oblique** (gate b): uniform 20-degree freestream through an empty box,
  3 inlets + 1 outlet. Must be preserved EXACTLY (max deviation 0.0 after
  200 steps, interior divergence 0.0) — constants are in the null space of
  every consistent operator. RESULT 2026-07-12: PASS (all 0.0).
- **pois** (gate c): plane Poiseuille, periodic reference (`pois_ref.ini`,
  forcing 8*nu) vs the inflow/outflow twin (`pois_io.ini`: parabolic
  Dirichlet inlet via `x_min_u_profile = parabola`, outlet at x_max, no
  forcing), restarted from the reference through a solver-minted template
  (`make_freestream_ics.py` — the restart metadata carries periodicity, so
  a periodic-run file cannot restart an in/outflow ini directly).
  RESULT 2026-07-12: profile vs reference 1.6e-3 (= O(h^2)) at x/lx = 0.5
  and 0.9, p(x) slope -0.0796 vs -0.08, nonlinearity 9.5e-5, last-cell
  p = 2.2e-3 (outlet-pinned level), drift 5.6e-16 over the second 5000
  steps. PASS.
  RE-RUN 2026-09-26 at niter 12 (was 6), CPU: profile dev 1.639e-3 / 1.539e-3
  (x/lx 0.5 / 0.9), slope -0.07958, nonlinearity 9.54e-5, last-cell p
  2.41e-3 (was 2.2e-3), drift 1.3e-15 (was 5.6e-16). PASS.
- **vortex** (gate d): Lamb-Oseen vortex (Gamma = 0.443, rc = 0.15) advected
  at U = 1 through the outlet. RESULT 2026-07-12: exits cleanly, perturbation
  energy at t = 2.5 is 5.2e-3 of initial (peak transient ~17% mid-exit,
  decaying); no convective-outlet fallback needed.
  RE-RUN 2026-09-26 at niter 12 (was 6), CPU: exits cleanly, but the
  reflected fraction (E(t=2.5)/E0) is **2.23e-2**, not 5.2e-3. Attributed by
  rerunning the current binary: niter 6 reproduces 5.34e-3, niter 24 gives
  3.23e-2 and niter 60 3.22e-2 — the 5e-3 was an artefact of the
  under-converged niter-6 projection damping the exit transient; the
  converged-projection value is ~3.2e-2. The metric is also noisy at the
  last snapshots (E/E0 at t = 1.75/2.0/2.25/2.5 = 0.049/0.228/0.059/0.022 at
  niter 12). Informational (no pass/fail threshold in the checker).
- **ranks** (gate e): oblique + pois_io, 1 == 4 ranks EXACT (max_abs 0);
  CPU vs GPU <= 1e-12 (measured 0.0). PASS.
  Re-run 2026-09-26 at niter 12 (pois_io; oblique is still niter 6), CPU:
  both 1 == 4 ranks EXACT. PASS (GPU not re-run).
- **config** (gate f): y walls declared `wall` == inferred (no declaration)
  bit-exact; `x_max_p_type = neumann` / `x_max_u_type = neumann` on the
  declared outlet error-stop with "contradicts the declared patch type".
  PASS. Re-run 2026-09-26 at niter 12, CPU: PASS (wall twin max_abs 0, both
  contradictions error-stop).

Notes:

- Field files store only the interior staggered low faces, so the outlet
  face u(nx+1) is not in the h5; the divergence check covers cells with all
  six stored faces and the exactness checks pin the rest.
- KEY DESIGN LESSON (found by the first pois run): the outlet face needs the
  zero-gradient PREDICTOR write (`apply_bc(..., outflow_copy=.true.)` after
  the momentum predictor and at init/restart). With the face only ever
  touched by the projection correction, its SHAPE never forgets its IC
  (corrections are smooth phi gradients): the run converged, drift-free, to
  a spurious steady state with a plug outlet profile and 0.2 crossflow.
  Inside the projection loop the copy stays OFF — the Dirichlet-p correction
  owns the face there.
