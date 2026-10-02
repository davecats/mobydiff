# A0 freestream gates — inlet/outlet patch types + Dirichlet-pressure outlet

Validation for phase A0 of the airfoil plan (`docs/next_session_airfoil.md`):
the `[boundary] <dir>_<side>_patch = wall | patch | inlet | outlet` face
concept (`resolve_face_bcs` derives the per-variable BC rows) and the
Dirichlet-pressure outlet in the projection (outlet face counted `2*d1f*mu`
in the Jacobi denominator, corrected with `d1f` against the MIRRORED phi
ghost). Since 2026-10-01 the outlet face is advanced by the momentum
predictor itself (`step.f90 predict_outlet_faces`); the zero-gradient copy
by `apply_bc` that the notes below describe is gone. **The current numbers
are in the last section.**

Run everything (one solver job at a time) with a built `build_cpu`/`build_gpu`:

    ./run_gates.sh              # or one group: oblique | pois | vortex | ranks | config
                                #             | low | refined | restart

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

## 2026-10-01: what the outlet does to these gates, and a prototype

`../cylinder/README.md` (last section) runs down the outlet pressure mode:
the predictor's zero-gradient copy RESETS the outlet face every substage, the
projection re-supplies the continuity correction each time, and the
incremental stored pressure integrates it in the last cell column. Two
numbers of this directory are that mechanism: the Lamb-Oseen exit (energy
regrowth 0.049 → 0.228 during the exit, reflected fraction 2.2e-2) and the
Poiseuille last-cell pressure (2.41e-3 where the held outlet pressure gives
G dx/2 = 1.25e-3). With the prototype of that section (branch
`proto/outflow-incremental`, not on main), same commands, CPU, niter 12:

| gate | main (`~/step9b_ref_binaries`) | prototype |
|---|---|---|
| oblique | exact (0.0) | exact (0.0) |
| pois: profile dev at x/lx 0.5 / 0.9 | 1.639e-3 / 1.539e-3 | 1.637e-3 / 1.521e-3 |
| pois: slope / nonlinearity | −0.079579 / 9.54e-5 | −0.079572 / 9.37e-5 |
| pois: last-cell p | 2.412e-3 | **1.247e-3** |
| pois: drift | 1.3e-15 | 5.6e-16 |
| vortex: E/E0 at t = 0.5 … 2.5 | 0.856 1.083 1.050 0.750 0.346 0.049 0.228 0.059 0.022 | 0.958 0.889 0.749 0.502 0.240 0.077 0.0079 0.0007 0.0001 |
| vortex: reflected fraction | 2.23e-2 | **9.8e-5** |

## 2026-10-01: the outlet face is predicted (`docs/next_session_outlet.md`)

The prototype's arithmetic is now the scheme, owned by the predictor: the
outlet face advances by its neighbour's non-pressure increment plus its own
pressure gradient against the held outlet pressure, on a low or a high face
alike. The reset row and its flag are gone. All of the below with the final
tree, `./run_gates.sh`, CPU (GPU where the driver uses it), niter 12.

| gate | result |
|---|---|
| oblique | exact (0.0), interior divergence 0.0 |
| pois: profile dev at x/lx 0.5 / 0.9 | 1.637e-3 / 1.521e-3 |
| pois: slope / nonlinearity | −0.079572 / 9.37e-5 |
| pois: last-cell p (exact G dx/2 = 1.25e-3) | **1.247e-3** |
| pois: drift | 5.6e-16 |
| vortex: E/E0 at t = 0.5 … 2.5 | 0.958 0.889 0.749 0.502 0.240 0.077 0.0079 0.0007 0.0001 |
| vortex: reflected fraction | **9.80e-5**, monotone |
| ranks | oblique and pois_io 1 == 4 EXACT; CPU vs GPU 0.0 / 5.6e-16 |
| config | wall twin exact; both contradictions error-stop |

Three new groups:

- **low** — the outlet on a LOW face, which had never run before.
  `outlet_high.ini` / `outlet_low.ini` are a mirrored pair with a
  cross-flow inlet, so the outflow is not parallel (the two schemes differ
  by 0.1 on it): high against mirrored low **3.4e-15** velocity, 3.1e-14
  pressure. The Poiseuille and the vortex cases run in −x
  (`make_freestream_ics.py --mirror`, the checkers' `--mirror`) print the
  SAME numbers as above to every digit, and their fields equal the x_max
  runs to 5.6e-16 / 7.7e-15 (Poiseuille) and 4.7e-15 / 8.1e-15 (vortex).
- **refined** — uniform oblique flow recovered FROM REST (`initial_u =
  initial_v = 0`) through a level-1 patch, 2000 steps, `dtmax 2.5e-3`.
  The exactness gates start from the uniform field and cannot see a halo
  that nothing writes; this one can.

  | patch | max deviation from uniform | stored p |
  |---|---|---|
  | touching nothing | 1.8e-13 | ≤ 9e-14 |
  | touching the x_max outlet | **2.9e-13** | ≤ 8e-14 |
  | touching the x_min outlet (mirrored) | 2.2e-13 | ≤ 8e-14 |

  On main (`154e48f`) the second leg froze at 2.28e-2 with the stored
  pressure growing linearly (the level-jump finding of the handout:
  cross-level exchange entries did not reach index `nb+1` on a physical
  face; fixed in `comm.f90 candidate_boxes`). 1 == 4 ranks EXACT on it.
  **Mind the time step in a variant of these legs**: `pecletmax = 0.5` is
  beyond RK3's limit on an isotropic fine level (0.31 in 2D, 0.21 in 3D);
  at `dtmax 1e-2` the patch sits in a sustained O(0.65) limit cycle.
- **restart** (`run_restart_outlet.sh`) — a run stopped and restarted at
  mid-length equals the continuous one at `max_abs 0`, five cases
  (`outlet_high`, `outlet_low`, `pois_io`, `lamboseen`, `../blasius`),
  CPU and GPU. The outlet face is state: on a low face it is index 1 of
  the velocity dataset, on a high face the snapshot carries it as
  `un_xmax` / `vn_ymax` / `wn_zmax`. A file without the plane (an older
  snapshot) is announced and falls back to the zero-gradient value; the
  run then differs by 1.9e-4 / 1.5e-7 / 4.4e-5 / 7.3e-6 in velocity on the
  four high-side cases after the second half.

**An initial-condition tool must supply the plane.** A template minted by
the solver already carries `un_xmax` from its own one-step run; a generator
that copies the template and overwrites only `un/vn/wn/pn` restarts the
outlet face from that unrelated field. `make_freestream_ics.py` and
`../blasius/make_blasius_ic.py` write the analytic plane; a tool that
cannot should delete the dataset.

Bit-exactness on the outlet cases, which the 7- and 9-case suites lack:
`run_bitexact_outlet.sh` (nofma binaries, absolute paths). Used for every
increment: the refactor of the final copy and the removal of the
post-predictor `apply_bc` at `max_abs 0` against the tree before them, and
the predictor row at `max_abs 0` against the PROTOTYPE binary on all seven
cases, CPU (1 and 4 ranks) and GPU. A reference binary older than the
outlet planes ignores them: compare against one only from ICs without the
plane.
