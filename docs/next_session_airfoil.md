# Next session(s) — airfoil flow case: freestream in/outflow, lift/drag, quasi-2D RANS+transition

## STATUS 2026-07-12 — A0 DONE (all gates), A1/A2 implemented, cylinder gates running

- **A0 COMPLETE.** Face concept (`wall|patch|inlet|outlet`, `resolve_face_bcs`,
  internal `outflow` type, `_profile = parabola` value profiles, explicit-key
  contradiction = hard error) + the Dirichlet-pressure outlet in the
  projection (`face_grad_denom` 2*d1f / `face_grad_corr` d1f-vs-mirrored-ghost
  pair, phi SCALAR_BC_MIRROR after every per-iteration exchange, jacobi_apply
  high- AND low-face outlet corrections). Gates in `validation/freestream/`
  (all PASS, results in its README): standard suite bit-exact (nofma, max_abs
  0, CPU+GPU); oblique freestream EXACT; in/outflow Poiseuille = periodic
  reference to O(h^2), p linear + outlet-pinned, drift 5.6e-16/5k steps;
  Lamb-Oseen exits (reflected fraction 5e-3); 1==4 ranks EXACT, CPU==GPU;
  wall-twin exact + contradiction error-stops.
- **AMENDMENT to §THE-gap item 5 (load-bearing, found by gate c):** skipping
  the outflow normal write EVERYWHERE leaves the outlet face SHAPE frozen at
  its IC (the projection correction is a smooth phi gradient) — the run
  converges drift-free to a spurious plug-profile steady state with O(0.2)
  crossflow. The correct split: `apply_bc(..., outflow_copy=.true.)` writes
  the ZERO-GRADIENT copy at the predictor stage (post-momentum + init/
  restart, which also fixes the restart hole — the face is not in the h5),
  and the copy stays OFF inside the projection loop, whose Dirichlet-p
  correction owns the face. This is the classical fractional-step outlet.
- **ALSO:** explicitly-set `[boundary]` `_type`/`_value` ini rows now beat the
  restart file's stored rows (config-is-authority, the sor precedent; needed
  so the contradiction check sees the ini's values).
- **A1/A2 COMPLETE (2026-07-13):** `src/modules/flow/airfoil/airfoil_flow.f90`
  ([case.airfoil] aoa/u_inf/chord/force_sample_interval/runtime_file; patch-
  type freestream composition; uniform init; penalization force with
  rank-count-independent reduction: per-block sums -> global-id scatter ->
  exact allreduce -> ordered sum), `tools/make_airfoil_stl.py` (cylinder +
  NACA 4-digit, trimesh+shapely+mapbox_earcut in the ibmc venv),
  `validation/cylinder/` (setup + Re 40/100 + empty + CV cross-check;
  results in its README). Gates: empty domain C_L = C_D = 0.0 EXACTLY +
  aoa = 5 freestream exact; Re 40 steady C_D = 1.6924 +- 2.7e-5, |C_L| 4e-4
  (1.5-1.6 unbounded band + ~6% 16D-Dirichlet blockage + first-order
  penalization D_eff ~ D+h); Re 100 St = 0.168 (spectral fundamental; the
  confined C_L carries a comparable-power 3rd harmonic — zero-crossing
  St-counting is wrong by 3x), mean C_D 1.448, mean C_L 2e-4; forces 1 vs 4
  ranks BYTE-IDENTICAL, CPU vs GPU 0.0.
- **Gauss/CV cross-check PASS (gate e):** border-flux C_D on the clean-p
  steady snapshot = 1.6906 / 1.733 / 1.803 for boxes [4,9]x[6,10] /
  [3,11]x[4,12] / [2,14]x[2,14] vs penalization 1.6924 — 0.1% / 2.4% /
  6.5%, all within discretization error. The independent validation of
  both the force statistic and the in/outflow faces.
- **pn-snapshot caveat (found by the CV gate):** long IBM runs at production
  niter = 6 accumulate a large velocity-neutral oscillating mode in stored
  pn (std ~4e2 here; the channel pn-drift family). Forces/dynamics immune
  (u-only), but instantaneous pn is useless for border fluxes. Clean-p
  recipe: copy the converged restart, ZERO pn, restart with niter = 60 for
  ~300 steps (no kick; the converged projection rebuilds the physical p in
  a few substages, forces hold steady throughout). Do NOT restart the
  POLLUTED p at niter = 60 (violent transient, C_L ~ O(100), dt collapse —
  the spurious grad-p loses its self-consistent sloppy-projection
  compensation) and do NOT run whole cases at niter = 60 (~15x cost).

Branch `claude/jacobi-interface`. Goal: a second flow type `airfoil`
(`src/modules/flow/airfoil/`) for quasi-2D (few cells in z, periodic)
immersed-boundary airfoil runs, culminating in a RANS + γ–Re_θt transition
validation. User-specified peculiarities (2026-07-10):

- (i) an **angle of attack** config computes the Dirichlet inflow velocity
  (u = U∞ cos α, v = U∞ sin α, w = 0);
- (ii) the **runtime statistics are the lift and drag coefficients**,
  computed by integrating the mean momentum equation in the volume (for
  volume penalization this is the penalization-term integral — see A2);
- (iii) a **freestream boundary condition**: Dirichlet velocity + Neumann
  pressure on inlet-type faces, Neumann velocity + Dirichlet pressure on
  the outlet.

This is a multi-session feature. Phase it A0→A3, each gated, exactly like
the block-refinement and IDDES work. A0 is solver core (the projection);
everything after is case-level composition.

## What already exists (inventory — read these first)

- **Per-face, per-variable BC types** (`boundary.f90`): `faceBcType(var,
  face)` Dirichlet(0)/Neumann(1) for u, v, w AND p, with values, parsed
  from `[boundary] <dir>_<side>_<var>_<type|value>` (config.f90
  `apply_boundary_value`). `apply_bc` (boundary.f90:311) already handles:
  normal velocity ON the face (Dirichlet value / Neumann from interior),
  tangential + pressure ghosts (Dirichlet midpoint mirror / Neumann
  derivative). The **inlet of (iii) is expressible today**: velocity
  Dirichlet + pressure Neumann is the default type combination.
- **Patch types** (T5 STEP 0): `[boundary] <dir>_<side>_patch = wall |
  patch`; `domain_face_is_wall` reads it. Declaring the freestream faces
  `patch` stops RANS from reading the Dirichlet inlet as a no-slip wall
  (no ω pinning, no dwall min-in — gated in STEP 0).
- **Generic scalar BC applicator** `apply_scalar_bc` (boundary.f90:417)
  with `SCALAR_BC_VALUE` — the hook built for exactly the scalar-inlet
  values A3 needs. rans.f90 already wraps it with per-scalar mode tables.
- **Flow-case registry** (`flow_case.f90`): select on `[case] name`;
  `channel_flow.f90` + `channel_stats.f90` are the template (config
  section `[case.airfoil]`, `apply_defaults` = BC composition,
  `initialise_fields`, `after_step` = runtime sampling, `finalize`,
  `runtime_file` append pattern in channel_stats.f90:463).
- **IBM penalization** is IMPLICIT: `mu = 1/(1 + dt_γ·coef)`
  (ibm.f90:553); the predictor writes `qs := (q + dt·rhs − dt_γ∇p)·mu`
  (step.f90:175-179) and every later write to a solid face (projection
  correction, jacobi_apply) is also ×mu. Consequence: at end of substage
  `q_final = mu·q_unmasked` exactly, so `coef·q_final =
  (q_unmasked − q_final)/dt_γ` — the penalization force density is exact
  in the discrete bookkeeping. A2 builds on this identity.
- **The T4/T5 inlet gaps, already itemized** (docs/next_session_iddes.md
  T4 entry): patch classification exists; the RANS scalars lack inlet
  ghost VALUES; first-order upwind was kept with channel-front evidence
  and the TVD revisit was explicitly deferred to "the flat-plate/inlet
  increment" — A3 is that increment's airfoil form.

## The face concept: ONE patch type drives everything (simplification,
## user request 2026-07-12)

The first draft of this plan scattered one physical fact — "this face is
an outlet" — across four independently-set axes: three velocity
`faceBcType` rows, the pressure row, the STEP-0 patch type, and the RANS
scalar mode tables, all of which the ini (or the case) had to keep
mutually consistent by hand. Simplify: the PATCH TYPE becomes the single
user-facing face concept, and everything else is DERIVED from it.

- `facePatchType` gains two values: `wall | patch | inlet | outlet`
  (config key unchanged: `[boundary] <dir>_<side>_patch`).
- One new `resolve_face_bcs(bc)` in boundary.f90 (runs where
  validate_patch_types runs — after parsing, when periodic_* is final,
  BEFORE the first update_boundary_values) seeds the per-variable rows of
  any DECLARED face, set-if-unset:
  - `wall`:   u,v,w Dirichlet (default 0), p Neumann — the classic no-slip
              composition, now writable as one key;
  - `inlet`:  u,v,w Dirichlet (values from `_value` keys or the case),
              p Neumann;
  - `outlet`: normal velocity → the INTERNAL `outflow` type (apply_bc
              skips the normal-component write; the projection owns the
              face), tangential Neumann 0, p Dirichlet 0;
  - `patch` / UNSET: exactly today's behaviour (raw per-variable keys,
    the tangential-Dirichlet wall inference) — existing inis bit-exact
    by construction.
- Explicit `<var>_type` keys still win, but a DIRECT contradiction with a
  declared patch type (e.g. `x_max_p_type = neumann` on an `outlet`) is a
  hard config error, not a silent override — the consistency constraints
  are enforced by construction instead of by ini discipline.
- The `outflow` velocity type (read_bc_type value 2) becomes an INTERNAL
  representation only — no ini ever writes it; declaring the face
  `outlet` is the interface.
- Consumers become pure functions of the patch type: the projection's
  outlet flags (below) = `facePatchType == PATCH_OUTLET`;
  `domain_face_is_wall` gains `inlet/outlet → .false.` (one line each);
  the A3 RANS scalar modes key off PATCH_INLET/PATCH_OUTLET instead of
  re-inferring "velocity-Dirichlet non-wall" from the type rows.

## THE gap: no Dirichlet-pressure outlet in the projection (A0's job)

`pressure_solver.f90` treats EVERY `FACE_PHYS` face as pinned/Neumann:
`face_grad` (line 366) returns 0 for FACE_PHYS in both the Jacobi
denominator (jacobi_compute_phi:239) and the velocity-face correction
(jacobi_apply:318). So a physical face's normal velocity is never
corrected — right for walls and for the Dirichlet inlet (prescribed flux),
fatal for an outlet: with net inflow prescribed and every face pinned, the
all-Neumann Poisson problem is INCOMPATIBLE (∫div dV ≠ 0 can never be
projected out; the residual stagnates and pressure drifts). The outlet of
(iii) — Neumann velocity, Dirichlet pressure — is what restores
solvability: the projection must own the outlet face flux.

Design (keep the denominator/correction pair consistent — the interface
work showed this is where instability hides):

1. **Outlet flag into the kernels.** Derive at init per-face flags from
   `facePatchType == PATCH_OUTLET`, pass them as small
   `pdirLow(3)/pdirHigh(3)` logicals into
   `jacobi_compute_phi`/`jacobi_apply` (map(to:), guarded by the block's
   own `physLow/High == FACE_PHYS` — the flags are global per domain
   face; a block only acts when its face IS that domain face).
2. **phi ghost = MIRROR at the outlet, every iteration.** After the
   per-iteration `exchange_scalar_halos(phi, ...)` (pressure_solver.f90:
   129) call `apply_scalar_bc(blk, bc, phi, mode)` with SCALAR_BC_MIRROR
   on outlet faces, NONE elsewhere — phi(face) = 0 is the increment of a
   HELD outlet pressure. Reuses STEP-0 machinery verbatim; do NOT rely on
   the phys ghost staying zero (the exchange's tangential extension can
   write physical halos).
3. **Metrics.** With the mirrored ghost the face-normal correction uses
   the REGULAR metric: `q(nx+1) += (phi(nx) − phi_ghost)·d1·mu =
   2·phi(nx)·d1·mu` — the half-cell Dirichlet gradient. The Jacobi
   DENOMINATOR must count that same sensitivity: the outlet face
   contributes `2·d1f·mu` (the ghost is −phi_i, so ∂flux/∂phi_i doubles).
   face_grad gets the new branch (or a thin wrapper) — the denominator
   caller and the correction caller need DIFFERENT values at this one
   face kind (2·d1f vs d1f-with-mirrored-ghost); keep both explicit and
   commented.
4. **jacobi_apply high-face branch.** The outlet x_max face `u(nx+1)` is
   a high face — today corrected only when 2:1-interface-owned
   (jacobi_apply:329). Extend the guard: also fire when the face is
   FACE_PHYS ∧ outlet. (The momentum predictor never writes high faces,
   so the outlet face evolves ONLY through this correction + its initial
   value — the standard "do-nothing" pressure outlet.)
5. **apply_bc must not stomp the corrected face.** apply_bc runs inside
   the projection loop (pressure_solver.f90:131) and would overwrite
   `u(nx+1)` with the Neumann copy each iteration. The internal `outflow`
   normal-velocity type (derived by resolve_face_bcs, see the face-concept
   section) makes apply_bc SKIP the normal-component write (tangential
   ghosts keep Neumann copy; pressure ghost keeps its Dirichlet mirror).
   The outlet face is then owned by the projection, which is the only
   consistent owner: you cannot impose zero-gradient velocity AND
   Dirichlet pressure on the same corrected field — the pressure
   condition wins; zero-gradient is only the predictor-stage stance.
6. **Solvability bonus:** one Dirichlet-p face makes the operator
   nonsingular — the known pressure null-mode drift (memory:
   pressure-volume-average-drift) cannot occur in this case. Chebyshev
   bounds: lmax = 2 (Gershgorin) still holds; lmin auto-derivation is
   unchanged (the pinned mode only helps).

Bit-exactness argument for A0: no existing ini declares a Dirichlet
pressure face, so every new branch is dormant → the standard suite must be
bit-exact BY CONSTRUCTION (still gate it, nofma CPU+GPU).

## A1 — the airfoil flow case (`src/modules/flow/airfoil/`)

`airfoil_flow.f90` registered in flow_case.f90 (`[case] name = airfoil`).
Host-side orchestration only (the type-bound flow-case pattern is fine
here; the no-polymorphism rule is for device code).

- Config `[case.airfoil]`: `aoa` (degrees), `u_inf` (default 1), `chord`
  (default 1), `force_sample_interval`, `runtime_file` (default
  forces.txt), optional `outlet_face = x_max` (default).
- `apply_defaults` composes the freestream (iii) from (i) in ONE
  vocabulary — patch types: x_min, y_min, y_max declared `inlet` with
  velocity values (U∞cosα, U∞sinα, 0); x_max declared `outlet`. That is
  the whole composition — resolve_face_bcs derives the per-variable rows
  (the face-concept section above). Explicit user keys in [boundary] win
  (set-if-unset, the channel apply_defaults idiom). z stays periodic
  (quasi-2D). Dirichlet freestream on the y faces is the standard
  penalization far-field; document that large |α| needs a taller domain.
- `initialise_fields`: uniform freestream everywhere (the IBM damps the
  interior of the body in the first steps; impulsive start). Optional
  small noise amplitude for the turbulent cases (block-origin-indexed
  like the channel fix — see the T2 initialise_channel_fields lesson).
- Geometry comes from the EXISTING file-based IBM path (mobygeom
  block-table + coef_blocks + dwall_blocks): nothing case-specific in the
  solver. A small `tools/` helper generates extruded STL from an airfoil
  coordinate file / NACA 4-digit formula (trimesh extrude; float32 STL
  vertices are part of the as-built geometry — T1 lesson) and a cylinder
  STL for the A2 gates.
- Quasi-2D: nz = nb (one block deep, e.g. 8), periodic z; refine_body
  refines z too (octree) — fine, the fine region just doubles nz locally.

## A2 — lift/drag runtime statistics (the momentum-balance force)

Two equivalent readings of "integrate the mean momentum equation in the
volume" (user question 2026-07-12 settled here). Summing the DISCRETE
momentum equation over all fluid cells telescopes every interior flux,
leaving exactly

  d/dt ∫u dV = Φ_border − F_pen,   Φ_border = ∮_outer (−uu·n − p n
               + (ν+ν_t)(∇u+∇uᵀ)·n) dS  (discrete face fluxes),

so the Gauss/control-volume form (outer boundary = domain border, inner
boundary = the airfoil, whose flux resultant is the sought force) and the
penalization integral are IDENTICAL in the time mean — the choice is
about conditioning and what is available instantaneously.

**Turbulent Reynolds stresses are automatically inside the penalization
integral.** Every momentum write — predictor, projection correction, the
eddy-viscosity SGS/RANS correction (step.f90:457, `dt_alpha*sgs_u*mu_u`),
the body force — is ×mu, so `coef·u_final` is the TOTAL momentum sink of
the model as solved: pressure, viscous, modeled turbulent stress and
resolved fluctuations all included, nothing to assemble. In the
outer-border form they must appear EXPLICITLY: the modeled deviatoric
ν_t stress on the border AND, when time-averaging an unsteady flow, the
resolved fluctuation flux −⟨u'u'⟩·n; and the border pressure is the
MODIFIED pressure (p absorbs (2/3)k in RANS runs — consistent as long as
the same p is used everywhere, but one more thing to get right.)

**Decision: the penalization integral stays the RUNTIME statistic**
(instantaneous C_L(t)/C_D(t) — the shedding-St gates need it; localized
at the body, no cancellation — border drag is a small difference of
large in/out momentum fluxes; no turbulence terms to assemble). **The
Gauss/CV border-flux balance is ADOPTED as the independent consistency
check**: (a) an A2 GATE on the steady cylinder (Re = 40: a single
converged snapshot, no averaging, no fluctuation terms — border fluxes
must reproduce the penalization C_D to discretization error; this also
independently validates the new in/outflow faces, since the inlet/outlet
fluxes dominate Φ_border); (b) an optional python postprocessing
diagnostic (`tools/` or the case directory) that decomposes the MEAN
force into border momentum / pressure / viscous / turbulent (modeled +
resolved) contributions from accumulated mean fields — physical insight,
not a runtime cost.

Mechanics of the runtime statistic: **F_on_body = ∫ coef·u dV**,
component-wise on the staggered grid; by the mu-masking identity above
this is EXACT discrete bookkeeping evaluated on the end-of-step field
(coef·u_final = momentum removed / dt_γ). No surface reconstruction, no
quadrature on the body.

- `after_step`, every `force_sample_interval` steps: one device reduction
  per component `F_d = Σ_blocks Σ_{i,j,k interior} coef(i,j,k,d,b) ·
  q(i,j,k,d,b) · V_d(i,j,k,b)` with the staggered cell volume from the
  d1 metrics (1/d1x(i,var)·1/d1y(j,var)·1/d1z(k,var)); `comm_allreduce_sum`
  across ranks. Interior staggered low faces (1..nb) count every face
  exactly once globally (the block-redundant layer stores the neighbour's
  face in the halo, not in 1..nb; domain-boundary high faces carry
  coef = 0 as long as the body does not touch the domain boundary).
- Coefficients: drag/lift unit vectors (cos α, sin α, 0) and
  (−sin α, cos α, 0); `C_D = 2 F·e_D/(ρ U∞² c L_z)`, `C_L` likewise
  (ρ = 1). Write `step, t, C_L, C_D` to `runtime_file` (channel_stats
  append pattern, channel_stats.f90:463). Moment coefficient is a later
  trivial extension (lever-arm weight in the same reduction) — out of
  scope now.
- Note in the case README: penalization force is TOTAL (pressure +
  friction combined); the split needs surface integration — not planned.

## A3 — RANS(+transition) inlet-awareness + the airfoil physics case

- **Scalar inlet values** (the T4 gap, hook ready): rans.f90's mode
  tables gain SCALAR_BC_VALUE at `PATCH_INLET` faces with the freestream
  values already computed from `[rans] tu` / `nut_ratio` at init
  (k∞ = 1.5(tu/100·U∞)², ω∞ from nut_ratio; γ∞ = 1; R̃e_θt,∞ = the tu
  correlation value T4 uses for the IC), and SCALAR_BC_COPY at
  `PATCH_OUTLET` — pure functions of the patch type, no re-inference
  from the velocity BC rows. nut wall-ghost handling is untouched.
- **Upwind stance**: KEEP first-order upwind initially, but the channel
  front-alignment argument does NOT carry over (an airfoil transition
  front is crossed by chordwise convection). Measure the smearing
  (γ-front width in cells on the SD7003 case) and report; the TVD/van
  Leer + second-scalar-halo increment stays separate and is triggered by
  that measurement, not assumed.
- **Physics sequence** (each its own run):
  1. Full-turbulent SST first (transition = false), moderate Re, vs
     XFOIL/literature C_L(α) slope and C_D magnitude — sanity, not a
     tight gate.
  2. `transition = true` on the **SD7003 at Re = 6e4, α = 4°** — the
     standard low-Re LSB/transition benchmark (Galbraith & Visbal LES,
     Windte et al. RANS-LM): gate on bubble presence, transition
     location ballpark (x_t/c ~ 0.5 ± 0.1 at tu ≈ 0.1%… use published
     γ–Re_θt results, not LES, as the bar), C_L/C_D within the scatter
     of published γ–Re_θt implementations (~10-15%).
  3. Resolved walls only (transition ∧ wall_function stays a hard
     error); y+₁ ≲ 1 at the surface via refine_body multi-level.
- **Grid reality check** (do this arithmetic before running): resolved
  y+ ≲ 1 at Re = 6e4 needs Δ ≈ 1-2e-3 c at the wall; from a far-field
  base grid that is 3-4 refine_body levels. Levels > 1 are lightly
  exercised (les_ibm used l1) — treat multi-level refine_body + mobygeom
  block-table as its OWN gate (uniform-flow preservation + dwall
  cross-check per level) before the physics run.

## Phase gates

- **A0** (solver core, everything else dormant):
  (a) standard suite bit-exact (nofma, max_abs 0, CPU AND GPU):
  min_channel, les_ibm ± refine_body, Beltrami y-slab 5 steps, turb180,
  wf180_y30, lam30t — the new branches are unreachable without a
  Dirichlet-p face;
  (b) uniform oblique flow (any α) through an EMPTY box with
  inlet/outlet: preserved EXACTLY (machine zero — constants are in the
  null space of every consistent operator), global mass residual
  round-off, projection divergence round-off;
  (c) inflow/outflow plane Poiseuille (Dirichlet parabola in, outflow
  out) matches the periodic-forcing reference profile
  (validation/poiseuille machinery) to discretization error, pressure
  gradient linear, level pinned at the outlet (no drift over 10k steps);
  (d) a Lamb–Oseen vortex advected through the outlet: no blow-up,
  bounded reflection (report the reflected fraction; a convective outlet
  is the documented fallback increment if it is ugly);
  (e) 1 == 4 ranks EXACT with in/outflow; CPU vs GPU at the usual level;
  (f) config validation: a patch declaration derives the expected rows
  (assert via a wall-declared twin == the per-variable-key twin,
  bit-exact) and a contradicting explicit key (`p_type = neumann` on an
  `outlet`) error-stops.
- **A1/A2** (case + forces): cylinder STL, quasi-2D:
  (a) Re = 40 steady: C_D ≈ 1.5–1.6 (literature band for penalization
  codes), C_L → 0;
  (b) Re = 100: vortex shedding, St ≈ 0.16–0.17, mean C_D ≈ 1.3–1.4,
  C_L oscillation symmetric about 0;
  (c) α ≠ 0 empty-domain freestream: C_L = C_D = 0 exactly (coef = 0
  everywhere);
  (d) force reduction 1 == 4 ranks exact, CPU == GPU;
  (e) the Gauss/CV cross-check: on the converged Re = 40 snapshot the
  outer-border flux balance (momentum + pressure + viscous; steady, so
  no fluctuation terms) reproduces the penalization C_D to
  discretization error — the independent validation of BOTH the force
  statistic and the in/outflow faces.
- **A3** (RANS + transition): the sequence above + the standard
  regression list stays bit-exact (RANS scalar BC changes are
  inlet-face-gated, dormant in channels).

## Watch for

- The denominator/correction consistency at the outlet (2·d1f vs
  mirrored-ghost d1f) — get it wrong and the projection loses SPD/
  convergence at exactly one face row, which reads as a slow mysterious
  drift, not a bang.
- apply_bc ordering inside the projection loop: the `outflow` type must
  skip ONLY the normal-component write; tangential and pressure ghosts
  still need their update each iteration.
- Do not touch `face_grad`'s existing branches — the 2:1 interface
  metrics are validated and locked; add the outlet branch strictly
  additively.
- The force identity needs the END-OF-STEP field (after the last
  exchange); sampling mid-substage breaks the mu bookkeeping.
- resolve_face_bcs ordering: after config parsing AND after periodic_*
  is final (the validate_patch_types slot), BEFORE the first
  update_boundary_values — the derived rows must seed pointBcValue.
- The CV cross-check uses the MODIFIED pressure consistently (RANS p
  absorbs (2/3)k) and, on unsteady averages, needs the resolved
  −⟨u'u'⟩·n flux as well as the modeled ν_t stress — the two things the
  penalization form gets for free (why it stays the runtime statistic).
- Freestream y-faces with Dirichlet v ≠ 0 inject mass through y — the
  outlet correction absorbs it; the empty-box gate (b) catches sign
  errors immediately.
- mobygeom STL: float32 vertices, length-periodic in z (extrude the full
  Lz and let the periodic images handle it — same assumption the
  coefficient machinery already makes).
- RANS at the inlet: the patch declaration (NOT the inference) is what
  keeps ω unpinned there — the case's apply_defaults must set it, and
  `report_patch_types` should show every freestream face as `patch
  (declared)`.
- Quasi-2D ≠ 2D: the solver is 3D; RANS solutions stay z-uniform, but
  a URANS shedding case is fine too. Keep nz = nb.
- `[flow] forcing_*` must stay 0 in airfoil inis (no volume forcing);
  Re in [flow] is the chord Reynolds number (U∞ = 1, c = 1 canonical).
- git: stage explicit paths; STL/h5/png outputs stay untracked.

## NEXT-SESSION PROMPT — A0 (+A1/A2 if it goes fast)

> Read `docs/next_session_airfoil.md` and CLAUDE.md. Branch
> `claude/jacobi-interface`. Implement phase A0 in two increments, gated
> separately. INCREMENT 1 — the face concept: extend
> `[boundary] <dir>_<side>_patch` to `wall | patch | inlet | outlet`;
> ONE `resolve_face_bcs` in boundary.f90 (the validate_patch_types slot,
> before the first update_boundary_values) derives the per-variable BC
> rows set-if-unset (wall: velocity Dirichlet 0 + p Neumann; inlet:
> velocity Dirichlet from the `_value` keys + p Neumann; outlet: normal
> velocity = the INTERNAL `outflow` type, tangential Neumann 0,
> p Dirichlet 0); explicit keys win but a direct contradiction with the
> declared patch type is a hard config error; `outflow` never appears in
> an ini. Consumers key off the patch type: domain_face_is_wall gains
> inlet/outlet → false; the projection flags and the later RANS scalar
> modes read PATCH_INLET/PATCH_OUTLET. INCREMENT 2 — the
> Dirichlet-PRESSURE outlet in the projection (pressure_solver.f90),
> the one real solver piece: PATCH_OUTLET flags into the kernels; phi
> ghost MIRRORED at outlet faces via apply_scalar_bc after each
> per-iteration phi exchange (do not rely on phys ghosts staying zero);
> the Jacobi DENOMINATOR counts 2·d1f·mu at the outlet face while the
> CORRECTION uses d1f against the mirrored ghost — keep the pair
> consistent and commented, SPD lives there; jacobi_apply's high-face
> branch extended to correct the outlet face (FACE_PHYS ∧ outlet),
> mirroring the 2:1 owned-face precedent; apply_bc skips the normal
> component of `outflow` faces only. All new branches are dormant
> without a declared inlet/outlet: gate (a) bit-exact (nofma, max_abs 0,
> CPU AND GPU) on min_channel / les_ibm ± refine_body / Beltrami y-slab
> (5 steps) / turb180 / wf180_y30 / lam30t. Physics gates in a new
> validation/freestream/: (b) uniform oblique flow through an empty box
> with inlet/outlet faces preserved exactly, mass + divergence
> round-off; (c) inflow/outflow Poiseuille vs the periodic reference,
> outlet-pinned pressure level, 10k-step drift-free; (d) Lamb–Oseen
> vortex exits, report the reflected fraction; (e) 1==4 ranks EXACT,
> CPU vs GPU usual level; (f) a wall-declared face == the equivalent
> per-variable-key ini bit-exact, and a contradicting explicit key
> error-stops. THEN, if A0 gates clean: A1 (the airfoil case module:
> [case.airfoil] aoa/u_inf/chord/force_sample_interval/runtime_file;
> apply_defaults = x_min/y_min/y_max inlet with (U∞cosα, U∞sinα, 0) +
> x_max outlet, set-if-unset; uniform-freestream init; quasi-2D nz = nb,
> periodic z) and A2 (the penalization-integral force F = ∫coef·u dV per
> component — the EXACT total force incl. the mu-masked eddy-viscosity
> stress, see the A2 section; device reduction + allreduce in
> after_step, C_L/C_D to the runtime file; gates: cylinder Re = 40
> C_D ≈ 1.5, Re = 100 St ≈ 0.165, empty domain C_L = C_D = 0 exactly,
> ranks/GPU exact, AND the Gauss/CV cross-check — outer-border momentum
> + pressure + viscous fluxes on the converged Re = 40 snapshot
> reproduce the penalization C_D to discretization error). Write the
> cylinder/airfoil STL generator in tools/ (trimesh, extruded, float32).
> WORKFLOW: one solver job at a time; sed-shortened inis for
> bit-exactness (memory: bit-exact-gates-short). Deferred, do NOT start:
> A3 (RANS scalar inlet values + SD7003 transition validation — its own
> session), the TVD upwind revisit (triggered by the measured γ-front
> smearing, not assumed), convective outlet (only if gate (d) demands
> it), the mean-force border-decomposition diagnostic (post-A3 nicety),
> moment coefficient, the GPU profiling task
> (docs/next_session_profiling.md).
