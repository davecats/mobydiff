# Next session(s) — IDDES: k-ω SST(+transition) RANS, WALE LES, hybrid blend

Branch `claude/jacobi-interface`. Goal: extend the solver with **RANS
(k-ω SST, optionally the 4-equation γ–Re_θt transition variant)** and an
**IDDES hybrid** that blends the existing WALE LES (fine regions) with SST
RANS (coarse/near-wall regions) — with ONE architecture serving all of:
no model, LES alone, RANS alone, IDDES. This is a large multi-session
feature; it is phased, each phase gated, exactly like the block-refinement
work. Terminology: "MAC" in this doc means the staggered marker-and-cell
cell arrangement (u/v/w on faces, scalars at cell centres) — NOT any
particular geometry. Validation runs on the existing case fleet (plane
channel, the `les_ibm` off-grid IBM channel, wavychannel); no new
geometries are required.

## The architectural key: one coupling variable

The solver already has a clean producer/consumer split around the eddy
viscosity, and EVERYTHING downstream of `nut` is model-agnostic and
validated (including across the 2:1 interface):

- consumer chain (DO NOT TOUCH — this is what keeps LES/no-model bit-exact):
  cell-centred `nut` → `exchange_scalar_halos(c, nut, blk)` (comm.f90:1115,
  2:1 restrict/prolong included) → `add_les_momentum_correction`
  (step.f90:358) → `update_timestep_limits` (step.f90:655, reads nut for the
  diffusive dt limit) → `fdm_h5_append_nut` (snapshot output).
- producer (model-specific): today `update_les_viscosity` (les.f90) —
  WALE / Smagorinsky algebraic kernels.

RANS and IDDES are just **different producers of the same cell-centred
`nut`**. The momentum coupling, halo exchange, GPU dt limit and io cost
zero new code, and the LES-only / no-model paths stay bit-exact by
construction because the consumer chain is moved, never edited.

Physics note (document in the code): the deviatoric Boussinesq stress is
what the nut correction applies; the `−(2/3)k δij` part is absorbed into
pressure, so in RANS/IDDES runs the output `p` is a modified pressure.
No code, but say so in a comment and in the config docs.

## Module layout (composition, no duplication)

```
src/modules/turbulence.f90   NEW  turb_type: model enum, nut, nut_sgs scratch,
                                  fd field, the shared metric tables (hoisted
                                  from les_type), orchestration + IDDES blend.
src/modules/les.f90          KEPT algebraic SGS viscosity only: q -> nut target
                                  array (WALE/Smag kernels + the shared
                                  velocity_gradient_tensor device helper).
src/modules/rans.f90         NEW  sst_type: k, omega (+ gamma, rethetat when
                                  transition on), their oldrhs, dwall/yeff,
                                  wall-cell marker; transport kernels,
                                  wall treatment, nut_rans assembly.
```

Rules that make it sleek:
- `les.f90` and `rans.f90` never reference each other. The IDDES blend
  lives ENTIRELY in `turbulence.f90` and touches exactly two places:
  (1) the k-destruction length `l_hyb` passed into the RANS kernel,
  (2) the final `nut = fd*nut_rans + (1-fd)*nut_sgs` assembly.
- Dispatch is a `select case` on an enum at kernel-launch level (the
  existing `les%model` idiom). NO Fortran polymorphism / type-bound
  procedures — they do not survive OpenMP target offload and fight the
  flat-arrays-mapped-once convention (CLAUDE.md).
- Derived types own flat contiguous allocatables, device-mapped once in
  `enter_*_data` / `exit_*_data` (mirror `bodyforce.f90`, the cleanest
  recent example).

Per-substage orchestration (in turbulence.f90, called from main.f90 where
`update_les_viscosity` + the nut exchange sit today, main.f90:172-179):

```fortran
select case (turb%model)
case (TURB_NONE);  return
case (TURB_LES);   call update_sgs_viscosity(les, ..., turb%nut)   ! today's kernel verbatim
case (TURB_RANS, TURB_IDDES)
    if (iddes) call update_sgs_viscosity(les, ..., turb%nut_sgs)   ! scratch target
    if (iddes) call compute_fd(rans, turb, ...)                    ! else fd == 1, l_hyb = l_RANS
    call rans_apply_wall_treatment(rans, ...)   ! pin wall-cell omega/nut BEFORE transport
    call rans_advance_scalars(rans, ..., dt_alpha, dt_beta, dt_gamma)  ! fused kernel
    call rans_scalar_halos_and_bcs(rans, c, blk, bc)
    call assemble_nut(turb, rans, ...)          ! nut_rans; blended with nut_sgs when iddes
end select
call exchange_scalar_halos(c, turb%nut, blk)                       ! unchanged
```

## What gets REUSED (the no-duplication inventory)

- `velocity_gradient_tensor` (les.f90:248, `!$omp declare target`): the SST
  production needs S and the γ-equation needs Ω — contractions of the same
  tensor WALE computes. Reuse the helper as-is (it was already factored for
  exactly this kind of sharing).
- The metric tables in `les_type` (`d1?m/0/p`, `p_from_*`, `filter_*`,
  `inv_d*`, les.f90:27-33): grid metrics, not LES physics. Hoist into
  `turb_type`; both producers read them. `filter_*` gives Δ_mesh for the
  IDDES length scale for free.
- `exchange_scalar_halos`: k, ω, γ, Re_θt are cell-centred scalars — same
  layout as nut/phi, so the validated 2:1 restrict/prolong path applies
  verbatim. Extra calls per substage, no new comm code. (If profiling later
  shows the 4 separate exchanges hurt, batching multiple scalars into one
  message is a Phase-4-style optimisation — do NOT design for it now.)
- The RK3 low-storage idiom: give the scalars their own `scal/oldrhs`
  pair and the same `dt_alpha/beta/gamma` treatment as the momentum
  predictor (step.f90 `momentum`). No new time-integration concept.
- The `ibm_aware` solid-cell test (staggered `ibm%coef > threshold`,
  les.f90:441-449): reused for (a) zeroing nut in solid (already done),
  (b) classifying IBM WALL CELLS (fluid cell with ≥1 solid staggered
  face), (c) masking scalar diffusive fluxes through solid faces.
- The correction-kernel pattern (LES SGS / bodyforce): `TURB_NONE`
  allocates nothing, touches nothing → disabled is bit-exact BY
  CONSTRUCTION, not by a `+0.0` argument.
- io: `fdm_h5_append_nut` unchanged; add the RANS scalars to snapshots the
  same way (cell-centred appends). Restart of the scalars rides the same
  block-table read/write machinery as nut/velocity.

## Formulation (what to implement)

Incompressible, ρ=1 throughout (all equations in kinematic form).

### k-ω SST (Menter 2003 constants)
- ν_t = a1 k / max(a1 ω, S F2), a1 = 0.31.
- k-eq: ∂k/∂t + conv = P̃_k − D_k + ∇·[(ν + σ_k ν_t)∇k]
  with P_k = min(ν_t S², 10 β* k ω); **D_k = k^{3/2}/l_hyb** is where
  IDDES enters (pure RANS: l_hyb = l_RANS = √k/(β* ω), giving the standard
  β* k ω).
- ω-eq unchanged from standard SST: α S² (α = γ-blend of 5/9, 0.44)
  − β ω² + diffusion + (1−F1)·CD_kω cross-diffusion.
- Set-1/2 constants blended by F1: σ_k = (0.85, 1.0), σ_ω = (0.5, 0.856),
  β = (0.075, 0.0828), β* = 0.09, κ = 0.41.
- F1, F2: standard Menter definitions on (k, ω, ν, dwall). Use **y_eff**
  (below) as the wall distance everywhere.

### γ–Re_θt transition variant (Langtry & Menter 2009, = OpenFOAM kOmegaSSTLM)
Two extra transported scalars per the user-provided spec in the session
notes: γ (production ca1 F_length S √(γ F_onset) (1 − ce1 γ), destruction
ca2 F_turb Ω γ (ce2 γ − 1), σ_γ = 1) and Re_θt (relaxation source
c_θt (Re_θt,corr − R̃e_θt)/t_scale, t_scale = 500ν/U², σ_θt = 2; note the
session spec's "50ν/U²" is a typo — the published model uses 500).
Coupling into SST: P̃_k = γ_eff P_k, D̃_k scaled by min(max(γ_eff,0.1),1).
Take the empirical correlations (F_length, Re_θc, F_onset, Re_θt,corr(Tu,λ_θ))
verbatim from Langtry & Menter (2009) / OpenFOAM `kOmegaSSTLM` — do not
re-derive; they are piecewise fits. Transition is a config sub-option and
REQUIRES resolved walls (below).

### IDDES blend (SST-IDDES, Gritskevich et al. 2012 calibration)
- l_RANS = √k/(β* ω); l_LES = C_DES Δ with C_DES = F1-blend of
  (C_DES1=0.78, C_DES2=0.61); Δ = min(C_w·max(Δx,Δy,Δz)-style IDDES mesh
  length — start with the simpler Δ = (ΔxΔyΔz)^{1/3} from `filter_*` and
  the IDDES Δ = min(max(0.15 d_w, 0.15 h_max, h_wn), h_max) as a step-2
  refinement; make Δ a small selectable).
- Blend function: implement the DDES shielding first —
  f_d = 1 − tanh((8 r_d)³), r_d = (ν_t + ν)/(κ² y_eff² √(Σ g_ij g_ij)) —
  reusing velocity_gradient_tensor for Σ g_ij². Full IDDES adds the
  elevating/WMLES branch (f_B, f_e, f_dt); phase it AFTER the DDES form is
  validated (the DDES f_d is a strict subset and the standard first gate).
- l_hyb = f_d l_RANS + (1 − f_d) l_LES  → replaces l_RANS in D_k only.
- ν_t assembly: nut = f_d·nut_rans + (1 − f_d)·nut_wale. Pure RANS is
  f_d ≡ 1 (same kernels, no LES eval — guard the two extra kernels with
  `model == TURB_IDDES`); pure LES never allocates the RANS state at all.

### Discretization of the 4 scalar equations (staggered grid)
- All scalars at CELL CENTRES (with p), halo layer 0:nb+1 like nut.
- ONE fused per-substage kernel updates all scalars: per cell compute the
  velocity-gradient tensor, S, Ω, F1, F2 and the correlations ONCE, then
  the sources for every transported scalar — they share all the expensive
  intermediates. Loop structure = the standard
  `do b / collapse` volume-kernel idiom.
- Convection: 2nd-order TVD (van Leer limiter) on the cell-centred scalars
  — one small `!$omp declare target` helper used by all equations.
  Central differences ring at transition fronts; momentum stays central.
- Diffusion: face-averaged effective diffusivity
  (ν_eff,face = arithmetic mean of neighbours) × central gradient; the
  per-equation σ is just a scalar factor.
- STIFFNESS (important for explicit RK3): treat the destruction terms
  point-implicitly (Patankar): sinks linear in the own variable are
  integrated as division, e.g. ω* = (ω + dt·P_ω)/(1 + dt·β ω), not as an
  explicit subtraction — near walls ω ~ 6ν/(β1 y²) is huge and a fully
  explicit sink forces dt → 0 or negative ω. Same treatment for the k sink
  (β* ω dt in the denominator) and the γ sinks. Enforce floors after each
  update: k ≥ 0, ω ≥ ω_min (small positive), 0 ≤ γ ≤ 1.
- BCs at FACE_PHYS: the scalars need a small cell-centred-scalar BC
  routine (apply_bc in boundary.f90 serves the staggered q only): ghost =
  mirror for Dirichlet (k_wall = 0), copy for Neumann (γ, Re_θt), and ω
  gets the Menter first-cell condition (below) — which is a WALL-CELL
  pinning, not a ghost value. FACE_CLOSED faces: zero-flux, exactly as the
  velocity machinery treats them (ghosts zeroed at init, no exchange).

## IBM wall treatment (from the Weber thesis §4.3, IBMkOmegaSST, translated)

Reference: Master_Thesis_Jannik_Weber.pdf pp. 35–38. Their collocated-FV/
matrix mechanics translate to our explicit staggered penalization solver
as follows.

1. **Wall distance.** `dwall = min(distance to domain wall patches,
   distance to immersed surface)`, used in F1/F2/f_d and the ω condition.
   Regularize (their Eq. 4.37–4.38): `y_eff = max(dwall, ½·min(Δx,Δy,Δz))`
   — without the floor, cell centres grazing the staircase surface get
   dwall → 0, ω ~ y⁻² explodes cell-to-cell along the boundary and the
   noise propagates through the ω transport. Sources:
   - file-based IBM: `mobygeom.py block-table` writes a per-leaf
     `dwall_blocks` dataset (evaluated at each leaf's level, exactly like
     `coef_blocks`); mobygeom already has the surface distance machinery.
   - analytic IBM: compute in the flow-case init (plane/wavy walls are
     closed-form).
   - domain-wall part: computed in the solver at init from the node lines
     (trivial), min'ed in.
2. **IBM wall cells** = fluid cells with ≥1 solid staggered face (the
   existing `ibm_aware` test). Classify ONCE at init into a per-cell byte
   marker in `sst_type`. These are the cell-centred stand-in for a wall
   patch.
3. **ω in wall cells: pinned, not transported.** OpenFOAM constrains the
   matrix BEFORE solving; the explicit-RK3 equivalent is to SET wall-cell
   ω at the START of the substage, before the transport kernel reads
   neighbours — then adjacent cells' diffusion sees the constrained value.
   Value (their Eqs. 4.39–4.42): y⁺ = C_μ^¼ y_eff √k/ν;
   ω_vis = 6ν/(β1 y_eff²) (β1 = 0.075), ω_log = √k/(C_μ^¼ κ y_eff);
   stepwise switch at y⁺_lam (≈ 11, the OpenFOAM nutWallFunction value).
   Solid cells: a constant benign ω (keeps halos/restart finite; the
   physical solution never reads it).
4. **k in wall cells: zero-gradient by discretisation, NO wall function**
   (their kqRWallFunction is just zeroGradient). For us: mask the k (and
   γ, Re_θt) diffusive flux through any solid staggered face — per-cell
   the same thing FACE_CLOSED does per-block — and k = 0 in solid cells.
5. **ν_t in wall cells: overwritten** after assembly: 0 on the viscous
   branch, ν(y⁺κ/ln(E y⁺) − 1) on the log branch (E = 9.8); 0 in solid
   (already the ibm_aware behaviour).
6. **Production in wall cells:** 0 on the viscous branch; on the log
   branch evaluate G from the tangential velocity relative to the local
   IB normal at y_eff. Get the normal locally as ∇dwall (one-sided
   differences of the dwall field) — no extra mobygeom output.
7. **Not applicable to us:** their fluid-masked residual norm (§4.2.3)
   exists for residual-based linear-solver convergence; our projection
   runs fixed Chebyshev-Jacobi iterations. Note it only if residual-based
   stopping is ever added. Their AMR §4.4 is our `refine_body`; the rule
   "all geometry-derived fields re-evaluated per level" is already our
   per-leaf convention — `dwall_blocks` follows it.

## Wall-treatment MODE switch (design fork — get this into config)

`[turbulence] wall_treatment = resolved | wall_function`:
- `resolved` (y⁺₁ ≲ 1): viscous branch only (Menter ω wall condition =
  Eq. 4.42's viscous limb, so it is the SAME code path with the switch
  pinned), no ν_t/G overwrite beyond ν_t = 0 in wall cells' viscous limb.
  Transition model ALLOWED.
- `wall_function` (coarse near-wall grids): full viscous/log blend +
  ν_t/G log corrections. Transition model FORBIDDEN (γ–Re_θt requires
  y⁺ ≲ 1; running it through log wall functions is meaningless) — hard
  config error, not a warning.

## Config (hierarchical: family in [turbulence], sub-models in their own sections)

The sections mirror the module layout one-to-one — [turbulence]→turb_type,
[les]→les_type, [rans]→sst_type:

```
[turbulence]
model = none | les | rans | iddes    (FAMILY; default none. When the key is
                                      absent, a configured [les] model implies
                                      les, so pre-[turbulence] inis run
                                      unchanged; an explicit none wins.)
# T5: the IDDES blend options (delta choice, f_e branch toggle) live here.

[les]           # SGS sub-model; REQUIRED for model = les and iddes
model = smagorinsky | wale
# cs, cw, delta_scale, ibm_aware (as today)

[rans]          # RANS sub-model; REQUIRED for model = rans and iddes (T2)
model          = sst
transition     = true|false        (default false)
wall_treatment = resolved | wall_function   (default resolved)
# SST constants only if overriding published defaults
# inlet/initial turbulence: tu (%), nut_ratio  -> k, omega initial+BC values
```

`[les]` is NOT deprecated — it is the canonical SGS section (T0 briefly
made it an alias; the hierarchical form replaced that same-day). rans /
iddes are recognized and rejected with "not implemented yet" until their
phases land; `model = wale` under [turbulence] is a hard error pointing at
[les]. Validation in `validate_turbulence_values` (family les with no
[les] model → error; later: transition ∧ wall_function → error; sst ∧ no
dwall source → error; etc.).

## Phase plan (each gated; never start N+1 with N unverified)

- **T0 — the hoist (pure refactor, bit-exact). DONE (validated 2026-07-05).**
  `turbulence.f90`/`turb_type` (TURB_NONE/TURB_LES live, TURB_RANS/TURB_IDDES
  reserved) owns `nut` + the metric tables; `les.f90` keeps only the algebraic
  SGS kernels (`update_sgs_viscosity(les, turb, ..., nut)` writes a
  caller-supplied nut target; `velocity_gradient_tensor` reads turb metrics);
  main.f90/step.f90 call sites switched to `turb`;
  `add_les_momentum_correction` renamed `add_eddy_viscosity_correction`
  (reads `turb%nut`, body verbatim). Config is the hierarchical scheme of
  the Config section above: `[turbulence] model = none|les|rans|iddes`
  (family; rans/iddes rejected until implemented; absent key + configured
  `[les] model` implies les) with `[les]` as the canonical SGS section
  (model/cs/cw/delta_scale/ibm_aware, not deprecated). Gate PASSED: bit-exact
  (`-Mnofma`/`-gpu=nofma`, max_abs 0 on un/vn/wn/pn + nut) vs the
  pre-refactor binary (4cd7c97) on min_channel (blocks + 2:1 + chebyshev,
  4-rank CPU), les_ibm channel + refine_body (file IBM + WALE ± 2:1),
  Beltrami y-slab (4-rank CPU) — CPU AND GPU; both builds green;
  `[turbulence]`-section run byte-identical to the `[les]`-alias run.
  Note: the SGS timing profiler moved with the hoist
  (`init_turbulence_profiler`, tag `turb_timing`, TURB_PROF_*) — console
  output name only.
- **T1 — dwall + wall cells. DONE (validated 2026-07-06).**
  `src/modules/rans.f90` / `sst_type` holds the geometry state only:
  `dwall` (raw cell-centred wall distance, ghost-inclusive — every value
  evaluated pointwise from geometry, no exchange), `yeff = max(dwall,
  ½·min(Δx,Δy,Δz))` (use THIS in the model), and the interior byte marker
  `wallcell` (0 fluid / 1 wall = ≥1 of the 6 staggered faces solid, the
  ibm_aware threshold test / 2 solid = all 6), plus
  `enter_/exit_rans_data`. dwall sources: file IBM reads the new per-leaf
  `dwall_blocks` tiles (`mobygeom.py block-table` writes them by default,
  trimesh unsigned surface distance at each leaf's level, `--no-dwall`
  opts out; the solver read cross-checks the file's blocks table like
  coef_blocks and hard-errors on legacy files without the dataset);
  analytic IBM uses `body_surface_distance` (ibm.f90: true Euclidean
  distance to the wavy wall by coarse scan + golden section, sharing the
  extracted `wavy_wall_height` with isInBody); the domain-wall part
  (non-periodic faces with Dirichlet tangential velocities) is min'ed in
  from the node-line ends. HOOK until [turbulence] model = rans exists:
  a `[rans]` section's presence builds the state at init;
  `[rans] dump_geometry = true` writes `<prefix>_ransgeom.h5`
  (blocks table + interior dwall/yeff/wallcell + per-block cell-centre
  coords; new self-contained parallel writer — the field-output path is
  untouched). Gates (validation/rans_geometry/, all PASS): flat les_ibm
  walls y_eff vs exact point-to-slab-box closed form max|err| = 0.0
  (float32 STL vertex quantization is part of the as-built geometry),
  single-level AND per-level under refine_body, wallcell exact; analytic
  wavy section vs independent scipy minimization 1.1e-16; 4-rank dump ==
  1-rank; GPU == CPU; regenerated block-table file's coef/masks/blocks
  byte-identical to the committed les_ibm file (dwall_blocks purely
  additive); full T0 case list (min_channel 4-rank, les_ibm ± refine_body,
  Beltrami y-slab, wavy section) bit-exact (nofma, max_abs 0) CPU AND GPU,
  and a [rans]-on run's fields bit-exact vs [rans]-off (init-only state).
- **T1b — geometry-agnostic analytic dwall. DONE (validated 2026-07-07).**
  New `src/modules/walldist.f90` computes the analytic-IBM dwall from the
  isInBody indicator ALONE (host-only init code; the indicator is a
  procedure argument, so gates drive the same machinery with other
  geometries): (1) surface point cloud by bisecting every finest-level
  cell-centre segment whose endpoints straddle the indicator
  (deterministic two-pass plane scan — thread/rank-count independent;
  every rank builds the same global cloud); (2) nearest cloud point via a
  kd-tree (median split, bbox pruning; chosen over the uniform-bin ring
  search, which degenerates to O((d/bin)^3) bin visits for cells far from
  the surface) with +-L image queries folded to the minimum image along
  periodic directions; (3) POLISH: re-sample the surface on a shrinking
  3x3x3 lattice around the current nearest point, halving the spacing
  until the error bound min(2s, 2s^2/d) reaches `[rans] dwall_tol`
  (default 1e-10) — the raw cloud distance overestimates by O(s^2/d),
  worst at the near-wall cells the model cares about. The wavy-specific
  `body_surface_distance` is DELETED from ibm.f90 (the scipy minimization
  in check_rans_geometry.py is the surviving specialized reference);
  production has no geometry-specific distance code. In periodic
  directions the indicator must be length-periodic (the assumption the
  coefficient machinery already makes at halo coordinates). Gates
  (validation/rans_geometry/, all PASS): (a) wavy.ini through the generic
  path vs scipy 2.28e-11; dwall_tol sweep 1e-2..1e-8 gives
  9.8e-4/2.5e-5/2.6e-7/2.6e-9 — monotone convergence with polish depth;
  (b) sphere straddling the periodic x boundary through the SAME
  machinery (`walldist_test`, src/test_walldist.f90) vs |r - R|: errors
  track tol (1.5e-4/1.6e-7/1.6e-10 at 1e-4/1e-7/1e-10); (c) analytic
  refine_body (wavy_refine.ini): level 0 = 2.35e-11, level 1 = 2.44e-11;
  (d) file-IBM path untouched — flat_l1/flat_refine still exactly 0.0;
  (e) T0/T1 case list bit-exact (nofma, max_abs 0) vs 2b0ff97, CPU AND
  GPU (min_channel 4-rank + GPU, les_ibm +- refine_body, Beltrami
  y-slab, wavy fields), [rans]-on == [rans]-off, ransgeom dump 4-rank ==
  1-rank == GPU.
- **T2 — SST, no transition, resolved mode. DONE (validated 2026-07-08;
  runs on a remote machine via `validation/rans_sst/run_gates.sh`, checks
  local).** rans.f90 owns the transport state (k/omg + oldrhs pairs +
  scratch, RK3 low-storage like the momentum predictor) and the fused
  per-substage kernel (`rans_substage`: constrained-cell pinning BEFORE
  the kernel reads neighbours → scalar ghosts → k/ω halos → transport →
  copyback → nut assembly); `[turbulence] model = rans` + `[rans] model =
  sst` (tu/nut_ratio initial state; transition/wall_function still
  rejected). k/ω snapshots + restart ride a generalized named-scalar io
  (`fdm_h5_append_scalar`/`fdm_h5_read_scalar`; absent datasets →
  reinitialize + warn). DEVIATIONS from this doc, both commented in
  rans.f90: (1) scalar convection is FIRST-ORDER UPWIND, not TVD van
  Leer — the limiter needs a second upwind cell that the single halo
  layer does not carry, and a block-edge fallback would break the
  nb/rank-independence invariant; revisit before the T4 transition
  fronts. (2) the ω cross-diffusion needed explicit-RK hardening: a flat
  freestream ω against the pinned wall rows drives ω through zero, and a
  floored ω flips F1→0 through the CD_kω branch, re-enabling the term
  with 1/ω amplification (observed ω→1e150 in 50 steps). Fix =
  wall-consistent ω IC (max of freestream and 6ν/(β1 y_eff²)) + Patankar
  sign-split (negative part into the point-implicit denominator) + the
  explicit positive part rate-limited to ≤ ω/dt_sub. Gates
  (validation/rans_sst/, all PASS): (a) laminar Re_τ 10, tu 1%: parabola
  to 2.4e-4, k → 3e-47 (at Re_τ 30/tu 5% the no-transition SST correctly
  finds its weakly-turbulent branch — that is model physics, not a bug);
  (b) Re_τ 180/395: U+ centreline 18.16/20.16 vs DNS 18.20/20.13
  (0.2%/0.15%), u_τ 1.001/1.002, log line to 4.9%/6.5% (the pure κ/B
  line deviates several % from real profiles — the DNS centreline anchor
  is the sharper criterion); (c) les_ibm off-grid IBM channel through
  the wall-cell ω pinning: log line to 4.3%, u_τ 0.966 — THE key IBM
  gate; (d) wall-band-refined channel: NO interface band
  (jump/local-variation ratios 0.58–1.11 in k/ω/nut/u), core matches the
  RESOLVED turb180 reference to ≤2.8% (the uniform coarse twin is an
  informational control only — its y+₁≈2.8 sublayer feeds a spurious
  core-k plateau); (e) LES/no-model case list bit-exact (nofma, max_abs
  0, incl. nut) vs 5851c2f CPU AND GPU; RANS 20-step run bit-exact
  across 1-rank == 4-rank == GPU (all fields incl. k/ω/nut). Found and
  fixed while gating (e): `initialise_channel_fields` filled only block
  slot 1 (Phase-0 relic) — any cold-start channel with `[blocks] nb` set
  got a rank-count-dependent, mostly-zero IC; fix loops all blocks with
  block-origin noise indexing (bit-exact for the default
  one-block-per-rank layout, where origin == the old rank-box offset).
- **T3 — wall_function mode. DONE (validated 2026-07-08).**
  `[rans] wall_treatment = wall_function` (config accepted; transition ∧
  wall_function stays a hard error, guard placed ahead of the T4
  rejection). On top of resolved T2, all branch-gated on the mode so
  resolved stays bit-exact: (1) constrained-cell ω (IBM wall cells AND
  domwall rows, both pinning sites host+device) becomes the stepwise
  viscous/log blend `omega_wall_blend` on the k-based y⁺ =
  C_μ^¼ √k y_eff/ν with y⁺_lam = 11.5301 (the fixed point of
  y⁺ = ln(E y⁺)/κ; κ = 0.41, E = 9.8, the OpenFOAM value); (2) wall-cell
  ν_t overwritten after assembly (`nut_wall_value`: 0 viscous,
  ν(y⁺κ/ln(E y⁺) − 1) log, continuous at the switch), and COPIED into
  the no-slip physical-face ghosts (`rans_apply_nut_wall_ghosts`) so the
  face-interpolated ν_t the momentum correction uses at the wall face IS
  the wall value — without the ghost copy the face sees ν_t,w/2 and the
  delivered wall shear (ν + ν_t,face) U₁/y₁ is wrong; (3) log-branch
  wall-cell k production G = (ν + ν_t,w)(|U_t|/y_eff) C_μ^¼√k/(κ y_eff)
  from the tangential velocity relative to the local wall normal
  (`sst%wnorm`, host-precomputed at init as normalized ∇dwall — RAW
  dwall, not the floored yeff — one-sided AWAY from solid staggered
  faces and no-slip physical faces, where dwall is V-shaped/mirrored and
  a centred difference sees spurious zero slope); the viscous branch
  keeps the resolved rules (IBM wall pk = 0, domwall pk normal). The ω
  cross-diffusion hardening and first-order upwind are untouched.
  Gates (validation/rans_sst/, all PASS): (a) wf180_y30/y45 (uniform
  ny = 6/4; ny = 6 is not nb-divisible → [blocks] nb unset, which also
  exercised rank-box blocks): implied U+ centreline vs DNS 18.20 to
  1.2%/0.7%, u_τ (delivered wall stress (ν+ν_t,1)U₁/y₁) = 1.0000;
  (b) sweep y⁺₁ = 5/15/22.5/30/45: centrelines −3.1%/+2.8%/+3.0%/+1.2%/
  +0.7% — a mild +3% buffer overshoot decaying into the log layer, NO
  double-counting dip (that would be a deficit); first cells below
  y⁺ 30 sit 10–19% above the resolved profile = the log-line error at
  the anchor cell, reported informationally (`--mode wallfn` in
  rans_channel_check.py gates rows with y⁺ ≥ 30 and near-centre rows
  against the RESOLVED turb180 profile, which carries the DNS anchor —
  the checker now also handles non-cubic rank-box blocks); (c) ibm180wf
  (the IBM channel at y⁺₁ ~ 2–3 through the wall-function blend, 200k
  steps on the local GPU at 18.6 ms/step): the steady profile matches
  the T2 resolved ibm180 field to ROUND-OFF (u 4.5e-16, nut 1.1e-16) —
  the k-based switch keeps every wall cell on the viscous branch, whose
  arithmetic is exactly the resolved treatment, and the RANS fixed point
  is hardware-independent;
  (d) resolved mode bit-exact vs T2 8991192 (nofma, max_abs 0 incl.
  k/ω/nut) on min_channel 4-rank, les_ibm ± refine_body, Beltrami
  y-slab, turb180 — CPU AND GPU; (e) wf180_y15 20-step: 1-rank == 4-rank
  EXACTLY; CPU vs GPU ≤ 2e-13 (pn) — NOT exact, unlike resolved (which
  still is): the `log()` intrinsic in the wall-function branch differs
  by an ulp between host and device libm (same class as the accepted
  les_ibm masking-branch 4.6e-14 precedent).
- **T4 — γ–Re_θt transition (resolved only).** The two extra scalars +
  correlations + the P̃_k/D̃_k coupling. Gate caveat: the canonical T3-series
  flat-plate cases need an inflow/outflow capability — CHECK first whether
  the bc machinery supports it; if not, gate on what the solver can run:
  (a) laminar channel at subcritical Re with γ active stays laminar
  (γ → 0 in the BL, P_k suppressed) where plain SST would transition;
  (b) developed turbulent channel with transition on reproduces the T2
  result (γ → 1); defer the flat plate to when inflow BCs exist. This
  phase is SEPARABLE — IDDES (T5) does not depend on it.
- **T5 — IDDES blend.** f_d (DDES shielding form first), l_hyb into D_k,
  the nut blend, Δ selectable. Gates: (a) f_d field sane on the developed
  channel (→1 at the wall through the RANS layer, →0 in the core) at a
  WMLES-style grid; (b) channel mean profile: no gross log-layer
  mismatch vs the T2 RANS and the pure-WALE LES references; (c) fd ≡ 1
  recovers T2 RANS answers, fd forced 0 recovers WALE behaviour in the
  core (consistency limits); (d) the les_ibm IBM channel runs IDDES
  stably with the wall treatment; (e) full-suite bit-exactness for
  model ≠ sst-iddes unchanged. THEN (separate increment) the full IDDES
  f_B/f_e elevating branch, gated on the same channel (log-layer-mismatch
  reduction is its entire purpose).

## Watch for

- The consumer chain (nut exchange, `add_eddy_viscosity_correction`,
  dt limits, `fdm_h5_append_nut`) is MOVED in T0 and then NEVER edited —
  that is the whole bit-exactness argument. If a later phase seems to need
  to touch it, stop and reconsider.
- No polymorphism, no allocatable components in arrays of derived types;
  enum dispatch + flat arrays + `enter_/exit_*_data` only (CLAUDE.md).
- Wall-cell ω pinning must happen BEFORE the transport kernel reads
  neighbours in the substage (the explicit analogue of OpenFOAM's matrix
  constraint) — pinning after the update leaks unconstrained values into
  neighbouring cells' diffusion for one substage.
- Point-implicit sinks + floors from day one; a naive explicit ω
  destruction with wall values ω ~ 6ν/(β1 y₁²) will NaN immediately at
  RANS-grid dt.
- The scalars ride `exchange_scalar_halos`, which is per-array — 4-5
  exchanges per substage on top of nut. Accept it now (correctness
  first); batching is a profiling-session item, and the profiling doc
  (`docs/next_session_profiling.md`) should gain a note that RANS adds
  scalar exchanges to the halo-bound budget.
- `restart`: k/ω(/γ/Re_θt) must round-trip through restart files; absent
  datasets (old restarts) → initialize from `tu`/`nut_ratio` and warn,
  mirroring how legacy layouts are handled elsewhere.
- Transition correlations: transcribe from Langtry & Menter (2009) /
  OpenFOAM `kOmegaSSTLM` and unit-test the pure functions on tabulated
  values host-side before they ever run in a kernel; note the session
  spec's t_scale typo (500ν/U², not 50).
- Solid-cell values (benign ω, k=0, nut=0) must be set AT INIT and after
  every restart read, so halos/io never carry garbage — same discipline
  as the pinned-face zeroing in Phase 2 of the block work.
- `git add` explicit paths only (see memory: git-staging-discipline);
  validation output .h5/.png are gitignored on purpose.

## NEXT-SESSION PROMPT (T3 done 2026-07-08; T5 — the IDDES blend — is next;
## T4 transition is separable and can come before or after)

> Read `docs/next_session_iddes.md` and CLAUDE.md. Branch
> `claude/jacobi-interface`. Execute phase T5 ONLY: the IDDES blend
> (SST-IDDES, Gritskevich et al. 2012), DDES shielding form FIRST —
> f_d = 1 − tanh((8 r_d)³), r_d = (ν_t + ν)/(κ² y_eff² √(Σ g_ij g_ij))
> reusing velocity_gradient_tensor; l_RANS = √k/(β* ω), l_LES = C_DES Δ
> (C_DES = F1-blend of 0.78/0.61, Δ = (ΔxΔyΔz)^{1/3} from filter_* to
> start, Δ selectable); l_hyb replaces l_RANS in the k-destruction ONLY;
> nut = f_d nut_rans + (1 − f_d) nut_wale with the blend living entirely
> in turbulence.f90 (les.f90/rans.f90 stay decoupled; pure RANS keeps
> f_d ≡ 1 and never evaluates the LES kernel; pure LES never allocates
> RANS state). Config: [turbulence] model = iddes (requires BOTH [les]
> model and [rans] model = sst; works with either wall_treatment).
> Respect the T2/T3 landmines (memory iddes-plan): the ω cross-diffusion
> hardening is load-bearing; scalar convection is first-order upwind by
> design; the wall-function paths are branch-gated — keep the blend out
> of the branch-gated wall-cell code. Gates per the doc's T5 list:
> (a) f_d sane on the developed channel at a WMLES-style grid (→1 at the
> wall, →0 in the core); (b) channel mean profile: no gross log-layer
> mismatch vs T2 RANS and pure-WALE references; (c) consistency limits:
> f_d forced 1 reproduces the T2/T3 RANS answers, f_d forced 0 recovers
> WALE in the core; (d) the les_ibm IBM channel runs IDDES stably with
> the wall treatment; (e) model ≠ iddes bit-exact vs T3 (nofma, max_abs
> 0, short case list incl. a RANS resolved AND a wall_function case) —
> the elevating f_B/f_e/f_dt branch is a SEPARATE later increment.
> WORKFLOW: long physics runs go through validation/rans_sst/
> run_gates.sh (or a new validation dir) on the big machine — the user
> rsyncs there and back; check scripts run locally; at most one solver
> job at a time on the local machine (the local GPU is fine for
> hour-scale runs, e.g. ibm180wf took ~200k steps there); bit-exactness
> comparisons use sed-shortened ini copies (~20 steps) per memory
> bit-exact-gates-short, never full-length runs. If doing T4 instead:
> first revisit the first-order-upwind front-sharpness question
> (rans.f90 deviation note) — γ–Re_θt fronts may need the TVD limiter,
> which needs a second upwind halo cell; resolve that BEFORE the
> correlations, and gate per the doc's T4 list (transition stays
> resolved-walls-only, the config guard already enforces it).
