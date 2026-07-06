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
- **T2 — SST, no transition, resolved mode.** The scalar-transport
  infrastructure (fused kernel, TVD convection helper, point-implicit
  sinks, floors, scalar BCs, halos, restart/io of k, ω) + nut_rans
  assembly + the ω wall pinning (viscous limb). Gates: (a) laminar
  channel stays laminar and matches the DNS parabola (`tools/
  check_parabolic_channel.py`) with k seeded small → decays; (b) developed
  turbulent channel Re_τ 180 (and one higher-Re, e.g. 395/590) vs the log
  law + DNS mean profile — RANS-quality agreement (a few %), using the
  developed-stats harness from `validation/channel_interface/`; (c) the
  les_ibm off-grid IBM plane channel reproduces (b) through the IBM wall
  treatment — THE key IBM gate; (d) 2:1 interface: band-refined RANS
  channel vs single-level — no interface band in k/ω/nut (the scalar
  exchange is the validated path, but gate it anyway); (e) CPU==GPU to
  round-off; LES/no-model cases still bit-exact.
- **T3 — wall_function mode.** Log-branch ω/ν_t/G in wall cells. Gate:
  the same channel deliberately coarsened to y⁺₁ ≈ 30–50 still recovers
  the log law and the correct bulk U; behaviour degrades gracefully as
  y⁺₁ varies (no double-counting dip).
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

## NEXT-SESSION PROMPT (T1 done 2026-07-06; T2 is next)

> Read `docs/next_session_iddes.md` and CLAUDE.md. Branch
> `claude/jacobi-interface`. Execute phase T2 ONLY: k-ω SST transport, no
> transition, resolved wall mode. Build in `rans.f90` on the T1 geometry
> state (`sst_type` already holds dwall/yeff/wallcell): the k and ω
> arrays + their oldrhs (cell-centred, halo 0:nb+1, RK3 low-storage like
> the momentum predictor), ONE fused per-substage transport kernel
> (velocity-gradient tensor / S / F1 / F2 computed once per cell — reuse
> `velocity_gradient_tensor`), TVD (van Leer) convection helper,
> face-averaged effective diffusivity, POINT-IMPLICIT sinks + floors
> (k ≥ 0, ω ≥ ω_min — a naive explicit ω destruction NaNs immediately),
> the Menter ω wall-cell pinning (viscous limb, BEFORE the transport
> kernel reads neighbours), a small cell-centred-scalar BC routine
> (Dirichlet mirror / Neumann copy at FACE_PHYS), scalar halos via
> `exchange_scalar_halos`, nut_rans = a1 k/max(a1 ω, S F2) assembled into
> `turb%nut` via the turbulence.f90 select-case (enable
> `[turbulence] model = rans` + `[rans] model = sst`; solid cells keep
> benign ω / k = 0 / nut = 0 set at init and after restart), restart/io
> of k and ω riding the block-table machinery (absent datasets →
> initialize from `[rans] tu`/`nut_ratio` and warn), and the k/γ/Re_θt
> diffusive-flux masking through solid staggered faces. Gates (doc T2
> list): (a) laminar channel stays laminar, parabola exact, seeded k
> decays; (b) developed turbulent channel Re_τ 180 + one higher Re vs the
> log law / DNS mean — RANS-quality (a few %); (c) the les_ibm off-grid
> IBM plane channel reproduces (b) through the IBM wall treatment — THE
> key gate; (d) band-refined RANS channel vs single-level: no interface
> band in k/ω/nut; (e) CPU==GPU to round-off; LES/no-model cases still
> bit-exact (nofma, max_abs 0) — the nut consumer chain must not be
> edited. Make a plan first; execute after. T3 (wall functions) → T5 per
> the doc, one phase per gate.
