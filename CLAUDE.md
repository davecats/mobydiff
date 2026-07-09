# CLAUDE.md

Guidance for Claude when working in this repository.

## What this is

mobydiff: an incompressible Navier-Stokes solver. Second-order finite
differences on a staggered Cartesian grid (uniform or stretched per
direction), RK3 time stepping, a segregated **damped-Jacobi pressure
projection** with optional **Chebyshev-Jacobi** acceleration (`[pressure]
accel = chebyshev`; on the `claude/jacobi-interface` branch — it replaced the
old coupled red-black SOR projection, which could not make the 2:1 interface
operators consistent), volume-penalization immersed boundary method (IBM),
optional LES. Fortran + MPI (3D Cartesian decomposition, 26-neighbour halos) with
OpenMP target offload for GPU. Entry points: `src/main.f90` (solver) and
`src/mobygrid.f90` (serial grid/preprocessing tool). Geometry
classification tooling lives in `tools/mobygeom*`.

## Build and run

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
./compile.sh cpu && ./compile.sh gpu   # builds build_cpu/ and build_gpu/
mpirun -n 1 ./build_gpu/main path/to/input.ini
```

- ALWAYS run executables through `mpirun`, even single rank.
- Build both CPU and GPU paths; the CPU path is the reference for debugging.
- The Ubuntu WSL host can be unstable after an update: if commands hang or
  fail oddly, retry before assuming a code problem.

## Coding conventions

- Performance matters, but the code must stay very easy for humans to read.
  Exception: `io.f90` / `field_hdf5.c` may be ugly if needed.
- Avoid duplicated code and complex interfaces that are not strictly needed.
- Comment code: explain intent and non-obvious choices, not syntax.
- GPU programming model is OpenMP target offload (`!$omp target teams
  distribute parallel do`, guarded by `#ifdef USE_OPENMP_OFFLOAD`). Switch a
  kernel to OpenACC only if OpenMP is genuinely missing a feature or leaves
  significant performance behind — and say so explicitly when you do.
- Derived types own flat contiguous allocatable arrays; map them to the
  device once in `enter_*_data`/`exit_*_data` routines (see `gpu_runtime.f90`,
  `blocks.f90`). No allocatable components inside arrays of derived types.

## Active work: block refinement + 2:1 interface (branch `claude/jacobi-interface`)

The block refactor (Phases 0–3) is complete and lives on `claude/jacobi-interface`
(forked from `claude/blocks` to rebuild the projection on a damped-Jacobi /
Chebyshev smoother). The CURRENT state and next steps are in
`docs/next_session_edges_les.md` (read it first); the master design is
`docs/block_refinement_strategy.md`. Goal:
BCM-style equal-size blocks (Nakahashi & Kim 2004; Jansson et al. 2019) to
enable 2:1 local refinement and removal of blocks buried inside the
immersed boundary. Phased, each phase verified before the next:

- Phase 0 (complete, validated 2026-06-11): the solver state lives in a
  `block_set_type`; `field_type` and the grid metric arrays are gone
  (`grid_type` keeps only generation parameters and the node lines;
  blocks slice from them via `slice_grid_direction`). Volume kernels
  loop `do b = 1, blk%nBlocks` folded into their collapse; `ibm%coef/mu`
  and `les%nut` carry the trailing block index. Bit-exact vs `6c03b67`;
  block-index GPU cost: none measurable.
  `tutorials/sailplane/sailplane_ibm_coeff.h5` was corrupt in git since
  at least `2b5e517`; regenerated per the tutorial README, recommitted.
- Phase 1 (complete, validated 2026-06-11): many same-level blocks.
  `[blocks] nb` (cubic, even, ≥4, must divide the global grid) makes the
  grid a uniform block lattice numbered along a Z-order Morton curve and
  split linearly over the ranks (`zorder_owner/start/count`, closed
  form); default = one block per rank box. Everything is per block:
  `physLow/physHigh` face masks drive the momentum starts, the SOR sweep
  window/Neumann terms and `colorOffset = modulo(sum(origin(:,b)),2)`
  inside the kernels; `apply_bc` points carry a block slot; LES tables,
  channel stats likewise. comm.f90 holds block-pair exchange entries
  (one per destination block × 26 directions, tangential extension into
  physical halos as the old rank boxes): same-rank entries are one flat
  device copy kernel overlapped with the MPI messages; off-rank entries
  form one message per peer rank in a canonical order both ends derive
  independently. Each block sweeps its open halo layer redundantly with
  the owner (the rank-level red-black trick one level down), which makes
  results EXACTLY independent of nb and rank count: channel/wavychannel/
  sailplane all bit-exact vs Phase 0 for nb=default and small nb, on
  1/2/8 ranks (channel nb=4 = 176k blocks). io writes one hyperslab per
  block (independent transfers) into unchanged global datasets, plus the
  Z-ordered `blocks` table (origin+level per global id, replacing
  `rank_local_range`); restart works on any rank count.
  GPU time/step (256x128x256 channel): nb=default at Phase-0 parity;
  nb=32 +19% (the (34/32)^3-1 halo-layer overhead), nb=16 +49%; tiny nb
  on big IBM cases is much worse (sailplane nb=10: ~25x) — choose nb=32+
  until Phase 4 tackles exchange overlap.
- Phase 2 (complete, validated 2026-06-12): removal of blocks buried
  inside the immersed boundary. Removable ⇔ block dilated by one halo
  cell solid at cell centres + all three staggered locations; analytic
  IBM classifies at init (`classify_active_blocks`), file-based IBM
  reads the `block_active` table written by `mobygeom.py block-active`
  into the coefficient file (absent table ⇒ keep all, warn).
  `[blocks] remove_solid = false` disables removal. Z-order ids are
  compacted to survivors (`zidOf = -1` for removed); `physLow/physHigh`
  are face kinds (`FACE_OPEN/FACE_PHYS/FACE_CLOSED`): closed faces are
  exact zero-flux via the wall mask machinery (momentum skip, sweep
  window, denom + face-correction merge(), halos and pinned faces
  zeroed once at init, no exchange entries toward removed blocks, the
  tangential extension keyed to "combined edge/corner neighbour
  absent"). apply_bc serves FACE_PHYS only. See strategy doc §7
  (updated: face masks instead of mu=0 at the face).
  Gates: wavychannel nb=4 (1150/125000 removed) and a sphere case
  (file path, 8/1728): fluid cells EXACTLY equal to no-removal,
  velocities ≤ ~1e-26 everywhere (SOLID·mu residual; only the
  decoupled solid-cell pressure differs O(1)); global mass residual
  1e-22; nb-unset / fully-fluid / remove_solid=false / no-flags all
  bit-exact vs Phase 1 (nofma). GPU s/step gain ≈ removed fraction
  (wavychannel −0.9%); tutorial-resolution bodies bury few blocks
  (sailplane at nb=10: zero) — payoff grows with volumetric bodies and
  finer grids. Beware: face-kind consumers must test `/= 0` (no-flux)
  vs `== FACE_PHYS` (BCs) — never treat the kind as arithmetic 0/1.
- Phase 3 (in progress): 2:1 refinement, gated sub-steps.
  - 3a (complete, validated 2026-06-12): multi-level infrastructure.
    Per-level node lines (midpoint subdivision); leaf block table
    (level, origin in level-l cells) replacing the single lattice, ids
    along the finest-lattice Morton curve; per-level `lidOf` lookup;
    `[blocks] refine = x0 x1 y0 y1 z0 z1` + `refine_levels` box
    refinement with 26-neighbour 2:1 smoothing; face kinds gain
    FACE_COARSE/FACE_FINE (interfaces error out until 3b); field io is
    the block-table layout ((nBlocksGlobal, nb^3) datasets, one
    independent row-range write per rank; legacy 3D restarts still
    read; no XDMF — reassemble via compare_fields.py --export-global).
    Gates (nofma, vs Phase 2 c87e1b0): channel nb=4 level 0 bit-exact;
    all-refined 64^3 (4096 level-1 leaves) bit-exact vs 128^3 — exact
    because midpoint subdivision of dyadic uniform lines is bitwise the
    doubled-resolution line.
  - 3b (complete, validated 2026-06-12): 2:1 interface transfer in the
    exchange entries. Entries carry an op (COPY/RESTRICT/PROLONG),
    direction and fine-quarter parity; sampling happens on the SOURCE
    side (src_samples: per dim, RESTRICT averages the 2 cell-centred or
    1 matching face-staggered fine samples — 8/4/4 totals for p /
    tangential / normal velocity — PROLONG injects the covering coarse
    value), so the wire always carries destination-point values. A
    coarse face is fed by up to 4 fine sub-entries (2 per edge, 1 per
    corner), enumerated in fixed child order for the canonical wire
    format. Corner extension generalizes to "combined neighbour
    occupied at ANY level". Momentum predicts the shared face on BOTH
    sides (pinning only PHYS/CLOSED — masking FACE_FINE froze qs=0 and
    the copy kernel zeroed the face: instant blow-up); the sweep
    denominator/correction masks (noflux) cover PHYS/CLOSED/FINE but
    not COARSE, per the fine-owns-face split that 3c completes.
    Gates: uniform (1, 0.5, 0.25) flow through a 216-block refined
    patch preserved EXACTLY (max dev 0.0 after 10 steps; dyadic grid +
    constants make even diffusion round-off vanish); channel-with-patch
    50 steps stable and bounded; channel nb=4 without refinement still
    bit-exact vs Phase 2.
  - 3c (complete, validated 2026-06-12): staggered interface treatment.
    DEVIATION from the original doc 6 (doc updated first): the LOW-side
    block owns the 2:1 shared face — unconditional fine-owns-face is
    unrealizable for fine-west orientations in this storage convention
    (fine face DOFs in unpredictable halos; restriction would write the
    coarse INTERIOR u(1) plane against the prolong reading it). Owner
    predicts and corrects the face normally; the other side's halo copy
    comes from the existing exchange (RESTRICT 4-sub-face average /
    PROLONG injection, both exactly conservative). Sweep masks:
    noflux_low = {PHYS, CLOSED}; noflux_high = any non-open kind (the
    high halo face only feeds the divergence).
    Gates: uniform flow through the patch still EXACT (0.0); global
    mass residual with the patch −3.5e-20 (velocity scale 7e-2);
    laminar channel patch vs uniform-fine reference converges at order
    2.45 (interface band) / 2.70 (away), error mildly localized at the
    interface; channel nb=4 without refinement bit-exact vs Phase 2.
  - 3d (complete for the analytic IBM, validated 2026-06-12):
    geometry-driven refinement. `[blocks] refine_body = true` classifies
    per level in ibmm (`classify_block_geometry`: touch = dilated block
    straddles the surface, buried = fully solid, on lines built by
    midpoint subdivision), and the leaf builder refines touched blocks
    plus a one-block 26-neighbour buffer to the finest level, applies
    2:1 smoothing, and removes buried leaves at every level.
    Gate (64^3 wavy-wall channel, nb=8, GPU): 1408 leaves (1024 fine at
    the wall+buffer) vs the 4096-leaf uniform-fine reference: mean-u
    profile in the refined region matches to 0.015% of peak; coarse far
    field differs by the expected truncation (~1.3% pointwise,
    coarse-averaged). Savings: 2.9x fewer cells, 2.7x faster GPU
    time/step (0.042 vs 0.114 s). Channel nb=4 without refinement
    remains bit-exact vs Phase 2.
  - 3d-file (complete, validated 2026-06-12): file-based IBM path +
    interface-relaxation rework. `mobygeom.py block-table` writes
    per-level `block_touch_l{l}`/`block_buried_l{l}` rasters, the
    `blocks` leaf table (its Python leaf builder mirrors the Fortran
    one; the solver cross-checks row-by-row at read) and per-leaf
    ghost-inclusive coefficient tiles `coef_blocks` evaluated at each
    leaf's level; the solver's `refine_body` accepts coefficient files
    (reads the masks), and legacy global-grid coefficient files stay
    readable for single-level runs (bit-exact vs the block-table format
    on the sailplane). Debugging the sailplane refine_body blow-up
    exposed that the 3c interface relaxation was unconditionally
    unstable (round-off-seeded pressure-jump mode at interfaces,
    per-step gain independent of dt/viscosity/sor — earlier gates were
    blind: uniform flow is exact under any consistent transfer, and
    band-refined channels keep pn = 0). Fixes, in docs §6: (1)
    tq-aware covering-cell source rows for edge/corner PROLONG entries;
    (2) blended pressure ghosts at PROLONG faces (ghost = (2 p_C +
    p_f)/3 uniform, weights from node lines, comm `entry_blend`); (3)
    BCM-style symmetric relaxation — every non-pinned face in the
    denominator and corrected by both adjacent cells (each side its own
    copy), per-colour exchanges apply same-level copies only
    (`exchange_halos(..., interp=.false.)`), the final full exchange
    reconciles copies conservatively to the owner's value. Gates:
    channel nb=4 bit-exact vs Phase 2 c87e1b0 (nofma, CPU 8 ranks +
    GPU); sailplane legacy vs block-table bit-exact; 3D-patch and
    x-band channels stable 1000 steps (formerly NaN by ~200); chanp 1
    vs 8 CPU ranks bit-exact (multi-level MPI path); uniform flow
    through a 3D refined patch exact (spread 0.0); refine_body
    sailplane stable, impulsive transient decaying, refined-region
    pressure 7x closer to uniform-fine than the unrefined run
    (pointwise velocities decohere at Re=1e5 — not a usable gate).
    A manufactured-linear-field halo audit (the `MOBY_HALO_AUDIT` hook,
    since removed in the cleanup) checked every exchange-written halo
    cell on the real layout (1.2M + 7.9M cells, 0 bad). Follow-up refactor (validated
    2026-06-12, bit-exact on the full case list): the exchange is one
    weighted gather (per-dim affine maps from `entry_gather_map`, ghost
    blend = destination-completion weights, no op branches in the
    kernels), and entries are ordered same-level-copies-first with
    prefix counts so the per-colour copy-only exchange is a prefix of
    the full one (shorter messages, no runtime filtering).
- 2:1 interface — turbulent validation (resolved 2026-06-29, branch
  `claude/jacobi-interface`). The energy-conserving **constant-1/2 interface is
  the DEFAULT** (`[blocks] interface_constant_half`): inject the velocity prolong
  + skip the cubic deep-halo reconstruction (the truncation-optimal cubic/metric
  weights break interface energy conservation and DESTABILIZE; const-1/2 is V&V
  2003's order-for-energy trade). The pressure projection is symmetric/SPD at the
  interface (`face_grad` composite stencil + conservative copy reconciliation),
  which is what keeps Chebyshev-Jacobi stable there.
  The u'/v' interface BANDS were a **momentum-reflux artifact**: the reflux
  replaces the coarse interface flux F_coarse by avg(F_fine); for the normal flux
  (q)^2, avg(of squares) >> (avg)^2 (Jensen), so it injects the fine-side resolved
  Reynolds-stress flux into the under-resolved coarse cell. It is conservation-
  CORRECT (verified: vanishes for uniform flow) but pumps the fluctuating stress
  onto the coarse cell — the textbook AMR coarsening band. **`momentum_reflux`
  should be OFF** for the generic 2:1 interface: developed-channel stats (Re_tau
  180, t=5..25) with reflux OFF remove the u' spike (excess 1.56->1.00) and v'
  step (kink 0.31->0.04) at ZERO cost to stability/divergence, and match the
  uniform-256 reference: -<u'v'> to 0.2%, no spurious band (x-y/z-y cross-sections
  clean across the interface). The only residual is a SMALL (~5%, isolated against
  a matched uniform-128 control) under-transmission of small-scale v'/w'
  fluctuation energy into the coarse core — the const-1/2 restriction is mildly
  dissipative (a loss, not a band; the mean transport is exact). Validation in
  `validation/channel_interface/` (developed/, interface_benchmark/, reference.ini
  / uniform128.ini); see `docs/next_session_edges_les.md` for the open items
  (edge/corner + LES validation; less-dissipative interface transfer).
- 2:1 interface — edge/corner + LES validation (DONE 2026-06-30, branch
  `claude/jacobi-interface`). EDGE/CORNER (no-LES): embedded core patch, all 6
  faces/12 edges/8 corners; const-1/2 + reflux-off clean (no band, developed
  stats). LES (WALE) across block refinement + the 2:1 interface VALIDATED in
  developed turbulence (coarse 64x48x64 Re_tau 180, vs a 128^3 no-LES filtered-DNS
  reference): mean U (log law) + Reynolds shear stress match to ~1%; the standard
  coarse-LES/WALE normal-stress bias (u' +5%, v'/w' -10%, partly the unfiltered-
  reference gap); WALE nut->0 at the wall; nut steps by the physical filter-width
  ratio (delta^2) across flat AND edge/corner interfaces with NO spurious band
  (velocity+nut band ratios 0.98-1.03). nut is written to field snapshots via
  `fdm_h5_append_nut` (no-LES output byte-identical). Case + figures + analysis in
  `validation/channel_interface/les/` (run_les.py, les_stats.py, fig_interface_rms.py).
  (LES<->IBM coupling now VALIDATED — see next bullet — so this CAVEAT is lifted.)
  RESIDUAL (the next no-LES task): a small v'-only spike at COARSE-OWNS y-faces
  (the Phase-3c low-block-owns-face orientation asymmetry; the fine cells get their
  interface-NORMAL velocity by prolong-injection of the under-resolved coarse face
  value). Reflux-off shrinks it to ~17% localized v' excess; fine-owns faces are
  clean. See `docs/next_session_interface_normal.md`.
- 2:1 interface — interface-NORMAL velocity asymmetry (RESOLVED 2026-06-30,
  ACCEPTED as the const-1/2 price; no code change). The coarse-owns y-face's ~9%
  v' excess (no LES, `vface_asym.py`: lower/upper v'_rms 0.85/0.78) is a **~2.7x
  spurious sub-coarse-cell fine-structure spike** (v_fine 0.235 vs the ~0.086
  interior/fine-owns baseline). MECHANISM: at a coarse-owns face the fine block
  PREDICTS its interface-normal face `q(1)` but its flux `vv_m=(q(0)+q(1))^2` reads
  the const-1/2 prolong-INJECTED coarse deep halo `q(0)` -> Jensen variance
  injection. Every lever is blocked: (1) the deep halo is the SOURCE but is on a
  stability knife-edge — a gentle `q(0)=q(1)` blew the patch up exponentially
  (~10^4/snapshot, <120 steps); the const-1/2 injection is load-bearing; (2) the
  shared-FACE reconciliation can't reach a predictor-sourced excess AND is
  storage-blocked (single halo layer `q(0:nb+1)` — the fine can't predict its
  high-halo face `q(ny+1)`, which needs `q(ny+2)`; a 2-layer halo would be needed
  and even then no fine data exists across the fine-owns face). So const-1/2's
  stability and this residual are the same trade. Full writeup +
  ruled-out levers: `docs/next_session_interface_normal.md` (RESOLVED header).
- LES<->IBM coupling — VALIDATED (2026-06-30, branch `claude/jacobi-interface`).
  The `ibm_aware` solid-cell nut masking (`les.f90`) exercised on an off-grid IBM
  plane-wall channel (file-based IBM from two wall-slab STLs via mobygeom, uniform
  y, walls mid-cell at y=0.259375/2.259375, fluid gap exactly 2.0, Re_tau~180).
  **NO solver code change**: the existing mask + WALE `sd2` already give the
  physical nut->0 into the IBM wall with NO spurious band, on a single grid AND
  across the 2:1 interface. The prime suspect (a spurious nut spike at the band
  cells, SGS reading the IBM velocity drop as resolved strain) does NOT
  materialise: band/core nut ratio 0.05 (gate 2). Gate 1 solid-cell nut==0 exact;
  gate 3 law of the wall recovered (U+~y+, log 2.44 ln y+ +5, bulk U=15.5); gate 4
  refine_body triple nut(y) smooth across the 2:1 interface, no band; gate 5
  stable (case a 4000 + case c 400 steps, no NaN); gate 7 CPU==GPU to 4.6e-14
  (masking branch). The IBM is IMPLICIT (`mu=1/(1+dt*coef)`, ibm.f90:506) so no dt
  restriction. CONVERGED developed-stats campaign (t=5..25, 51 snapshots, all 3
  cases ~25600 steps) confirms it: gate 2 band/core nut 0.05; gate 3 a_wale +
  b_none + grid-aligned `../les/` collapse on `2.44 ln y+ +5` (WALE bulk U=15.09 >
  no-LES 14.69 — SGS raises the log-layer U toward the reference); gate 4 nut steps
  2.03x up into the coarse core across the 2:1 interface, a smooth step with NO
  band (fine wall bands carry lower nut, the physical filter-width step). Case +
  driver + analysis + committed prereqs + figure (`ibm_les_profiles.png`) in
  `validation/channel_interface/les_ibm/` (README + RESUME_STATUS). DONE.
- Code cleanup (DONE 2026-06-30, branch `claude/jacobi-interface`). Removed the
  19 `MOBY_*` testing/diagnostic hooks (63 refs). Category A (pure diagnostics:
  PROJONLY/PREDONLY/DIVDUMP/RHSDUMP/TERMDUMP/MANUF/KEBAL/KESKEW-env/PHASETIME/
  HALO_AUDIT/RESLOG/STEPDIV) deleted outright with their buffers, slot-parking and
  the whole `main.f90` `contains` block; Category B (algorithmic toggles) collapsed
  to the validated production branch and the losing branch deleted: CHEB* env →
  config `accel = chebyshev`; PHIINTERP → inject (the dead `doInterp` two-pass
  scalar-exchange path removed); VELINJECT → `[blocks] interface_constant_half`;
  IFFILT → filter removed (production α=0); NORECON → the `const_half` guard alone.
  Pure refactor: bit-exact (max_abs 0, un/vn/wn/pn) vs the
  pre-cleanup `-Mnofma`/`-gpu=nofma` binary, CPU AND GPU, on min_channel (blocks +
  2:1 interface + Chebyshev), les_ibm channel + refine_body (file IBM + WALE LES ±
  2:1) and the Beltrami y-slab interface regression. Retired diagnostic drivers
  (`momentum_interface/run_gate.sh`, `interface_benchmark/run_benchmark.py`) carry a
  RETIRED header.
- Production-config lockdown (DONE 2026-07-01, branch `claude/jacobi-interface`).
  The solver now carries ONLY the validated 2:1-interface configuration; the
  interface config toggles are gone (not options -- the validated behaviour is a
  must): `interface_constant_half` removed and the const-1/2 transfer hardwired
  ON (the cubic/metric reconstruction path + `reconstruct_interface_halos` /
  `lim_extrap` / the velocity-prolong tangential-interp branch deleted);
  `momentum_reflux` removed (it was the u'/v' coarse-cell band artifact -- the
  reflux machinery `reflux_*` + `refluxF`/`refluxCorr` deleted); `interface_skew`
  removed (experimental, rode the reflux). The four orphaned dump post-processors
  (`tools/{rhsband,rhsterms,divsum,momsum}.py`) deleted. Bit-exact CPU+GPU on the
  full suite (every production case already ran const-1/2 / reflux-off). Recover
  the removed numerics from history: 4149aa0 (reflux), 1428641 (const-1/2 default
  + interface_skew), 9343a3c / 902e30a (deep-halo reconstruction), df697d8 (corner
  cubic), 61499af (the reflux-band finding); the MOBY_* hooks from 5fcdd0c.
- Spatially-varying volumetric body force (DONE 2026-07-01, branch
  `claude/jacobi-interface`). Config-gated `f(x)` added to the momentum predictor
  ON TOP of the constant `[flow] forcing_*`. Own module `src/modules/bodyforce.f90`
  + `bodyforce_type` owning a flat `f(1:nb,1:nb,1:nb,NVEL,nBlocks)` device-mapped
  array (`enter_/exit_bodyforce_data`). New `[force]` config on `dns` (`enabled`
  default false; `type = profile|file|custom`; `profile = constant|sine` with
  `amp_{x,y,z}`, `k_{x,y,z}`/`dir`; `file`). `profile` fills `f` at each
  component's staggered coord at init; `file` reads fx/fy/fz from an HDF5 velocity-
  layout field (`io.f90 read_force_file`, reuses `fdm_h5_read_field`); `custom`
  = the user edits the `update_bodyforce(bf, blk, dns, g, t)` hook to fill `bf%f`
  in the RK loop (public `bf%f` + `bodyforce_zero`/`_update_to_device`/`_from_device`).
  KEY DESIGN: the force is a SEPARATE correction kernel `add_bodyforce_correction`
  (step.f90, parallel to the LES SGS pass: `qs += dt_alpha*f*mu`, `oldrhs += f`),
  so the fused predictor kernel is byte-for-byte untouched -> disabled is bit-exact
  BY CONSTRUCTION (not a `+0.0`/FMA argument). The `*ibm%mu` mask zeroes the force
  in solid cells (intended; not re-masked). Gates: disabled bit-exact (max_abs 0)
  vs the pre-feature binary on min_channel (blocks+2:1+cheb), les_ibm channel (file
  IBM+WALE) and refine_body, Beltrami y-slab — CPU AND GPU; enabled: a constant `f`
  reproduces the `forcing_x` trajectory to round-off (~5e-15 u), sine profile +
  file source both have the expected effect and are CPU==GPU bit-exact. Design in
  `docs/next_session_bodyforce.md`.
- IDDES phase T0 — turbulence-module hoist (DONE 2026-07-05, branch
  `claude/jacobi-interface`). Pure refactor opening the k-ω SST / IDDES plan
  (`docs/next_session_iddes.md`): new `src/modules/turbulence.f90` with
  `turb_type` (model enum TURB_NONE/TURB_LES live, TURB_RANS/TURB_IDDES
  reserved; the `nut` array; the grid-metric tables hoisted from `les_type`).
  `les.f90` keeps only the algebraic SGS kernels — `update_sgs_viscosity(les,
  turb, blk, dns, ibm, nut)` writes a caller-supplied nut target,
  `velocity_gradient_tensor` reads the turb metrics. step.f90/main.f90
  rewired to `turb`; `add_les_momentum_correction` →
  `add_eddy_viscosity_correction` (reads `turb%nut`; the nut consumer chain
  moved verbatim — never edit it, that is the whole bit-exactness argument
  for the later phases). Config is HIERARCHICAL, mirroring the modules:
  `[turbulence] model = none|les|rans|iddes` selects the FAMILY (rans/iddes
  rejected until implemented; an SGS name here is a hard error pointing at
  [les]); `[les]` is the canonical SGS section (model =
  none|smagorinsky|wale + cs/cw/delta_scale/ibm_aware, NOT deprecated);
  `[rans]` arrives in T2. When the [turbulence] key is absent, a configured
  `[les] model` implies the les family (explicit none wins), so every
  existing ini runs unchanged. The SGS timing profiler moved to
  turbulence.f90 (`turb_timing` tag, TURB_PROF_*). Gate: bit-exact
  (nofma, max_abs 0 incl. nut) vs pre-refactor 4cd7c97 on min_channel
  (blocks + 2:1 + chebyshev, 4-rank CPU), les_ibm channel + refine_body
  (file IBM + WALE ± 2:1), Beltrami y-slab — CPU AND GPU; a run with
  explicit `[turbulence] model = les` is byte-identical to the implied-
  family run.
- IDDES phase T1 — wall distance + IBM wall cells (DONE 2026-07-06, branch
  `claude/jacobi-interface`). New `src/modules/rans.f90` / `sst_type` holding
  ONLY the SST geometry state: cell-centred `dwall` (ghost-inclusive, every
  value pointwise from geometry — no exchange), `yeff = max(dwall,
  ½·min(Δx,Δy,Δz))` (the model must use yeff), interior byte `wallcell`
  (0 fluid / 1 wall = ≥1 of the 6 staggered faces solid, the ibm_aware
  threshold test / 2 solid = all 6) + `enter_/exit_rans_data`. dwall sources:
  file IBM = per-leaf `dwall_blocks` tiles (`mobygeom.py block-table` writes
  them by default at each leaf's level like coef_blocks, `--no-dwall` opts
  out; the solver read cross-checks the file's blocks table and hard-errors
  on legacy files without the dataset); analytic IBM = `body_surface_distance`
  (ibm.f90, coarse scan + golden section to the wavy wall; `wavy_wall_height`
  extracted, shared with isInBody); domain walls (non-periodic faces with
  Dirichlet tangential BCs) min'ed in from the node-line ends. HOOK until
  `[turbulence] model = rans` exists (still rejected): a `[rans]` section's
  presence builds the state at init; `[rans] dump_geometry = true` writes
  `<prefix>_ransgeom.h5` (blocks table + interior dwall/yeff/wallcell +
  per-block cell centres; self-contained parallel writer — the field-output
  path untouched). Gates in `validation/rans_geometry/` (all PASS): flat
  les_ibm walls vs the exact slab-box closed form, max|err| 0.0 single-level
  AND per-level under refine_body, wallcell exact (the STL float32 vertex
  quantization is part of the as-built geometry — reference must round the
  planes through float32); analytic wavy section vs an independent scipy
  minimization 1.1e-16; 4-rank dump == 1-rank; GPU == CPU; the regenerated
  block-table file is byte-identical to the committed les_ibm one apart from
  the added dwall_blocks; T0 case list (min_channel 4-rank CPU + GPU, les_ibm
  ± refine_body, Beltrami y-slab, wavy section) bit-exact (nofma, max_abs 0)
  CPU AND GPU, and [rans]-on fields bit-exact vs [rans]-off (init-only).
- IDDES phase T1b — geometry-agnostic analytic dwall (DONE 2026-07-07,
  branch `claude/jacobi-interface`). New `src/modules/walldist.f90`: the
  analytic-IBM wall distance is computed from the isInBody indicator ALONE
  (host-only init; the indicator is a procedure argument): surface point
  cloud by bisecting finest-level cell-centre segments that straddle the
  indicator (deterministic two-pass plane scan), kd-tree nearest point
  (bbox pruning; ±L image queries folded to the minimum image on periodic
  dims), then a shrinking-3x3x3-lattice POLISH to `[rans] dwall_tol`
  (default 1e-10; the raw cloud distance overestimates by O(s²/d)). The
  wavy-specific `body_surface_distance` is DELETED (the scipy minimization
  in check_rans_geometry.py is the surviving specialized reference); in
  periodic dims the indicator must be length-periodic. Gates all PASS
  (validation/rans_geometry/): wavy generic-vs-scipy 2.3e-11 with monotone
  dwall_tol-sweep convergence; a sphere across the periodic boundary
  through the SAME machinery (`walldist_test`, src/test_walldist.f90)
  tracks tol down to 1.6e-10; refine_body per-level 2.4e-11; flat file-IBM
  gates still 0.0; T0/T1 case list bit-exact (nofma, max_abs 0) CPU AND
  GPU; ransgeom dump 1-rank == 4-rank == GPU.
- IDDES phase T2 — k-ω SST transport, resolved wall mode (DONE 2026-07-08,
  branch `claude/jacobi-interface`). rans.f90 owns k/ω (+oldrhs, scratch;
  RK3 low-storage like the momentum predictor) and the fused per-substage
  kernel: constrained-cell ω pinning BEFORE the kernel reads neighbours
  (IBM wall cells via wallcell, domain no-slip rows via a new domwall
  byte; viscous limb 6ν/(β1 y_eff²)), cell-centred scalar ghosts (k
  mirror-0 at walls, ω copy), k/ω halos via exchange_scalar_halos,
  point-implicit sinks + floors, solid-face diffusive-flux masking, nut =
  a1 k/max(a1 ω, S F2) with wall/solid-cell nut = 0. Config
  `[turbulence] model = rans` + `[rans] model = sst` (tu, nut_ratio;
  transition→T4, wall_function→T3 rejected). k/ω snapshots + restart via
  the new named-scalar io (`fdm_h5_append_scalar`/`fdm_h5_read_scalar`;
  absent → reinit + warn). TWO documented deviations (comments in
  rans.f90): scalar convection is FIRST-ORDER UPWIND (van Leer needs a
  2nd upwind halo cell; a block-edge fallback would break
  nb/rank-independence — revisit before T4 fronts); the ω cross-diffusion
  needed explicit-RK hardening (wall-consistent ω IC + Patankar
  sign-split + rate-limited positive part — a floored ω flips F1→0 via
  the CD_kω branch and the 1/ω term cascades to 1e150 otherwise). Gates
  all PASS (validation/rans_sst/, long runs remote via run_gates.sh):
  laminar Re_τ10 parabola 2.4e-4 + k→3e-47; Re_τ 180/395 U+ centreline
  0.2%/0.15% vs DNS, u_τ 1.001/1.002; les_ibm IBM channel log line 4.3%
  (THE key IBM gate); wall-band-refined channel: no interface band
  (jump ratios ≤1.11), core ≤2.8% vs the resolved reference; LES/
  no-model bit-exact (nofma, max_abs 0) vs 5851c2f CPU+GPU; RANS 1-rank
  == 4-rank == GPU exactly. ALSO FIXED (pre-existing, found by the rank
  gate): initialise_channel_fields filled only block slot 1 — cold-start
  channels with [blocks] nb set got a rank-dependent mostly-zero IC; now
  loops all blocks with block-origin noise indexing (bit-exact for the
  default layout).
- IDDES phase T3 — wall functions (DONE 2026-07-08, branch
  `claude/jacobi-interface`). `[rans] wall_treatment = wall_function`
  (Weber Eqs. 4.39-4.42 / the OpenFOAM omega+nutk wall functions), all
  branch-gated on the mode so resolved stays bit-exact: constrained-cell
  ω (IBM wall + domwall rows) = stepwise viscous/log blend on the k-based
  y⁺ (switch y⁺_lam = 11.5301 = the ln(E y⁺)/κ fixed point; κ = 0.41,
  E = 9.8); wall-cell ν_t = ν(y⁺κ/ln(E y⁺) − 1) on the log branch AND
  copied into the no-slip physical-face ghosts (the momentum correction
  interpolates ν_t to faces — without the ghost copy the wall face sees
  ν_t,w/2 and the delivered wall shear is wrong); log-branch wall-cell
  P_k = (ν+ν_t,w)(|U_t|/y_eff)C_μ^¼√k/(κ y_eff) with U_t tangential to
  the precomputed `sst%wnorm` = normalized ∇dwall (RAW dwall, not the
  floored yeff; one-sided away from solid staggered faces and no-slip
  physical faces, where dwall is V-shaped/mirrored). transition ∧
  wall_function is a hard config error placed ahead of the T4 rejection.
  Gates (validation/rans_sst/, all PASS): y⁺₁ = 30/45 channels hit the
  DNS centreline anchor 18.20 to 1.2%/0.7%, u_τ = 1.0000 (delivered
  wall stress (ν+ν_t,1)U₁/y₁); the y⁺₁ = 5/15/22.5/30/45 sweep degrades
  gracefully (mild +3% buffer overshoot, NO double-counting dip; first
  cells below y⁺ 30 carry the textbook 10-19% log-line error,
  informational); ibm180wf (IBM channel y⁺₁ ~ 2-3 → viscous branch, 200k
  steps on the local GPU) matches T2 resolved ibm180 to ROUND-OFF
  (u 4.5e-16 — the viscous branch IS the resolved arithmetic and the
  RANS fixed point is hardware-independent); resolved mode
  bit-exact vs 8991192 (nofma, max_abs 0 incl. k/ω/nut, CPU AND GPU) on
  min_channel/les_ibm ± refine_body/Beltrami y-slab/turb180;
  wall-function 20-step: 1-rank == 4-rank EXACT, CPU vs GPU ≤ 2e-13 (the
  wall-function `log()` differs an ulp between host/device libm;
  resolved stays exactly CPU==GPU). rans_channel_check.py gained
  `--mode wallfn` (gates y⁺ ≥ 30 + near-centre rows against the RESOLVED
  turb180 profile — the DNS anchor is transitive) and non-cubic rank-box
  block support (the ny = 6 case runs with [blocks] nb unset).
- IDDES phase T4 — γ–Re_θt transition, resolved walls only (DONE
  2026-07-09, branch `claude/jacobi-interface`). `[rans] transition =
  true` (∧ wall_function stays a hard config error): γ and R̃e_θt ride the
  fused substage kernel (shared gradients/S/Ω/F1/F2), RK3 oldrhs pairs,
  Patankar sinks in the OpenFOAM split (+P − ce1·P·γ + E − ce2·E·γ —
  implicit coefficient nonnegative on both sides of the γ = 1/ce2
  destruction sign flip), zero-gradient ghosts, exchange_scalar_halos,
  "gamma"/"rethetat" named-scalar io + restart (absent → reinit + warn;
  arrays are 1-cell dummies when off, the wnorm uniform-map idiom).
  Correlations transcribed VERBATIM from OpenFOAM kOmegaSSTLM.C as pure
  declare-target functions, unit-tested by `src/test_transition.f90` (26
  tabulated values, every branch). Coupling: P̃_k = γ·P_k, k-destruction
  ×min(max(γ,0.1),1) in the point-implicit denominator (×1.0 exact when
  off), F1 = max(F1, F3). γ_eff = γ; γ_sep is a marked later increment.
  STEP-0 decision (deviation comment in rans.f90): first-order upwind
  KEPT — flat plate deferred (an inlet IS composable from Dirichlet
  velocity + Neumann pressure faces — user note 2026-07-09 — but the
  RANS layer is not inlet-aware: domain_face_is_wall reads a Dirichlet
  inlet as a no-slip wall, and the scalars lack inlet ghost values);
  channel fronts are wall-normal with ~zero cross-front velocity
  (measured D_num/D_phys 7.3e-5, `t4_front_check.py`). FOUND WHILE GATING: R̃e_θt's
  diffusivity σ_θt(ν+ν_t) = 2(ν+ν_t) exceeds the Peclet dt budget (2×) —
  explicit diffusion checkerboards to 1e6 in ~40 steps; FIX = its
  diffusion DIAGONAL is point-implicit (same steady state; load-bearing,
  do not simplify). Gates (validation/rans_sst/ `t4` group, all PASS):
  lam30t (Re_τ 30 / tu 5%, where plain SST self-sustains — control lam30:
  parabola off 12.2%, k 0.37) laminarizes: parabola 1.6e-3, wall γ 0.024,
  mean-k 9.1e-3 = the γ-floor residual (stationary t=150→300; check runs
  --k-max 0.02); laminart parabola 2.4e-4, k → 8.6e-16; turb180t γ ≥
  0.999 (y⁺ ≥ 30), U+ centreline 18.44 vs DNS 18.20 (1.3%), u_τ 1.0009;
  transition=false bit-exact vs T3 25ef6ed (nofma, max_abs 0 incl.
  k/ω/nut, CPU AND GPU) on min_channel / les_ibm ± refine_body / Beltrami
  y-slab / turb180 / wf180_y30; transition-on 1==4 ranks EXACT, CPU vs
  GPU ≤ 2.8e-14 (exp/pow ulps, the T3 log() class); restart round-trip
  proven with a changed tu (read ≈122 vs reinit 584), legacy warns +
  reinits. NEXT IDDES phase: T5 (IDDES blend, DDES shielding first) —
  prompt at the end of `docs/next_session_iddes.md`; reject transition
  under model = iddes until validated there. T5 STEP 0 (decided
  2026-07-10): `facePatchType` (wall|patch, non-periodic faces only,
  default = today's tangential-Dirichlet inference so existing inis stay
  bit-exact; `periodic_*` untouched) replaces the `domain_face_is_wall`
  inference — a Dirichlet INLET currently misclassifies as a no-slip
  wall — plus ONE generic cell-centred scalar BC applicator in
  boundary.f90 (mechanics only; rans passes mode tables). Deferred by
  user decision: augmented-q scalar batching (profiling phase, see
  docs/next_session_profiling.md) and the flat-plate inlet increment.
- ALSO PENDING: **Profile + optimise** the GPU step for the 2:1-refined channel.
  The last hard profile is STALE (the reflux that was 23% is removed; the
  `MOBY_PHASETIME` timer is deleted): re-profile first with a minimal removable
  phase timer, then attack the dominant cost (likely the projection's
  per-Jacobi-iteration halo exchanges). Every change is a scheduling refactor and
  must stay bit-exact. Full plan + next-session prompt in
  `docs/next_session_profiling.md` (Phase-4 overlap sketch in
  `docs/nonblocking_overlap_strategy.md`, which predates the Chebyshev-Jacobi
  solver and needs updating).

## Verification

- Pure refactors must be bit-exact vs. the pre-refactor code: compare
  output fields with `tools/compare_fields.py`. Build BOTH sides with
  `-Mnofma` (CPU) / `-Mnofma -gpu=nofma` (GPU) for the comparison —
  default FMA contraction makes the compiler introduce 1-2 ulp
  differences for arithmetically identical source.
- Channel sanity: `tools/check_parabolic_channel.py`. IBM cases:
  `tutorials/sailplane/` (coefficient file path),
  `tutorials/wavychannel/` (analytic `set_ibm_coeff`).
- Refinement phases: uniform-flow preservation across interfaces and global
  mass conservation to round-off (see strategy doc §11 for the full list).
- Never declare a phase done with failing builds or unverified results.
