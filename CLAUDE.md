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
OpenMP target offload for GPU. Entry points: `src/moby_solve.f90` (the
solver; the build keeps a `main` symlink for the older scripts) and
`src/moby_prepare.f90` (the MPI-parallel preprocessor writing the case
file, `docs/prepare_solve_strategy.md`). `tools/mobygeom*` is the RETIRED
Python preprocessor, kept as the cross-implementation validation
reference.

## Build and run

```bash
module load toolkits/nvhpc/25.9
./compile.sh cpu && ./compile.sh gpu   # builds build_cpu/ and build_gpu/
mpirun -n 4 ./build_cpu/moby_prepare case.ini case.h5   # file-geometry cases
mpirun -n 1 ./build_gpu/moby_solve path/to/input.ini    # main is a symlink
```

- ALWAYS run executables through `mpirun`, even single rank.
- Build both CPU and GPU paths; the CPU path is the reference for debugging.
- Additional GPU hosts over SSH (shared filesystem, same paths): istmcetus
  (2x A6000 cc86 — the local build_gpu binary runs; pin CUDA_VISIBLE_DEVICES
  and check nvidia-smi first, shared machine), istmcorax (RTX 5090 cc120 —
  own `build_gpu_corax/`; no 25.9 modulefile, export
  PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.9/{compilers/bin,comm_libs/12.9/hpcx/latest/ompi/bin}).

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
- **After `enter_*_data`, the host and device copies are independent.** Host
  code that READS a mapped array sees whatever the HOST last wrote (device
  kernels do not update it), and a host WRITE is invisible to the device
  until an explicit update. The two escapes are `!$omp target update
  from(x)` (device → host, before a host read) and `... to(x)` (host →
  device, after a host write) — both `#ifdef USE_OPENMP_OFFLOAD`-guarded.
  This is a real defect class: the 2026-08-05 RANS cold-start IC fix was
  correct on CPU and a pure no-op on GPU for exactly this reason. Never
  argue a site is safe — TEST it: make the host-side change and check the
  GPU output MOVES. Full audit (every call site, verdict + probe) in
  `docs/next_session_verification.md` §2.

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
  reinits.
- IDDES T5 STEP 0 — domain-face patch types + generic scalar BCs (DONE
  2026-07-10, commit 91f129f). `[boundary] <dir>_<side>_patch = wall|patch`
  (`facePatchType`; non-periodic faces only — a periodic declaration is a
  config error checked in init_boundary_faces; absent = the tangential-
  Dirichlet inference, so existing inis are bit-exact by construction;
  resolved types print at RANS init). `domain_face_is_wall` lives in
  boundary.f90 and reads the declaration first — declare a Dirichlet
  velocity INLET `patch` so RANS stops classifying it as a no-slip wall.
  ONE generic cell-centred scalar ghost applicator `apply_scalar_bc`
  (per-face SCALAR_BC_NONE/COPY/MIRROR/VALUE over the bc point lists;
  VALUE = the future scalar-inlet hook) replaced the duplicated rans ghost
  kernels (now thin mode-table wrappers). Gates: bit-exact vs T4 e227e68
  (nofma, max_abs 0 incl. all RANS scalars, CPU AND GPU) on min_channel /
  les_ibm ± refine_body / Beltrami y-slab / turb180 / wf180_y30 / lam30t;
  declared wall == inferred exactly; y_min=patch removes the dwall min-in
  (ransgeom dwall == 2−y exactly) + the omega pinning.
- IDDES T5 — DDES-shielding blend (DONE 2026-07-10; first of the two T5
  increments). `[turbulence] model = iddes` (needs [les] SGS + [rans] sst;
  transition/wall_function under iddes are hard errors until validated).
  The blend lives in turbulence.f90 and touches exactly two things: the
  POINT-IMPLICIT k-destruction coefficient √k/l_hyb (`iddes_k_sink_coeff`;
  l_hyb = fd·l_RANS + (1−fd)·C_DES·Δ, C_DES = F1-blend of 0.78/0.61, Δ =
  raw (ΔxΔyΔz)^⅓ without delta_scale) on the rans kernel's iddes branch —
  the pure-RANS branch keeps the T2 arithmetic VERBATIM, which is the
  bit-exactness argument — and nut = fd·nut_rans + (1−fd)·nut_sgs
  (`blend_iddes_nut`). CONVENTION LANDMINE: the stored fd is the
  RANS-RETENTION weight tanh((8 r_d)³) = 1 − f_d^Spalart; implementing
  Spalart's formula verbatim with this blend hands the WALL layer to WALE
  (measured fd(wall)=0, +16% log-layer error before the flip).
  `velocity_gradient_tensor` HOISTED les.f90 → turbulence.f90 (the module
  graph runs turbulence → les/rans; the shielding needs it below both
  producers). turb_type gains nut_sgs/fd (1-cell dummies off-iddes) +
  `[turbulence] fd_force` (validation hook; 0 = SGS limit, 1 = RANS
  limit); fd rides the named-scalar io ("fd"). Gates (validation/iddes/,
  all PASS): fd = 1.000 below y⁺ 5 → <0.001 in the core (handover EARLY,
  y⁺ 5–60 — the known DDES-with-resolved-content behaviour; the f_B/f_e
  elevating branch is increment 2's whole purpose); developed-channel
  log-layer mean U within 3.0%/2.8% of the pure-WALE / T2-RANS
  references (full t=5..25 average); fd_force=0 BIT-EXACT vs pure WALE (the IEEE blend
  identity); fd_force=1 holds converged turb180 to 9.8e-13 over 2000
  steps; iddes_ibm stable 2000 steps; model ≠ iddes bit-exact vs
  post-STEP-0 (nofma, max_abs 0, CPU AND GPU) on the full standard list;
  iddes 1==4 ranks EXACT; CPU vs GPU ≤ 2e-14 (tanh ulps). Deferred:
  augmented-q scalar batching (profiling), the flat-plate inlet increment.
- IDDES T5 increment 2 — full IDDES elevating/WMLES branch (DONE
  2026-07-10; Gritskevich et al. 2012 SST-IDDES). In our RANS-RETENTION
  convention IDDES's f̃_d = max(1 − f_dt, f_B) is simply fd =
  max(fd_dt, f_B): fd_dt = tanh((C_dt1 r_dt)³) with r_dt = ν_t/(κ² y_eff²
  |∇u|) — ν_t ALONE and C_dt1 = 20 (both Gritskevich differences from the
  DDES r_d); l_hyb = fd (1 + f_e) l_RANS + (1 − fd) l_LES
  (iddes_k_sink_coeff, still point-implicit; fd = 1 ∧ f_e = 0 reduces to
  β* ω — the pure-RANS branch stays verbatim). f_e = max(f_e1 − 1, 0)·
  (1 − max(f_t, f_l)) (f_t/f_l on r_dt/r_dl, C_t = 1.87, C_l = 5.0,
  Ψ = 1: WALE needs no low-Re correction); the geometric pieces (f_B,
  f_e1 on α = 0.25 − d_w/h_max, and the mesh length Δ) are STATIC,
  host-precomputed at init (`init_iddes_geometry`, called from main.f90
  between init_turbulence and the device maps; d_w = sst%yeff as a plain
  array). `[turbulence] iddes_delta = iddes` (default: min(max(0.15 d_w,
  0.15 h_max, h_wn), h_max), h_wn = spacing along the dominant |∇dwall|
  axis) | `cbrt`; evaluation toggles `iddes_cdt1` (20; 8 = DDES) and
  `iddes_clip` (Spalart max(0, l_RANS − l_LES) ≡ l_LES → min(l_RANS,
  l_LES)). The WALE blend nut = fd nut_rans + (1 − fd) nut_sgs is KEPT
  with the SAME fd (textbook SST-IDDES has no separate SGS model; ours is
  the validated variant); fd_force now zeroes f_e (force = 1 ⇒ l_hyb =
  l_RANS exactly). Gates (validation/iddes/, all PASS): fd = 1.000
  through y+ ≤ 15 (f_B guarantee, d_w ≈ 0.53 h_max = y+ 18.6), band means
  0.872/0.034 at y+ 5–25/25–60 (DDES: 0.67/0.10) — the handover moved
  outward; log-layer mean U (full t = 5..25) 0.5%/0.7% vs pure-WALE /
  T2-RANS (DDES increment: 3.0%/2.8% — a ~5x improvement); toggles: cdt1=8
  marginally worse, clip never binds, cbrt Δ 3x worse vs WALE ⇒ Gritskevich
  defaults stand; fd_force limits unchanged (0 = bit-exact WALE, 1 = holds
  turb180 to 9.8e-13/2000 steps); iddes_ibm stable; model ≠ iddes
  bit-exact vs pre-increment (nofma, max_abs 0 incl. all RANS scalars,
  CPU AND GPU, 7-case list); iddes 1==4 ranks EXACT; CPU vs GPU 20 steps
  EXACT (max_abs 0). NEXT (IDDES track): flat-plate inlet increment,
  transition/wall_function under iddes (hard errors until validated).
- Airfoil phases A0–A2 (DONE 2026-07-12/13, branch `claude/jacobi-interface`).
  A0, the face concept: `[boundary] <dir>_<side>_patch = wall|patch|inlet|
  outlet` is the ONE user-facing face axis; `resolve_face_bcs`
  (boundary.f90, the validate_patch_types slot) derives the per-variable BC
  rows set-if-unset (explicit `_type` keys win; a contradiction with the
  declaration is a hard config error; explicit ini `_type`/`_value` rows
  now also beat the restart file's — config-is-authority). Outlet ⇒
  internal BC_OUTFLOW normal velocity + Dirichlet p; `_profile = parabola`
  gives per-point boundary values (Poiseuille inlet). A0, the projection
  outlet (pressure_solver.f90): outlet FACE_PHYS faces enter the Jacobi
  DENOMINATOR as 2·d1f·mu (`face_grad_denom`) while the velocity CORRECTION
  uses d1f against the MIRRORED phi ghost (`face_grad_corr`;
  apply_scalar_bc SCALAR_BC_MIRROR after every per-iteration phi exchange)
  — the SPD pair; jacobi_apply corrects outlet high AND low faces.
  LOAD-BEARING (found by the Poiseuille gate): the outflow face needs the
  zero-gradient PREDICTOR write — `apply_bc(..., outflow_copy=.true.)`
  post-momentum and at init/restart, OFF inside the projection loop; a
  face touched only by phi corrections keeps its IC shape forever and the
  run converges drift-free to a WRONG steady state (plug outlet profile,
  O(0.2) crossflow). A1: `[case] name = airfoil`
  (`src/modules/flow/airfoil/airfoil_flow.f90`; [case.airfoil] aoa/u_inf/
  chord/force_sample_interval/runtime_file; x_min,y_min,y_max inlet at
  (U∞cosα, U∞sinα, 0) + x_max outlet via patch types, set-if-unset;
  uniform-freestream init; the case `after_step` interface gained the
  `ibm` arg). A2: C_L/C_D = the penalization integral ∫coef·u dV per
  component on the end-of-step field; per-block sums → global-id scatter →
  EXACT allreduce (one contributor per entry) → ordered final sum ⇒ forces
  BYTE-IDENTICAL across rank counts (measured; CPU vs GPU 0.0).
  `tools/make_airfoil_stl.py` writes extruded cylinder/NACA-4-digit STLs
  (ibmc venv gained shapely + mapbox_earcut). Gates all PASS: the standard
  7-case suite bit-exact (nofma, max_abs 0, CPU AND GPU — every A0/A1/A2
  branch is dormant without a declared inlet/outlet); `validation/
  freestream/`: oblique box EXACT (0.0), in/outflow Poiseuille = periodic
  reference to O(h²) with p linear + outlet-pinned and 5.6e-16 drift,
  Lamb–Oseen exits (reflected fraction 5e-3), 1==4 ranks EXACT, wall-twin
  bit-exact + contradicting key error-stops; `validation/cylinder/`:
  empty domain C_L=C_D=0.0 EXACTLY, Re 40 C_D 1.6924±3e-5 (unbounded band
  1.5–1.6 + ~6% 16D-Dirichlet blockage + first-order penalization
  D_eff≈D+h), Re 100 St 0.168 from the spectral FUNDAMENTAL (the confined
  C_L carries a comparable-power 3rd harmonic — zero-crossing counting
  reads 3×St), mean C_D 1.448, mean C_L 2e-4; the Gauss/CV border-flux
  cross-check reproduces the penalization C_D to 0.1% (tight box) / 2.4% /
  6.5% (largest box) — the independent validation of the force statistic
  AND the in/outflow faces. CAVEAT: niter=6 IBM runs accumulate a large
  VELOCITY-NEUTRAL oscillating mode in stored pn (the channel pn-drift
  family) — dynamics and the u-only penalization force are immune, but
  border fluxes need a clean-p snapshot: ZERO pn in a copy of the
  converged restart and rerun ~300 steps at niter=60 (no kick, forces hold
  steady; restarting the POLLUTED p at niter=60 instead transients
  violently — the accumulated spurious ∇p loses its self-consistent
  sloppy-projection compensation — and whole-case niter=60 is ~15x cost).
  NEXT: A3 (RANS scalar inlet values keyed on PATCH_INLET + SD7003 Re 6e4
  transition benchmark; TVD upwind revisit only if the measured γ-front
  smearing demands it) — `docs/next_session_airfoil.md`.
- Airfoil phase A3 (DONE 2026-07-14, increments 0-3 all gated, branch
  `claude/jacobi-interface`; results + findings in
  `docs/next_session_airfoil.md` STATUS and the per-case READMEs).
  INCREMENT 0 (081f387): 3-level refine_body gates
  (`validation/multilevel_body/`, cylinder, 224/192/496 leaves) — uniform
  oblique u,v,w flow EXACT (0.0 incl. pn) across every interface level on
  CPU 1/4 ranks + GPU via a zero-force twin; per-level dwall 3.6e-15 vs an
  exact prism reference. FIXED: trimesh proximity returned dwall O(1e-6)
  high near quantized surfaces → mobygeom dwall tiles use igl
  point_mesh_squared_distance (exact; flat T1 gates still 0.0).
  INCREMENT 1 (25dbffe): RANS scalar inlets — rans.f90 ghost-mode tables
  are pure functions of the DECLARED patch type (SCALAR_BC_VALUE at
  PATCH_INLET with init-computed freestream values: k∞ = 1.5(tu/100 U∞)²
  from the face velocity magnitude, ω∞ = k∞/(nut_ratio ν) unblended,
  γ∞ = 1, R̃e_θt∞ = the T4 λ=0 correlation; COPY at PATCH_OUTLET); 7-case
  suite bit-exact CPU+GPU (dormant in channels); inlet channel holds
  k∞/ω∞ (0.2 %/2.3 %), 1==4 EXACT (`validation/rans_inlet/`).
  INCREMENT 2 (…afdcda0): NACA 0012 SST sanity PASS
  (`validation/naca0012/`: 12c×12c×0.1875c, base 512²×8, refine_body 5
  levels, Δ=1.465e-3c, dtmax 4e-4 — explicit eddy diffusion vs the
  molecular-only Peclet limiter): C_L(0/4/8) = −0.001/0.384/0.745, slope
  0.093/deg (85 % of 2π, 12c-blockage class), C_D(0) 0.0186. LANDMINES
  FOUND: (1) **penalization forces REQUIRE keep-buried files** (cf68225):
  refine_body's removed core absorbs pressure loading through FACE_CLOSED
  faces outside the coef bookkeeping — first run read C_L 0.018 while the
  flow carried Γ ⇒ C_L 0.37; `mobygeom block-table --keep-buried`; the
  cylinder was immune only via its legacy no-block_active file. (2) File
  IBM is **second-order already**: mobygeom always writes the graded
  sharp-interface coefficients (Σ((d0−d)/d)/d0², = the analytic
  USE_IBM_SECONDORDER formula); LE "chequerboard" = its staircase residual
  (~1 % rms u fan at the nose, ω 40 %@2 cells→0@40, no (−1)^(i+j) mode) +
  slice replication at level interfaces; escalation = calibrated
  smoothed-mask/Brinkman (post-A3). (3) mobygeom's block-mask integral
  image now builds in 2 GB chunks (window_solid_counts; the airfoil L4
  lattice needed ~70 GB monolithic).
  INCREMENT 3 (2d514ca, 50f03bd): SD7003 Re 6e4 / α=4 / tu 0.1 %
  transition benchmark PASS — after implementing **γ_sep (LM Eq. 18)**,
  which the first run proved REQUIRED (bubble separated at x/c 0.223 but
  the shear layer stayed laminar, k ≤ 5e-5: separation-induced transition
  IS the SD7003 mechanism; integral C_L/C_D sat in the published band
  anyway — field-level k/γ gates caught it). γ_eff = max(γ, γ_sep) into
  P_k/dkfac only; unit-tested; 7-case suite BIT-EXACT incl. lam30t
  (γ_sep exactly 0 sub-critical). Results: x_s 0.22 (pub 0.22-0.30),
  x_t(k-onset) 0.427 (pub RANS-LM 0.53-0.58, early edge), reattachment
  ~0.56-0.68 (pub 0.65-0.70), turbulent k 3.7e-2, C_L 0.562 / C_D 0.0267.
  **Measured transition-front smearing: 104 level-4 cells (0.152c) across
  k = 10→1000 k∞ — the TVD/van-Leer + second-scalar-halo increment is now
  measurement-justified as follow-up** (do not start without a session).
  γ = 1 in the FREESTREAM (LM convention): transition metrics must use
  near-wall k onset, not γ bands (check_sd7003.py). Selig STLs:
  `make_airfoil_stl.py selig` (+ --resample; SD7003 committed in
  validation/sd7003/). Remote GPU hosts for parallel runs: see "Build and
  run" (istmcetus/istmcorax; memory: remote-hosts).
- 2D block refinement R2D-0..3 (DONE 2026-07-15, commits becbda3/6384aac/
  f0ce7a7/9dd33c2 + follow-up, branch `claude/jacobi-interface`;
  docs/next_session_refine2d.md STATUS header, gates in
  validation/refine2d/ + naca0012/sd7003 READMEs). `[blocks] refine_dims =
  xyz (default octree) | xz` (quadtree: blocks refine in x,z only; y keeps
  ONE global — possibly stretched — node line at all levels). Mechanics: a
  per-direction mask (dns%block_refine_mask / blk%refMask) makes every
  level scaling `2**(l*mask(d))`; 2x1x2 children; CANONICAL mixed Morton
  ids in xz mode (y tile in key bits 42+ over the 2D x,z Morton curve —
  mobygeom/make_channel_restart mirror it); `refine_dims` file attribute
  (xz only, xyz files byte-identical; restart cross-checks BLOCK-layout
  files only); exchange entries generated per-direction (x/z faces 2 fine
  sub-entries; y faces the NEW conforming-normal type: 4 in-plane
  sub-entries, copy form in every unrefined dim, ghost blend degenerates
  to identity) with ZERO kernel changes — the per-dim affine gather maps
  cover it; projection face_grad gains a per-dim `refined` flag (conforming
  y faces use the uniform d1f in BOTH denominator and correction — the SPD
  pair; refined dims keep the locked 2/3-4/3 composite). [case.airfoil]
  span = z|y + make_airfoil_stl --span y re-orient the quasi-2D airfoil
  (chord x, LIFT z, span y periodic never refined). Every phase: 7-case
  nofma suite bit-exact CPU+GPU. Gates: mobygeom==Fortran xz tables
  row-by-row; all-refined-xz == doubled-x,z twin bit-exact (subdivided
  lines); uniform oblique flow EXACT (0.0 incl. pn) through 3-level xz
  patches/body twins, 1==4 ranks, CPU==GPU; Beltrami patch order 2.83 with
  interface error ~30% BELOW the validated octree at equal base; Re_tau 180
  xz wall-band channel (shared stretched y line): NO interface band (jump
  ratios <=1.04), core == the validated reflux-off signature vs uniform128;
  analytic refine_body dwall in xz 2.3e-11; L5-xz NACA fan bench REPRODUCES
  the R1a collapse identically at 6748 vs 65094 leaves = 0.046 s/step (31%
  of L4-3D, 12% of L5-3D on the 5090); AoA sweep 0/4/8 at L5-xz passes all
  bands (slope 88% of 2pi, drag ~0.002 below L4-3D — the halved D_eff~D+h
  bias); SD7003 L4-xz == L4-3D to EVERY checker digit (x_t 0.427, smear
  104.0 cells, C_L 0.5617/C_D 0.0267). FINDING (records the next TVD
  motivation): at L5-xz the SD7003 separation-induced transition does NOT
  fire — k stays at freestream, C_D reads laminar-high — the first-order
  upwind gamma/Re_thetat front (104 cells = 0.152c) no longer triggers
  gamma_sep at the finer bubble; the deferred TVD/van-Leer scalar
  increment is now doubly measurement-justified (own session; then rerun
  the L5-xz SD7003).
- Prepare/solve split P0 (DONE 2026-07-16, branch `claude/jacobi-interface`).
  Strategy in `docs/prepare_solve_strategy.md` (AMPHIBIOUS-style two-step,
  gated P0-P3; read it first): consolidate mobygrid+mobygeom into ONE
  MPI-parallel Fortran preprocessor reusing the solver's own init kernels;
  the case file IS the block-table coefficient file (`[ibm] coeff_file`).
  P0 = analytic geometry, zero new geometry code: `src/moby_prepare.f90`
  (CMake `moby_prepare`; `moby_prepare case.ini case.h5`) runs classify_* /
  init_block_set (Z-order world split) / set_ibm_coeff /
  fill_body_distance_analytic (now public, array-arg) and writes via
  `write_case_file` (io.f90) + `fdm_h5_case_*` (field_hdf5.c), each writer
  the exact inverse of its reader (blocks + coef_blocks + masks +
  block_active + dwall_blocks). Solver solve path UNCHANGED. Gates all PASS
  (`validation/prepare/`): wavy / wavy_refine / wavysolid (1150/125000 =
  the Phase-2 count) solve-from-file bit-exact (max_abs 0) vs inline,
  ransgeom dumps identical, prepare 1==4 ranks identical files, GPU==CPU
  from the same file at tol 0, 7-case suite bit-exact nofma CPU+GPU.
  Prepare with the CPU build (canonical; GPU build computes coef on device,
  libm ulps).
- Prepare/solve split P1 (DONE 2026-07-16, branch `claude/jacobi-interface`).
  STL geometry in moby_prepare with zero mobygeom involvement:
  `src/modules/geometry_stl.f90` (binary STL reader, BVH, majority-vote
  ray-parity inside test with deterministic large-rotation degenerate
  retries, minimum-image periodic queries, EXACT BVH point-triangle dwall)
  behind the analytic indicator signature — `body_indicator_i` moved to
  ibmm, `classify_*` / `fill_body_distance_analytic` take the indicator as
  an argument (solver passes isInBody), `set_ibm_coeff_host` = host twin
  of the device kernel (KEEP IN LOCKSTEP comment). `[ibm] stl_file` is a
  prepare-only input (solver hard-errors without coeff_file). CONVENTIONS
  (comments in geometry_stl.f90): distance images a periodic dim only when
  the mesh is narrower than the cell (full-span padded slabs are their own
  periodic continuation; their skin is not a wall); dwall ghosts beyond a
  periodic boundary carry the minimum-image distance (walldist convention;
  mobygeom stores base-mesh there — interiors agree to round-off). Gates
  all PASS (`validation/prepare/run_gates_stl.sh`): flat/flat_refine vs
  the COMMITTED mobygeom refs (blocks + every mask IDENTICAL, coef <=
  2.1e-8 rel, interior dwall <= 2.7e-11), 1-step solve prepared-vs-
  committed <= 1e-10; sphere vs generated mobygeom ref (masks identical,
  coef 4e-7, dwall 1.7e-16); sphere float32-EXACTLY translated onto the
  x-periodic boundary -> rolled masks + BIT-IDENTICAL tiles (0.0); P0
  analytic gates 22/22 through the refactor; 7-case suite bit-exact vs
  P0 binaries CPU+GPU.
- Prepare/solve split P1b (DONE 2026-07-17, branch `claude/jacobi-interface`).
  The big committed geometries from prepare-built files + mobygeom
  retirement (gates `validation/prepare/run_gates_big.sh`, status table in
  its README). Features the cases demanded: repeatable `[ibm] stl_file`
  (paths with spaces), `stl_scale`/`stl_translate` (float64 v*scale+t,
  mobygeom order), ASCII STL parsed straight to float64 (trimesh
  rounding), `[blocks] keep_buried` (zeroes buried masks in the analytic
  classify branch; LOAD-BEARING for penalization forces), solid-possible
  bbox cull (`stl_cull_box` -> optional cullLo/cullHi in classify_*) +
  OpenMP lattice loops (L5 airfoil level-4 lattice = 1.7e7 blocks).
  Gates all PASS: NACA-0012 L5 / SD7003 L5 (keep_buried+dwall) /
  sailplane (ASCII 55k-tri CAD, transform, nb=10) vs regenerated mobygeom
  refs - leaf tables + every mask level IDENTICAL, ZERO classification
  flips (76M/72M/93M points), graded coef in the near-grazing envelope
  (<=1.2e-4/6.2e-5/1.2e-3 rel; 16320/11680/2 outliers - ((d0-d)/d)/d0^2
  blows up relatively as d->d0, mu unaffected), interior dwall <=1.3e-10;
  200-step GPU solves prepared-vs-mobygeom <=1.5e-7 fields / 8-digit
  C_L,C_D; sailplane 1-step vs the COMMITTED legacy file BIT-EXACT.
  Prepare outruns mobygeom (SD7003 L5: 5m18s vs 6m49s at --jobs 16).
  mobygeom geometry subcommands RETIRED (README_mobygeom + docstring),
  kept as the cross-implementation reference. LANDMINE:
  tools/compare_fields.py reassembles block-table snapshots onto the
  FINEST lattice (69 GB at L5) -> OOM-killed (can take the whole session
  with it); use the chunked validation/prepare/compare_snapshots.py for
  deep-refinement snapshots.
- Prepare/solve split P2 (DONE 2026-07-17, branch `claude/jacobi-interface`).
  Rank-split geometry classification, the last redundant heavy prepare
  stage: classify_active_blocks / classify_block_geometry take optional
  (nsplit, isplit) and classify a contiguous flattened-raster range
  (closed-form raster decode, unchanged cell arithmetic); the wrappers
  merge with an elementwise integer MAX (`comm_allreduce_max_int`,
  MPI_COMM_WORLD -- world-rank indexed, prepare has no Cartesian
  topology). Only the owner can write a 1 => the merge is EXACT and masks
  are identical on any rank count. moby_prepare AND the solver's inline
  analytic path (main.f90) pass their communicator, so analytic
  refine_body init parallelizes too. Gates all PASS: NACA L5 P2 case file
  dataset-IDENTICAL (h5same) to P1b's; flat_refine identical to the
  committed mobygeom reference; 7-case suite bit-exact vs P1b binaries
  CPU+GPU; P0 22/22 + P1 16/16 re-runs (their 1==4-rank identity gates
  exercise the split+merge directly). Timing at 4 ranks x 4 threads:
  NACA L5 prepare 7m13s -> 3m56s, flat_refine 38s -> 22.6s. Still
  redundant per rank (cheap, documented): leaf-table build, STL load/BVH,
  the analytic walldist cloud.
- Prepare/solve split P3 (DONE 2026-07-17, branch `claude/jacobi-interface`).
  Retire and rename, completing the split: `main.f90` -> `moby_solve.f90`
  (CMake target `moby_solve`; the build keeps a `main` SYMLINK -- committed
  scripts/tutorials launch "main"). `mobygrid.f90` DELETED, absorbed into
  the case file: write_case_file appends the node lines + mobygrid-format
  attrs (`fdm_h5_case_append_grid`), so a prepared case file serves
  directly as the retired mobygeom's --grid-file (sphere gate generates
  its reference that way; legacy mobygrid grid files stay readable; the
  dead write_grid_export/fdm_h5_write_grid deleted). Production
  preprocessing consolidated: naca0012/sd7003 setup.sh run moby_prepare
  (venv only for STL generation), run_gates_big.sh references read the
  grid from case files, sailplane README documents the prepare flow.
  Gates all PASS: 7-case suite bit-exact (nofma, via the main symlink) vs
  the P2 binaries CPU AND GPU; P0 22/22 + P1 16/16 re-runs. The
  prepare/solve split (P0-P3) is COMPLETE: one Fortran executable pair,
  one case-file contract, mobygeom/mobygrid retired to validation
  references.
- Boundary layer B0 — laminar Blasius precursor (DONE 2026-07-19, branch
  `boundaryLayer` off claude/jacobi-interface; docs/next_session_boundary_layer.md).
  `[boundary] <x-face>_<u|v>_profile = blasius` (+`blasius_theta`; value key =
  U_inf, shooting-solved table in boundary.f90, outer asymptote beyond eta=12)
  and the `[grid.<d>] one_sided` key for the natural distribution (flag existed,
  key did not; NOT in snapshot attrs — restart inis keep [grid.y]).
  `tutorials/turbulentBoundaryLayer/`: Re_theta,in=100, 400x100x4 theta units,
  inlet(Blasius u+v)/outlet/wall/outlet(top) faces; flow = template.ini mint ->
  make_blasius_ic.py analytic IC -> blasius2d.ini -> compare_blasius.py
  (independent ODE; stations <= 0.7, the last 15% is the outlet zone). Gates
  PASS: theta <= 1.13%, H <= 0.34%, du/Ue <= 2.4e-3, in-layer dv <= 6.6e-2,
  steady to 6e-4; top entrainment ~9-19% below Blasius aloft (p=0 top, info).
  Dormant-path bit-exact vs 94a9249 (pois_io 200 steps, max_abs 0).
  FINDING (open, solver-level): **chebyshev + niter=6 + Dirichlet-p outlets +
  dt~0.5 is unstable on long steady runs** — 2-dx pressure mode from
  round-off, e-fold ~36 t.u., IC-independent, saturates at O(1) divergence;
  STABLE controls: plain Jacobi niter=6, chebyshev niter=60, chebyshev
  dtmax=0.25 (free-slip top still rings => the x outlet alone suffices).
  A per-step resonance: the sign-oscillating 6-term Chebyshev residual
  polynomial against the momentum/BC map (the A2 velocity-neutral pn family
  gone velocity-active through outlet faces); dt detunes it. The case ships
  accel-off. Latent risk for long chebyshev+outlet runs (cylinder/naca
  horizons were too short to show it).
- Passive scalars S0-S5a (DONE 2026-08-03/04, branch `scalar` off
  `boundaryLayer`; plan + every deviation, finding and gate number in
  `docs/next_session_scalar.md` STATUS header, gate machinery in
  `validation/scalar/README.md`). N user-selectable transported scalars, each
  with its own Pr/Sc. THE ONE DESIGN DECISION: **scalars are extra variables of
  `blk%q`** (`dns%nVar = NVAR + dns%nScalar`, scalar `is` at `VAR_S0+is`;
  `SCR_S0+is` in `qs`/`oldrhs` — p has no scratch, so THE TWO INDICES DIFFER),
  which makes the metric arrays (a scalar sits exactly at the `VAR_P` position),
  the 2:1 halo exchange (cell-centred ⇒ the PRESSURE transfer: 8-cell restrict +
  blended `(2 p_C + p_f)/3` prolong ghost — deliberately unlike the RANS
  scalars' plain injection) and the io (one file, one collective write) come for
  free. `count = 0` is bit-exact BY CONSTRUCTION and gated at max_abs 0 (nofma,
  CPU AND GPU) on the standard 7-case suite at every increment.
  `src/modules/scalar.f90` owns `scalar_type` (allocatable per-scalar arrays —
  there is NO `MAX_SCALARS`), the fused `scalar_transport` (one GPU launch for
  any N: collapse(4) + inner scalar loop) and `scalar_finish` (qs→q, ghosts, ONE
  batched exchange), both OUTSIDE the projection.
  - S0/S1: config `[scalar]`/`[scalar.N]`, `apply_scalar_bc_q` (a var-indexed
    twin of apply_bc's cell-centred branch — do NOT extend apply_bc, it runs
    `nIter`x inside the projection), patch-derived BC defaults; 2nd-order
    central divergence-form convection on the p-cell's own face velocities,
    `convection = skew` honoured, face-flux diffusion, momentum RK3 verbatim.
    **NO upwind option** (project stance; the sharp-front over/undershoot is a
    known limitation until the shared TVD increment lands). Called BEFORE
    `momentum` — after it the velocity is the non-solenoidal predictor with
    stale halos.
  - S2: `D_face = 1/(Re Pr) + ½(nut_L+nut_R)/Pr_t(face)` reading `turb%nut`, so
    ONE path covers LES/RANS/IDDES; `prt_model = kays` (Kays-Crawford, pure
    declare-target, small-x series branch or it loses all precision at large
    Pe_t) unit-tested in `src/test_scalar.f90`. Peclet limiter scaled by
    `nu_eff = ire/Pr_min + nut/Pr_t,min`.
  - S3: `ibm%coef` gains its `VAR_P` column ONLY when `[scalar]` is configured
    (**`ibm%mu` does NOT** — `mu_s = 1/(1+dt_gamma coef_p/Pr)` is Pr-dependent
    and formed inline); `ibm_wall = dirichlet` (penalization, solid cell == the
    body value to the last bit) | `adiabatic` (six-face convective AND diffusive
    masking, symmetric ⇒ `∫s dV` conserved exactly). File path reads the NEW
    OPTIONAL `coef_p_blocks` case-file dataset — scalar-free case files stay
    byte-identical, and scalars + a file without it is a hard error naming the
    fix (**any pre-S3 coefficient file must be re-prepared**; the generated
    zero-force twins are fixed by re-running their generator, not moby_prepare).
    **FINDING: the A2 penalization integral does NOT transpose to a Dirichlet
    scalar** — `coef_p (s_body − s)` is `1e28 x 0 = 0` in every solid cell and
    the integral sees only 63 % of the heat; the body heat release is measured
    as staircase-interface flux + graded-cell penalization (validated by a
    discrete energy budget to 3.9e-4).
  - S4: `src/modules/scalar_stats.f90` — seven columns per scalar per row
    (`<s>`, `<s²>`, `<u_c s>`, and on the LOW and HIGH y face the convective and
    the TOTAL flux) built with the TRANSPORT KERNEL's own face diffusivity, so a
    wall row's `J` IS the exact discrete wall flux and `theta_tau`/Nusselt are
    exact rather than reconstructed; `stats_layout = profile|plane`, per-level
    files, restart-continued. DEVIATION: a solver-level facility called from
    `moby_solve.f90`, NOT a case component — the case `after_step` interface
    carries neither `sc` nor `turb%nut`. `tools/scalar_stats.py` reads it;
    `compare_fields.py` with no dataset arguments now discovers datasets.
  - S5a: the Kader/Jayatilleke THERMAL WALL FUNCTION for `[rans]
    wall_treatment = wall_function` (was S2's hard config error), delivered
    exactly as T3 delivers the wall shear — **as a wall-cell eddy
    DIFFUSIVITY**, so the ordinary face-flux discretisation reproduces the
    wall-function flux with no special-cased flux anywhere (measured identity
    1e-15). THE T3 LESSON TRANSPOSED: the wall-cell value must also be copied
    into the no-slip ghosts or the face-interpolated diffusivity, and the
    delivered flux, are halved. `kappa`/`E` are USE-ASSOCIATED from rans.f90 —
    one definition of the log law. The wall cell uses the CONSTANT `prt` even
    under `prt_model = kays` (deliberate: P and the log branch are defined with
    a constant Pr_t); gate (x2) measures what that costs.
  - **THE PASSIVE-SCALAR PLAN IS CONCLUDED (decided 2026-08-07).** S0–S5a are
    the shipped feature set, and the verification debt behind them is closed
    (`docs/next_session_verification.md`: every gate group re-measured, plus
    a host/device staleness audit that came back clean). The two remaining
    plan items are reclassified and are BOTH LOW PRIORITY:
    **S5b** (TVD/van-Leer convection) moves to the comm/halo track — it is a
    halo-DEPTH change in comm.f90 (second upwind cell, the per-dim affine
    gather maps, every 2:1 transfer), shared with the RANS transition
    scalars, and its measured motivation is the SD7003 gamma front, not
    anything a passive scalar failed; its gates are `validation/
    interface_suite/` + `validation/refine2d/`. **S5c** (Boussinesq) is
    parked: the `[force] type = custom` hook is in and gated, so it is user
    code away. WHAT CONCLUDING RATIFIES: scalar convection stays 2nd-order
    CENTRAL with no upwind option, so sharp fronts over/undershoot — the
    documented stance, now the shipped behaviour. The live scalar work is
    conjugate heat transfer, `docs/next_session_conjugate.md` C1, which gets
    its own session (it is a third `ibm_wall` mode on the same cut cells the
    penalization pins).
  - IBM THERMAL WALL FUNCTION — case built + PARTIALLY gated 2026-08-05
    (S5a's open item; every S5a gate was a domain wall).
    `validation/scalar/ibmwf180.ini`: the les_ibm wall slabs on a COARSE grid
    (ly 2.5 / ny 8), where the classified wall cells are CUT cells (centre
    inside the solid, one fluid staggered face) that the Dirichlet
    penalization pins (u 2e-27, theta 1.6e-29) while the wall function reads
    them. ONE body value serves both walls, so the case drives the scalar
    with isothermal walls + a constant volumetric `source`: the steady budget
    is then CLOSED-FORM (fluxes telescope; heat into the body =
    `source*V_fluid`) and NO reference run is needed -- gated to 1.4e-15.
    LOG BRANCH -- CLOSED at Re_tau 1000 (`ibmwf1000.ini`, GPU on istmcetus).
    KEY INSIGHT: **coarsening the grid is not the y+ lever at an immersed
    wall.** The cut cell's velocity is penalized, so its k stays small
    however coarse the grid; y+ only picks up the bounded growth of
    y_eff ~ dwall <= dy/2 (ibmwf180 already has 6x ibm180wf's y_eff and still
    converges to y+ 5.8, the conduction branch). The lever is `nu`: at fixed
    u_tau = 1, y+ ∝ Re. `ibmwf1000.ini` is the same ini with re = 1000 and
    nothing else (the case file is re-prepared because the coefficients carry
    the 1/Re scaling; dwall/yeff/wallcell come out BIT-IDENTICAL, so y+ moves
    only through nu and k). Converged: y+_k **38.3-41.0** (mean 39.7), the
    log branch fires on **128/128** wall cells for BOTH wall functions, the
    wall-cell nut matches the independent transcription **exactly (0.0, on
    GPU)**, and the closed-form budget holds to **2.3e-15**. So the
    penalization and the thermal wall function coexist correctly on the same
    cut cell in the log branch. BONUS control: the same case run with the
    PRE-FIX GPU binary read 132 % high, reproducing the double count at a
    second Re and on a second device.
  - FOUND BY THAT GATE, FIXED 2026-08-05: **the S4 body-heat diagnostic
    double counted whenever `ibm_value = 0`.** The staircase/penalization
    split assumes a solid cell contributes exactly 0 to the penalization sum
    (it holds the body value to the last bit -- the S3 FINDING), so its heat
    is carried by the staircase term. That cancellation is an ARTEFACT OF THE
    VALUE: it is bitwise only when `ibm_value` is large enough to swallow the
    O(1e-29) penalization residual under its own ulp. With `ibm_value = 0` the
    residual survives, `coef_p*(0 - 1.6e-29)` is O(0.1) per cell, and the heat
    was counted twice. MEASURED on the same physics with the two conventions
    (theta differs by a constant; the staircase column is bit-identical to 13
    digits): `ibm_value = 1` gave the exact closed-form answer 0.46263770630106
    (rel dev 1.7e-14), `ibm_value = 0` gave 1.0556 -- 128 % high. FIX: exclude
    solid cells from the penalization accumulator explicitly
    (`scalar_stats.f90`), which is a NO-OP at `ibm_value = 1` -- re-gated, the
    S4 heat gates read their recorded numbers (solver-vs-Python 1.3e-15,
    energy budget 4.49e-04, adiabatic exactly 0) and the fixed binary at
    `ibm_value = 0` now reproduces the `ibm_value = 1` numbers to 1e-15. The
    diagnostic is now value-independent, which is the invariance it must have.
    LESSON: never rely on a floating-point cancellation as a classification.
  - FIXED 2026-08-05 (pre-existing, found while gating S5a): a cold-started
    RANS run was **not rank/nb-independent**. `init_rans_transport`'s
    `k = 1.5 (tu/100 |u|)^2` IC interpolates the cell velocity from the two
    staggered faces, so at a block's LAST interior cell it read the halo
    `q(nb+1)` — which `moby_solve.f90` did not fill until after the whole
    init block: `k` came out a factor **4** low on the last plane of EVERY
    block (`omega` following it down to the viscous limb), one bad plane per
    block, so the answer depended on the decomposition. Visible in x only on
    the channels because their `v`/`w` vanish in the IC; a nonzero-`v`/`w`
    cold start was wrong on the high face of every block in all three
    directions. FIX: `apply_bc` + `exchange_halos` now run BEFORE the
    `[rans]` init block (the later calls are idempotent — ghosts and halos
    only, from interior data nothing modifies in between), so non-RANS cases
    stay bit-exact. **PLUS `!$omp target update from(blk%q)` after it**: those
    calls write the DEVICE copy while `init_rans_transport` is HOST code, so
    without the update the fix was a pure NO-OP on GPU (measured: the fixed
    GPU binary reproduced the pre-fix result bit-for-bit and kept the
    x-dependent k, x-spread 1.49e-02 vs 0.0 on CPU). This is the general
    trap: **any host-side code reading a device-mapped array after its
    `enter_*_data` sees a stale host copy** (and any host-side WRITE is
    invisible to the device until a `target update to`) — the
    `target update from(ibm%coef)` a few lines below is the same pattern.
    Only the GPU bit-exactness suite could catch it; the CPU one looked
    complete. **Exactly one instance of this class has been found; the audit
    for the others is `docs/next_session_verification.md`**, which also
    carries the gate groups not re-run after these two fixes. Only the k/omega INITIAL CONDITION of a cold start
    changes; the converged answer does not (`wf180_y30`/`wf180_y45`
    reproduce the T3 gate to every printed digit), but cold-started RANS
    snapshots written before this date no longer reproduce bit-for-bit.
    `run_bitexact.sh` vs the S5a binaries now reads the INTENDED signature on
    CPU **and GPU**: min_channel / les_ibm / les_ibm_refine / beltrami_yslab
    max_abs 0, and turb180 / wf180_y30 / lam30t moved (20-step transients
    from the corrected IC: un ~2e-2, k ~3e-1; turb180's GPU deviation is
    IDENTICAL to its CPU one to every digit). Same for `run_bitexact_s3.sh`:
    8 of 9 scalar cases max_abs 0, `turbsst` moved (it is cold-started RANS).
    The S2/S5a PHYSICS is unmoved — re-measured, not assumed: the S5a
    `theta_tau` sweep still reads +1.39/+2.28/+5.64/+7.22 %, the resolved
    reference still 0.053413, Kader still 17.1/0.1/0.7 %, and the S5a
    determinism gate is now max_abs 0 for the X SPLIT too (CPU vs GPU 5.6e-17).
    **THE REFERENCE SET TO USE IS `~/s5c_ref_binaries/`** (CPU+GPU nofma
    solve/prepare, commit-pinned to `8f60944`, PROVENANCE inside).
    `~/s5a_ref_binaries/` is stale for cold-started RANS cases, and so is the
    **GPU** binary of `~/s5b_ref_binaries/` — found 2026-08-07: it was
    archived in the four minutes BEFORE the `target update from(blk%q)` was
    added, so it is the "correct on CPU, no-op on GPU" intermediate state and
    gives output IDENTICAL to s5a's GPU binary on turb180 (its own
    PROVENANCE claimed both fixes; that holds for the CPU binary only). It
    silently failed a GPU gate before it was caught. LESSON, now in the s5c
    provenance: **cut a reference set from a COMMIT, never from a working
    tree mid-edit, and record the hash.** LESSON: `q`'s halos are invalid
    throughout init — any init-time consumer reading a neighbour must
    exchange first. Write-up in `validation/rans_sst/README.md`.
  - LANDMINE (cost an hour): a local `nVar` in blocks.f90 SHADOWED the
    use-associated `NVAR` parameter (Fortran is case-insensitive) and silently
    allocated a zero-size dimension. The local is now `nQ`.
  - Conjugate heat transfer at the immersed interface is the planned third
    `ibm_wall` mode: `docs/next_session_conjugate.md` (PROPOSAL, not started;
    derivation in `docs/conjugate/conjugate_ibm.tex`).
- ALSO PENDING: **Profile + optimise** the GPU step for the 2:1-refined channel.
  **STEP 1 (re-profile) IS DONE 2026-08-07; step 2 (optimise) is not started.**
  The phase timer is back as a CONFIG key, not an env hook: `[output] profile
  = true` builds chron.f90's `step_timing` profiler (six `STEP_PROF_*`
  buckets: momentum, syncface_exchange, projection_jacobi,
  projection_exchange, projection_bc, scalar_transport);
  `pressure_projection` takes it as an OPTIONAL argument so the production
  path never reads the clock. Fresh split (A6000, 100 steps, min_channel
  refined 128x64x8, niter=6 cheb): loop 31.3 ms/step = projection_exchange
  44.5% + projection_jacobi 33.9% + syncface 7.9% + momentum 6.1% +
  projection_bc 4.2%; **all halo exchange 52.4%**, and on the heavier
  128x64x128 refined case 42.4% with the projection's exchange (79.1 ms)
  costing as much as its Jacobi kernels (82.1 ms). So the target is the
  projection's per-iteration exchange, i.e. the doc's hypothesis 1
  (core/shell nonblocking overlap) — ceiling ~40% of the step. The nb lever
  is measured too: unrefined 128x64x128, nb 8/16/32 = 48.6/39.3/36.4
  ms/step (**-25%**) with the exchange share 44.6%/35.2%/35.1% — but the
  refined channel is PINNED to nb=8 (its 24-cell wall band is 3 block rows;
  nb=16 does not divide it), so bigger blocks are not the fix there.
  LANDMINE: `nb` must divide the global grid in every direction, so an nb
  sweep is impossible on the nz=8 case. Every step-2 change is a scheduling
  refactor and must stay bit-exact. Full plan, the measured table and the
  next-session prompt in
  `docs/next_session_profiling.md` (Phase-4 overlap sketch in
  `docs/nonblocking_overlap_strategy.md`, which predates the Chebyshev-Jacobi
  solver and needs updating).
- Verification debt + host/device staleness audit (DONE 2026-08-07, branch
  `scalar`, `docs/next_session_verification.md`). Every gate group left
  un-measured by the 2026-08-05 fixes was re-run and reproduces; the audit of
  host-side consumers of device-mapped arrays came back CLEAN (verdict table
  + the GPU probe behind each row in §A; the class is now a coding convention
  above). TWO OPEN ITEMS it produced, neither a scalar feature:
  - **A `t_final`-terminated run takes one extra step whose `dt` is the
    accumulated round-off in `t_current`, and that step's `pn` is amplified
    by `1/dt`** — measured |pn| 1.5e6 at step 60001 vs 9.1 at 60000 with the
    velocity identical to 8.5e-6 (10401 steps instead of 10400;
    `t_current = 29.999999999975433` against a 1e-12 stopping tolerance).
    **The FINAL snapshot of any `t_final` run is therefore a bad restart** —
    this is the real mechanism behind the "never restart from the campaign's
    final `*_50001.h5`" landmine, which had been blamed on the niter = 6
    pn-drift mode. Fixed at the gate level (`run_gates_s2.sh last_periodic`);
    and FIXED in the solver (`trim_dt_for_final_time` suppresses a
    `remaining` below `FINAL_STEP_FRACTION*dns%dt`, 1e-6). The test must be
    RELATIVE: an absolute snap at `run_should_continue`'s tolerance would NOT
    catch this, since the round-off in `t_current` grows like `N eps t_final`
    (6.9e-11 here) and outruns the fixed 1e-12 floor, while
    `remaining/dt = 4.9e-08` identifies it scale-free. Gated three ways:
    inert on all 32 bit-exactness case-runs, effective on the leg that
    produced it (10400 steps not 10401, `max|pn|` 9.09 not 1.5e6, last real
    step bit-identical), and a GENUINE final partial step still taken.
    Details in §B1.
    **CORRECTED 2026-08-28.** The first form of that fix signalled the end by
    SETTING `dns%dt = 0` and letting the loop's `dt <= 0` exit fire — but
    `dns%dt` is written into every snapshot's metadata and read back by the
    restart, so the final snapshot of a `t_final` run became unusable in a
    NEW way (`config.f90: time step must be positive`). It traded one bad
    final restart for another. `trim_dt_for_final_time` is now a LOGICAL
    FUNCTION that reports "no step left" and never leaves `dt` non-positive;
    the trajectory is unchanged (same exit step, a genuine partial step still
    gets `dt = remaining`), only the recorded metadata differs.
    **And the claim "every suite case is nsteps-terminated" was FALSE**:
    `run_gates_s2.sh`'s `les_legs()` REWRITES `t_final` in a generated
    variant, so its `les` and `band` legs are `t_final`-terminated. That is
    what caught this — `RT_turbles.h5` came out with `dt = 0.0` at
    `t = 29.999999999975433`, the S2 `les` leg died on the restart, and
    because that leg deletes `turbles_*.h5` before it fails it took S4 down
    with it. LESSON: a claim of the form "no case exercises this path" must
    be checked against the GENERATED inis, not only the committed ones.
  - **`check_scalar_turb.py cmd_band` compares RAW `theta'_rms`**, so the two
    campaigns' wall-flux difference enters the ratio as a global factor:
    measured core ratio 1.0365 == the flux ratio 1.0360, and normalising by
    each run's own `theta_tau` leaves a localized −2.6 % DEFICIT at the
    interface (the const-1/2 restriction's known dissipation — the opposite
    sign from a spurious band). Normalising inside the checker is suggested,
    not applied. Write-up in `validation/scalar/README.md`.

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
