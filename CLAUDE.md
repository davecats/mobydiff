# CLAUDE.md

Guidance for Claude when working in this repository.

## What this is

mobydiff: an incompressible Navier-Stokes solver. Second-order finite
differences on a staggered Cartesian grid (uniform or stretched per
direction), RK3 time stepping, a segregated **damped-Jacobi pressure
projection** with optional **Chebyshev-Jacobi** acceleration (`[pressure]
accel = chebyshev`; it replaced the old coupled red-black SOR projection, which
could not make the 2:1 interface operators consistent on the OPERATOR of the
day; red-black is selectable again
via `[pressure] solver = redblack` and, since R1, runs across a 2:1 interface
too — see `validation/redblack_interface/`), volume-penalization immersed
boundary method (IBM),
optional LES, RANS (k-omega SST, transition, wall functions) and IDDES, and
**passive scalars** with conjugate heat transfer at the immersed interface
(`src/modules/scalar.f90`, `[scalar]`/`[scalar.N]`). Fortran + MPI (3D
Cartesian decomposition, 26-neighbour halos) with OpenMP target offload for
GPU. Entry points: `src/moby_solve.f90` (the solver; the build keeps a `main`
symlink for the older scripts) and
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

## Active work

**THE BRANCHES ARE CONSOLIDATED AND THE WORK IS ON `main` (2026-09-25).**
`boundaryLayer` (the CaNS/SIMSON trip) and `scalar` (passive scalars +
conjugate heat transfer) were merged into
`optimiseBlockRefinement_parentBoundaryLayer`, itself a descendant of
`claude/jacobi-interface`; `main` was then fast-forwarded to it and that branch
retired, so **`main` is the working branch** -- the long feature-branch phase is
over. `multiGPU` and `claude/blocks` are deleted (the latter archived as the
tag `archive/claude-blocks`, since it was not contained), and `boundaryLayer`
and `scalar` are deleted too once `main` carried them (no tag needed -- they
are contained). The one branch still holding unmerged work is
`claude/jacobi-interface` — six RANS/airfoil features (`[rans]`
kpin_box/ktrip_box/kpin_dwall/boostconv, `[case.airfoil] steady_tol`,
refine_body_box/levels) plus the 2026-08 `naca/rans` tutorial state; keep it
until those are wanted or explicitly dropped. The log below is one list, in
rough chronological order: block refinement, then RANS/IDDES, then the
airfoil, then passive scalars, then the performance campaign.

The block refactor (Phases 0–3) is complete (the projection was rebuilt on a
damped-Jacobi / Chebyshev smoother, which is why `claude/jacobi-interface` was
forked from `claude/blocks`). The CURRENT state and next steps are in
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
  ids in xz mode (2D x,z Morton key in the high bits over the y tile in the
  low 21 — y sat in the HIGH bits until 2026-08-28, which made the order
  y-major and so put the whole 2:1 interface on a rank boundary at EVERY rank
  count; peer traffic 20.5x lower after the swap, `docs/next_session_2to1_penalty.md`;
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
- CHT validation campaigns — CONSOLIDATED into `tutorials/cht/` (2026-09-21).
  `tutorials/cht/channel/` = the Flageul-matched turbulent channel (flat,
  grid-aligned interface, so the cut-face coefficient is EXACT and what is
  measured is the conjugate physics); `tutorials/cht/pipe/` = the Neuhauser
  NekRS pipe (the first CURVED conjugate interface). Each directory holds what
  is needed to RUN the case; its `asset/` holds the reference data, the
  comparison scripts, the figures and the report, and runs without any large
  file. `tutorials/cht/pipe/asset/neuhauser_profiles.npz` is a 250 kB
  reduction of the 11.8 GB published archive — `extract_neuhauser.py` rebuilds
  it and figures 1-5 come out BYTE-IDENTICAL either way. The superseded
  channel campaigns (1-4) are gone from the tree and live in git history +
  `tutorials/cht/channel/asset/CAMPAIGN_NOTES.md`. `make_geometry_stl.py` and
  `check_annulus.py` moved to `tools/` (shared by the C2/F5 gates and the pipe
  tutorial). The FEATURE gates stay in `validation/conjugate/`.
- Conjugate heat transfer at the immersed interface — increments **C1, C2 and
  C3 DONE** (C1 2026-08-27/28, C2 and C3 2026-08-28, branch `scalar`; plan +
  every deviation and gate number in the STATUS header of
  `docs/next_session_conjugate.md`, commands and measurements in
  `validation/conjugate/README.md`).
  **C2 — the tangential term the C1 baseline drops: measured, and SHIPPED
  DISABLED.** `[scalar] indicator_interval` reduces `e_face = h_d s_t/(T_R −
  T_L)` over the cut faces (max, rms, count) — §3's own closed form, so it IS
  the relative error the baseline makes, and it ships ENABLED as a free
  diagnostic; `[scalar.N] tangential_correction` adds `s_t (k_loc − k_face)`
  as its OWN flux divergence after the C1 term (inert by construction when
  off, not by a `+0.0` argument) and is DEFAULT OFF because three independent
  measurements say so: a crossover ratio `r* ≈ 0.027` that the DNS-like
  `|∇_tT|/|∂_nT| ~ 10⁻²` sits below; the `r = ∞` BVP where the corrected FACE
  flux is exact to 1e-11 and the corrected FIELD is 2.6× WORSE (a
  finite-difference divergence wants the face AVERAGE, not the midpoint
  value); and — decisive — a cylinder at κ_s = 10³ where the correction is
  40–65× worse than dropping the term, because `|∇_tT| ∝ 2/(1+κ_s) → 0` in the
  isothermal limit while the discrete estimate of it does not. The measured
  fix (`k_area`, the face AREA-weighted mean, which removes the paradox
  exactly on a plane) is recorded in the README and NOT implemented: it cannot
  rescue the feature, because the binding constraint is the `s_t` estimate on
  curved interfaces. **Anyone reopening this starts from `s_t`, not from the
  multiplier.**
  **C3 — the fluid-fraction-weighted capacity, the Nusselt diagnostic, and the
  TIME STEP.** (1) `C_cell = f + (1−f)C_s` with `f` from the closed
  plane-in-box form (`plane_box_fraction`; `clip(½ + φ_c/h)` is not an
  approximation of it but its degenerate limit, the one a GRID-ALIGNED wall
  lands in). `f` is pure geometry, built once at init beside `φ` (`sc%vfrac`)
  and read by the transport kernel, the time-step limiter AND the checkers —
  it rides the snapshot as `vfrac` (write only), which is what lets a gate
  state the invariant the SOLVER conserves. Gated: field error order **1.99,
  2.00** and decay rate order **1.99, 2.00** on a two-material eigenmode,
  against **1.71, 1.09** and **1.28, 0.54** from the archived C2 binary on the
  same inis — the pointwise capacity costs a full order, measured with a
  control, not argued. (2) The conjugate branch of the interface-heat
  diagnostic (`scalar_stats.f90`) = the Nusselt number: the sum over cut faces
  of the scheme's own flux, penalization column EXPLICITLY zero (there is none
  in this mode, and `coef_p` is finite in a graded fluid cell — the 2026-08-05
  lesson applied before it could bite). Gated three ways: a closed form
  (**4.572e-13**), an independent Python transcription (**3.553e-16** / 
  **2.709e-16**), and a control-volume budget with the flow on that falls like
  `dt²` and does not move with `niter`. FOUND: the PRE-C3 diagnostic on a
  conjugate case reported a MEANINGLESS number, so C1's "smoke-gated" columns
  were also wrong. (3) **The cut-cell share in `scalar_conjugate_peclet_rate`
  is now 2, not 3** — C2 measured the old convention at 96 % of the Gershgorin
  bound and found the case that attains it; `share = 2` gives a 1.56× margin
  and that case runs. Cost: `dt × 2/3` AT CUT CELLS ONLY. Every C1/C2 number
  was therefore re-measured, not assumed.
  **C1 — the baseline** (derivation in `docs/conjugate/conjugate_ibm.tex`).
  `[scalar.N] ibm_wall = conjugate` is a
  THIRD mode beside `dirichlet`/`adiabatic` (re-expressing those two through
  the new arithmetic could not be bit-exact, so they are untouched): the solid
  stops being a boundary condition and becomes a REAL unknown carrying
  `solid_k` = κ_s and `solid_rhocp` = C_s (both ≡ 1 in the fluid), plus
  `solid_init` / `solid_source` / `contact_resistance`. **THE WHOLE SCHEME IS
  ONE FACE COEFFICIENT**: at a face whose two cell centres straddle the
  interface, the face diffusivity becomes the distance-weighted harmonic mean
  `dm/(w/κ_L + (1−w)/κ_R + R_c dm/h)` on the level-set fraction
  `w = φ_L/(φ_L − φ_R)`; every other face keeps today's kernel line verbatim.
  The OBLIQUITY LEMMA is why this needs no new data: φ stores the
  PERPENDICULAR distance, so the direction cosine cancels in the ratio and `w`
  is the true arm fraction at any orientation, with no normal ever computed —
  **no new dataset, no case-file format change**. `sc%phi = ±dwall` is
  ghost-inclusive, signed by the cell-centred IBM marker
  (`|coef(VAR_P)| > 1e20` is exactly "this centre is in the solid", because
  `set_ibm_coeff` writes SOLID/Re there and grades only fluid-centred cells),
  and comes from the two EXISTING dwall producers. Around the coefficient:
  convection HARD-masked on solid and cut faces with the skew term's
  divergence built from the SAME masked velocities (so a uniform scalar is
  still preserved exactly, in the solid too); the flux divergence divided by
  the local capacity; ν_t in neither the solid nor a cut face; no
  penalization. HARD CONFIG ERRORS (each silently produces a wrong answer):
  `remove_solid = true`, `refine_body` without `keep_buried`, `ibm_value`,
  a `solid_*` key on a non-conjugate scalar, no immersed body, wall functions.
  The 2:1 precondition is CHECKED at init, not assumed (a cut face may not sit
  on a coarse/fine block face — the coefficient is a same-level arm).
  Gates (all PASS, `validation/conjugate/run_gates_c1.sh`): 1D two-material
  slab with the cut swept through a full cell × κ_s over five decades —
  max|θ − exact| **≤ 9.7e-15** over 28 pairs, and `w` itself **2.2e-14** vs
  the analytic cut, rebuilt from the case file without the solver;
  capacity-independent steady state; contact resistance exact (R_c = 0.25
  gives **0.0**); `Σ C θ dV` drift **1.2e-16 relative** with the flow on;
  1 == 4 ranks and CPU == GPU **max_abs 0** on both geometry paths; 5/5 config
  guards rejected; `[scalar] count = 0` and every `ibm_wall /= conjugate` run
  **max_abs 0** vs `~/s5c_ref_binaries` on the 7-case and 9-case suites, CPU
  AND GPU. TWO TIME-STEP FINDINGS, both from gate 1 and both recorded in
  `scalar_conjugate_peclet_rate`: (1) the explicit limit is NOT the max over
  materials of `α = κ/C` — a cut face's `k_face` reaches `max(κ_L,κ_R)` but
  feeds the cell on the OTHER side, whose capacity is the other material's, so
  a fluid cell against a κ_s = 1000 solid sees 1000× the fluid rate even at
  α_s = α_f; the rate is built from the ACTUAL face coefficients instead
  (also 500× less conservative at w = ½); (2) a cut cell attains the
  Gershgorin factor `ρ ≤ 2A_ii/C_i` that the uniform interior never excites,
  so its rate is raised — measured, `pecletmax` 0.3 blows up / 0.2 is stable.
  C1 doubled it (share 6 → 3); C2 measured that even so the default sat at
  96 % of the bound and found an oblique high-contrast case that went to NaN
  there, and **C3 went to share 2** (a 3× rate at cut cells, a 1.56× margin).
  NEXT (optional, its own session): **C4, conducting sharp corners** — the
  COCO wedge model. It rests on the same pointwise-flux premise C2 falsified,
  so read C2's verdict and the README's "the way out" BEFORE starting it; the
  first thing to settle is the `s_t` estimate on a curved interface, not the
  corner eigensolution. Everything else in the conjugate plan is done.
- Profile + optimise the GPU step for the 2:1-refined channel — **STEP 1
  (re-profile) DONE 2026-08-07, and step 2 has since been done many times
  over: see the performance-campaign bullets below**, which supersede every
  number in this one. Kept for the two things it established. (1) The fresh
  split that set the target (A6000, 100 steps, min_channel refined
  128x64x8, niter=6 cheb): loop 31.3 ms/step = projection_exchange 44.5% +
  projection_jacobi 33.9% + syncface 7.9% + momentum 6.1% + projection_bc
  4.2%; **all halo exchange 52.4%**, and on the heavier 128x64x128 refined
  case 42.4% with the projection's exchange (79.1 ms) costing as much as its
  Jacobi kernels (82.1 ms). That reading was right about WHERE and wrong
  about WHY — the 2026-09-11 kernel-timeline work showed the exchange's cost
  was per-LAUNCH overhead, not transport. (2) The nb lever: unrefined
  128x64x128, nb 8/16/32 = 48.6/39.3/36.4 ms/step (**-25%**) with the
  exchange share 44.6%/35.2%/35.1% — but the refined channel is PINNED to
  nb=8 (its 24-cell wall band is 3 block rows; nb=16 does not divide it), so
  bigger blocks are not the fix there. LANDMINE, still live: `nb` must divide
  the global grid in every direction, so an nb sweep is impossible on the
  nz=8 case. **The profiler it describes is GONE**: the `scalar` branch's six
  `STEP_PROF_*` buckets in chron.f90 and the optional `prof` argument to
  `pressure_projection` were superseded at the merge by `profiling.f90`'s
  three nested profilers (step/proj/exch), which the same `[output] profile
  = true` key drives; chron.f90 keeps only the generic `profiler_type` both
  used. `docs/next_session_profiling.md` carries a SUPERSEDED header.
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
- Rank-to-GPU mapping + campaign re-measurement (DONE 2026-09-10, branch
  `optimiseBlockRefinement_parentBoundaryLayer`). The 708-753 us `mpi_wait` that
  three reports called "the blocked exchange's stall" was GPU AFFINITY: **the two
  ends of a cross-node link must sit in the same affinity class** (mixed 604 us,
  both-far-but-matched 81, both-near 75 — matching is worth ~8x, NIC-affinity a
  further ~1.3x). `comm.f90 select_target_device` now serves a node's first and
  last local ranks first, identically on every node (node-homogeneity is
  load-bearing: a per-node order cost base_jacobi 14%); `MOBY_GPU_ORDER` overrides
  the device order; single-node runs keep the identity mapping. Campaign
  re-measured, 46 runs in one allocation: 1/2/4 ranks and base_jacobi move ±0.3%
  (the control), while at 8/16 ranks rect +19.4/+24.5%, refined_yp82
  +21.4/+26.8%, red-black +25.3/+29.1%, refined_big +13.6/+20.7%. **Block tax
  1.310 → 1.059 (8 ranks), 1.499 → 1.130 (16) — there is no node-boundary block
  tax**, and `results_horeka_2026-09-07.md` carries a SUPERSEDED-IN-PART header.
  The 2:1 headline is untouched (0.98 coarse-cell-equivalents at 4 ranks).
  `tools/moby_tune.sh` finds the mapping by measurement on a new machine (it
  recovers the affinity classes with no topology input); `tools/h5maxdiff` does
  field comparison where h5py/h5diff are unavailable. The owed bit-exactness gate
  is DISCHARGED (2026-09-10, job 5139976): the 7-case suite runs on HoreKa
  (`overheadTest/horeka/exchange/run_mapgate.sh`) against the pre-change
  `build_gpu/moby_solve.ref`, all seven at max_abs 0 incl. every RANS scalar,
  deliberately WITHOUT nofma (the arithmetic source is byte-identical, so a
  production-flag comparison is strictly tighter). Not covered: `les_ibm` +
  `refine_body`, whose `IC_refine.h5` is generated by setup.sh, not committed. Diagnostics from this track:
  `exchange balance:` line and `[output] exchange_barrier` (skew vs transfer).
- The 2:1 exchange at scale — ATTRIBUTED (2026-09-10, jobs 5139461/5139581,
  `overheadTest/results_horeka_exchange_2026-09-10.md`). The "1.9x exchange per
  cell" is **not volume of any kind**: the refined case exchanges 0.954x the
  points PER CELL of its single-level twin, and cross-level points are
  **16.43% of the total at 1, 4 and 8 ranks alike** (16.38% on a case with 4x the
  leaves) — the Morton split does not cut the interface. Nor is it a per-point
  penalty: the new `exch_timing copy_cross` bucket (the cross-level kernel timed
  apart from `local_copy`, which now means SAME-LEVEL only and is 0 on a
  single-level grid) puts a cross-level point at 0.114–0.115 ns against a
  same-level point's 0.097–0.107 — **1.15x, not 2–4x**. The cost is **kernel
  LAUNCHES**: each carries a fixed cost before it moves anything (same-level copy
  56–67 us, cross-level copy 84.6–87.3, pack ~92–114, unpack ~83–100, against the
  projection's `jacobi_compute_phi` at 19), and there are ~141 per step (39 pack +
  39 unpack + 39 same-copy + **24** cross-copy) = 9.7 ms/step for `rect_jacobi`
  and 11.7 for `refined_yp82` — 80% of the refined case's device-local exchange
  and the entire difference between the twins, who pay the same launch bill on
  138 M and 61 M cells. `base_jacobi` with zero local copy points reads **0.1 us**,
  which is what makes it per-LAUNCH rather than per-round. The 2:1 interface's own
  exchange cost is 2.03 ms/step of launch against 0.36 ms of transfer — **85%
  launch, 15% data** — and a THIRD launch cost sits outside the exchange:
  `jacobi_compute_phi` launches at 24-26 us whatever the mesh (the cheapest launch
  in the solver, 43-46 ps/cell everywhere) while `jacobi_apply` is 46 us
  single-level and **137 us refined**, because it launches `interface_correct`
  whenever interfaces are present and that kernel does nb^2 work against the
  sweep's nb^3 — 91 us x 18 calls = **1.63 ms/step**. Refinement-specific total
  **3.66 ms/step = 6.1% of the refined 8-rank step, almost all launches: that is
  the 2:1 tax at scale.** CONFIRMED AT 16 RANKS (job 5139977): cross-level 16.43%
  at a fourth rank count, `rect` 13.18 and `refined_yp82` 13.23 ms/step of
  device-local exchange on 138.41 M and 60.56 M cells, and the 4/8-rank fits
  predict the 16-rank projection kernels to 0.3% and the exchange kernels to 3.4%
  without refitting. **31.8% of the refined 16-rank step is launch-fixed cost**
  against 1.7 ms/step of actual halo data movement. Also: **A0 re-run at 8 ranks fails again** (5.7 ms of
  compute between the posts and the Waitall leaves `mpi_wait` at 1.06–1.44x of
  baseline and adds exactly 39 x 5.7 ms to the step) — **P2/overlap is CLOSED**
  and the 2-rank probe's scope caveat is discharged. MECHANISM NAMED (2026-09-11, job
  5141872, `overheadTest/results_kernel_timeline_2026-09-11.md`): **every launch of
  an exchange kernel copies the ENTIRE `comm_type` object host-to-device — 7952 B,
  all 45 of its array descriptors — plus 14-23 separate descriptor/scalar copies**,
  whatever the kernel names. 26 rank-1 descriptors x 152 B + 19 rank-2 x 200 B =
  7752 B accounts for the block; the trace's other copy sizes (8/152/200 B) are
  scalars and the two descriptor ranks. No other kernel pays one — the projection
  kernels move 8-148 B and cost 5-22 us against the exchange kernels' 62-88. The
  arrays are ALREADY resident from `init_block_exchange`'s `target enter data`, so
  none of the ~10 kB per launch is data the device needs. **Recoverable
  7.8-9.4 ms/step = 19-23% of the refined 16-rank step.** FIXED (2026-09-11, job
  5142027, same report §7): **`comm_type` was the ONLY derived type in the solver
  that did not map its PARENT object** — `blk`, `bf`, `bc`, `g`, `ibm`, `sst`,
  `turb` all do, it is the convention this file states, and comm.f90 was the one
  place it was missed, so there was no device copy of `c` for the components to
  attach into. One line each way (`map(to: c)` before the component maps,
  `map(delete: c)` after), plus a local copy of `activeVars` — the ONE thing a
  kernel reads from `c` that changes per call, which a resident `c` would
  otherwise leave stale. Measured: H2D per exchange-kernel launch **14–23 copies
  / ~10 kB → 0–1 copies / 0–16 B**; per step 3484 copies and 1.13 MB → 456 and
  ~0; per-launch fixed cost pack **88.9→18.3** us, unpack 78.5→18.0, same-level
  copy 62.8→14.6, cross-level 81.6→18.3 — every exchange kernel now at the floor
  the projection kernels already had, with the two `c`-free controls moving 0.8 %.
  Step at 4 ranks: refined **−8.36 %**, single-level **−3.37 %**, against 8.51 and
  6.62 ms/step predicted from the per-launch table. BIT-EXACT (max_abs 0) on both
  production cases and all 7 suite cases incl. every RANS scalar. CONFIRMED AT SCALE
  (job 5142047, 4 nodes, both binaries in one allocation): refined **−20.14 %** at
  16 ranks and −13.84 % at 8; single-level −9.21 % and −6.35 %. **The absolute
  saving is rank-independent** — 8.64/8.34/8.34 ms/step refined at 4/8/16 ranks,
  6.97/6.94/6.21 single-level — a per-launch constant times a fixed launch count,
  landing on whatever the step is. Device-local exchange at 16 ranks **13.15 →
  4.82 ms/step (31.8 % → 14.6 % of the step)**, and the 2:1 per-cell tax at 16
  ranks falls **1.403 → 1.234** as a side effect (removing a per-LAUNCH cost helps
  the case with fewer cells per launch more). MATRIX RE-RUN (2026-09-14,
  job 5142973, `results_horeka_2026-09-14.md`; `results_horeka_2026-09-10.md` §9
  now carries a SUPERSEDED-IN-PART header). Its `ref` column reproduces the
  2026-09-10 tables on different nodes, so the `new` column is the fix alone.
  Worth at 16 ranks: base +8.7%, rect +11.2%, refined +18.9%, red-black +23.1%,
  refined_big +9.2%; absolute saving a constant 6.5-9.9 ms/step from 2 ranks up
  and 1.6-3.6 ms at ONE rank (no peers ⇒ only the 39+24 copy launches ever paid
  the blob, not all 141). Strong scaling at 16 ranks base 71→78, rect 66→74,
  refined 49→60, red-black 44→56, refined_big 80→86%; block tax 1.130→1.100.
  **TWO PUBLISHED CONCLUSIONS REVISED: the 2:1 machinery costs 1.004
  coarse-cell-equivalents at 4 ranks, not 0.98** (the added fine cells cost what
  the coarse cells beside them cost, to 0.4%; the old value was below 1 only
  because the single-level twin carried more of the per-launch bill), **and
  red-black erodes to 0.808 at 16 ranks, not 0.852** (0.76→0.81, half the erosion
  published). `mpi_wait` did not move but its share grew to 10-12% of the smaller
  step — the largest single exchange item again, and the one neither partitioning
  nor overlap can touch. Largest launch count left is the
  projection's (`interface_correct` 3 kernels × 18 calls + `jacobi_apply` 36),
  already at the per-launch floor — a kernel-count question worth ~1 ms/step.
  **`jacobi_apply` is now the biggest single item at 33.1% of the post-fix refined
  16-rank step** (34.4% single-level), and ncu says why
  (`results_ncu_apply_2026-09-14.md`, job 5144931): it is **occupancy-limited, not
  wasteful.** Traffic is near-minimal everywhere (compute_phi 1.10x, apply k1
  1.07x, k2 1.08x of the source-counted minimum -- published as 1.35x until
  2026-09-14, when the minimum was found to count ibm%mu as one array rather than
  its three staggered components) and load/store sectors-per-request
  match the control, so there is NO coalescing defect and no wasted bytes. apply
  costs 2.87x compute_phi because it moves 1.98x the bytes it genuinely needs AND
  runs at 0.83x the efficiency, and the efficiency gap is registers: k1 at **59
  regs → 45.3% occupancy → 64.8% of peak DRAM**, k2 at **88 regs → 29.6% → 44.0%**,
  same grid/block, zero local memory in both. Levers, sized: cut k2's registers to
  ≤64 (~2.5 ms/step at 16 ranks, 7-8%) and fuse k1 into k2 (~1.0 ms/step, 3%).
  CHECK ANY REGISTER CHANGE against `launch__registers_per_thread` and
  `sm__warps_active`, not a step time — cutting registers by spilling looks like
  progress and loses on the clock.
  **THE REGISTER LEVER IS TAKEN (2026-09-14, jobs 5145099/5145100/5145114,
  `results_apply_registers_2026-09-14.md`): k2 is 88 -> 64 registers and the time
  followed.** `face_grad_corr` depends on the block and the FACE-NORMAL index
  alone, so evaluating it inside the `collapse(4)` body recomputed it nb^2 times
  and -- the part that cost something -- kept ten array bases live for every
  thread. It is now precomputed into the module arrays `cfLow(idx,d,b)` (low
  faces, carrying the interface zeroing) and `cfHigh(d,b)` (the outlet high
  face), which are STATIC unlike `rdenom` (face kinds from the leaf table,
  metrics from the node lines) and so are formed ONCE on the host and mapped
  once. `STACK`/`LOCAL` stay 0. ncu, both binaries on one node: occupancy
  **29.6 -> 46.0%**, DRAM **42.5 -> 60.5% of peak**, k2 **7565 -> 5309 us
  (-29.8%)** with the traffic unchanged to the last digit printed (10.31
  doubles/cell) -- the same bytes, issued by 1.56x the warps. `proj_timing:
  apply` **-21.8% / -20.1%** and the STEP **-8.6% / -7.8%** at 4 and 8 ranks on
  both production cases, with `sweep` an unmoved control (+1.0 to -0.3%) and two
  separate allocations agreeing to 0.1%. Unlike the `map(to: c)` fix this saving
  is PER-CELL, not per-launch: it halves from 4 to 8 ranks (17.5 -> 8.8 ms
  `rect`, 7.7 -> 4.0 `refined`) while staying ~8% of the step. Bit-exact by
  construction and measured so: Pass G on both production cases (138M / 60M
  points) and the 7-case suite at max_abs 0 incl. every RANS scalar, CPU and GPU,
  deliberately WITHOUT nofma. L1 (splitting the high-face planes into their own
  kernel) was NOT taken and should not be -- it costs a launch and the threshold
  it was for is already crossed. AT 16 RANKS (job 5145120, 4 nodes, `new`
  built from a PINNED worktree so later work could not contaminate it): step
  **-8.2% / -6.5%**, apply **-21.3% / -18.7%**, its share **34.7 -> 29.7%** and
  **33.9 -> 29.5%** -- and the forecast written before that job ran was right to
  0.3 ms and 0.7 percentage points, with its stated caveat firing in the stated
  direction (the refined case's fractional gain softens from 20.9% at 4 ranks to
  18.7% at 16, as 3.8 Mcell/GPU stops filling the machine). The same logs settle
  what `results_horeka_2026-09-14.md` section 6 said was owed: apply at 16 ranks
  IS 34.7% / 33.9% of the step, as inferred. **THE SAME LEVER, TAKEN ON THE OTHER TWO
  PROJECTION KERNELS (2026-09-14, jobs 5145507/5145518,
  `results_rdenom_registers_2026-09-14.md`): `compute_rdenom` 110 -> 80,
  `jacobi_compute_phi` 94 -> 80.** `face_grad_denom` and the divergence metric
  `d1?(idx,VAR_P,b)` are static in exactly the same way, so they join the tables
  as `dnLow`/`dnHigh`/`d1P` (d1P MUST stay a separate factor: the diagonal is
  `(dnLow*mu + dnHigh*mu)*d1P` and folding it in would distribute the multiply).
  ncu: rdenom **-27.5%** (DRAM 28.9 -> 39.8% of peak, occupancy 23.7 -> 35.0),
  phi **-14.1%** (53.8 -> 62.3, 29.1 -> 35.0), the two apply kernels unmoved to
  0.1% as controls. `sweep` **-12%**, `setup` **-19%**, STEP **-2.2 to -2.7%** at
  4 and 8 ranks (16 not measured; scaling gives ~-2%). FOUND BY THE REGISTER
  TABLE, NOT THE PHASE TABLE: **`compute_rdenom` had never been profiled** and
  was running at 28.9% of peak DRAM -- the worst in the solver -- because at 3
  calls/step it never ranked in a bucket. It is now COMPUTE-limited (SM 49.6% >
  DRAM 39.8%): the remaining cost is very likely its per-cell fp64 divide, the
  one the reciprocal-rdenom change removed from compute_phi for -50.5% of that
  kernel. That is the next thing to look at, with ncu before any code.
  **GATE LANDMINE, and it bit here:** this change MOVES expressions, so the
  compiler fuses `a*b+c` differently on the two sides and the PRODUCTION-flag
  comparison fails at 1e-15..1e-12 with nothing wrong. It needs
  `-Mnofma -gpu=nofma` on BOTH sides -- which this machine had no build for:
  `compile.sh` gained `cpu_nofma`/`gpu_nofma` modes (production modes byte-
  identical to before) and `horeka/exchange/submit_nofma_gate.sh` runs Pass G +
  the 7-case suite with them: all 9 comparisons **max_abs 0**. Use the
  production-flag gate ONLY when the arithmetic source is byte-identical between
  the two binaries (a device index, a launch parameter); use the nofma gate for
  anything that moves an expression.
  **WORK THAT DID NOT NEED DOING (2026-09-14, jobs 5145541/5145549/5145554,
  `results_stepwork_2026-09-14.md`): a further -8.0 to -8.2% of the step, larger
  than either register increment, because these mostly make kernels NOT RUN.**
  (1) `update_ibm_mu` ran an fp64 DIVIDE per ghost-inclusive cell x 3 components
  x 3 substages to compute `mu = 1/(1+dt*0) = 1` on every body-free case --
  channels, boundary layers, Beltrami. One device reduction on the first call now
  decides whether any coef is non-zero ON THIS RANK (local is correct: mu is
  pointwise, so the fields are bit-identical either way) and the kernel never runs
  again: **-99.5%** on the bucket. (2) The trip force refreshed the WHOLE domain
  every substage although its own `ex < -50` cutoff confines the envelope to under
  1% of it: now the v component only (f_u/f_w are always zero and already zero
  from allocation), over a block list built once at init (**40 of 256 blocks** on
  the rank that owns the trip -- the solver prints it, treat it as a gate), and
  the spanwise Fourier sums g_k(z) -- which depend on z ALONE and change only when
  the random walk redraws, once per trip_ts -- are tabulated per (k, block)
  instead of evaluated per cell: **-88.6%** on the bucket. Both bit-exact at
  PRODUCTION flags (no expression moves), 9 comparisons per job at max_abs 0.
  (3) **`step_momentum` IS occupancy-limited -- and the register lever does NOT
  work on it.** First ncu on the solver's largest kernel: 128 registers, 23.96%
  occupancy, 32.50% of peak DRAM, traffic 1.05x its source-counted minimum,
  against its own sibling (same grid/block, 66 regs) at 41.09% and 74.22%.
  Packing its NINE Laplacian coefficient arrays into three (`lapX(LAP_M/LAP_0/
  LAP_P,i,var,b)`, stencil first) changed the register count **by nothing** --
  for this kernel the array bases are not the binding constraint, the ~30 live
  values of three fused component blocks are. And `-gpu=maxregcount:102` reaches
  102 only with `STACK:16`, while GLOBALLY it acts as a TARGET, not a cap, raising
  every kernel toward it (`jacobi_apply` k2 64 -> 96, `compute_rdenom` 80 -> 100,
  `interface_correct` 54 -> 80) and undoing the day's work. The pack was kept only
  because it was measured separately: momentum **-3.0%** from CONTIGUITY alone
  (registers and occupancy unchanged, DRAM 32.50 -> 34.52). NEXT for momentum is a
  different idea -- splitting the fused predictor into three per-component kernels
  -- not another register trick.
  **RE-RUN THE CAMPAIGN MATRIX before quoting any ratio: ~15% has come off the
  step across 2026-09-14 and every published number in
  `results_horeka_2026-09-14.md` has a moved denominator.**
  Probe: `horeka/exchange/run_ncu.sh` (ncu needs
  `--bind-to none`; `--page raw` is WIDE; the solver's stdout lands in ncu.csv). `horeka/exchange/analyse_launch_traffic.py`
  reproduces the per-launch traffic table from an nsys sqlite export. Cheap concrete win meanwhile: fold the local copy into the pack
  kernel (independent, both before the `Waitall`, 2.3 ms/step, bit-exact).
  Recommended AGAINST: fusing the cross-level kernel into the same-level one — it
  saves 2.0 ms of launch but puts the interface gather's ~128 registers on all 39
  rounds. `docs/next_session_profiling.md` and the Phase-4 sketch in
  `docs/nonblocking_overlap_strategy.md` predate all of this.
- **The between-iteration velocity exchange, reduced to what the divergence
  reads (DONE 2026-09-15, jobs 5145805/5145806,
  `results_divhalo_2026-09-15.md`): -3.0 to -4.9% of the step.** 15 of the 18
  velocity halo rounds per step sit BETWEEN projection iterations, and the only
  velocity halo anything reads before the next divergence is `q(nb+1)` of the
  component NORMAL to that face. The solver already exploited that at ONE rank;
  with peers it fell back to the full 26-direction three-component shell because
  the entry list was not partitioned for it. `entry_round` now gives the
  enumeration a THIRD round -- pure `+axis` same-level face COPY entries first,
  then the remaining same-level copies, then the cross-level ones -- so a
  divergence round is a per-peer **prefix of the copy prefix** (the `copyOnly`
  trick one level down), derived identically on both ends with no negotiation,
  carrying ONE variable per entry (`lDivVar`/`sDivVar`/`rDivVar`, the face
  normal). `dsSlot` and its three kernels are GONE: the same-rank half is that
  same prefix, one kernel, peers or not. COMPLETENESS IS FROM THE ENTRY LIST,
  not from a gate: `entry_boxes` gives `off(d)=+1` the plane `nb(d)+1`,
  `off(d)=-1` the plane 0, and `off(d)=0` reaches `nb(d)+1` only through the
  tangential extension -- i.e. only at HALO indices of the other dims -- so
  inside the range the divergence reads, the pure `+axis` FACE entries are the
  only writers; and the cross-level ones write nothing there anyway
  (`interface_normal_dim` marks that face the low-side block's own and unpack
  skips it). Measured: wire volume **6.46x/6.36x** smaller, device-local
  **6.20x/6.14x** (the entry arithmetic predicts 6.24x at `nb = 64 44 48`);
  `proj vel_exchange` **-52.8/-53.5%** (`rect` 4/8 ranks) and **-46.0/-46.7%**
  (`refined`); STEP -3.21/-4.55% and -2.93/-4.93%. **The launch count does NOT
  fall** (39 pack + 39 unpack + 39 local copy + 24 cross, both sides) -- the
  saving is bytes and wait, and the biggest single bucket is `local_copy`
  (-27%), not the message, which is why it is worth something at one rank too.
  FOUND WHILE GATING, unexplained: the REORDERING moves `phi_exchange`, whose
  volume is untouched, by **+7.3/+9.6% at 4 ranks and -1.7/-2.9% at 8** -- the
  scalar exchange's scatter locality changed (the copies used to run all 26
  directions of block 1 then all 26 of block 2; they now run three faces of every
  block then the other 23 of every block). It eats 7-15% of the 4-rank gain and
  is why those two fell short of the pre-registered band. Gates at PRODUCTION
  flags (nothing moves, values are copied not recomputed): Pass G on both
  production cases (138M/60M points) and the 7-case suite, all `max_abs 0` incl.
  every RANS scalar; `min_channel` 1 == 2 == 3 == 4 ranks EXACTLY (four peer
  topologies -- the load-bearing gate, since it is what the single-rank reduced
  path was validated against); `validation/redblack_interface` 1 and 4 ranks
  `max_abs 0` (red-black does not use this path, but the reordering changes its
  wire layout). LANDMINE: `tools/h5maxdiff` is a BUILD PRODUCT, not tracked, so a
  pinned worktree lacks it -- and `run_exchange.sh` deletes Pass G's snapshots
  whether or not the comparison ran, so a missing comparator LOSES the gate
  rather than merely failing to report it (`submit_divhalo.sh` now builds it).
  **ITS REAL SHAPE, from the 3-column matrix (2026-09-17, job 5147466,
  `results_horeka_2026-09-17.md` §4): the gain GROWS with rank count --
  -1.0% at 1 rank, +2.7/+3.3% at 4, +4.5/+5.1% at 8, +8.6 to +12.0% at 16 --
  and the change is a NET LOSS at one rank and for RED-BLACK at every rank
  count (-3.8 to -0.8%).** Two ranks counts were not enough to see that. The
  cause is the half of the change that is not the saving: making the divergence
  set a PREFIX reorders the entry list, and that costs every OTHER exchange.
  Measured cleanly at 1 rank, where no MPI confounds it: `phi_exchange`
  **+11.9% (`rect`) / +24.2% (`refined`)** against **-1.7% on `base_jacobi`,
  which has ONE block and so nothing to reorder** -- same rank count, same
  binaries, entry count the only variable. Red-black never runs the reduced
  round (its `mpi_wait` is unmoved, 96.2 -> 97.4 us at 16 ranks), so it pays the
  reordering and collects nothing: `proj vel_exchange` +23.1% and `local_copy`
  +30.7% at 4 ranks. THE FIX, measurement-justified and not yet taken: keep the
  two-round enumeration and drive the divergence round from explicit index lists
  (`lDivEnt`/`sDivEnt`/`rDivEnt` + point prefixes) instead of a prefix -- ~6
  extra integer arrays, the kernels otherwise unchanged. That design was
  considered first and rejected on "more state" grounds; the measurement
  overturns them. **TAKEN (2026-09-17, job 5149889,
  `results_divlist_2026-09-17.md`), and it works.** The enumeration is restored
  EXACTLY as it was and the subset is compacted out of the finished lists
  (`?DivEnt`/`?DivVar`/`?DivOff`/`?DivPt`, `peerSend/RecvDivOff` now compaction
  outputs); both ends still select the same entries because the predicate is a
  pure function of the op and the direction. It is also SIMPLER than the prefix
  where it counts: the lists inherit peer-major order, so a peer's points are
  one contiguous range and `find_entry` is gone from all three divergence
  kernels. Measured, three columns in one allocation at 1/4/8 ranks: **red-black
  back to `3c2903a` within -0.07/+0.04/-0.10%** (its `proj vel_exchange` 15.982
  -> 15.984 ms at 1 rank against the prefix's 21.636), **Jacobi +4.2/+5.1%
  (`rect`) and +3.8/+5.2% (`refined`) at 4/8 ranks vs `3c2903a`**, i.e.
  +0.6 to +0.9% over the prefix; `phi_exchange` at 1 rank back to +0.1/+0.2% of
  base against the prefix's +12.6/+24.4%. Bit-exact both ways: 9 comparisons at
  `max_abs 0` vs the prefix (Pass G + the 7-case suite incl. every RANS scalar)
  and `max_abs 0` vs `3c2903a` on a RED-BLACK case, which is what proves the
  enumeration is restored rather than merely similar. ONE PRE-REGISTERED
  PREDICTION MISSED: the 1-rank step was predicted to become a gain and came out
  at PARITY (-0.02%/+0.21%) -- at one rank the old `dsSlot` path was already
  good, and three specialised plane-copy kernels against one generic index-list
  kernel is a wash. **The campaign ratios in `results_horeka_2026-09-17.md` were
  measured with the PREFIX and are now slightly pessimistic** (Jacobi) to 1-4%
  pessimistic (red-black); 16 ranks was not re-measured.
  **FULL SWEEP MEASURED (2026-09-23, job 5150030, `results_horeka_2026-09-23.md`),
  and it CORRECTS one claim.** Whole divergence work `3c2903a -> b9414bd`:
  `rect` +0.0/+3.9/+4.2/+5.2/**+13.5%** at 1/2/4/8/16, `refined` +8.2% at 16,
  `refined_big` +9.3%, and **red-black -0.0/-0.1/-0.0/-0.4/+0.5% -- neutral at
  EVERY rank count, measured rather than chained**, which is the whole point of
  the index list. CORRECTION: the index list's OWN contribution at 16 ranks is
  NOT measurable (-1.9 to +0.6% across five configs) -- the claim that it "adds
  a little at 4-16 ranks" holds for 1-8 (+0.8 to +1.3%) and not for 16. The
  control that says so: the `ref -> mid` column repeats job 5147466's comparison
  in a different allocation and lands within **3 percentage points** at 16 ranks,
  so **16-rank differences below ~3% are not resolvable in one allocation of
  this matrix**. REVISED: block tax **1.015/1.034/1.024/1.000/1.007** (was
  1.036/1.053/1.038/1.010/1.033) -- the blocked single-level case now costs
  essentially nothing at every rank count; strong scaling at 16 ranks base 80%,
  rect 80, refined 59, redblack 52, big 88; the 2:1 machinery 1.001x
  coarse-cell-equivalents at 4 ranks and 0.879x at 16.
  **THAT PASS IS DONE (2026-09-25, job 5159149, `results_horeka_2026-09-25.md`):
  three columns, six configs (rough_jacobi joins run_matrix.sh), 79 runs.**
  The rdenom work is worth **+3.0 to +3.5%** on body-free cases and **+1.9 to
  +2.8%** on the rough-wall body case, **essentially FLAT in rank count**. Its
  CONTROL column (`ref -> mid` = HEAD with the narrowing disabled) earned its
  place: it is **not zero but a systematic -0.4 to -0.6%** on every config that
  calls `compute_rdenom`, and **exactly 0.0 on red-black, which does not call
  it** -- so the roughness feature is dormant as claimed AND the `rdenomBlocks`
  indirection costs ~0.5% while the kernel runs. Quote `ref -> new`, not
  `mid -> new`, which that 0.5% flatters. The indirection is not a residual to
  fix: on a body-free rank in HEAD the list is empty and the kernel never runs.
  `rough_jacobi` classifies **exactly 25% body blocks at every rank count** and
  its gain tracks the pre-registered `1 - (body/all)` model to three decimals at
  2 and 4 ranks (ratio 0.739/0.750), drifting to 0.51 at 16 -- **do not
  extrapolate the model to 16 ranks**. The roughness itself costs **+6.2 to
  +6.9%** of the step against its body-free twin (more in HEAD than in `mid`
  because the narrowing speeds the twin up more -- same overhead, smaller
  denominator). **METHOD CORRECTION to the 09-23 claim that "16-rank differences
  below ~3 points are not resolvable": that is true of ONE pair across TWO
  allocations. A pattern consistent across configs and rank counts INSIDE one
  allocation resolves well under 1% -- the control above pins 0.5%.** LANDMINE
  designed around: `config.f90` has NO `case default`, so an unknown key in a
  known section is SILENTLY IGNORED -- a pre-`7b2bc2a` binary runs
  `rough_jacobi` with `wall_shape` discarded, falls back to the 2D wavy wall and
  reports success, so the old column must be gated with `CONFIGS`.
  NOT ATTEMPTED, and the only idea left here: a reduced round for RED-BLACK
  itself. `redblack_sweep` iterates `0..hi`, sweeping the lower halo layer
  redundantly, so a cell at `i=0` reads the low halo plane of ALL THREE
  components where Jacobi's divergence reads the high plane of one -- its
  reduced set is ~12 of 18 face-component units against Jacobi's 3, a ~1.4x cut
  needing its own read-set proof (a cell at `(0,hi2,k)` reads an EDGE,
  `q(0,hi2+1,k,V)`).
- **`rdenom` is STATIC on a body-free rank (DONE 2026-09-17, job 5150041,
  `results_rdenom_static_2026-09-17.md`): +2.8 to +3.4% of the step measured at
  100 steps, 4.1-4.9% asymptotically.** The `rdenom` comment said the metric
  tables are static "unlike rdenom: rdenom follows `ibm%mu`, which
  `update_ibm_mu` rewrites every substage". **That expired on 2026-09-14**, in a
  DIFFERENT increment of the same day: `update_ibm_mu` returns immediately on a
  rank with no body, so `mu = 1.0` for the whole run and `rdenom` is as static as
  the tables it is contrasted with -- it was being recomputed 3x a step to get
  the same answer. New `ibm_mu_is_unit(ibm)` exposes the cached flag and the
  projection forms `rdenom` once (`rdenomStatic`); a rank WITH a body keeps
  recomputing, because there `mu` really does follow `dt_gamma`. Bit-exact by
  construction (same inputs, same expression); gates `max_abs 0` at production
  flags incl. Pass G and `les_ibm`, which is the case WITH a body and so the
  correctness half. **This is task 2 of the handout but NOT the change it
  proposes** -- that one attacks the per-cell fp64 divide, on the expired
  premise. The divide is untouched.
  **EXTENDED THE SAME DAY TO CASES WITH A BODY, which are COMMON in production
  (user, 2026-09-17 -- the first version helped only body-free cases and was
  aimed at the wrong half).** `mu = 1/(1+dt*coef)` is EXACTLY 1.0 wherever coef
  is zero whatever dt does, so the dt-dependence is confined to blocks that
  actually hold coefficients: `ibm_body_blocks` (one device reduction per block,
  once per run) replaces `ibm_mu_is_unit`, the first projection fills every block
  and then narrows `rdenomBlocks` to the body ones. The win is geometric, not
  binary -- `1 - (body blocks / all blocks)` of a 4-5% bucket: body-free 0/N
  (never runs again); `les_ibm`, plane walls spanning the domain, **256/640**
  (setup 8.756 -> 4.003 ms/step, **-54%**); `sailplane` at nb=10, a compact body
  in a large domain (the airfoil shape), **48/4500**. CAVEAT, not a defect: with
  `nb` UNSET there is one block per rank, the body touches it, and there is
  nothing to narrow (sailplane reads 1/1) -- the gain needs block granularity,
  which production cases set anyway. The fraction is PRINTED at init like the
  trip force's block list, because a silent 100% looks identical to a silent 0%.
  The handout's divide question now survives only for cells inside body blocks.
  **MEASURED ON A PRODUCTION-SHAPED BODY CASE (2026-09-19, job 5151978,
  `results_rough_2026-09-19.md`)**, which needed a new benchmark because every
  overheadTest config is body-free: `[ibm] wall_shape = eggcarton` adds the 3D
  sinusoidal roughness `h*sin(kx x)*sin(kz z)` (MacDonald, Chung, Hutchins, Ooi
  & Sandberg JFM 2017) and `configs/rough_jacobi.ini` is `rect_jacobi`'s grid,
  blocks and flow plus that wall. The `wavy` branch is textually unchanged so
  every analytic case stays bit-exact; `set_ibm_geometry` is applied by
  moby_solve AND moby_prepare (LOAD-BEARING -- the indicator drives coefficients,
  classification and wall distance, so one-sided application would make a
  prepared case file describe a different wall). Results at 200 steps: body-free
  `rect` step **-3.5%** (setup -82%), rough `rough_jacobi` step **-2.4%** (setup
  -60%) on **64/256 = 25% body blocks**. The pre-registered model -- saving
  proportional to `1 - (body blocks / all blocks)` -- predicted a 61.5% setup
  drop and 0.75x the absolute saving; measured 60.5% and 0.736. **The one-time
  residual is now MEASURED, not inferred**: `setup` after the change reads
  1.218 ms/step at 100 steps and 0.615 at 200, i.e. a constant 122 ms of
  first-call allocation, so the bucket fraction keeps improving with run length
  (-67% at 100 steps, -82% at 200) and production runs are thousands of steps.
  The case is a BENCHMARK, not validated physics -- height and wavelengths are
  chosen to be resolved and to sit in the first y-block, nothing physical should
  be quoted from it. **CORRECTED AT THE 2026-09-25 MERGE:** when this case was
  designed I wrote that the solver had no generic scalar transport and that the
  MacDonald forced-convection configuration would need its own track. That was
  true of THIS branch and false of the project — the `scalar` branch, now merged,
  carries `src/modules/scalar.f90` (passive scalars, `[scalar]`/`[scalar.N]`,
  conjugate heat transfer at the immersed interface). **The MacDonald &
  Hutchins rough-wall forced-convection benchmark therefore needs NO new
  physics**, only a case file combining `[ibm] wall_shape = eggcarton` with a
  `[scalar.N]` section: the natural follow-up to `configs/rough_jacobi.ini`.
  MEASUREMENT LANDMINE, and it bit: at 100 steps the bucket falls only 67%, NOT
  because the skip half-works (a CPU 10-vs-40-step run drops `setup`/step 3.90x
  with the TOTAL constant -- it runs ONCE PER RUN) but because what is left is
  **one-time allocation, zeroing and device mapping of `phi`/`delta`/`rdenom`**
  (~3.7 GB, 867 ms, on `rect` at 1 rank) that a 100-step benchmark charges to
  every step. Netted out, the kernel is 6.52 ms/call on 138.4 M cells = ~0.68
  TB/s = 44% of peak, matching the independently measured 39.8%. **A short
  benchmark UNDERSTATES any change that removes per-step work from a bracket
  that also holds first-call setup.**
- **The campaign matrix, three columns (DONE 2026-09-17, job 5147466,
  `results_horeka_2026-09-17.md`). SUPERSEDES `results_horeka_2026-09-14.md`
  sections 2-4 and every ratio quoted from it here.** 23 runs at each of
  `55bee89` (control, reproduces the 09-14 table), `3c2903a` (+ registers and
  step work) and `95312d7` (+ divergence halo), one allocation. Total off the
  step since `map(to: c)`: base +18.0 to +21.3%, rect +19.0 to +25.2%, refined
  +18.3 to +23.4%, refined_big +21.2 to +24.0%, **red-black only +5.4 to +7.8%**.
  The two increments run in OPPOSITE directions with rank count: per-cell work
  removed shrinks as the per-GPU problem shrinks (rect +19.8% at 1 rank ->
  +15.3% at 16), message volume removed grows (-1.0% -> +11.8%). REVISED:
  **block tax 1.036/1.053/1.038/1.010/1.033** at 1/2/4/8/16 (was
  1.053/1.095/1.070/1.043/1.100); strong scaling at 16 ranks base 82%, rect 82,
  refined 65, redblack 57, big 89; the 2:1 machinery **1.002x** coarse-cell-
  equivalents at 4 ranks, 0.986x at 8, **0.900x** at 16 (was 0.826).
  **THE LARGEST REVISION: red-black against Jacobi went 0.804 -> 0.993 at 16
  ranks -- red-black was ~20% faster and is now at PARITY.** Nothing was done to
  red-black; Jacobi got faster and it did not.

- **Branch consolidation (DONE 2026-09-25, job 5163314,
  `docs/next_session_merge_to_main.md` STATUS header).** `boundaryLayer` and
  `scalar` are merged into this branch; `multiGPU` and `claude/blocks` hold
  nothing it needs; `claude/jacobi-interface` still does (six RANS/airfoil
  features, verified absent here key by key) and must NOT be deleted.
  TWO THINGS WORTH CARRYING FORWARD.
  (1) **The trip maths changed and the optimisation survived it.** The
  `boundaryLayer` port replaces the Schlatter & Orlu unit-rms Fourier
  coefficients with the exact CaNS/SIMSON form (random PHASES on a flat
  spectrum, mode 0 + nmodes/2 harmonics, amp/nmodes scaling, plus a steady
  realisation `trip_amp_s`) — and it keeps the same separable Gaussian envelope
  and the same `ex < -50` cutoff, which is the ONLY reason the 2026-09-14
  block-list optimisation still applies. `fill_trip_span` now tabulates three
  spanwise signals instead of two. The printed active-block fraction moved
  (4/16 on the gate case vs the old 40/256) because the default `trip_x0` moved
  15 -> 10, NOT because the envelope changed — check that print, do not assume
  it. Gate: `max_abs 0` vs their bodyforce.f90 verbatim, 1 and 4 ranks, with
  `trip_ts` small enough to redraw the walk ten times inside the run.
  (2) **A MERGE HAS TWO PARENTS AND THEY CAN DISAGREE ON PURPOSE.** `scalar`
  carries the 2026-08-05 fix for the cold-started RANS initial condition (`k` a
  factor 4 low on the last plane of every block, because `init_rans_transport`
  read a halo nothing had filled yet). So the merged head is `max_abs 0`
  against `origin/scalar` on turb180 / wf180_y30 and DIFFERS from the pre-merge
  head by 0.18 / 0.88 / 2.09 on the three RANS cases — a `max_abs 0` there
  would mean the fix was LOST, and `submit_merge_gate.sh` asserts the
  difference rather than tolerating it. Everything neither branch's physics
  touches is `max_abs 0` against the pre-merge head: min_channel 1 and 4 ranks,
  beltrami_slaby, les_ibm, and Pass G's two production cases at 138 M / 60 M
  cells (run with `trip_amp = 0`, since the trip maths changed).
  ALSO: the two branches had each built a per-phase step timer on the SAME
  `[output] profile` key; `profiling.f90`'s three nested profilers won and
  chron.f90's six `STEP_PROF_*` buckets are gone, with a new `PROF_SCALAR`
  bucket carrying over the one thing they had that it lacked.
  The 1e-14 velocity / 1e-13 pressure gap between this branch and `origin/scalar`
  on Beltrami-flow cases is NOT the merge — it is the 2026-09-02
  reciprocal-`rdenom` change, which gave up projection bit-exactness at exactly
  that magnitude; the two PARENTS compared directly reproduce it to the last
  digit.
- **Features ported from `claude/jacobi-interface` (2026-09-25, branch
  `port/jacobi-interface-features`).** Four of the six that branch held alone,
  plus the one the feature list missed. **`boostconv` and `kpin_dwall` were
  deliberately DROPPED** -- boostconv on its own dossier (V1 negative on
  turb180, best config 2.7x SLOWER than plain marching; V2's win INVALIDATED
  because the recombination suppressed the ktrip strip, k 7.6e-3 -> 7.9e-6),
  kpin_dwall because the instability it patched was addressed at its root by
  skew convection and no validated ini uses it.
  - **The skew lockdown was the one that mattered, and it was not on the list.**
    main carried a PRE-lockdown snapshot: the `[flow] convection` key still
    existed and still DEFAULTED TO DIVERGENCE, so any ini omitting it ran the
    form that branch found unstable at 2:1 interfaces. Skew is now hardwired and
    the key error-stops. **The snag: scalar.f90 read the SAME key for a
    DIFFERENT operator** -- its "skew" is the FULL subtraction, i.e. the
    ADVECTIVE form u.grad s, which preserves a uniform scalar for any advecting
    field but gives up the exact conservation `validation/scalar/conserve.ini`
    gates; the momentum kernel subtracts a HALF (energy neutrality). The
    scalar's choice therefore moved to its own `[scalar] convection =
    divergence (default) | advective`. Gate: 9/9 max_abs 0 against the toggle
    binary running `convection = skew`, with ibmwavy as the non-vacuous leg (the
    new key moves theta by 4.1e-10 between its settings, so max_abs 0 there is a
    real equivalence; conduction/prsweep/wave have no advection and gate the
    momentum half only).
  - `[rans] kpin_box` / `ktrip_box` (the OpenFOAM fvOptions forced-transition
    pair). Dormant bit-exact; the pin proven EXACT by pinning the whole domain
    and changing `tu` 1 -> 10 (un/vn/wn/pn/nut/k all max_abs 0; omega does
    differ, correctly, since kpin pins k alone and nut = 0 makes omega inert),
    against a control where the same tu change without the pin moves k by 2.72.
  - `[blocks] refine_body_levels` + `refine_body_box`: a CAP on body-driven
    refinement, raised locally by boxes that lift the cap rather than filling
    their volume, so the refinement follows the surface BAND (~1000 leaves on
    the NACA nose against ~12700 volumetric). Dormant identical (case-file
    datasets and the inline solve path both max_abs 0); cap 0 gives 0 refined,
    a half-domain box gives exactly half.
  - Runtime CONTROL-VOLUME forces replacing the penalization integral, plus
    `[case.airfoil] steady_tol`. **DELIBERATE DEVIATION:** the branch makes a
    missing `cv_box` an `error stop`, which leaves EIGHTEEN inis
    (validation/naca0012, validation/sd7003, tutorials/naca) unable to start --
    they are broken on that branch today for this reason. Here it is a loud
    per-run WARNING that disables force sampling instead. Boxes were NOT
    invented for them: the budget is sensitive to the per-face `p_inf`
    subtraction and to borders crossing a 2:1 interface, so an unvalidated box
    yields numbers nobody has checked. **That is the top follow-up.**
  CPU dormancy vs the pre-port binary: min_channel (4 ranks), beltrami_slaby,
  turb180, lam30t and conduction all `max_abs 0`, plus conduction's `s1` named
  explicitly -- **`tools/h5maxdiff` with no dataset arguments compares
  un/vn/wn/pn and the RANS scalars ONLY, so a passive-scalar case gated that way
  never looks at the scalar.** `les_ibm` is the one case not yet run.
  Also brought over: `tutorials/naca/rans` (the converged OpenFOAM comparison,
  C_L 0.5199 vs 0.5142, Cp_min matching to four digits) and the cv_forces /
  skew / naca docs. `claude/jacobi-interface` is now down to the naca LES
  kickoff, the boostconv module and the C10/C11 analysis history.

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
