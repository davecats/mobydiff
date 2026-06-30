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
    `MOBY_HALO_AUDIT=1` (hook in main.f90) audits every
    exchange-written halo cell against manufactured linear fields on
    the real layout (1.2M + 7.9M cells, 0 bad) — run it FIRST when an
    interface case misbehaves. Follow-up refactor (validated
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
  CAVEAT: LES is still **IBM-UNAWARE in practice** — validated for channel flow
  only; the LES<->IBM coupling (the `ibm_aware` solid-cell nut masking) was NOT
  exercised and is untested with block refinement / across the 2:1 interface.
  RESIDUAL (the next no-LES task): a small v'-only spike at COARSE-OWNS y-faces
  (the Phase-3c low-block-owns-face orientation asymmetry; the fine cells get their
  interface-NORMAL velocity by prolong-injection of the under-resolved coarse face
  value). Reflux-off shrinks it to ~17% localized v' excess; fine-owns faces are
  clean. See `docs/next_session_interface_normal.md`.
- NEXT (no LES): interface-NORMAL velocity treatment — replace the orientation-
  dependent coarse-owns(inject)/fine-owns(restrict) ownership with a two-sided
  symmetric conservative reconciliation so both faces behave like the clean
  fine-owns one. Handout: `docs/next_session_interface_normal.md`.
- Phase 4: performance (overlap, see `docs/nonblocking_overlap_strategy.md`).

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
