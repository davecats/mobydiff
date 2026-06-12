# CLAUDE.md

Guidance for Claude when working in this repository.

## What this is

mobydiff: an incompressible Navier-Stokes solver. Second-order finite
differences on a staggered Cartesian grid (uniform or stretched per
direction), RK3 time stepping, coupled velocity-pressure red-black SOR
projection, volume-penalization immersed boundary method (IBM), optional
LES. Fortran + MPI (3D Cartesian decomposition, 26-neighbour halos) with
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

## Active work: block refactor (branch `claude/blocks`)

The plan is `docs/block_refinement_strategy.md` — read it first. Goal:
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
  - 3c: staggered flux matching, FACE_COARSE/FACE_FINE sweep masks.
  - 3d: geometry-driven refinement in mobygrid/mobygeom.
- Phase 4: performance (overlap, see `docs/nonblocking_overlap_strategy.md`).

If the branch `claude/blocks` does not exist yet, create it from the
current checkout (`codex/les`) before committing: `git switch -c claude/blocks`.

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
