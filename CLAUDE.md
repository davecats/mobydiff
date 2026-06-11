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
  `block_set_type` with one block per rank; `field_type` and the grid
  metric arrays are gone (`grid_type` keeps only generation parameters
  and the node lines; blocks slice from them via `slice_grid_direction`).
  Volume kernels (`step.f90`, `pressure_solver.f90`, `ibm.f90`, `les.f90`)
  loop `do b = 1, blk%nBlocks` folded into their collapse; `ibm%coef/mu`
  and `les%nut` carry the trailing block index. Still rank-shaped and
  serving block slot 1 until Phase 1: comm.f90 send/recv boxes, io.f90
  rank-box datasets, the `apply_bc` point list, channel statistics, the
  LES 1D metric tables, the `uStartX`-style face masks and the scalar
  `colorOffset`.
  - Validated bit-exact vs `6c03b67` with both trees built `-Mnofma`
    (CPU) / `-Mnofma -gpu=nofma` (GPU): channel_kmm180 restart 11 steps
    (CPU 8 ranks + GPU 1 rank), sailplane 1 step (GPU), wavychannel
    5 steps (CPU 8 ranks). With default FMA contraction the binaries
    differ by 1-2 ulps/step (compiler instruction selection on the
    re-indexed source, not arithmetic) — for "pure refactor" gates,
    compare `-Mnofma` builds of both sides.
  - GPU cost of the block index: none measurable (channel 50 steps,
    interleaved x5: 0.201 vs 0.204 s/step, noise ±8% on this WSL box).
  - `tutorials/sailplane/sailplane_ibm_coeff.h5` was corrupt in git
    (unreadable by HDF5 since at least `2b5e517`); regenerated per the
    tutorial README and recommitted.
- Phase 1 (next): many same-level blocks per rank, block-pair halo
  exchange entries (replace every "slot 1" above), Z-order distribution,
  per-block face masks and red-black `colorOffset` from
  `modulo(sum(blk%origin(:,b)), 2)`.
- Phase 2: removal of solid blocks (`FACE_CLOSED` zero-flux faces).
- Phase 3: 2:1 refinement (restrict/prolong in pack/unpack, fine-owns-face,
  finest-level buffer at the wall, per-level node lines via midpoint
  subdivision).
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
