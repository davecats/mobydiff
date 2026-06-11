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
./compile.sh all          # builds build_cpu/ and build_gpu/
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

- Phase 0 (in progress): `block_set_type` container, one block per rank,
  bit-identical to current behaviour.
  Done so far (may be uncommitted in the working tree):
  - `init.f90`: `slice_grid_direction` extracted from `init_grid_direction`
    (shared metric construction, pure refactor).
  - `src/modules/blocks.f90`: `block_set_type`, Phase-0 builder, device
    mapping, `block_set_matches_grid` bitwise check, `subdivide_node_line`.
  - `CMakeLists.txt`: module added.
  - `main.f90`: temporary scaffold building + verifying + freeing the block
    set after `init_grid` (remove when the solver actually uses blocks).
  Next: commit, build, run a case to validate; then migrate solver kernels
  (`step.f90`, `pressure_solver.f90`, `ibm.f90`, `comm.f90`, `io.f90`) to
  read `blk%` arrays with an outer block loop (`collapse(4)`).
- Phase 1: many same-level blocks per rank, block-pair halo exchange,
  Z-order distribution.
- Phase 2: removal of solid blocks (`FACE_CLOSED` zero-flux faces).
- Phase 3: 2:1 refinement (restrict/prolong in pack/unpack, fine-owns-face,
  finest-level buffer at the wall, per-level node lines via midpoint
  subdivision).
- Phase 4: performance (overlap, see `docs/nonblocking_overlap_strategy.md`).

If the branch `claude/blocks` does not exist yet, create it from the
current checkout (`codex/les`) before committing: `git switch -c claude/blocks`.

## Verification

- Phase 0 must be bit-exact vs. the pre-refactor code: compare output
  fields with `tools/compare_fields.py`.
- Channel sanity: `tools/check_parabolic_channel.py`. IBM case:
  `tutorials/sailplane/`.
- Refinement phases: uniform-flow preservation across interfaces and global
  mass conservation to round-off (see strategy doc §11 for the full list).
- Never declare a phase done with failing builds or unverified results.
