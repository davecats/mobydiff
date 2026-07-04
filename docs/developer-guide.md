# Developer guide

This page orients a developer to the source layout, the core data model, and the coding and
GPU conventions. The design records for the block-refinement work are the
[design notes](index.md#design-notes--internal-history); the top-level `CLAUDE.md` is a
detailed engineering log of the refactor phases.

## Source layout

```
src/
  main.f90                 solver entry point (RK3 time loop)
  mobygrid.f90             serial grid / preprocessing tool (exports node lines)
  modules/
    config.f90             .ini parsing and validation
    init.f90               runtime/domain state shared across modules
    blocks.f90             block set: leaf blocks, Z-order ids, 2:1 refinement
    comm.f90               MPI halo exchange + 2:1 interface transfer entries
    boundary.f90           physical boundary-condition application
    step.f90               momentum predictor, corrector, body-force/LES kernels
    pressure_solver.f90    damped-Jacobi / Chebyshev–Jacobi projection
    ibm.f90                immersed boundary (volume penalization)
    les.f90                subgrid model (Smagorinsky / WALE)
    bodyforce.f90          spatially varying volumetric force f(x)
    io.f90                 HDF5 field / stats / restart I/O
    field_hdf5.c           low-level HDF5 hyperslab helpers
    gpu_runtime.f90        OpenMP target-offload data mapping
    chron.f90              timers
    flow_case.f90 / flow_case_base.f90   flow-case dispatch
    flow/
      generic_flow.f90     the "generic" case (IBM / external flows)
      case_config_helpers.f90
      channel/             channel case: initializer, profiles, statistics
```

Two executables come out of the build: `main` (the solver) and `mobygrid` (a serial tool
that writes the exact node coordinates an IBM preprocessor needs).

## Data model

The solver state lives in a **block set** (`block_set_type`). There is no monolithic field
type or global metric array: the grid (`grid_type`) keeps only the generation parameters and
the per-direction node lines, and blocks slice their local geometry from them.

- A block is an `nb × nb × nb` box of cells with a one-cell halo. Blocks are numbered along a
  **Z-order (Morton) curve** and split contiguously over MPI ranks.
- Volume kernels loop `do b = 1, blk%nBlocks` folded into their `collapse`; field arrays carry
  a **trailing block index** (e.g. `ibm%coef(...,b)`, `les%nut(...,b)`).
- Per-block face descriptors drive everything geometry-dependent: momentum start indices, the
  projection sweep window and Neumann terms, and the red/black colour offset are all derived
  from the block's origin and face kinds.
- **Face kinds** classify each block face: `FACE_OPEN`, `FACE_PHYS` (a physical boundary),
  `FACE_CLOSED` (a zero-flux wall from solid-block removal), and `FACE_COARSE` / `FACE_FINE`
  (a 2:1 interface). Consumers must test the *kind*, not treat it as an arithmetic 0/1 — a
  no-flux face is any non-open kind, while boundary conditions apply only to `FACE_PHYS`.

`comm.f90` builds one exchange entry per destination block × direction. Same-rank entries are
a single device copy overlapped with the MPI messages; off-rank entries form one message per
peer in a canonical order both ends derive independently. 2:1 entries carry an operation
(copy / restrict / prolong) and sample on the source side so the wire always carries
destination-point values.

## GPU programming model

The GPU path is **OpenMP target offload**, from the same source as the CPU build, guarded by
`#ifdef USE_OPENMP_OFFLOAD`:

```fortran
!$omp target teams distribute parallel do collapse(...)
```

- Derived types own **flat, contiguous allocatable arrays**, mapped to the device once in
  `enter_*_data` / `exit_*_data` routines (see `gpu_runtime.f90`, `blocks.f90`). Do not put
  allocatable components inside arrays of derived types.
- Prefer OpenMP. Switch a kernel to OpenACC only if OpenMP genuinely lacks a needed feature
  or leaves significant performance behind — and say so explicitly when you do.
- One MPI rank drives one GPU; GPU builds require a GPU-aware MPI so device halo buffers can
  be sent directly.

## Coding conventions

- **Readability first.** Performance matters, but the code must stay easy for humans to read.
  The I/O layer (`io.f90`, `field_hdf5.c`) is the only place allowed to be ugly if it must be.
- Avoid duplicated code and interfaces more complex than strictly needed.
- Comment **intent and non-obvious choices**, not syntax.
- Match the surrounding code's naming, idiom, and comment density.

## Verification discipline

- A pure refactor must be **bit-exact** vs. the pre-refactor code — build both sides with
  `-Mnofma` (CPU) / `-Mnofma -gpu=nofma` (GPU) and compare with `tools/compare_fields.py`.
- Build both the CPU and GPU paths; the CPU build is the reference for debugging.
- Never declare work done with failing builds or unverified results.

See [Validation](validation.md) for the full verification workflow.
