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
(`src/modules/scalar*.f90`, `[scalar]`/`[scalar.N]`). Fortran + MPI (3D
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
  `docs/history/next_session_verification.md` §2.
- **A variable that appears ONLY as the bound of a mapped array section
  (`map(to: a(1:n))`) must also be read by ordinary code** -- as a loop
  bound, typically. nvfortran 25.9 drops the assignment of a variable it
  sees used nowhere else, the section is then sized from an uninitialised
  value, and the run stops with "partially present on the device", or does
  not, at random (found 2026-10-01 in the first form of
  `predict_outlet_faces`: two runs in three). The boundary kernels loop
  over all `npts` points and skip for this reason.

## State of the solver (2026-10-24)

The full work log (block refactor, 2:1 refinement, RANS/IDDES, the airfoil,
passive scalars and conjugate heat transfer, the prepare/solve split, the
performance campaign, the outlet rework, the cleaning passes) is
**`docs/history/worklog_2026.md`** -- one bullet per finished item with its
deviations, findings and gate numbers; **append new work there**. The closed
handouts are indexed in `docs/history/README.md`; the master design is
`docs/block_refinement_strategy.md`; the two reviews are
`docs/numerics_review_2026-09-26.md` (findings F1-F9, numerics) and
`docs/cleanup_review_2026-10-06.md` (F10-F16, duplication, simplifications,
the profile). **The current handout (2026-10-25) is `docs/next_session_sailplane_junction.md`**
(the wing-fuselage junction study; its section 2 is the ledger, step 1 = the
XFOIL + AVL lifting line). The two before it are closed and
stay in `docs/` until the next cleaning pass moves them to `docs/history/`:
`docs/next_session_code_simplification.md` (items 1-12, closed 2026-10-21)
and `docs/next_session_revalidation.md` (the one-by-one numerics
re-validation of `validation/README.md`, the index of every surviving case;
closed 2026-10-24 -- its section 7 ledger reads 36 DONE, 1 RETIRED, the
boundary layer STANDS, no PENDING row). What is still open is listed at the
end of this section.
`main` is the working branch (consolidated 2026-09-25); the only unmerged
history is the tag `archive/claude-jacobi-interface` (`boostconv`,
`kpin_dwall`, both deliberately dropped).

**What the solver is, in one pass through the time step.** The state lives
in a `block_set_type` of equal-size blocks (`[blocks] nb`, optional: unset =
one block per rank of the prepare run, the nb rule) on a Z-order curve whose
bit order `moby_prepare` chooses from the geometry (`min_surface_key_order`,
recorded as `block_key_order`; `report_partition` prints the cross-rank /
cross-node face cells at init -- read it on every multi-node run). 2:1
refinement by boxes or by the body (`refine`, `refine_body` +
`refine_body_levels`/`_box`, `refine_dims = xyz | xz`), buried blocks
removed (`keep_buried` is LOAD-BEARING for penalization forces and
conjugate scalars). The exchange (`comm.f90`) is one weighted gather over
entries enumerated once (`collect_pieces`), same-level copies first, a
divergence subset as index lists, every transfer extended tangentially into
physical halos; the communicator is a dup of WORLD and ranks are world
ranks. The momentum predictor is one fused RK3 kernel (skew convection
hardwired; `[flow] convection = divergence` re-admitted for the CaNS
boundary-layer tutorial only; the `if (skew)` guard is a REGISTER lever, do
not fold it) + separate correction passes (eddy viscosity, body force /
trip over a block list, the penalization state correction with the RATIONAL
factor `1/P3` -- third order, no libm); the declared outlet face is PREDICTED
(`predict_outlet_faces`: the neighbour's non-pressure right-hand side and
its own pressure gradient against the held datum, the TOTAL pressure where
fluid enters, its own penalization on a body; low and high faces one
statement) and is restart STATE (`un_xmax` planes). The projection is
damped Jacobi / Chebyshev (`[pressure] accel = chebyshev`, `niter = 12`,
`sor` ignored) or red-black SOR (`solver = redblack`, `sor = 1.5`), both on
the static metric tables (`dnLow/dnHigh/d1P`, `cfLow/cfHigh`, index 0 for
the sweep's redundant halo layer; `d1P` stays a separate factor), `rdenom`
recomputed on body blocks only, `interface_correct` one launch over a static
face list, the outlet's phi ghost held at zero as data. Turbulence:
`[turbulence] model = none | les | rans | iddes` (WALE/Smagorinsky; k-omega
SST with wall functions, gamma-Re_theta transition, inlet freestream
values; SST-IDDES with the WALE blend in the RANS-retention convention);
`nut` is read BACK on a RANS/IDDES restart (2026-10-21). Passive scalars
are extra variables of `blk%q` (`[scalar]`/`[scalar.N]`, `ibm_wall =
dirichlet | adiabatic | conjugate`, central convection, the conjugate face
coefficient on the level-set fraction, `scalar_*` modules). Time step:
`cflmax` and `pecletmax` bound the SUM over directions (RK3 limits sqrt(3)
and 0.628), both reductions over the INTERIOR (2026-10-21), `dtmax` the cap;
do not push both limits at once. The runtime text traces (the scalar heat
file, the airfoil forces, the boundary-layer line) RESUME across a restart
since 2026-10-23 (`io.f90` `resume_runtime_log`: rows up to the resumed
run's first sample kept, later ones dropped, then append) and an airfoil
restart primes the control-volume momentum at setup, so a chunked chain
writes one continuous trace with no dmdt = 0 row per chunk.

**The prepare/solve contract.** Every run reads a CASE FILE written by
`mpirun -n N moby_prepare input.ini` (`[case] file`, default
`<field_prefix>.case.h5`; `moby_prepare --case-file-name` prints the rule,
`--check-config` validates an ini): grid, block size, leaf table, face
kinds, IBM coefficients (`[ibm] stl_file` / analytic `wall_shape`), wall
distance. The solver never builds one; a stale file is caught by the input
echo (`case_inputs`), the first differing key named. Prepare with the CPU
build on the rank count you solve with (an unset nb is one block per
prepare rank); drivers call `tools/prepare_if_missing.sh`. `config.f90` has
a `case default`: an unknown section, key or value STOPS the run. A restart
takes state (`step`, `t_current`, `cfl`, `dt`) from the field file, run
length from the ini if set, and for everything else the ini is authority:
the grid keys must agree with the file, a run-control key left at its
default must agree too, else the run stops naming both values
(`docs/configuration.md` "[restart]", 2026-10-21).

**Gate rules (what every change is measured against).** Reference binaries
are cut from a COMMIT, with the hash in a PROVENANCE file:
**`~/pass2_ref_binaries` (`d821341`; CPU + GPU, production and nofma, solve
+ prepare) is the current set** -- every change since it is `max_abs 0` on
the suites at nofma (the 2026-10-21 interior-only Courant reduction moves
only a Courant-bound cold start from rest with an inlet, which no suite
holds; last confirmed 2026-10-23 after `52ee737`: 9 / 11 / 7 cases, CPU and
GPU). Use the NOFMA builds (`./compile.sh cpu_nofma | gpu_nofma`) for any change
that moves an expression and the production flags only when the arithmetic
source is byte-identical; the GPU compiler contracts a moved statement
differently (the "FMA class", 1e-15..1e-12) and a production-flag difference
of that size on the GPU with an exact CPU is not a defect. `./compile.sh`
takes ONE mode per call. The suites: `validation/scalar/run_bitexact.sh`
(9 cases: channels, 2:1, LES, RANS, red-black, Beltrami), `run_bitexact_s3.sh`
(11 scalar cases incl. the mixed-mode ones), `validation/freestream/
run_bitexact_outlet.sh` (7 outlet cases; the other suites have NO outlet)
and `run_gates.sh` groups `low refined restart`, `validation/body_outlet/
run_gates.sh` (body on an inlet/outlet plane, `pocket`), `validation/
conjugate/run_gates_c1.sh`, `validation/prepare/run_gates*.sh` (P0/P1/order),
`validation/courant_limit` and `diffusion_limit`, `validation/redblack_interface`,
`validation/block_nb`; drivers take absolute `REF`/`NEW` paths and a `MODE`
label, and `tools/h5maxdiff` (a build product; rebuild it with the tree)
compares `un vn wn pn` + the RANS scalars by default -- name a passive
scalar or an outlet plane explicitly. Rank-count identity (1 == 2 == 4) is
a gate in its own right for anything touching the exchange or the halo
layer. `OMPI_MCA_hwloc_base_binding_policy=none`; always `mpirun`, even for
one rank (MPI singleton init fails here); `module load` in the same shell
as the run. Never edit a driver while it runs (bash reads it incrementally);
a gate group that reuses another group's output compares whatever binary
THAT group ran; a gate that passes on one build is not a gate (the CPU
builds survived a stack overflow the GPU host code did not). GPU
statistics/heat files are not bit-reproducible run to run (atomic order);
the fields are. Short benchmarks understate a per-step saving in a bucket
that also holds first-call setup; name the device of every timing
(RTX 3060 here, A6000 on istmcetus, RTX 5090 on istmcorax, A100 on HoreKa);
ncu needs permissions the workstation lacks -- registers from `cuobjdump
-res-usage`, times from nsys; `-gpu=maxregcount` is a TARGET, not a cap.

**Landmines still live.** `nb` must divide the global grid in every
direction. The committed `validation/blasius/blasius.case.h5` is refused by
nofma builds (its stretched line is one bit apart; the node-line check is
exact): strip `[case] file` and prepare per run, as the suites do. A
solver-minted template carries the outlet plane of its own run: an IC tool
must write the plane or delete it. `tools/compare_fields.py` reassembles
block-table snapshots on the FINEST lattice (69 GB at L5; use
`validation/prepare/compare_snapshots.py`); comparing block-layout files by
row POSITION is wrong across key orders (the tools match rows). A floating-
point cancellation is never a classification (the `ibm_value = 0` double
count). `x/6.0d0` is a second fp64 divide in a per-DOF kernel: write
reciprocals of non-dyadic constants as multiplies. `set_ibm_coeff_host` and
the device kernel are KEPT IN LOCKSTEP (comment). Penalization forces need
`keep_buried`; the A2 penalization integral does not transpose to a
Dirichlet scalar. The 2:1 const-1/2 transfer and the coarse-owns-face v'
residual are one trade (accepted); the SD7003 transition front is the
measured motivation for the deferred TVD/second-halo increment (S5b). The
`t_final` final-step trim, the cold-start RANS IC, the absent-body restart,
the level jump on a physical face, the body on an outlet plane are all FIXED
and gated -- the worklog has each mechanism.

**Open, in order.** The re-validation pass is DONE (2026-10-21 -> 10-24):
every long number on file is now measured on a NAMED binary at the pass's
conventions (Chebyshev `niter = 12` where the ini is Chebyshev, `cflmax =
1.0`; the NACA keeps its `niter = 18`, the pipe ran red-black 6, the
boundary layer stands at its red-black 6) -- the `RE-VALIDATED` blocks in
the tutorial and `validation/*/` READMEs are the records; the pass found and
fixed the Q4 Neumann ghost row (`a991a1c`) and the three restart defects of
the text traces (`52ee737`). (1) The performance ideas measured but not
taken: the momentum predictor split by component, the red-black reduced
exchange round, 7-5 (a weighted Morton split), the `apply_bc` face-row skip,
the `rans_transport_kernel` register cut (the twin of the scalar one,
`docs/next_session_code_simplification.md` item 1), the fp64 divide in
`compute_rdenom` on body blocks, the unexplained `phi_exchange` locality
shift (`docs/cleanup_review_2026-10-06.md` section 8). (2) The feature tail:
the IDDES flat-plate inlet, transition/wall functions under IDDES. (3) The
sailplane wing-fuselage junction (decided 2026-10-24): a DNS / wall-resolved
LES of the junction region -- half the fuselage, a few chords of wing -- with
IDDES/RANS elsewhere, the laminar-turbulent transition at the junction the
question; strategy proposed 2026-10-24 (worklog), the PLAN is
`docs/next_session_sailplane_junction.md` (2026-10-25: geometry measured --
AR 30, root chord 0.831 m, flapped section, no horizontal tail, no fillet;
steps 1-7 with their decisions D1-D5; `tutorials/sailplane/junction/
stl_sections.py` slices, measures and pitches the STL), no case yet. Dropped
2026-10-24 by decision: S5b (TVD + second halo), C4 (conducting corners), the
MacDonald rough-wall forced-convection case.

## Verification

- Pure refactors must be bit-exact vs. the pre-refactor code: compare
  output fields with `tools/compare_fields.py`. Build BOTH sides with
  `-Mnofma` (CPU) / `-Mnofma -gpu=nofma` (GPU) for the comparison —
  default FMA contraction makes the compiler introduce 1-2 ulp
  differences for arithmetically identical source.
- Channel sanity: `tools/check_parabolic_channel.py`. IBM cases:
  `tutorials/sailplane/` (coefficient file path),
  `validation/prepare/wavy*.ini` (analytic `set_ibm_coeff`). The index of
  every surviving test case is `validation/README.md`.
- Refinement phases: uniform-flow preservation across interfaces and global
  mass conservation to round-off (see strategy doc §11 for the full list).
- Anything touching an inlet/outlet face or the exchange next to a physical
  boundary: `validation/freestream/run_gates.sh` (groups `low`, `refined`,
  `restart` since 2026-10-01) and `run_bitexact_outlet.sh` -- the 7- and
  9-case suites contain no outlet. A uniform-flow gate started from the
  uniform field cannot see a halo that nothing writes: the `refined` group
  starts from rest. An immersed body on an inlet or outlet plane:
  `validation/body_outlet/run_gates.sh` (since 2026-10-05).
- Anything touching the time-step limiter, the convection or the diffusion
  operator: `validation/courant_limit/run_gate.sh` and
  `validation/diffusion_limit/run_gate.sh` (the RK3 limits, bracketed;
  `cflmax` and `pecletmax` bound the sum over directions). A run that prints
  `WARNING: dt x sum_d(nu/h_d^2) ... exceeds the RK3 limit` is beyond the
  limit of its grid: lower the fixed dt or `pecletmax` unless the field is
  exactly 1D/2D.
- Never declare a phase done with failing builds or unverified results.
