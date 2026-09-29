# Next session: numerics review step 7 — `moby_prepare` does ALL preprocessing

STATUS (2026-09-29): **7-0 … 7-4 DONE and gated in one session** (review
section 10, step 7 MEASURED: every number, and the control that explains the
only non-zero GPU comparison). 7-5 is the plan entry below and is NOT
started. Deviations from the text below, all recorded in the MEASURED
block: the solver keeps `build_level_lines` (a pure function of the level-0
line) instead of reading per-level lines from the file — the file carries
level 0 only, as P3 wrote it; the "restart cross-check that STAYS" did not
exist (only dataset extents were checked) and was ADDED; les_ibm's
committed `ibm_coeff.h5` was a SECOND legacy-layout user and was converted
once (`ibm_coeff_case.h5`); the analytic-IBM inis needed no driver edits
because the solver prepares in-process (their explicit prepare+solve legs
live in `run_gates_step7.sh`); a derived nb carries NO parity guard (the
red-black colour offset is per block origin); `[mpi] dims` in a prepare
input describes the solve, not the prepare run. Handout written 2026-09-29; the `nb` rule
below was settled the same day — `nb` stays optional, the rank-box LAYOUT
goes. Steps 0–6 of the
review's execution sequence (`docs/numerics_review_2026-09-26.md` section
10) are DONE and gated on `main`, and step 2d (the pipe re-run) is closed —
see `docs/next_session_numerics_steps.md` STATUS. This is the next item.
Read `CLAUDE.md`, the review's section 1 (the design), section 10 step 7
(the sequence) and `docs/prepare_solve_strategy.md` (what P0–P3 already
built) before touching anything. Then execute the increments below, in
order, each closed by its gate before the next starts. Do not start steps
8–12. Record every measured number in the review under step 7 (a MEASURED
block) and update this file's STATUS when you stop.

## Why (one paragraph)

The solver still has three init paths (`moby_solve.f90:94-120`: analytic
`refine_body`, analytic `remove_solid`, plain), and the file path rebuilds
the leaf table from the masks and then cross-checks the file row by row
(`read_ibm_coeff_file`, the `ierr == 2` branches; `check_block_table` in
`field_hdf5.c`). `moby_prepare` refuses a body-free case and an unset
`[blocks] nb` (`moby_prepare.f90:90-95`). After this step the case file is
the single source of truth for grid, leaf table and face kinds; the solver
reads it and never rebuilds it; there is ONE builder, invoked by
`moby_prepare` or by `moby_solve` itself when the file is missing.
Expected: about −300 lines across `moby_solve.f90` / `ibm.f90` /
`blocks.f90` / `comm.f90`, one fewer concept, and the natural home for a
per-leaf weight column (load balancing, later).

## What is already true (do not rebuild it)

- The case-file CONTRACT exists (P0–P3): header attrs, `blocks` leaf table
  (origin + level per global id, Morton order), node lines per level
  (`fdm_h5_case_append_grid`), `coef_blocks`, optional `coef_p_blocks`,
  masks, `block_active`, `dwall_blocks`. `write_case_file` (io.f90:742) and
  `fdm_h5_case_*` are each the exact inverse of a reader.
- A reader of the `blocks` table exists in `field_hdf5.c` (`H5Dopen2(file,
  "blocks")` at ~1206, used by the cross-check and the restart layout
  check at ~844). The restart reader already reads block-table snapshots.
- `[blocks] nb` is PER-DIRECTION (`nb = 64 44 48`), each ≥ 4, even, dividing
  the grid (config.f90:551-560). `validation/block_nb` proves the tiling is
  bookkeeping: `nb = 8` == `nb` unset (rank box) == `nb = 32 16 8`, all
  `max_abs 0`, 1 == 4 ranks. That property is what makes the rank-box
  removal bit-exact by construction.
- Prepare with the CPU build is canonical: the GPU build computes the
  analytic coefficients on the device (libm ulps). P0 measured "GPU solve ==
  CPU solve from the SAME file at tol 0".

## Increments

**7-0 — reference set.** Cut `~/step7_ref_binaries/` (CPU + GPU, production
AND nofma, `moby_solve` + `moby_prepare`, PROVENANCE with the hash) from the
CLEAN commit that carries this handout. Everything below compares against
it. Half an hour. (`~/numrev_ref_binaries` is `6db60ef`, pre-step-6; its
`apply_bc` differs, so it is not a reference for this step.)

**7-1 — prepare accepts body-free cases.** Drop the `[ibm] enabled` refusal;
a body-free case writes attrs + node lines + `blocks` + face kinds, and NO
coefficient/mask/dwall datasets. The solver's readers tolerate their
absence (they already return `found = 0` for optional datasets — check
each: `coef_blocks` is currently mandatory once `coeff_file` is set).
Drop the `nb` refusal too, replacing it by THE DEFAULT RULE (decided
2026-09-29, see "The `nb` rule" below): an unset `nb` is chosen by the
builder from its own rank count and stored in the file.
Gate: `min_channel` and `beltrami/slab_y` prepared then solved == inline
solve, `max_abs 0` at production flags, CPU 1 + 4 ranks and GPU.

**7-2 — the solver READS the leaf table.** On the file path,
`init_block_set` takes the `blocks` table (and the node lines) from the
case file and builds only the derived lookups (`lidOf` per level,
`levelOffset`, `nBlocksGlobal`, the Z-order split `idStart/nBlocks`);
`build_leaf_table` / `build_level_lines` / `classify_*` are no longer
called on the solve path. The row-by-row cross-check of coefficient and
dwall tiles against an ini-derived rebuild (`ierr == 2`) then compares the
file with itself and goes; the RESTART cross-check (a snapshot's `blocks`
table against the case file's) STAYS — it guards a stale restart.
Gate: the 7-case and 9-case suites from prepared files, `max_abs 0` at
production flags vs 7-0 (CPU 1 + 4 ranks, GPU); `validation/prepare`
P0 22/22, P1 16/16, `run_gates_big.sh` sailplane leg; `validation/block_nb`.

**7-3 — one builder, two entry points.** Move `moby_prepare.f90`'s body into
a module routine (`prepare_case(input_file, case_file, c)`); `moby_prepare`
calls it; `moby_solve` calls it IN-PROCESS when the case file named by the
ini is absent (`moby_solve --prepare` forces it), so every tutorial stays
one command. The case-file key becomes the ONE user-facing name — proposal
`[case] file = x.h5`, with `[ibm] coeff_file` accepted as an alias for one
release and a deprecation line printed. STALENESS is attribute-by-attribute,
not a hash: the readers already take `leng`, `re`, `nb`; add the refine
keys, `periodic_*`, `refine_dims`, `remove_solid`, `keep_buried`, the STL
list + transform, and `[scalar] count > 0` (the `coef_p_blocks` column) to
the stored attrs and hard-error on a mismatch naming the key. `re` IS a
case-file input (the coefficients carry `1/Re`).
Gate: a tutorial run with the case file deleted regenerates it and
reproduces the 7-2 result `max_abs 0`; a deliberately stale file (change
`re` in the ini) error-stops with the key named; `config` gate: alias
accepted, both keys set = error.

**7-4 — delete.** Remove: the inline classify dispatch
in `moby_solve.f90` and the file branches of `classify_active_mask` /
`classify_refinement_masks` (ibm.f90 ~895-975); the legacy global-layout
coefficient reader (`fdm_h5_read_ibm_coeff` + the fallback in
`read_ibm_coeff_file`, ibm.f90:420-475); `DIST_RANKBOX` and every branch
keyed on `distMode` (`init_block_set` else-branch blocks.f90:214-220,
`face_kind` ~999, `build_block_metrics`, comm.f90 `rank_box`/`local_range`
and the rank-box peer enumeration in `init_block_exchange`,
`set_serial_local_size`); the device `set_ibm_coeff` kernel (the host twin
becomes THE kernel; keep the `KEEP IN LOCKSTEP` comment's partner honest by
deleting it). The `[blocks] nb`-less live inis (list below) stay AS THEY
ARE: under the `nb` rule they run one block per rank exactly as today.
`validation/block_nb`'s "`nb` unset" gate stays and now states the rule
(unset == the rank-derived lattice == `nb = 8`, `max_abs 0`).
Gate: 7-case + 9-case suites `max_abs 0` at production flags vs 7-0 on the
UNCHANGED inis (the nb-independence property is what makes this exact —
any non-zero here is a real defect); `validation/freestream`, `redblack_interface`,
`rans_inlet` 1 + 4 ranks; the analytic-IBM cases (list below) as
prepare(CPU)+solve pairs on CPU AND GPU vs the old CPU inline binary
(`max_abs 0` — the P0 property); `min_channel` 1 == 2 == 3 == 4 ranks EXACT.

**7-5 (not this session; write the plan entry only).** A per-leaf weight
column written by prepare (body block, cut cells × the conjugate rate
share, trip blocks, refined blocks) and a weighted prefix sum replacing
the closed-form `zorder_start/count`; results stay rank-count independent
because the split moves work, never arithmetic.

## The `nb` rule (decided 2026-09-29)

The block size is a property of the CASE FILE, never of the run: that is
the one thing the rank-box layout cannot satisfy (its block is
`globalSize/dims`, which changes with the rank count), and it is why the
layout goes while the KEY stays optional.

1. Explicit `[blocks] nb` wins, with today's guards (per direction even,
   ≥ 4, divides the grid; ≥ 4 is a typo guard from Phase 1a with no
   derivation, "even" is the red-black colour continuity across block
   faces, `colorOffset = modulo(sum(origin), 2)`).
2. Unset: the builder — `moby_prepare`, or `moby_solve` preparing
   in-process (7-3), which then knows the SOLVE's rank count — forms the
   Cartesian decomposition the solver would (`MPI_Dims_create`) and sets
   `nb = globalSize/dims` per direction: ONE BLOCK PER RANK, i.e. exactly
   today's rank box (same halo footprint, same ownership, bit-identical
   fields), but written into the file. On 1 rank that is the whole grid as
   one block. No floor applies; the parity check applies only where it
   matters — a direction split into > 1 block, or periodic — and the
   session must CHECK what the current rank-box path does with an ODD
   periodic direction (a single wrapped block of odd length breaks the
   checkerboard today too) and keep that behaviour.
3. If the decomposition does not divide the grid evenly (today's rank
   boxes may be unequal, a lattice cannot express that): fall back to the
   grid-only rule — per direction the largest even divisor ≥ 4 up to a
   ceiling (64) — and PRINT it.
4. Refinement and buried-block removal keep requiring an explicit `nb`
   (they do today).
5. Whatever was chosen is printed at init and stored in the case file; a
   file prepared on R ranks solves on any rank count that still owns ≥ 1
   block per rank, else the existing "rank owns no blocks" stop.

Gate consequence: the P0 leg "prepare 1 == 4 ranks identical files" holds
for explicit `nb` only; add a leg for the rank-derived default — DIFFERENT
`blocks` tables, SAME fields at tolerance 0.

## The lists the increments need

Live inis with NO `[blocks] nb` (archived `overheadTest/horeka/results_*`
copies excluded). Under the `nb` rule they need NO edit; the list is here
so the session can run a sample of them as the "unset == today" gate:
`tutorials/channel_kmm180/input.ini`, `tutorials/sailplane/input.ini` (also
the last user of the LEGACY coefficient file `sailplane_ibm_coeff.h5`: move
the tutorial to its prepare flow, which `run_gates_big.sh` already runs at
`nb = 10`, and retire the committed legacy file to a reference),
`tutorials/turbulentBoundaryLayer/{cold_start,production,production_stats}.ini`,
`.../overheadTest/horeka/configs/base_jacobi.ini`,
`.../overheadTest/singleLevel/{base_jacobi,base_redblack}.ini`,
`validation/beltrami/uniform.ini`, `validation/block_nb/base.ini`,
`validation/channel_interface/{reference,uniform128}.ini`,
`validation/conjugate/{kasagi_channel,kasagi_uniform}.ini`,
`validation/rans_sst/wf180_y30.ini` (grid 8 × 6 × 8 — a good test of rule
2 on 4 ranks), `validation/scalar/{conduction,prsweep,smoke,wferr,
wfs180_y05,wfs180_y15,wfs180_y30,wfs180_y45}.ini` (`conduction` has nx = 4).

Live inis on the INLINE analytic IBM path (`[ibm] enabled`, no
`coeff_file`, no `stl_file`) — they become prepare+solve pairs and their
drivers gain a prepare step:
`validation/conjugate/wavy.ini`, `validation/prepare/{wavy,wavy_refine,
wavysolid}.ini` (already the P0 gates), `validation/rans_geometry/{wavy,
wavy_refine}.ini`, `validation/scalar/{ibmwavy,ibmwavyr}.ini`,
`tutorials/turbulentBoundaryLayer/overheadTest/horeka/configs/rough_jacobi.ini`
(`wall_shape = eggcarton`; `set_ibm_geometry` is applied by BOTH binaries
today — keep that, it is load-bearing).

## Landmines, inherited

- Build environment on this workstation: `source /etc/profile.d/lmod.sh &&
  module load toolkits/nvhpc/25.9` in the SAME shell as `./compile.sh`;
  without the source the build silently picks the system gfortran MPI
  module file. `build_gpu_corax` on istmcorax is free again (the pipe run
  ended 2026-09-28); the 15 snapshots in `~/pipe_rerun_f2/` (85 GB) can go.
- The generated suite inputs (`IC_turbles.h5`, `IC_turbslab.h5`,
  `IC_refine.h5`, `ibm_coeff_blocks.h5`) are symlinks into
  `mobydiff.scalar` / `mobydiff.bl`, not committed.
- `run_bitexact*.sh` key their output prefix on `MODE`: two concurrent
  suites with the same label in one tree delete each other's snapshots.
- `tools/h5maxdiff` with no dataset arguments never compares a passive
  scalar; the suite drivers name their datasets explicitly.
- `config.f90` has no `case default`: an unknown key in a known section is
  silently ignored, so a gate that changes ini keys must be run with a
  binary that KNOWS the key on both sides, or the old side runs something
  else and reports success.
- Production-flag comparison is valid only when the arithmetic source is
  byte-identical (this step copies/reads values, so it is); anything that
  moves an expression needs `cpu_nofma`/`gpu_nofma` on both sides.

## Stop conditions

A `max_abs` that is not 0 and not explained by the ONE expected move (the
analytic coefficients now always come from the CPU host kernel, so a GPU
analytic-IBM case compared against the OLD GPU inline binary differs by
libm ulps — compare against the old CPU binary instead); a live case whose
grid admits no `nb` under rules 2 AND 3 (report it, do not force one); a
restart file
that the new reader rejects for a reason other than a genuinely stale
table. Commit after each increment with "step 7-N" in the subject.
