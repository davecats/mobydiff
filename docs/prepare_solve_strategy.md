# Strategy: splitting mobydiff into `moby_prepare` + `moby_solve`

Status: PROPOSAL (2026-07-16). No code has been changed. Modeled on the
AMPHIBIOUS (CPL) preproc/solver split, adapted to mobydiff's block/leaf
architecture.

## 1. The AMPHIBIOUS model, in one paragraph

`preproc.cpl` is itself an MPI program (one rank per candidate cube). Each
rank independently builds the global grid and FD metric tables, loads the
STL files directly (raw binary read, ~100 lines), builds a BVH over the
triangles, classifies every staggered point by majority-vote ray-parity
casting, computes graded sharp-interface IBM coefficients from ray/triangle
intersection distances along the coordinate axes, and finally all ranks vote
on which cubes contain zero fluid points; those are pruned. Output: one
self-contained `cube-data.N.bin` per surviving cube plus `general-info.dat`.
`solver.cpl` only reads cube files, maps cubes to however many ranks it was
given (static balancing + `MPI_Dist_graph_create_adjacent`), and time-steps.
Everything is one language, the geometry work is embarrassingly parallel,
and nothing geometric is recomputed at solve time.

## 2. Where mobydiff stands today

We already HAVE a two-step architecture — it is just spread over three
tools, two languages, and one redundant recomputation:

| Step | Tool | Language | Parallel? |
|---|---|---|---|
| Grid export | `src/mobygrid.f90` | Fortran, serial (hard-errors on >1 rank) | no |
| Geometry: masks, leaf table, `coef_blocks`, `dwall_blocks` | `tools/mobygeom.py` (~3200 lines) | Python (trimesh + libigl) | `multiprocessing`, single node |
| Leaf table + exchange + (analytic) coefficients | solver init (`main.f90`) | Fortran + MPI (+GPU) | yes |

What is actually wrong with it:

1. **The leaf builder exists twice.** `build_leaf_table` (blocks.f90:537)
   and its rule-for-rule Python mirror `build_leaf_table_py`
   (mobygeom.py:1732), plus mirrored Morton curves, level windows and
   subdivided lines. The solver cross-checks the file row-by-row
   (`check_block_table`, field_hdf5.c:1244) precisely because the mirror
   can drift. This is the single largest maintenance hazard in the
   preprocessing chain.
2. **The slow parts are Python.** Winding-number classification over the
   padded-bbox window of the finest lattice, per-leaf STL tile evaluation,
   and the 2 GB-chunked integral-image reductions are numpy + libigl under
   `multiprocessing` — single node, no GPU, and the deep-refinement cases
   (B11 airfoil L5/L6) already needed windowing heroics (`mobygeom.py:1439`,
   `:1518`) to fit at all.
3. **Two handshakes.** mobygrid.h5 exists only so mobygeom agrees with the
   solver about node lines; the solver then re-derives everything from the
   ini anyway and re-verifies the file.
4. **A venv dependency** (trimesh, libigl, scipy, rtree) for any
   STL-geometry run.

What is already RIGHT about it (and better than AMPHIBIOUS — keep):

- The prepared file is **rank-count independent**: the solver splits the
  global leaf table by Z-order (`zorder_start/count`) at load. AMPHIBIOUS
  bakes cubes=ranks at preproc time and needs a balancing pass in the
  solver. We should NOT copy that; our Morton split is the cleaner design.
- The dataset contract already exists and the solver already reads all of
  it: `blocks`, `block_touch_l{l}`/`block_buried_l{l}` (+windows),
  `coef_blocks`, `dwall_blocks`, `block_active`, `refine_dims`
  (io.f90 / field_hdf5.c readers, ibm.f90:313, rans.f90:370).
- The io layer already does independent per-block/row-range parallel-HDF5
  transfers in both directions — exactly what a parallel prepare needs.

## 3. Feasibility verdict

**Yes — and cheaply, because the decisive abstraction already exists.**
The analytic pipeline computes everything downstream of geometry from a
single point-membership indicator:

- masks: `classify_active_blocks` / `classify_block_geometry`
  (ibm.f90:440/520) scan `isInBody`;
- coefficients: `set_ibm_coeff` (ibm.f90:719) builds the graded
  sharp-interface formula `Σ((d0−d)/d)/d0²` via `add_neighbor_coeff` +
  `bisection` (ibm.f90:391–432) — from `isInBody` alone;
- wall distance: `walldist.f90` (phase T1b) computes `dwall` to
  `dwall_tol ≈ 1e-10` from an indicator passed as a **procedure argument**
  (surface cloud + kd-tree + polish), periodic dims folded by minimum
  image.

mobygeom's outputs are by construction the same quantities (its coefficient
formula is documented as identical, mobygeom.py:987). So the entire STL
port reduces to: **one new Fortran module providing an STL-backed inside
test**, plugged into machinery that is already geometry-agnostic or one
procedure-argument away from it. AMPHIBIOUS proves the required geometry
kernels are small: binary STL reader (~100 lines), median-split BVH
(~120 lines), Möller–Trumbore ray/triangle + parity vote (~150 lines).

## 4. Target architecture

Two executables built from the same module library (compile.sh gains one
target next to `main` and the test drivers):

```
moby_prepare  (MPI, CPU)                     moby_solve  (MPI, CPU/GPU)
  read ini (config.f90)                        read ini (config.f90)
  node lines (init.f90)                        node lines (init.f90)
  geometry indicator:                          read case file:
    analytic (ibm.f90 isInBody)                  blocks table  -> cross-check
    or STL   (new geometry_stl.f90)              masks         -> leaf table
  masks  (ibm.f90 classify_*)                    coef_blocks   -> ibm%coef
  leaf table (blocks.f90 build_leaf_table)       dwall_blocks  -> sst%dwall
  coef tiles (ibm.f90 set_ibm_coeff)           exchange entries (comm.f90)
  dwall tiles (walldist.f90)                   device maps, time stepping
  write ONE case file (io.f90/field_hdf5.c)
```

Key decisions:

- **The case file IS the current block-table coefficient file**, extended
  minimally: add the mobygrid grid datasets (`x_nodes/y_nodes/z_nodes`,
  grid attrs — writer `write_grid_export` already exists, io.f90:404) and a
  config-echo attribute set into the same HDF5. No new format; every
  existing reader keeps working; mobygeom-produced files stay loadable
  during the transition. `mobygrid.h5` as a separate artifact disappears
  (`mobygeom.py` can read the same attrs from the case file for
  cross-validation runs).
- **The solver keeps rebuilding the leaf table from the masks and
  cross-checking the file row-by-row** (today's file-IBM refine_body path,
  unchanged). This preserves the strongest consistency guarantee in the
  code for free. Trusting the file's `blocks` table outright (skipping the
  rebuild) is a *later, optional* optimization for extreme lattices — keep
  the cross-check behind a flag if ever done.
- **The solver keeps its inline analytic path** (`set_ibm_coeff` on device,
  no prepare step needed). Channels, wavy walls, Beltrami and all
  bit-exactness gates stay zero-preprocessing. `moby_prepare` is required
  exactly where mobygeom is required today: file/STL geometry — and
  becomes *available* for heavy analytic refine_body cases (the redundant
  per-rank classification at init, hard-capped at 2e8 finest cells,
  ibm.f90:703, moves offline).
- **prepare is CPU-only MPI.** The kernels are per-leaf independent; MPI
  ranks (+ optional OpenMP threads inside a rank) are plenty and keep the
  tool runnable on login/pre-post nodes. No OpenMP-offload in prepare —
  simplicity wins; revisit only if a measured case demands it.

## 5. Reuse map (the no-duplication argument)

`moby_prepare.f90` is a ~200-line driver. Every stage is an existing
routine; the only genuinely new code is the STL module and two writers.

| prepare stage | implementation | status |
|---|---|---|
| config parse | `read_runtime_config` (config.f90) | exists |
| node lines | `init_grid`/`build_node_line` (init.f90), `build_level_lines` (blocks.f90) | exists |
| geometry indicator (analytic) | `isInBody` (ibm.f90) | exists |
| geometry indicator (STL) | `geometry_stl.f90`: reader + BVH + ray-parity inside test | **NEW (~500–700 lines)** |
| per-level touch/buried masks | `classify_block_geometry` (ibm.f90:520) with the indicator as a procedure argument (same generalization walldist already made in T1b) | exists, small refactor |
| solid-block pruning / `block_active` | `classify_active_blocks` (ibm.f90:440); `--keep-buried` equivalent = config flag writing zeroed buried masks (mobygeom.py:1959 semantics) | exists |
| leaf table | `build_leaf_table` (blocks.f90:537) — the ONE implementation; `build_leaf_table_py` retires | exists |
| coefficient tiles | `set_ibm_coeff` (ibm.f90:719) per local leaf, ghost-inclusive window; indicator-agnostic via `bisection` | exists, host loop |
| dwall tiles | `build_walldist` + queries (walldist.f90), indicator-driven | exists |
| rank split of the work | `zorder_start`/`zorder_count` (blocks.f90) — same split the solver uses | exists |
| file writes | `write_grid_export` (io.f90:404), `write_block_table` (field_hdf5.c:260), mask writers; row-range writers for `coef_blocks`/`dwall_blocks` mirroring the existing readers (field_hdf5.c:1284/1393) | mostly exists; **NEW: 2 C writer functions** (same hyperslab pattern as `fdm_h5_write_field`) |

The one deliberate interface change: `isInBody` today is a fixed
`declare target` function (analytic wavy wall). Host-side consumers
(classify_*, bisection path in prepare, walldist) take the indicator as a
procedure argument instead — walldist.f90 already works exactly this way,
so this extends a proven idiom rather than inventing dispatch. The GPU
in-solver analytic kernel keeps calling the concrete `isInBody` directly
(declare-target procedure pointers are not worth the pain; the device path
never needs STL).

## 6. The new STL module (`src/modules/geometry_stl.f90`)

Scope (deliberately minimal, AMPHIBIOUS-proven):

- **Reader**: binary STL (80-byte header, count, 50-byte records), float32
  vertices widened to float64 exactly — the float32 quantization is part of
  the as-built geometry (T1 gate finding; mobygeom convention). Multiple
  files concatenate. ASCII STL: detect and error with a pointer to a
  converter (mobygeom keeps that job).
- **BVH**: median-split over triangle centroids, leaf size ~4
  (find-inside.cpl is a direct template; walldist.f90's kd-tree shows the
  house style for flat-array trees).
- **Inside test**: ray-parity casting with the fixed direction table +
  degenerate-hit perturbation retry + majority vote (3 rays, escalate to
  11 on disagreement — AMPHIBIOUS's scheme). Parity is
  orientation-independent, so no normal-fixing/repair pass is needed for
  watertight meshes (which our convention already demands —
  `stress-stl-watertightness` exists in mobygeom for this). This choice
  avoids porting fast winding numbers (hierarchical dipole approximation)
  entirely.
- **Periodicity**: fold query points by minimum image in periodic dims
  (walldist idiom); the indicator must be length-periodic exactly as the
  analytic contract states.
- Exposed as the same indicator signature the analytic path satisfies, so
  masks/coefficients/dwall flow through unchanged shared code.

Deferred, behind the same interface (P2, only if measured):
- direct BVH segment-intersection crossing query for `add_neighbor_coeff`
  (one traversal instead of ~40 indicator calls per crossing bisection);
- exact point-triangle nearest distance for dwall (replaces the
  cloud+polish path; mobygeom moved to libigl's exact query for the same
  reason — but note the indicator+polish path is already validated to
  1e-10, so this is purely performance).

## 7. Parallelization of prepare

- **Leaf-parallel everywhere it matters.** After the (cheap, redundant,
  every-rank — same as the solver) leaf-table build, split leaves by the
  existing Z-order closed form. Each rank computes coefficient + dwall
  tiles for its leaf range and writes its contiguous row range with
  independent parallel-HDF5 transfers — identical pattern to field output.
  No communication except the final metadata.
- **Masks before the leaf table** are lattice-window-parallel: split the
  per-level padded-bbox window rows across ranks, `MPI_Allreduce` the
  (windowed, small) byte rasters. This replaces mobygeom's integral-image
  chunking machinery with plain distributed loops.
- **Far-field shortcut**: keep mobygeom's optimization — leaves whose
  padded window misses the geometry bbox get zero coef tiles without any
  indicator call, and dwall via plain nearest-point queries
  (mobygeom.py:2024 semantics). This was the "hours to minutes" win; port
  the *idea*, not the code.
- Every STL query is read-only after BVH build: one BVH per rank,
  embarrassingly parallel over points. This is the part that is Python-slow
  today and mechanically fast in Fortran.

## 8. Phased plan (each phase gated before the next)

**P0 — prepare for analytic geometry (no new geometry code). DONE
2026-07-16 (gates in `validation/prepare/`).**
`src/moby_prepare.f90` (CMake target `moby_prepare`) runs the solver's own
init pipeline — config, node lines, `classify_refinement_masks` /
`classify_active_mask` (analytic branch), `init_block_set` split over the
WORLD ranks by the solver's Z-order closed form, `init_ibm` +
`set_ibm_coeff`, `fill_body_distance_analytic` (now public in rans.f90,
takes the dwall array instead of `sst`) — and writes the case file through
`write_case_file` (io.f90) / the `fdm_h5_case_*` C writers (field_hdf5.c),
each the exact inverse of its reader: attrs + `blocks` +
`coef_blocks` always, per-level `block_touch_l{l}`/`block_buried_l{l}`
(refine_body, full rasters = the legacy no-window convention),
`block_active` (remove_solid, the default), `dwall_blocks` (`[rans]`).
The solver is UNCHANGED on its solve path (the file IS a coefficient file;
`[ibm] coeff_file = case.h5` is the whole user-facing switch).
SCOPE DEVIATIONS from the plan above: the grid datasets are NOT yet in the
case file (mobygrid still serves mobygeom; absorb in P3); prepare is
required to have `[blocks] nb`, `[ibm] enabled` and no `coeff_file` input.
*Gates (all PASS, `validation/prepare/run_gates.sh`)*: wavy (single-level
+ [rans]), wavy_refine (refine_body 2 levels + [rans]), wavysolid (200^3
nb=4 removal, exactly the Phase-2 1150/125000 buried count) — solve from
the prepared file bit-exact (max_abs 0, un/vn/wn/pn) vs the inline
analytic run, ransgeom dumps identical; prepare 1==4 ranks identical
files; GPU solve from the CPU-prepared file == CPU solve at tolerance 0
(and GPU inline == GPU from-file exactly on wavy_refine); 7-case suite
main_ref vs main bit-exact nofma CPU+GPU. Canonical prepare = the CPU
build (a GPU-built prepare computes coefficients on the device, libm-ulp
caveat).

**P1 — STL module. DONE 2026-07-16 (gates in `validation/prepare/`,
`run_gates_stl.sh`).**
`src/modules/geometry_stl.f90`: watertight binary-STL bodies behind the
analytic indicator signature — `body_indicator_i` MOVED to ibmm (walldist
re-exports it), `classify_*`/`fill_body_distance_analytic` take the
indicator as an argument (the solver passes `isInBody` explicitly — same
function, now via argument), and `set_ibm_coeff_host` is the host twin of
the device coefficient kernel over any indicator (declare-target procedure
arguments are not portable; the twins are marked KEEP IN LOCKSTEP). The
module: BVH (median split), majority-vote ray-parity inside test (3→11
fixed directions, deterministic large-rotation retry on degenerate casts —
small perturbations cannot escape the edge zone when a near-surface origin
hits at tiny t), minimum-image periodic queries, and the EXACT BVH
point-triangle distance for dwall (the indicator-driven walldist polish
was pathologically slow on STL — millions of near-surface parity casts —
and the exact query is what mobygeom/igl computes anyway: interior
agreement to round-off). `[ibm] stl_file` (whitespace-separated file list)
is a moby_prepare input; the solver hard-errors on it without a
`coeff_file`.
TWO CONVENTIONS DISCOVERED while gating (documented in geometry_stl.f90):
(1) distance queries image a periodic dim ONLY when the mesh is narrower
than the cell there — an STL spanning the full cell (padded wall slabs) is
its own periodic continuation and its overhanging skin is interior to the
periodic union, not a wall; membership (parity) always images. (2) dwall
ghost cells beyond a periodic boundary carry the periodic minimum-image
distance (the solver's analytic-walldist convention); mobygeom stores the
base-mesh distance there — interior cells agree to round-off.
*Gates (all PASS)*: flat + flat_refine vs the COMMITTED mobygeom
block-table references — blocks and every per-level mask IDENTICAL (incl.
768/56 buried leaves), coef 5.6e-9/2.1e-8 rel, interior dwall
7.5e-12/2.7e-11; 1-step solve prepared-vs-committed file ≤ 1e-10 fields,
ransgeom ≤ 1e-9, wallcell identical; sphere (curved, buried leaves,
all-periodic) vs a generated mobygeom reference — masks identical, coef
4.0e-7, interior dwall 1.7e-16; sphere exactly float32-translated onto the
x-periodic boundary — masks/blocks the exactly rolled copy and coef/dwall
tiles BIT-IDENTICAL (0.0), the zero-tolerance minimum-image gate; the P0
analytic gates all still PASS through the refactor (22 incl. the GPU
cross-check); 7-case suite bit-exact vs the P0 binaries CPU AND GPU.
DEFERRED to P1b: re-gating the big committed mobygeom cases (sailplane,
NACA, SD7003 — physical gates from prepare-built files) and the mobygeom
retirement; prepare's per-rank-redundant classification is the P2
parallelization target (flat_refine: 38 s at 4 ranks, dominated by
level-1 lattice parity casts).

**P2 — parallel + fast.**
MPI split as §7; far-field shortcut; optional direct segment-intersection
and exact-distance queries; OpenMP threading inside ranks if profiles ask.
*Gates*: 1==N ranks byte-identical output; timing vs mobygeom on the
NACA/B11 deep-refinement cases (target: minutes on one node where mobygeom
needed hours serial-Python).

**P3 — retire and rename.**
`main.f90` → `moby_solve.f90` (pure rename; compile.sh/module docs
updated). `mobygrid.f90` deleted (absorbed; the case file carries the
grid). mobygeom's geometry subcommands (`stl-ibm-coeff`, `block-active`,
`block-table` and `build_leaf_table_py` + Morton mirrors) get a RETIRED
header and survive only as the independent cross-implementation reference
for prepare's validation suite (the same role scipy kept for walldist);
STL generators/checkers (`make-*-stl`, watertightness tests) stay live.
*Gates*: 7-case suite bit-exact through the rename; tutorial READMEs
updated to the two-step flow.

## 9. Risks and open questions

- **Non-watertight / dirty STLs.** mobygeom silently repairs
  (trimesh `process=True`, `fix_normals`); prepare reads raw triangles.
  Parity casting tolerates unoriented normals but not holes or
  self-intersections. Policy: prepare validates (parity vote disagreement
  rate, AMPHIBIOUS's `contested` counter) and hard-errors pointing at
  mobygeom's repair/check tooling rather than repairing silently. The
  committed validation STLs are watertight; re-gate each.
- **Coefficient ulp differences vs mobygeom files.** Existing turbulent
  IBM results were produced from mobygeom files; prepare-built files will
  differ in last-bit coefficient values. The physical gates (not
  bit-compares) are the right acceptance criterion for P1(b)/(c) — plan
  the gate list accordingly and keep the mobygeom path readable until done.
- **`isInBody` device/host split.** The procedure-argument refactor must
  not disturb the GPU analytic kernel (declare-target constraint). The
  clean line: device kernel keeps the concrete function; only host-side
  callers generalize. Verify with the standard bit-exactness suite.
- **Prepare memory at L5/L6.** The windowed-mask machinery must be reused
  as-is (never materialize a finest-lattice raster); the leaf-parallel
  split bounds per-rank tile memory automatically. The known scale points
  (B11: ~7.2e9 finest-lattice blocks) are the sizing test.
- **Config drift between prepare and solve.** Both parse the same ini; the
  solver's existing cross-checks (blocks table row-by-row, `refine_dims`,
  `block_nb`, lx/ly/lz/re attrs) already catch stale files. Add the ini's
  geometry-relevant keys as case-file attributes echoed at solve init for
  friendlier diagnostics.

## 10. What this buys

- One language, one leaf builder, one grid definition — the Python mirror
  and the mobygrid handshake disappear.
- Geometry preprocessing becomes MPI-parallel (and trivially
  multi-node), replacing single-node Python multiprocessing.
- No venv/trimesh/libigl dependency for production runs.
- The solver init sheds nothing it needs and keeps every consistency
  cross-check; channels stay zero-preprocessing.
- mobygeom is retained as an independent reference implementation —
  exactly the validation asset a numerical-geometry port wants.
