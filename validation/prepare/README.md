# moby_prepare gates (docs/prepare_solve_strategy.md)

Two gate groups: `run_gates.sh` (P0, analytic geometry, bit-exactness) and
`run_gates_stl.sh` (P1, STL geometry vs mobygeom references + exact
shift-invariance). See the P1 section at the bottom.

## P0 — analytic geometry

P0 splits the analytic-geometry preprocessing out of the solver init into
the `moby_prepare` executable, with **zero new geometry code**: prepare runs
the solver's own classification/coefficient/wall-distance kernels and
writes one case file in the block-table coefficient-file format the solver
already reads (`[ibm] coeff_file`). The gate is therefore genuine
bit-exactness: the solve from the prepared file must equal the inline
analytic run to the last bit.

## Cases

| ini | exercises |
|---|---|
| `wavy.ini` | single-level blocks (nb 8), coef_blocks + block_active + dwall_blocks (`[rans]` geometry + ransgeom dump) |
| `wavy_refine.ini` | refine_body, 2 levels: per-level touch/buried masks + multi-level coef/dwall tiles |
| `wavysolid.ini` | solid-block removal at scale (200^3, nb 4): the Phase-2 analytic wavy-wall geometry, 1150/125000 blocks buried -> compacted blocks table |

## Run

```bash
./run_gates.sh                 # default build ../../build_cpu_nofma
./run_gates.sh ../../build_cpu
```

Per case: prepare on 1 and 4 ranks (case files identical — the Z-order row
split makes the file rank-count independent), inline analytic solve vs
solve-from-file (`tools/compare_fields.py --tolerance 0`), and for the
`[rans]` cases the `*_ransgeom.h5` dumps (dwall/yeff/wallcell) compared
exactly (`h5same.py`; h5diff is not installed everywhere).

## Status 2026-07-16 (P0, nofma builds)

All 20 CPU checks PASS: fields bit-exact (max_abs 0) on all three cases,
case files 1==4 prepare ranks identical, ransgeom dumps identical,
wavysolid removes exactly the Phase-2 count (1150/125000).

GPU (`GPU_BUILD=../../build_gpu_nofma ./run_gates.sh`): the GPU solve from
the CPU-prepared wavy_refine case file matches the CPU solve from the same
file at tolerance 0; on this case even the GPU INLINE run matches the
solve-from-file exactly (max_abs 0 incl. pn) and its ransgeom dump is
identical to the CPU inline one. The 7-case regression suite (main_ref vs
main) passes bit-exact on CPU AND GPU -- every solver-side P0 change
(writers, the fill_body_distance_analytic signature, set_serial_local_size)
is dormant in a normal solve.

Notes:
- prepare with the **CPU build** is canonical: the GPU build computes the
  coefficients on the device, whose libm ulps can differ from the host's.
- The solver still classifies/rebuilds its leaf table from the file masks
  and cross-checks the file's blocks table row-by-row — a stale or
  differently-built case file is a hard error, exactly as for mobygeom
  files.

## P1 — STL geometry (`run_gates_stl.sh`)

`[ibm] stl_file` (moby_prepare input only) loads watertight binary STLs
behind the analytic indicator signature (`geometry_stl.f90`: BVH +
majority-vote ray-parity inside test, exact BVH point-triangle wall
distance), so masks/coefficients/dwall flow through the SAME machinery as
analytic bodies.

| gate | reference | expectation |
|---|---|---|
| `flat.ini`, `flat_refine.ini` | mobygeom block-table files (`../rans_geometry/ibm_coeff_blocks_l{1,2}.h5`, les_ibm wall slabs; not in git — generate them with `../rans_geometry/setup.sh`) | blocks + all per-level masks IDENTICAL; coef ≤ 1e-6 rel (indicator bisection vs exact ray crossings); interior dwall ≤ 2e-9 |
| flat solve | 1-step solve, prepared file vs committed file | fields ≤ 1e-10, ransgeom dwall/yeff ≤ 1e-9, wallcell identical |
| `sphere.ini` | freshly generated mobygeom reference (needs the ibmc venv; skipped without it) | same identity/tolerance classes on a curved, buried-leaf body |
| `sphere_shift.ini` | the SAME mesh float32-EXACTLY translated onto the x-periodic boundary | masks/blocks the exactly rolled copy; coef/dwall tiles bit-identical (gates the minimum-image logic with zero tolerance) |

## Status 2026-07-16 (P1, CPU nofma build)

All gates PASS. Measured: flat/flat_refine/sphere blocks + masks identical
to mobygeom (incl. 768 and 56 buried leaves); coef 5.6e-9 / 2.1e-8 /
4.0e-7 relative; interior dwall 7.5e-12 / 2.7e-11 / 1.7e-16; the shift
gate is EXACTLY 0.0 on every dataset.

Conventions found while gating (documented in geometry_stl.f90):
- **dwall ghost cells beyond a periodic boundary**: prepare stores the
  periodic minimum-image distance (the solver's analytic-walldist
  convention); mobygeom stores the base-mesh distance. compare_case.py
  compares interior cells and reports the ghost gap informationally.
- **distance imaging rule**: a periodic dim is imaged only when the mesh
  is narrower than the cell there — an STL spanning the full cell (the
  padded wall slabs) is its own periodic continuation and its overhanging
  skin is interior to the periodic union, not a wall. Membership (parity)
  always images.
- STL dwall uses the exact BVH point-triangle query, NOT the
  indicator-driven walldist machinery (millions of near-surface parity
  casts; the exact query is also what mobygeom/igl computes — interior
  agreement to round-off).

## P1b — the big committed geometries (`run_gates_big.sh`)

> **Note (2026-09-26):** `validation/naca0012` and `validation/sd7003` were
> removed, so `run_gates_big.sh` now runs the sailplane leg only. The airfoil
> rows below are the recorded 2026-07-17 results; their drivers are in git
> history.

The production airfoil/sailplane cases re-gated from prepare-built files
(P1b additions: repeatable `[ibm] stl_file` for paths with spaces,
`stl_scale`/`stl_translate` = mobygeom's float64 `v*scale + translate`,
ASCII STL parsing straight to float64 like trimesh, `[blocks] keep_buried`
= mobygeom's `--keep-buried`, a solid-possible bbox cull + OpenMP in the
classify loops — the L5 airfoil lattices are 1.7e7 blocks).

Status 2026-07-17, all PASS (references regenerated with the retired
mobygeom, `--jobs 16`):

| case | identity gates | tolerance gates | solve gate |
|---|---|---|---|
| NACA 0012 L5 (`validation/naca0012`, keep_buried + dwall) | 25418 leaves + all 5 mask levels identical, 0 classification flips /76M | graded coef ≤ 1.2e-4 (16320 grazing cells of 76M), interior dwall 1.3e-10 | 200 GPU steps: fields ≤ 1.5e-7, forces 8-digit |
| SD7003 L5 (`validation/sd7003`) | 23836 leaves + masks identical, 0 flips /72M | coef ≤ 6.2e-5 (11680 grazing), dwall 1.1e-10 | 200 GPU steps: fields ≤ 1.2e-7 |
| sailplane (`tutorials/sailplane`, ASCII 55k-tri CAD mesh + scale/translate, nb=10) | 18000 leaves identical, 0 flips /93M | graded coef ≤ 1.2e-3 (2 grazing cells of 93M) | 1 step vs the COMMITTED legacy file: BIT-EXACT (max_abs 0) |

Timing: prepare beats mobygeom on the same case (SD7003 L5: 5m18s at 4
ranks x 4 threads vs 6m49s at `--jobs 16`); NACA L5 prepare 7m13s.

Notes:
- Graded-coefficient outliers are NEAR-GRAZING crossings: `((d0-d)/d)/d0²`
  has exploding relative sensitivity as d→d0 while the absolute value
  stays negligible (`mu` unaffected); the compare gates classification
  flips (must be 0) separately from the graded tolerance.
- Deep-refinement snapshots must be compared with `compare_snapshots.py`
  (chunked per-block): `tools/compare_fields.py` reassembles onto the
  finest lattice — 69 GB for the L5 airfoil grid — and gets OOM-killed.
- mobygeom's geometry subcommands are RETIRED for production
  (tools/README_mobygeom.md) and survive as the reference implementation
  these gates compare against.

## Step 7 — prepare does ALL preprocessing (`run_gates_step7.sh`)

Numerics review step 7 (`docs/next_session_prepare_everything.md`): a case
PREPARED by `moby_prepare` and SOLVED from the case file must reproduce the
reference binary's inline solve at `max_abs 0`, at PRODUCTION flags (the file
path copies values, nothing is recomputed). The driver covers body-free
cases with an explicit `[blocks] nb` (prepare on 1 == 4 ranks writes
identical files), body-free cases with NO `nb` (the nb rule: prepare on P
ranks stores `nb = grid/dims(P)` — different `blocks` tables, same fields;
the 1-rank file is one block and stops on 4 ranks with "rank owns no
blocks", rule 5) and an analytic-body case.

```bash
REF=~/step7_ref_binaries/moby_solve_cpu NEW=../../build_cpu ./run_gates_step7.sh
REF=~/step7_ref_binaries/moby_solve_gpu NEW=../../build_gpu PREP=../../build_cpu/moby_prepare \
    MODE=gpu SOLVE_RANKS=1 ./run_gates_step7.sh
```

Status 2026-09-29 (7-1/7-2, ref `d2249b1`): CPU 20/20 + wavy 3/3, GPU 16/16,
all `max_abs 0`. Derived sizes: wf180_y30 (8×6×8) → `8 6 8` on 1 rank,
`4 3 8` on 4 (an odd nb is accepted — parity only matters for refinement);
conduction (4×16×4) → `4 16 4` / `2 8 4`.

## The block order (`run_gates_order.sh`, 2026-10-01)

`docs/next_session_after_step9.md` item 2. The leaf table's row order — the
key the linear rank split cuts — is a bit permutation of the finest-lattice
block coordinates (`blocks.f90 leaf_key`). `moby_prepare` now chooses the
significance of the bits from the geometry (`min_surface_key_order`:
recursive bisection across the plane with the fewest cells, a periodic
direction's first cut counted twice) and records it in the case file as the
attribute `block_key_order` (directions 1..3, most to least significant). A
case file without the record — everything written before this date, and
every mobygeom file — is read with the legacy interleave (x lowest, z on
top), unchanged. `refine_dims = xz` keeps its own key and writes no record.
Field files are sliced by row in their writer's order; the restart reader
matches rows to blocks on (origin, level) (`field_hdf5.c block_row_map`), so
a snapshot or a generated IC restarts under any order.

```bash
OLD=~/step9b_ref_binaries/build_cpu NEW=../../build_cpu ./run_gates_order.sh
OLD=~/step9b_ref_binaries/build_gpu NEW=../../build_gpu ./run_gates_order.sh
```

`OLD` is a build of a tree that predates the record (it writes legacy-order
files). Three cases: the `min_channel` tutorial (xyz, two refine boxes, 800
leaves, order `x x y x y x y x y z`, 748 rows move), the analytic wavy wall
under `refine_body` (352 leaves, `y x x y x y x y z z`), and an nb-less
128 x 8 x 16 box prepared on 8 ranks (one block per rank, `x z y` — the
layout of the 2026-09-30 node-boundary finding; it carries the default flat
immersed wall, so it is a body case with non-cubic blocks).

| group | what | result (CPU and GPU builds) |
|---|---|---|
| `unit` | `keyorder_test`: the rule's bit orders on seven layouts, and the builder + `partition_surface` on a scaled `base_jacobi` (8 ranks on 2 nodes: 4224 cells across ranks, 128 across nodes) | PASS |
| `order` | legacy-order file vs re-prepared file: the same leaves; the builder's row order equals the Python mirror's (`tools/partition_analysis.py --order minsurface`), and the mirror's legacy order equals the old file's; fields `max_abs 0` between both files at every rank count (1 and 4; 8) and against the old solver; the body cases' coefficient and wall-distance tiles identical at tolerance 0 once rows are matched | PASS |
| `print` | the solver's `partition: face cells shared across ranks …` line equals `partition_analysis.py --order file` on its case file, incl. the coarse/fine faces of the two refined cases | PASS |
| `restart` | a snapshot written under either order continues identically under both (`max_abs 0`); the cross-order restart prints its note, an own-order one does not | PASS |
| `tool` | `tools/h5maxdiff` (rows matched on origin + level) reads `max_abs 0` across the two orders | PASS |

43 checks on CPU, 43 on GPU. The comparison helpers follow the same rule:
`compare_case.py` and `compare_snapshots.py` match rows through `h5rows.py`
(a prepared xyz file against a committed mobygeom reference differs in row
order now), and so does the `ransgeom` check of `run_gates_stl.sh`. Re-run
with this tree: P0 `run_gates.sh` 26/26, P1 `run_gates_stl.sh` 15/15
(flat / flat_refine vs the committed mobygeom files, sphere vs a generated
one), and the nofma 7-case + 9-case suites against `~/step9b_ref_binaries`
at `max_abs 0`, CPU (4 ranks / 1 rank) and GPU — the new side re-prepares
its case files in the new order and restarts the suites' legacy-order ICs
through the row map.

**`run_gates_step7.sh` is no longer a `max_abs 0` gate, with or without this
change.** Its reference is the step 7-0 binary, the only one that still
builds inline, and that binary predates F1 (step 8) and the penalization
factors (step 9): 18 checks pass and 26 fail with this tree AND, identically
check for check, with `~/step9b_ref_binaries` (control run 2026-10-01). What
it gated — a prepared case file solving like the inline build — is covered
since by the suites' prepare-per-side legs and by `run_gates_order.sh`.

The cross-NODE number of the `partition:` line cannot be exercised on one
workstation; it is unit-tested with a given node map and measured on HoreKa
(`tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange/PREREGISTERED_order.md`).
