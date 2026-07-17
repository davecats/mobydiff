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
| `wavysolid.ini` | solid-block removal at scale (200^3, nb 4): the Phase-2 wavychannel geometry, 1150/125000 blocks buried -> compacted blocks table |

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
| `flat.ini`, `flat_refine.ini` | committed mobygeom block-table files (`../rans_geometry/ibm_coeff_blocks_l{1,2}.h5`, les_ibm wall slabs) | blocks + all per-level masks IDENTICAL; coef ≤ 1e-6 rel (indicator bisection vs exact ray crossings); interior dwall ≤ 2e-9 |
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
