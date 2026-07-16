# moby_prepare P0 gates (docs/prepare_solve_strategy.md)

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
