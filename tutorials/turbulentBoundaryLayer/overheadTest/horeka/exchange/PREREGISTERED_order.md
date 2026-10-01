# Pre-registered: the minimum-surface block order (written before job submission, 2026-10-01)

**Change.** `blocks.f90 leaf_key` (xyz mode) is a bit permutation whose bit
significance is chosen once per case by `min_surface_key_order` — recursive
bisection across the plane that carries the fewest cells, a periodic
direction's first cut counted twice — and recorded in the case file
(`block_key_order`). A file without the record is read with the legacy
interleave (x lowest, z on top). xz mode keeps its own key. Restart readers
match rows on (origin, level). No kernel is touched: the arithmetic source of
the two binaries is byte-identical, only block OWNERSHIP moves.

**Why.** `results_horeka_2026-09-30.md` section 3: since step 7 an nb-less
case is one block per rank in Morton order, whose top bit is z whatever the
geometry. `base_jacobi` at 8 ranks on two nodes (2 x 2 x 2 blocks of
2048 x 88 x 96, z periodic) puts the 2048 x 88 z-faces on the node boundary:
1,441,792 face cells where an x cut takes 33,792; `mpi_wait` 3.3 → 58.8 ms,
step 77.7 → 132.0 ms; `--map-by node` on the same binary and file: 77.6 ms.

**Binaries.** ref = `5bdc5eb` (the head before the change; its case files
are legacy order), new = the order commit. Both re-prepare their own case
file in every run directory, so each column runs its own order. Shipped
configs on both sides, `--map-by numa --bind-to core` (the DEFAULT
placement), 200 steps.

**Predictions.**

1. The solver's init line, new column (the ref binary has no such line):
   `base_jacobi` n=8 `partition: face cells shared across ranks 2262016
   (8 ranks), across nodes 33792 (2 nodes)`; `rect_jacobi` n=8 across ranks
   236544, across nodes 33792; every n=4 run `across nodes 0 (1 nodes)`.
   (tools/partition_analysis.py: the legacy order reads 1,441,792 across
   nodes for `base_jacobi` and the same 33,792 for `rect_jacobi`.)
2. `base_jacobi` n=8: ref **132 ms ± 1 %**, `mpi_wait` ≈ 59 ms (the 09-30
   number reproduced); new **77.7 ms ± 1 %**, `mpi_wait` **≈ 3.3 ms**
   (3.0–3.6).
3. `base_jacobi` n=4 (one node, no cross-node link): unchanged within 0.5 %.
4. `rect_jacobi` n=4 and n=8: unchanged within 0.5 %. Its order changes
   (x x x x x y x z z y against the interleave) but the lattice's long
   direction already put x on the coarse cuts, so the cross-node count is the
   same 33,792.
5. `refined_yp82_rect_jacobi` n=4 and n=8 (`refine_dims = xz`): the key is
   untouched and the case file carries no record: unchanged within 0.5 %.
6. Field gate at PRODUCTION flags (nothing is recomputed, values are owned by
   different ranks): `base_jacobi` and `rect_jacobi`, n=8, 20 steps, final
   snapshot ref vs new through `tools/h5maxdiff`, which must report that the
   row orders differ and **max_abs 0** on un/vn/wn/pn.
7. 16 ranks (4 nodes, a separate submission on `accelerated`; never measured
   since step 7): `base_jacobi` is 4 x 2 x 2 blocks of 1024 x 88 x 96, order
   x x y z, across nodes 101,376 (legacy 1,475,584). New: the block tax
   `rect / base` back inside the published 16-rank range, 1.00–1.03
   (09-23: 1.007, 09-25: 1.030); ref: `base_jacobi` wait-dominated, at least
   30 % slower than new.

**What would falsify it.** `base_jacobi` n=8 new still near 132 ms with the
init line reading 33,792: the cross-node face count is not what sets
`mpi_wait` and the 09-30 placement test was explained by something else
(GPU affinity classes of the two ends, see the 2026-09-10 finding) — then
compare the `gpu binding` lines of the two runs. The init line reading
1,441,792 on the new column: the case file was not re-prepared (a stale
legacy file in the run directory) — check `prepare.log`.
