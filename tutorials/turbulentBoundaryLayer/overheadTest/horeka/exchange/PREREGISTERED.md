# Pre-registered readings — the 2:1 exchange at scale

Written **before** job `moby_exch` ran, from the committed logs of job 5139351
(`results_scaling/new/`) and job 5138715 (`results_balance/`). This campaign has
twice published a mechanism that the next measurement overturned, so what the
numbers are expected to say is written down first and the report is graded
against it.

## What is already derivable from the committed logs

These are not predictions about new data — they are arithmetic on data that is
already in the repository, and they are what motivated the volume sweep added to
Pass 1. They are recorded here so the report cannot be read as having discovered
them in the new runs.

**(A) The refined case's exchange POINTS per cell are the same as the
single-level twin's.** At 16 ranks, `exchange sizes` gives 6.42 M points
(5.735 M local + 0.690 M send) for 60.56 M cells, against 15.39 M for 138.41 M
cells: 1.061e-1 vs 1.112e-1 points per cell, i.e. **0.954x**. The handout's
"1.9x exchange per cell" is therefore **1.9x TIME PER POINT**, not extra points.
Pre-registered branch 1 of the handout ("cross-level points are a small fraction
and the 1.9x is COPY volume") cannot fire as stated: there is no extra volume of
either kind.

**(B) `pack` and `unpack` are almost volume-independent.** Per exchange round,
across five configurations at 8 and 16 ranks, send points per rank span 43 k to
575 k (13x) while `pack` spans 99 to 168 us and `unpack` 86 to 125 us. An affine
fit gives pack = **92 us + 0.13 ns/pt**, unpack = **86 us + 0.07 ns/pt**.

**(C) `local_copy` has a fixed part too**, though smaller: `rect` and `nb16` at
8 ranks (1.86 M and 7.05 M local points per rank) give **31 us + 0.111 ns/pt**.

Summed, that is **~210 us per round of cost that no volume change can touch**,
against 39 rounds per step — 8.2 ms/step, or **20 % of the refined 16-rank step
(41.2 ms)** and 12 % of the single-level twin's (67.8 ms). It is the same
absolute number in both, which is exactly why the refined case looks 1.9x more
expensive per cell: it carries 44 % of the cells and pays the identical fixed
bill.

## Predictions for the new runs

| # | prediction | what refutes it |
|---|---|---|
| 1 | Cross-level (`restrict` + `prolong`) points are **≤ 20 %** of the refined case's exchange points at 1, 4 and 16 ranks. | a share above 20 % — then the interface transfer really is the volume, and the handout's branch 2 fires. |
| 2 | The op split is essentially **rank-count independent** (the interface is a geometric property of the leaf table). | a share that grows with rank count — then the Morton split is cutting the interface. |
| 3 | `copy_cross` costs **2–4x per point** what `local_copy` (same level) costs, so it takes 25–45 % of the copy time while carrying ~15 % of the points. | ≈ 1x — then the interface kernel is not special and the cost is granularity. |
| 4 | The volume sweep confirms (B)/(C) **within one allocation**: pack and unpack per round vary by less than 2x while their points vary by more than 10x. | per-round time tracking points — then the exchange is bandwidth-bound after all and the fixed-cost reading is wrong. |
| 5 | **A0 fails again** at 8 and 16 ranks: inserting a multi-millisecond compute kernel between the posts and the `Waitall` leaves `mpi_wait` at ≥ 0.9x its baseline. | `mpi_wait` collapsing — then MPI does progress during a kernel at these message sizes, P2 (overlap) reopens, and the 2-rank A0 conclusion was scope-limited exactly as it warned. |

## What would follow

- 1 + 3 + 4 ⇒ the lever is the **per-round fixed cost** (kernel launch and map-clause
  marshalling: ~1300 host-to-device memcpys per step were already seen by nsys and
  never chased), not block size, not the entry structure, and not the transport.
  That is the handout's third branch, reached by a different route than it
  anticipated.
- 2 refuted ⇒ partitioning is back on the table for the refined case specifically.
- 5 refuted ⇒ P2 reopens and everything else waits.
