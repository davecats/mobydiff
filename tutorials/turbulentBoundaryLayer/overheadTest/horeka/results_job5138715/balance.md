# What a rank is waiting for — the be48d44 exchange diagnostics

## Pass A — timing and per-rank wait balance

`wait r0` is the rank-0 bucket every earlier report quoted; min/mean/max
are that same accumulated `mpi_wait` reduced across ALL ranks, per round.

| run | r x N | /node | s/step | wait r0 | min | mean | max | max/min | argmax | peers | copy pts |
|---|---|---|---|---|---|---|---|---|---|---|---|
| base_jacobi | 2x1 | 2 | 0.34930 | 53.0 | 53.0 | 62.4 | 71.9 | 1.35 | 1 | 1 | 1458888 |
| base_jacobi | 2x2 | 1 | 0.34942 | 78.2 | 78.2 | 94.2 | 110.3 | 1.41 | 1 | 1 | 1458888 |
| base_jacobi | 4x1 | 4 | 0.18656 | 129.2 | 124.8 | 153.3 | 189.3 | 1.52 | 2 | 3 | 1458888 |
| base_jacobi | 4x2 | 2 | 0.18575 | 123.5 | 123.5 | 144.8 | 176.5 | 1.43 | 2 | 3 | 1458888 |
| base_jacobi | 4x4 | 1 | 0.20599 | 669.3 | 669.3 | 698.8 | 722.9 | 1.08 | 2 | 3 | 1458888 |
| base_jacobi | 8x2 | 4 | 0.10261 | 120.1 | 120.1 | 130.6 | 143.6 | 1.20 | 4 | 7 | 0 |
| base_jacobi | 8x4 | 2 | 0.16060 | 1621.6 | 1621.6 | 1638.5 | 1647.9 | 1.02 | 5 | 7 | 0 |
| blk8_jacobi | 4x1 | 4 | 0.18985 | 190.6 | 170.7 | 192.5 | 218.2 | 1.28 | 2 | 3 | 136704 |
| blk8_jacobi | 4x2 | 2 | 0.18903 | 181.5 | 173.1 | 186.9 | 209.9 | 1.21 | 2 | 3 | 136704 |
| blk8_jacobi | 8x2 | 4 | 0.10015 | 123.4 | 111.8 | 127.3 | 139.9 | 1.25 | 5 | 5 | 0 |
| blk8_jacobi | 8x4 | 2 | 0.19898 | 2512.0 | 2512.0 | 2610.2 | 2696.1 | 1.07 | 5 | 5 | 0 |
| nb16_jacobi | 4x1 | 4 | 0.26289 | 59.2 | 59.2 | 74.9 | 85.7 | 1.45 | 2 | 2 | 56765088 |
| nb16_jacobi | 4x2 | 2 | 0.26358 | 76.9 | 76.9 | 103.1 | 124.0 | 1.61 | 2 | 2 | 56765088 |
| nb16_jacobi | 8x2 | 4 | 0.17323 | 858.5 | 855.2 | 890.7 | 925.9 | 1.08 | 4 | 2 | 56422944 |
| nb16_jacobi | 8x4 | 2 | 0.14460 | 175.3 | 175.3 | 200.8 | 206.8 | 1.18 | 7 | 2 | 56422944 |
| rect_jacobi | 2x1 | 2 | 0.38194 | 49.4 | 49.4 | 60.8 | 72.2 | 1.46 | 1 | 1 | 15316352 |
| rect_jacobi | 2x2 | 1 | 0.38259 | 77.7 | 77.7 | 94.7 | 111.7 | 1.44 | 1 | 1 | 15316352 |
| rect_jacobi | 4x1 | 4 | 0.19937 | 60.8 | 60.8 | 72.7 | 83.6 | 1.38 | 1 | 2 | 15169152 |
| rect_jacobi | 4x2 | 2 | 0.20005 | 78.2 | 78.2 | 94.9 | 107.9 | 1.38 | 2 | 2 | 15169152 |
| rect_jacobi | 4x4 | 1 | 0.20106 | 120.8 | 120.8 | 138.8 | 146.5 | 1.21 | 1 | 2 | 15169152 |
| rect_jacobi | 8x2 | 4 | 0.13622 | 738.5 | 734.8 | 768.0 | 802.3 | 1.09 | 4 | 2 | 14874752 |
| rect_jacobi | 8x4 | 2 | 0.11176 | 149.8 | 149.8 | 177.2 | 185.2 | 1.24 | 7 | 2 | 14874752 |
| refined_yp82_rect_jacobi | 2x1 | 2 | 0.17861 | 29.8 | 29.8 | 46.4 | 62.9 | 2.11 | 1 | 1 | 6378680 |
| refined_yp82_rect_jacobi | 2x2 | 1 | 0.17898 | 45.8 | 45.8 | 67.9 | 90.1 | 1.96 | 1 | 1 | 6378680 |
| refined_yp82_rect_jacobi | 4x1 | 4 | 0.09937 | 45.4 | 45.4 | 59.0 | 73.3 | 1.61 | 1 | 2 | 6286680 |
| refined_yp82_rect_jacobi | 4x2 | 2 | 0.09953 | 53.8 | 53.8 | 73.7 | 84.9 | 1.58 | 2 | 2 | 6286680 |

### Rank independence (`L2_div` must match across placements)

All configs: identical `L2_div` across every placement at each rank count.

## Pass B — barrier: arrival SKEW vs TRANSFER

`exchange_barrier = true` puts an MPI_Barrier before every Waitall, so
the barrier absorbs the arrival skew (`skew_barrier`) and what is left
in `mpi_wait` is transfer. **Serialising — its step times mean nothing.**

WHAT THIS DOES AND DOES NOT PIN DOWN. A barrier costing X us means the
ranks genuinely arrive ~X us apart (an 8-rank barrier on arrived ranks
costs tens of us, not hundreds) -- so the SKEW column is a direct
measurement. The transfer column is only a LOWER bound on wire time: a
transfer that finished during the barrier leaves ~0 behind, so a small
residual proves the transfer is not the critical path, not that it was
fast in isolation. A LARGE residual is the informative case -- that much
transfer remained after every rank had arrived.

| run | r x N | skew/round | transfer/round | skew share | pass-A wait r0 |
|---|---|---|---|---|---|
| base_jacobi | 8x2 | 64.8 us | 63.8 us | 50 % | 120.1 us |
| blk8_jacobi | 8x2 | 83.1 us | 49.1 us | 63 % | 123.4 us |
| nb16_jacobi | 8x2 | 878.6 us | 1.2 us | 100 % | 858.5 us |
| rect_jacobi | 4x1 | 63.8 us | 8.7 us | 88 % | 60.8 us |
| rect_jacobi | 8x2 | 773.0 us | 1.2 us | 100 % | 738.5 us |
| rect_jacobi | 8x4 | 160.0 us | 7.2 us | 96 % | 149.8 us |

## Verdict, against the pre-registered table

| run | max/min | skew vs transfer | reading |
|---|---|---|---|
| base_jacobi 2x1 | 1.35 (rank 1) | - | no barrier pass |
| base_jacobi 2x2 | 1.41 (rank 1) | - | no barrier pass |
| base_jacobi 4x1 | 1.52 (rank 2) | - | no barrier pass |
| base_jacobi 4x2 | 1.43 (rank 2) | - | no barrier pass |
| base_jacobi 4x4 | 1.08 (rank 2) | - | no barrier pass |
| base_jacobi 8x2 | 1.20 (rank 4) | 65 / 64 us | **transport** — equal totals and real transfer left after synchronisation; fix in the exchange, not the partitioning |
| base_jacobi 8x4 | 1.02 (rank 5) | - | no barrier pass |
| blk8_jacobi 4x1 | 1.28 (rank 2) | - | no barrier pass |
| blk8_jacobi 4x2 | 1.21 (rank 2) | - | no barrier pass |
| blk8_jacobi 8x2 | 1.25 (rank 5) | 83 / 49 us | **transport** — equal totals and real transfer left after synchronisation; fix in the exchange, not the partitioning |
| blk8_jacobi 8x4 | 1.07 (rank 5) | - | no barrier pass |
| nb16_jacobi 4x1 | 1.45 (rank 2) | - | no barrier pass |
| nb16_jacobi 4x2 | 1.61 (rank 2) | - | no barrier pass |
| nb16_jacobi 8x2 | 1.08 (rank 4) | 879 / 1 us | **rotating skew** — ~879 us of per-round arrival spread, but equal per-rank TOTALS (max/min 1.08): the late rank rotates. Not a fixed load imbalance, and not transport. |
| nb16_jacobi 8x4 | 1.18 (rank 7) | - | no barrier pass |
| rect_jacobi 2x1 | 1.46 (rank 1) | - | no barrier pass |
| rect_jacobi 2x2 | 1.44 (rank 1) | - | no barrier pass |
| rect_jacobi 4x1 | 1.38 (rank 1) | 64 / 9 us | **rotating skew** — ~64 us of per-round arrival spread, but equal per-rank TOTALS (max/min 1.38): the late rank rotates. Not a fixed load imbalance, and not transport. |
| rect_jacobi 4x2 | 1.38 (rank 2) | - | no barrier pass |
| rect_jacobi 4x4 | 1.21 (rank 1) | - | no barrier pass |
| rect_jacobi 8x2 | 1.09 (rank 4) | 773 / 1 us | **rotating skew** — ~773 us of per-round arrival spread, but equal per-rank TOTALS (max/min 1.09): the late rank rotates. Not a fixed load imbalance, and not transport. |
| rect_jacobi 8x4 | 1.24 (rank 7) | 160 / 7 us | **rotating skew** — ~160 us of per-round arrival spread, but equal per-rank TOTALS (max/min 1.24): the late rank rotates. Not a fixed load imbalance, and not transport. |
| refined_yp82_rect_jacobi 2x1 | 2.11 (rank 1) | - | no barrier pass |
| refined_yp82_rect_jacobi 2x2 | 1.96 (rank 1) | - | no barrier pass |
| refined_yp82_rect_jacobi 4x1 | 1.61 (rank 1) | - | no barrier pass |
| refined_yp82_rect_jacobi 4x2 | 1.58 (rank 2) | - | no barrier pass |
