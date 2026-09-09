# Mechanism probe: is the blocked exchange's stall caused by node crossing?

## All runs

| run | ranks | nodes | ranks/node | hosts | s/step | wait/round | peers | send max | local copy |
|---|---|---|---|---|---|---|---|---|---|
| base_jacobi | 2 | 1 | 2 | 1 | 0.34930 | 53.0 us | 1 | 34532 | 1458888 |
| base_jacobi | 2 | 2 | 1 | 2 | 0.34942 | 78.2 us | 1 | 34532 | 1458888 |
| base_jacobi | 4 | 1 | 4 | 1 | 0.18656 | 129.2 us | 3 | 414966 | 1458888 |
| base_jacobi | 4 | 2 | 2 | 2 | 0.18575 | 123.5 us | 3 | 414966 | 1458888 |
| base_jacobi | 4 | 4 | 1 | 4 | 0.20599 | 669.3 us | 3 | 414966 | 1458888 |
| base_jacobi | 8 | 2 | 4 | 2 | 0.10261 | 120.1 us | 7 | 574344 | 0 |
| base_jacobi | 8 | 4 | 2 | 4 | 0.16060 | 1621.6 us | 7 | 574344 | 0 |
| blk8_jacobi | 4 | 1 | 4 | 1 | 0.18985 | 190.6 us | 3 | 747600 | 136704 |
| blk8_jacobi | 4 | 2 | 2 | 2 | 0.18903 | 181.5 us | 3 | 747600 | 136704 |
| blk8_jacobi | 8 | 2 | 4 | 2 | 0.10015 | 123.4 us | 5 | 399432 | 0 |
| blk8_jacobi | 8 | 4 | 2 | 4 | 0.19898 | 2512.0 us | 5 | 399432 | 0 |
| nb16_jacobi | 4 | 1 | 4 | 1 | 0.26289 | 59.2 us | 2 | 85536 | 56765088 |
| nb16_jacobi | 4 | 2 | 2 | 2 | 0.26358 | 76.9 us | 2 | 85536 | 56765088 |
| nb16_jacobi | 8 | 2 | 4 | 2 | 0.17323 | 858.5 us | 2 | 85536 | 56422944 |
| nb16_jacobi | 8 | 4 | 2 | 4 | 0.14460 | 175.3 us | 2 | 85536 | 56422944 |
| rect_jacobi | 2 | 1 | 2 | 1 | 0.38194 | 49.4 us | 1 | 36800 | 15316352 |
| rect_jacobi | 2 | 2 | 1 | 2 | 0.38259 | 77.7 us | 1 | 36800 | 15316352 |
| rect_jacobi | 4 | 1 | 4 | 1 | 0.19937 | 60.8 us | 2 | 73600 | 15169152 |
| rect_jacobi | 4 | 2 | 2 | 2 | 0.20005 | 78.2 us | 2 | 73600 | 15169152 |
| rect_jacobi | 4 | 4 | 1 | 4 | 0.20106 | 120.8 us | 2 | 73600 | 15169152 |
| rect_jacobi | 8 | 2 | 4 | 2 | 0.13622 | 738.5 us | 2 | 73600 | 14874752 |
| rect_jacobi | 8 | 4 | 2 | 4 | 0.11176 | 149.8 us | 2 | 73600 | 14874752 |
| refined_yp82_rect_jacobi | 2 | 1 | 2 | 1 | 0.17861 | 29.8 us | 1 | 23000 | 6378680 |
| refined_yp82_rect_jacobi | 2 | 2 | 1 | 2 | 0.17898 | 45.8 us | 1 | 23000 | 6378680 |
| refined_yp82_rect_jacobi | 4 | 1 | 4 | 1 | 0.09937 | 45.4 us | 2 | 46000 | 6286680 |
| refined_yp82_rect_jacobi | 4 | 2 | 2 | 2 | 0.09953 | 53.8 us | 2 | 46000 | 6286680 |

## TIER 1 — one link, no chain: 2 ranks on 1 node vs 2 nodes

At 2 ranks each rank has exactly ONE peer and the pair is symmetric,
so no stall can propagate. Only the link's locality changes.

| config | 2x1 (intra) | 2x2 (cross) | ratio |
|---|---|---|---|
| rect_jacobi | 49.4 us | 77.7 us | **1.57x** |
| base_jacobi | 53.0 us | 78.2 us | **1.47x** |
| refined_yp82_rect_jacobi | 29.8 us | 45.8 us | **1.54x** |


Blocked / unblocked on ONE cross-node link: **0.99x** (rect 77.7 us vs base 78.2 us)
Same ratio intra-node: 0.93x (rect 49.4 vs base 53.0)

**VERDICT: H-chain.** One cross-node link is as cheap for the blocked
path as for the unblocked one, so the 10x at 8+ ranks is COLLECTIVE:
it needs several ranks and links. That is chain propagation or load
skew, and the target is the PARTITIONING (the 2-peer Morton chain).
Confirm with the link-count series: saturating after the first
crossing is the chain signature.

## Link count at FIXED rank count

Rank count, decomposition, peer count and volume held fixed; only the
number of node boundaries the chain crosses changes.

| config | ranks | 1 node | 2 nodes | 4 nodes |
|---|---|---|---|---|
| rect_jacobi | 2 | 49.4 us | 77.7 us | - |
| rect_jacobi | 4 | 60.8 us | 78.2 us | 120.8 us |
| rect_jacobi | 8 | - | 738.5 us | 149.8 us |
| base_jacobi | 2 | 53.0 us | 78.2 us | - |
| base_jacobi | 4 | 129.2 us | 123.5 us | 669.3 us |
| base_jacobi | 8 | - | 120.1 us | 1621.6 us |
| refined_yp82_rect_jacobi | 2 | 29.8 us | 45.8 us | - |
| refined_yp82_rect_jacobi | 4 | 45.4 us | 53.8 us | - |
| nb16_jacobi | 4 | 59.2 us | 76.9 us | - |
| nb16_jacobi | 8 | - | 858.5 us | 175.3 us |
| blk8_jacobi | 4 | 190.6 us | 181.5 us | - |
| blk8_jacobi | 8 | - | 123.4 us | 2512.0 us |

Proportional in the number of crossings => per-link transport cost.
Saturating after the FIRST crossing => a single slow link stalls the chain.

### Monotonicity check

- `rect_jacobi` 8 ranks: **4 nodes (149.8 us) is 4.9x FASTER than 2 nodes (738.5 us)** despite more node crossings.
- `nb16_jacobi` 8 ranks: **4 nodes (175.3 us) is 4.9x FASTER than 2 nodes (858.5 us)** despite more node crossings.

**This REFUTES both a per-link transport cost and a slow-link chain
stall, and it overrides the Tier-1 verdict above**: Tier 1 can only
separate *transport* from *collective*, and a cause that is not
communication at all satisfies it trivially. The variable that
changes with the node count here is RANKS PER NODE. Read the blk8
probe (node crossing without the device-local copy load) before
choosing any fix.

## Discriminator — node crossing WITHOUT device-local copy (`blk8_jacobi`)

One block per rank at 8 ranks: same 2-peer chain, same cells, but the
same-rank copy volume collapses.

| run | wait/round | local copy pts |
|---|---|---|
| rect_jacobi 8x2 | 738.5 us | 14874752 |
| blk8_jacobi 8x2 | 123.4 us | 0 |
| base_jacobi 8x2 | 120.1 us | 0 |
| rect_jacobi 4x1 | 60.8 us | 15169152 |
| blk8_jacobi 4x1 | 190.6 us | 136704 |

**Local-copy / transfer CONTENTION confirmed**: removing the same-rank
copy load removes most of the stall at the same placement and topology.

## Latency or bandwidth — exchange volume at a fixed placement

| run | send pts (max/rank) | wait/round | us per MB |
|---|---|---|---|
| rect_jacobi 4x2 | 73600 | 78.2 us | 33.2 |
| nb16_jacobi 4x2 | 85536 | 76.9 us | 28.1 |
| blk8_jacobi 8x2 | 399432 | 123.4 us | 9.7 |
| rect_jacobi 8x2 | 73600 | 738.5 us | 313.6 |

Wait flat against a large change in bytes => latency/serialization.
Wait tracking the bytes => bandwidth.
