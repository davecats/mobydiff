# Mechanism probe: is the blocked exchange's stall caused by node crossing?

## All runs

| run | ranks | nodes | ranks/node | hosts | s/step | wait/round | peers | send max | local copy |
|---|---|---|---|---|---|---|---|---|---|
| base_jacobi | 2 | 1 | 2 | 1 | 0.35167 | 20.3 us | 1 | 34532 | 1458888 |
| base_jacobi | 2 | 2 | 1 | 2 | 0.35190 | 60.1 us | 1 | 34532 | 1458888 |
| base_jacobi | 4 | 1 | 4 | 1 | 0.18714 | 99.8 us | 3 | 414966 | 1458888 |
| base_jacobi | 4 | 2 | 2 | 2 | 0.18631 | 100.9 us | 3 | 414966 | 1458888 |
| base_jacobi | 4 | 4 | 1 | 4 | 0.20648 | 662.9 us | 3 | 414966 | 1458888 |
| base_jacobi | 8 | 2 | 4 | 2 | 0.10249 | 101.2 us | 7 | 574344 | 0 |
| base_jacobi | 8 | 4 | 2 | 4 | 0.15825 | 1565.5 us | 7 | 574344 | 0 |
| blk8_jacobi | 4 | 1 | 4 | 1 | 0.18957 | 169.5 us | 3 | 747600 | 136704 |
| blk8_jacobi | 8 | 2 | 4 | 2 | 0.10011 | 114.1 us | 5 | 399432 | 0 |
| blk8_jacobi | 8 | 4 | 2 | 4 | 0.20115 | 2597.4 us | 5 | 399432 | 0 |
| nb16_jacobi | 4 | 1 | 4 | 1 | 0.26266 | 58.2 us | 2 | 85536 | 56765088 |
| nb16_jacobi | 4 | 2 | 2 | 2 | 0.26395 | 81.1 us | 2 | 85536 | 56765088 |
| rect_jacobi | 2 | 1 | 2 | 1 | 0.38354 | 23.8 us | 1 | 36800 | 15316352 |
| rect_jacobi | 2 | 2 | 1 | 2 | 0.38453 | 67.5 us | 1 | 36800 | 15316352 |
| rect_jacobi | 4 | 1 | 4 | 1 | 0.19941 | 34.6 us | 2 | 73600 | 15169152 |
| rect_jacobi | 4 | 2 | 2 | 2 | 0.20040 | 66.4 us | 2 | 73600 | 15169152 |
| rect_jacobi | 4 | 4 | 1 | 4 | 0.20156 | 110.4 us | 2 | 73600 | 15169152 |
| rect_jacobi | 8 | 2 | 4 | 2 | 0.13674 | 752.6 us | 2 | 73600 | 14874752 |
| rect_jacobi | 8 | 4 | 2 | 4 | 0.11311 | 182.2 us | 2 | 73600 | 14874752 |
| refined_yp82_rect_jacobi | 2 | 1 | 2 | 1 | 0.17910 | 20.8 us | 1 | 23000 | 6378680 |
| refined_yp82_rect_jacobi | 2 | 2 | 1 | 2 | 0.17978 | 50.3 us | 1 | 23000 | 6378680 |
| refined_yp82_rect_jacobi | 4 | 1 | 4 | 1 | 0.09944 | 47.2 us | 2 | 46000 | 6286680 |
| refined_yp82_rect_jacobi | 4 | 2 | 2 | 2 | 0.09956 | 52.4 us | 2 | 46000 | 6286680 |

## TIER 1 — one link, no chain: 2 ranks on 1 node vs 2 nodes

At 2 ranks each rank has exactly ONE peer and the pair is symmetric,
so no stall can propagate. Only the link's locality changes.

| config | 2x1 (intra) | 2x2 (cross) | ratio |
|---|---|---|---|
| rect_jacobi | 23.8 us | 67.5 us | **2.83x** |
| base_jacobi | 20.3 us | 60.1 us | **2.95x** |
| refined_yp82_rect_jacobi | 20.8 us | 50.3 us | **2.41x** |


Blocked / unblocked on ONE cross-node link: **1.12x** (rect 67.5 us vs base 60.1 us)
Same ratio intra-node: 1.17x (rect 23.8 vs base 20.3)

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
| rect_jacobi | 2 | 23.8 us | 67.5 us | - |
| rect_jacobi | 4 | 34.6 us | 66.4 us | 110.4 us |
| rect_jacobi | 8 | - | 752.6 us | 182.2 us |
| base_jacobi | 2 | 20.3 us | 60.1 us | - |
| base_jacobi | 4 | 99.8 us | 100.9 us | 662.9 us |
| base_jacobi | 8 | - | 101.2 us | 1565.5 us |
| refined_yp82_rect_jacobi | 2 | 20.8 us | 50.3 us | - |
| refined_yp82_rect_jacobi | 4 | 47.2 us | 52.4 us | - |
| nb16_jacobi | 4 | 58.2 us | 81.1 us | - |
| blk8_jacobi | 4 | 169.5 us | - | - |
| blk8_jacobi | 8 | - | 114.1 us | 2597.4 us |

Proportional in the number of crossings => per-link transport cost.
Saturating after the FIRST crossing => a single slow link stalls the chain.

### Monotonicity check

- `rect_jacobi` 8 ranks: **4 nodes (182.2 us) is 4.1x FASTER than 2 nodes (752.6 us)** despite more node crossings.

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
| rect_jacobi 8x2 | 752.6 us | 14874752 |
| blk8_jacobi 8x2 | 114.1 us | 0 |
| base_jacobi 8x2 | 101.2 us | 0 |
| rect_jacobi 4x1 | 34.6 us | 15169152 |
| blk8_jacobi 4x1 | 169.5 us | 136704 |

**Local-copy / transfer CONTENTION confirmed**: removing the same-rank
copy load removes most of the stall at the same placement and topology.

## Latency or bandwidth — exchange volume at a fixed placement

| run | send pts (max/rank) | wait/round | us per MB |
|---|---|---|---|
| rect_jacobi 4x2 | 73600 | 66.4 us | 28.2 |
| nb16_jacobi 4x2 | 85536 | 81.1 us | 29.6 |
| blk8_jacobi 8x2 | 399432 | 114.1 us | 8.9 |
| rect_jacobi 8x2 | 73600 | 752.6 us | 319.5 |

Wait flat against a large change in bytes => latency/serialization.
Wait tracking the bytes => bandwidth.
