# Mechanism probe: it is NOT the network, and the block tax is placement-dependent

Job 5135435, 4x HoreKa Green nodes (hkn[0529-0530,0533-0534]), commit
`a11e355`, clean tree, 200 steps/run. 18 of 23 runs completed — the `blk8` and
`nb16` tiers were **skipped because their configs were never staged** (the
package on the cluster predates the commit that added them; the driver reported
`MISSING CONFIG` and continued). Everything below rests on the 18 that ran.

## Validation

- **Placement honoured exactly** in all 18 runs: requested node count == nodes
  in `hosts.txt`, no `VOID` markers, no failures.
- **`L2_div` bit-identical** across every placement and rank count within a
  config (1.07282926E-05 for both single-level configs).
- **39 exchange rounds per step** in every run, as in the main matrix.
- **All runs stable**: marginal/cumulative rate 0.992–0.997, monotone flat.

## TIER 1 — the decisive pair: H-transport is dead

At 2 ranks each rank has exactly ONE peer, the pair is symmetric, and no stall
can propagate. Only the link's locality changes.

| config | 2x1 (intra) | 2x2 (cross) | own ratio |
|---|---|---|---|
| `rect_jacobi` (blocked) | 23.8 us | 67.5 us | 2.83x |
| `base_jacobi` (unblocked) | 20.3 us | 60.1 us | 2.95x |
| `refined_yp82` | 20.8 us | 50.3 us | 2.41x |

**Blocked / unblocked on the same single cross-node link: 1.12x.** Both pay the
same ~2.4–3x for crossing a node — ordinary fabric cost, identical for both. The
blocked path is **not** intrinsically slow cross-node.

## The link-count test kills the chain story too

`rect_jacobi`, 8 ranks, identical decomposition, identical volume, identical
rank count — only the layout differs:

| layout | cross-node links in the chain | `mpi_wait`/round | s/step |
|---|---|---|---|
| 4 ranks/node on **2** nodes | **1** | **752.6 us** | 0.13674 |
| 2 ranks/node on **4** nodes | **3** | **182.2 us** | 0.11311 |

**Three times the crossings, 4.1x faster.** No per-link transport cost and no
"one slow link stalls the chain" can produce that. Both hypotheses from the
probe's README are refuted.

## What does predict it

| case | GPUs/node | off-node traffic | local copy/rank | wait |
|---|---|---|---|---|
| `rect` 4 ranks, 4/node, 1 node | 4 | no | 3.79 M | 34.6 us |
| `rect` 8 ranks, 4/node, 2 nodes | 4 | **yes** | 1.86 M | **752.6 us** |
| `rect` 8 ranks, 2/node, 4 nodes | 2 | yes | 1.86 M | 182.2 us |
| `base` 8 ranks, 4/node, 2 nodes | 4 | yes | **0** | 101.2 us |

Slow only when **all three** are present: many ranks per node, off-node traffic,
and a large device-local copy load. Remove any one and it is fast. This is the
contention hypothesis, now with four supporting cells — and `blk8_jacobi`, which
removes *only* the local copy while holding everything else fixed, is exactly
the missing fifth. **It did not run. Until it does this remains inference, not
proof.**

## The block tax is largely a PLACEMENT artifact

Same 8 GPUs, same job, only the layout differs:

| layout | `base` s/step | `rect` s/step | **block tax** |
|---|---|---|---|
| 4 ranks/node, 2 nodes | 0.10249 | 0.13674 | **1.334** |
| 2 ranks/node, 4 nodes | 0.15825 | 0.11311 | **0.715** |

The main matrix ran 8 and 16 ranks at **4 ranks/node** — the bad regime for the
blocked case. Spread 2 per node, the blocked case is **faster than the
unblocked one**. The headline "block tax 1.315 at 8 ranks, 1.504 at 16" from
`results_horeka_2026-09-07.md` is therefore **specific to that layout**, not an
intrinsic property of blocking.

The 2:1 conclusion is unaffected: it was measured between two blocked configs
at identical placement, so the layout cancels.

## The two decompositions have OPPOSITE placement optima

| config | 4 ranks/node | 2 ranks/node | |
|---|---|---|---|
| `base` (3-D Cartesian, 7 peers) | 101.2 us | 1565.5 us | **15.5x worse spread** |
| `rect` (Morton chain, 2 peers) | 752.6 us | 182.2 us | **4.1x better spread** |

`base` is hurt by spreading — with 7 peers, thinning the ranks pushes nearly all
of them off-node. `rect` is hurt by packing — four GPUs per node each running a
large local-copy kernel while the node also drives off-node traffic. No single
layout serves both, which is itself a result: **placement must be chosen per
decomposition, and it was not a free parameter in any earlier measurement.**

## Consequences

1. **The partitioning rewrite is no longer the indicated fix.** It was proposed
   to break a 1-D chain believed to serialize on a slow link. That mechanism is
   refuted: more cross-node links made it faster.
2. **The cheapest real gain available today is a launch flag.** At 8 ranks,
   2/node is 17 % faster end-to-end than 4/node for the blocked case (0.1132 vs
   0.1369 s/step) on the same number of GPUs — at the cost of twice the nodes.
3. **The earlier scaling numbers need a placement-corrected re-run** before they
   are quoted again.

## Next, in order

1. **Stage the two missing configs and re-run** (`blk8`, `nb16`). The driver
   skips the 18 completed runs, so this is ~10 minutes. `blk8` at 8 ranks
   4/node is the confirmation: fast (~100 us) proves local-copy contention;
   slow (~750 us) refutes it and leaves peer-count/message structure.
2. **Re-run the main matrix at 2 ranks/node** to get placement-corrected block
   tax and scaling.
3. Only then choose the fix. If contention is confirmed, the target is
   scheduling — separating the local copy from the off-node transfer (different
   stream, or issued after the MPI posts) — which is a bit-exact rescheduling
   change, not a repartitioning.
