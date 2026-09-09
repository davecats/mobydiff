# Mechanism probe, round 2: all three hypotheses are dead, and placement dominates everything

Job 5137106 completed the five runs skipped last time (`blk8`, `nb16`); the
other 18 were carried over by the resume logic from job 5135435. Commit
`a11e355`, clean, placement honoured exactly in all 23, no failures.

**Caveat, stated up front:** the two groups ran on different nodes
(hkn[0529-0530,0533-0534] vs hkn[0512-0513,0520,0524]) two days apart. Node-to-
node variation measured across earlier jobs is ~6 %; every effect below is 4x
or larger, so the cross-job comparison is safe — but it is a cross-job
comparison, and the confirmation run proposed at the end should re-measure the
`rect` reference in the same allocation.

## Hypothesis 1 — per-link transport cost: REFUTED

`rect_jacobi`, 8 ranks, identical decomposition and volume:

| layout | cross-node links | wait/round |
|---|---|---|
| 4 ranks/node, 2 nodes | 1 | 752.6 us |
| 2 ranks/node, 4 nodes | 3 | 182.2 us |

Three times the crossings, 4.1x faster.

## Hypothesis 2 — a slow link stalling the 1-D chain: REFUTED

Tier 1, two symmetric ranks with one peer each and no chain to propagate
anything: blocked/unblocked on the same single cross-node link is **1.12x**
(67.5 vs 60.1 us). Both pay the same ~2.4-3x to cross a node. The blocked path
is not intrinsically slow off-node, and the link-count result above kills the
propagation story.

## Hypothesis 3 — device-local copy contending with off-node transfer: REFUTED

This was the one `blk8_jacobi` was built to confirm, and at first sight it does:
at 8 ranks / 4 per node — the anomaly's exact regime — `blk8` (zero local copy)
runs at **114.1 us** against `rect`'s **752.6 us**, while being *worse* on every
other axis (5 peers vs 2, 5.4x the send volume). Copy load looks decisive.

**`nb16_jacobi` refutes it.** At matched placement (4 ranks, 2 nodes, off-node
traffic present):

| run | local copy / rank | copy kernel | wait/round |
|---|---|---|---|
| `rect_jacobi` | 3.79 M pts | 16.67 ms | 66.4 us |
| `nb16_jacobi` | **14.19 M pts** | **61.56 ms** | **81.1 us** |

**3.7x the copy load buys 1.22x the wait.** A mechanism that turns copy volume
into stall cannot produce both this and the 6.6x gap above.

The honest reading of `blk8` is that it changes copy **and** peer count **and**
volume simultaneously, so it cannot attribute the 6.6x to any one of them —
and `nb16`, the one config that varies copy alone, was never run in the regime
where the anomaly lives.

## What the data actually shows: placement dominates, and the optimum inverts

| config | peers | 2 nodes | 4 nodes | effect of spreading |
|---|---|---|---|---|
| `base_jacobi` 8 ranks | 7 | 101.2 us | 1565.5 us | **15.5x worse** |
| `blk8_jacobi` 8 ranks | 5 | 114.1 us | 2597.4 us | **22.8x worse** |
| `rect_jacobi` 8 ranks | 2 | 752.6 us | 182.2 us | **4.1x better** |

Decompositions with many peers collapse when their ranks are spread thin —
a **15-23x** penalty, far larger than anything in the original report. The
2-peer Morton chain does the opposite. There is no placement that serves both.

Two further facts worth keeping:

- The `rect` 8-rank / 4-per-node value is **reproducible across two jobs, two
  node sets and two mapping policies** (708.3 us with `--map-by numa` in job
  5133554; 752.6 us with `ppr:4:node` here), so it is a real property of that
  configuration, not a scheduling accident.
- `blk8` at 4 ranks on one node is *slower* than `rect` there (169.5 vs
  34.6 us) despite zero copy — consistent with peer count mattering when there
  is no network at all.

## Consequence

**No exchange measurement in this project is interpretable without stating the
rank-to-node placement.** Placement moves `mpi_wait` by up to 23x — more than
blocking, more than refinement, more than the smoother, more than every
optimisation measured to date combined. The 2026-09-07 scaling numbers were all
taken at one arbitrary placement (4 ranks/node, the SLURM default from
`--ntasks-per-node=4`) and must be re-read as "at that placement". The 2:1
result survives unchanged, because it compared two blocked configs at identical
placement.

## The three runs that would close this — ~5 minutes

The missing cells are all in the anomaly's regime (8 ranks, 4 per node):

1. **`nb16_jacobi:8:2`** — the direct test. Same placement as the anomaly, 7.6x
   `rect`'s copy load per rank, same 2-peer chain. If copy drives the stall this
   must be *worse* than 752 us; if it lands near `rect` or below, copy is
   conclusively out and peer count is the remaining candidate.
2. **`blk8_jacobi:4:2`** — a zero-copy point in the mild regime, to complete the
   2x2 of (copy, placement).
3. **`rect_jacobi:8:2` re-run in the same allocation** — removes the cross-job
   caveat from every comparison above.

These are already added to `run_placement.sh`'s `SPECS`; the resume logic skips
the 23 completed runs, so a resubmission runs only these.

## What I would NOT do next

Another black-box placement sweep. Three mechanisms have now been proposed and
killed by data, which is progress, but the remaining question — *what is a rank
actually waiting for* — is not answerable from aggregate `mpi_wait`. The next
step after the three runs above is a small, bit-exact, env-gated diagnostic:
**per-rank wait min/max/argmax** (is one rank stalling everyone, or are all
ranks equally slow?) and a barrier immediately before the `Waitall` to split
"waiting for a slow peer" from "slow transfer". That distinction decides
between a partitioning fix and a transport fix, and no amount of placement
sweeping will reveal it.
