# HoreKa: the 2:1 machinery is free, the RANK SPLIT is the problem

4x HoreKa Green nodes = 16x A100-SXM4-40GB, job 5133554, commit `a11e355`
(clean tree), 23 runs of 400 steps, no failures. Package + driver:
`overheadTest/horeka/`; raw logs in its `moby-2to1-run/results/`.

## Validation before any number

- **Provenance**: commit `a11e355`, `dirty 0`, nvfortran 25.3, A100-SXM4-40GB.
- **Leaf gates EXACT**: 448 leaves (small refined) and 1792 (big refined) in
  every run — the refinement placement did not move.
- **Rank binding**: every run's Cartesian dims match its directory name.
- **Rank independence**: `L2_div` and `Linf` are **bit-identical across 1..16
  ranks** within each config, and `rect_jacobi` equals `base_jacobi` exactly —
  blocking does not change the answer.
- **Drift**: cumulative average moves < 0.6 % over the last 200 steps in all 23
  runs; marginal/cumulative rate 0.994–0.999. Dedicated nodes, no contention.

## Headline 1 — the 2:1 interface machinery costs nothing

`refined_big_rect_jacobi` shares level-0 grid, block shape, refine box and
solver with `rect_jacobi` exactly; it just refines y-tile 1, turning 1024
leaves into 1792 and 138.41 M cells into 242.22 M. **The 103.81 M added cells
cost the same as, or less than, the coarse cells they sit beside:**

| ranks | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell average |
|---|---|---|---|---|
| 4 | 0.19978 | 0.34642 | 1.4126 | **0.979x** |
| 8 | 0.13487 | 0.21466 | 0.7686 | **0.789x** |
| 16 | 0.08971 | 0.13520 | 0.4382 | **0.676x** |

1.75x the cells for 1.73x / 1.59x / 1.51x the time. The ratio falls below 1 at
8 and 16 ranks because the refined case keeps more work per GPU while the
single-level twin starves (8.65 M cells/GPU at 16 ranks), so read the 4-rank
number, 0.979, as the honest one: **within measurement error the 2:1 interface
adds no marginal cost.** This is the first like-for-like measurement of it —
every earlier figure compared two different level-0 grids.

## Headline 2 — the block tax explodes at the node boundary

`rect_jacobi` vs `base_jacobi` (identical grid, identical results, blocking is
the only difference), inverse per-GPU throughput `t_step * n / cells`:

| ranks | Mcell/GPU | base | rect | big | **block tax** |
|---|---|---|---|---|---|
| 1 | 138.41 | 4.933 | 5.188 | — | 1.052 |
| 2 | 69.20 | 5.048 | 5.518 | — | 1.093 |
| 4 | 34.60 | 5.383 | 5.774 | 5.721 | 1.073 |
| 8 | 17.30 | 5.927 | 7.795 | 7.090 | **1.315** |
| 16 | 8.65 | 6.893 | 10.370 | 8.931 | **1.504** |

Strong-scaling efficiency (vs each config's smallest rank count):

| config | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| `base_jacobi` (unblocked) | 100 % | 98 % | 92 % | 83 % | **72 %** |
| `rect_jacobi` (blocked) | 100 % | 94 % | 90 % | 67 % | **50 %** |
| `refined_big` (blocked + 2:1) | — | — | 100 % | 81 % | **64 %** |
| `refined_yp82` small | 100 % | 91 % | 81 % | 53 % | **36 %** |
| `refined_yp82` red-black | 100 % | 90 % | 79 % | 47 % | **31 %** |

## Where it goes: `mpi_wait`, and nothing else

Attribution of the block tax (rect − base), by phase:

| ranks | total loss | `exch/mpi_wait` | `proj/vel_exchange` | `proj/phi_exchange` | `exch/local_copy` |
|---|---|---|---|---|---|
| 4 | +13.5 ms (+7.3 %) | −1.9 ms | +4.9 | +2.9 | **+13.4 (99 %)** |
| 8 | +32.3 ms (+31.5 %) | **+22.8 (70 %)** | +16.6 | +10.2 | +9.3 |
| 16 | +30.1 ms (+50.4 %) | **+22.3 (74 %)** | +18.2 | +6.8 | +5.7 |

Two different regimes. Up to 4 ranks the block tax is small and is
**device-local copy** work. From 8 ranks it is **`mpi_wait` inside the
projection's two exchanges**, and it is large.

### It is not volume, not rounds, and not bandwidth

**Every configuration performs exactly 39 exchange rounds per step**, at every
rank count. Per-round `mpi_wait`:

| config | n=2 | n=4 | **n=8** | **n=16** |
|---|---|---|---|---|
| `base_jacobi` (unblocked) | 38.9 us | 121.0 | **124.8** | **143.4** |
| `rect_jacobi` (blocked) | 39.3 | 71.7 | **708.3** | **715.7** |
| `refined_big` | — | 75.4 | **781.2** | **804.1** |
| `refined_yp82` small | 41.2 | 55.3 | **449.7** | **437.7** |

Every blocked config jumps ~10x between 4 and 8 ranks — **exactly the node
boundary** (4 GPUs/node). The unblocked config does not move (121 -> 125 us).
And at 16 ranks the blocked case sends **4.3x FEWER points to 5.5x fewer
peers** than the unblocked one, yet waits **5x longer**:

| at 16 ranks | peers | total send pts | MB/round | `mpi_wait`/round | effective |
|---|---|---|---|---|---|
| `base_jacobi` | 11 | 4 735 872 | 151.5 | 143.4 us | **67.0 GB/s** |
| `rect_jacobi` | 2 | 1 104 000 | 35.3 | 715.7 us | **3.3 GB/s** |
| `refined_big` | 2 | 1 380 000 | 44.2 | 804.1 us | 3.7 GB/s |

Same hardware, same round count, comparable message sizes (base 656 kB/peer,
rect 258 kB/peer — the *bigger* messages are the fast ones). The blocked
exchange achieves **3 GB/s where the unblocked achieves 67 GB/s**. It is also
**flat from 8 to 16 ranks** while total volume doubles — the signature of a
latency/serialization limit, not a bandwidth one.

### The mechanism: the Z-order linear split is a 1-D chain

The send-size report gives it away. In every blocked run, at every rank count:

```
send pts/rank min = q,  max = 2q      (q = 36800 / 46000 / 23000)
peers/rank (max)  = 2
```

`min` is one peer and `max` is two, with an identical per-peer quantum. So the
closed-form linear split of the Morton curve gives each rank exactly its
**curve predecessor and successor** — a **1-D chain**, whose two end ranks have
a single peer. Peer count and per-peer volume do not change from 4 to 16 ranks;
only whether those two links cross a node does.

That explains every observation: a chain has no parallelism to hide latency (2
concurrent transfers against the unblocked decomposition's 11), and a single
slow link — the cross-node one — stalls its neighbours, which stall theirs,
39 times per step. It predicts the flatness from 8 to 16 ranks, since the cost
is set by the slowest link rather than by aggregate volume, and that is what is
measured.

**This is a hypothesis with strong circumstantial support, not a proof.** The
decisive experiment is cheap and is proposed below.

## Headline 3 — red-black does NOT reduce communication

Red-black at its calibrated `niter = 3` performs **the same 39 exchange rounds**
as Jacobi at `niter = 6` (3 iterations x 2 colours = 6 exchanges), with the same
per-round wait (449 vs 450 us at 8 ranks). Its advantage is entirely compute-side
and therefore **shrinks as the exchange takes over**:

| ranks | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| red-black / Jacobi | 0.763 | 0.763 | 0.783 | 0.852 | **0.892** |

24 % faster on one GPU, 11 % on sixteen. The earlier hope that red-black would
cut communication is **disproved**: it trades exchange rounds one-for-one.

## The prize

If the blocked exchange's per-round wait fell to the unblocked level — proven
achievable on this hardware, by `base_jacobi`, at 10x the byte volume:

| config | ranks | s/step now | if fixed | **gain** |
|---|---|---|---|---|
| `rect_jacobi` | 16 | 0.08971 | 0.06739 | **24.9 %** |
| `refined_big` | 16 | 0.13520 | 0.10943 | **19.1 %** |
| `refined_yp82` | 16 | 0.05615 | 0.04467 | **20.4 %** |
| `refined_yp82` red-black | 16 | 0.05010 | 0.03857 | **23.0 %** |
| `rect_jacobi` | 8 | 0.13487 | 0.11212 | 16.9 % |
| `refined_big` | 8 | 0.21466 | 0.18906 | 11.9 % |

**19–25 % of the step at 8–16 GPUs**, and it grows with rank count. For
comparison, the entire single-rank kernel-optimisation track of the last weeks
(diagonal hoist + reciprocal) was worth 9.5 %.

## What this changes

The optimisation track has been aimed at the wrong target. `jacobi_apply` and
`jacobi_compute_phi` are the right targets **for a single GPU**, and the
per-cell 2:1 cost — the thing this whole campaign was built to measure — turns
out to be **zero**. At production rank counts the cost is the **rank split**,
which is not a 2:1 problem at all: it hits the unrefined blocked case just as
hard (block tax 1.504 vs the refined case's 1.55 relative to unblocked).
