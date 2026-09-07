# Mechanism probe — why does the blocked exchange stall at the node boundary?

## The observation to explain

From `overheadTest/results_horeka_2026-09-07.md` (16x A100, commit `a11e355`).
Per-round `mpi_wait`, 39 rounds per step in **every** configuration:

| config | n=2 | n=4 | **n=8** | **n=16** |
|---|---|---|---|---|
| `base_jacobi` (unblocked) | 38.9 us | 121.0 | **124.8** | **143.4** |
| `rect_jacobi` (blocked) | 39.3 | 71.7 | **708.3** | **715.7** |

Every blocked config jumps ~10x between 4 and 8 ranks — exactly where the job
stops fitting on one node. The unblocked control does not move. At 16 ranks the
blocked case sends **4.3x fewer points to 5.5x fewer peers** yet waits **5x
longer**: 3.3 GB/s effective against 67 GB/s on the same hardware, with the
*bigger* messages being the fast ones. The wait is **flat from 8 to 16 ranks**
while volume doubles.

So it is not volume, not round count, and not bandwidth. Worth **19–25 % of the
step at 8–16 GPUs** if fixed.

## The four data points any explanation must fit

|  | local copy | crosses a node | wait/round |
|---|---|---|---|
| `base` n=4 | 1.46 M | no | 121 us |
| `base` n=8 | **0** | yes | 125 us |
| `rect` n=4 | 15.2 M | no | 72 us |
| `rect` n=8 | 14.9 M | **yes** | **708 us** |

Only the last cell is slow. Two hypotheses survive:

- **H-chain** — the linear Morton split gives every rank exactly its curve
  predecessor and successor (the size report shows `min = q`, `max = 2q` with
  an identical per-peer quantum at every rank count). Two concurrent transfers
  instead of the unblocked decomposition's eleven, and one slow cross-node link
  stalls its neighbours, 39 times per step. Predicts the 8->16 flatness.
- **H-contention** — the device-local copy kernel and the cross-node transfer
  contend. Intra-node peers go over NVLink/`cuda_ipc` and never touch the path
  the copy engine is saturating; cross-node peers do. This fits all four cells
  above, which H-chain alone does not (it says nothing about `base` n=8 being
  free).

## Pre-registered predictions

Written before the runs, so the verdict is read off the data.

### Tier 1 — the decisive pair: 2 ranks, 1 node vs 2 nodes

At **2 ranks** each rank has exactly **one** peer, the pair is symmetric
(36 800 points each), and **no stall can propagate** — there is no chain and no
third rank to be imbalanced against. Only the link's locality changes.

**The discriminator is `rect(2x2) / base(2x2)`** — blocked against unblocked on
the *same* single cross-node link. Each config's own 1-node-to-2-node jump is
*not* the test: the unblocked control is itself expected to rise when it starts
crossing (it goes 38.9 -> 121 us in the main matrix), so its jump carries no
information alone. Controlling against it removes the fabric from the question.

| `rect(2x2)/base(2x2)` | reading |
|---|---|
| **> 3x** | **H-transport.** A single cross-node link is intrinsically slow *in the blocked path* — two symmetric ranks, one peer each, no chain and no third rank to be imbalanced against — while the unblocked control over the same link is fine. Re-partitioning would **not** fix it; the cause is the blocked path's buffers or message structure. |
| **< 1.6x** | **H-chain.** One crossing is as cheap for the blocked path as for the unblocked one, so the 10x is *collective* — it needs several ranks and links. **The partitioning is the target.** |
| 1.6–3x | Intermediate; weight the link-count series and the `blk8` probe before choosing a fix. |

### Tier 2 — link count at fixed rank count

4 ranks on 1 / 2 / 4 nodes (0 / 1 / 3 crossings) and 8 ranks on 2 / 4.
Rank count, decomposition, peer count and volume are all held fixed.

- **Proportional** to the number of crossings -> per-link transport cost.
- **Saturating after the first crossing** -> one slow link stalls the chain
  (H-chain).

### Tier 3 — `blk8_jacobi`: node crossing WITHOUT the local-copy load

`nb = 1024 176 96` gives exactly 8 blocks (verified with `leaftable_test`), so
at 8 ranks there is **one block per rank** and no same-rank pairs: the copy
volume collapses while the 2-peer chain, the cell count and the placement stay
identical to `rect_jacobi`. This is the missing fourth cell of the table above.

- `blk8` 8x2 fast (~125 us) -> **H-contention confirmed.**
- `blk8` 8x2 slow (~700 us) -> local copy is innocent; H-chain or the message
  structure.

### Tier 4 — `nb16_jacobi`: volume at fixed placement

Changes the exchange volume several-fold without moving anything else.
Wait flat -> latency/serialization. Wait tracking bytes -> bandwidth (which the
main matrix already argues against, so this is a cross-check).

## Running it

Transfer the whole `horeka/` package (this directory is inside it), then on a
login node:

```bash
cd horeka/mechanism
WS=/hkfs/work/workspace/scratch/<user>-<ws> bash setup_mechanism.sh
```

It reuses the main campaign's clone and HDF5 if present — same commit
`a11e355`, deliberately: this probe exists to explain that matrix, so it must
not run different code. 23 runs of 200 steps inside one 4-node allocation,
roughly 25–35 min; 2 h requested.

Results: `$RUN_DIR/results_mechanism/mechanism.md`, written by
`collect_mechanism.py`, which prints the tables **and the verdict** against the
predictions above. Raw logs sit one directory per run alongside it.

## Reading it — check these before the verdict

1. **Placement is evidence, not assumption.** Every run captures
   `--display-map`; `hosts.txt` lists the nodes actually used. A run labelled
   `_N2` whose `hosts.txt` has one entry is void — OpenMPI ignored
   `--map-by ppr:N:node`. If that happens, override the launcher:
   `MPIRUN_EXTRA=... ` or replace `ppr` with an explicit `--host`/`--nodelist`.
2. **GPU binding.** With 2 ranks/node only 2 of the 4 GPUs are used — intended,
   since total rank count is what is held fixed. Confirm each rank still reports
   one target device.
3. **`nsteps = 200`** here, not 400. Every number is a ratio *within* this
   matrix; do not compare its absolute s/step against the main campaign's.
4. **Rank independence** still holds: `L2_div` must be identical across
   placements of the same config at the same rank count. If it is not,
   something is wrong with the run, not with the network.

## If the verdict is H-contention

The fix is scheduling, not partitioning: separate the device-local copy from
the cross-node transfer (issue local copies after the MPI posts, or on a
different stream), which is a bit-exact rescheduling change.

## If the verdict is H-chain

The fix is the rank split: replace the closed-form linear Morton split with a
partitioning giving each rank a compact 3-D neighbourhood — more peers, no
chain. `tools/partition_analysis.py` already scores candidates offline by
shared area and can be pointed at these leaf tables before any solver change.

Either way the prize is the 19–25 % of the step quantified in the main report.
