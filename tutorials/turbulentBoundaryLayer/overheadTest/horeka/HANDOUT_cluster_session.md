# Handout — run this session ON HoreKa

Everything needed is in the repo. No copy-back: run, analyse and write up on the
cluster, commit and push from there.

## Start here

```bash
cd $WS/moby-2to1-code            # the clone; $WS = your workspace
git fetch origin && git checkout optimiseBlockRefinement_parentBoundaryLayer
git pull --ff-only               # you need be48d44 or later (exchange diagnostics)
git log --oneline -1
```

`$WS` is `/hkfs/work/workspace/scratch/<user>-<ws>/...`. Do **not** work in
`$HOME`. If `$RUN_DIR` (`$WS/moby-2to1-run`) already exists, keep it: the
drivers resume, skipping completed runs.

**Read first**, in this order — they are short and each one closes a door you
would otherwise re-open:

1. `tutorials/turbulentBoundaryLayer/overheadTest/results_horeka_2026-09-07.md`
2. `.../results_horeka_mechanism_2026-09-08.md`
3. `.../results_horeka_mechanism2_2026-09-09.md`
4. `.../horeka/mechanism/README.md` (the pre-registered predictions)
5. `CLAUDE.md`

## Where the question stands

The 2:1 refinement machinery is **free** (added cells cost 0.979x a coarse cell
at 4 ranks, like-for-like). What is expensive is the exchange at scale, and
**three mechanisms for it have been proposed and killed by data**:

| hypothesis | killed by |
|---|---|
| per-link cross-node transport cost | `rect` 8 ranks: 3 crossings (182 us) is **4.1x faster** than 1 crossing (753 us) |
| a slow link stalling the 1-D Morton chain | Tier 1: blocked/unblocked on one cross-node link is **1.12x** |
| device-local copy contending with off-node transfer | `nb16` carries **3.7x** `rect`'s copy load for **1.22x** the wait |

What *is* established: **placement moves `mpi_wait` by up to 23x**, and its
optimum inverts with peer count (many-peer decompositions collapse when spread;
the 2-peer chain improves). No exchange number in this project is meaningful
without stating the rank-to-node layout.

The open question is not *where* the time goes but **what a rank is waiting
for**, and aggregate `mpi_wait` cannot answer it. Two diagnostics were added for
exactly this (commit `be48d44`, bit-exact, gated on the workstation):

- **`exchange balance:`** — one line at the end of every profiled run: `mpi_wait`
  min / mean / max across ranks and **which rank** holds the max.
  `min ~ max` = everyone waits equally (fabric/transport). Wide spread with a
  stable argmax = **everyone waits for one late rank** (skew/partitioning).
- **`[output] exchange_barrier = true`** — an `MPI_Barrier` before every
  `Waitall`. The barrier then absorbs the **arrival skew** into its own
  `skew_barrier` bucket, and what remains in `mpi_wait` is **transfer**.
  Serialising: use it for attribution only, never quote its step times.

## Phase 1 — close the copy question (~5 min of compute)

Three runs already in `SPECS`; resume skips the 23 that are done.

```bash
cd $RUN_DIR
bash mechanism/setup_mechanism.sh          # re-stages scripts+configs, submits
# or, if staging is already current:
sbatch --export=ALL,CODE_DIR=$CODE_DIR,RUN_DIR=$RUN_DIR,HDF5_ROOT=$HOME/hdf5 \
       mechanism/submit_mechanism.sh
```

- **`nb16_jacobi:8:2`** is the test: the anomaly's exact placement, the same
  2-peer chain, **7.6x** `rect`'s copy load per rank.
  Worse than 753 us ⇒ copy is back in. Near or below ⇒ copy is out for good.
- `nb16_jacobi:8:4` and `blk8_jacobi:4:2` complete the (copy × placement) 2x2.
- **Delete `results_mechanism/rect_jacobi_r8_N2/` first** so it re-runs in this
  allocation. Every comparison against it is currently cross-job (different
  nodes, two days apart, ~6 % node variation); this removes that caveat.

## Phase 2 — what is a rank waiting for (~10 min)

Two passes over the same small set, into **separate** results directories.

```bash
cd $RUN_DIR
S="rect_jacobi:8:2 rect_jacobi:8:4 base_jacobi:8:2 nb16_jacobi:8:2"

# (a) balance only -- the `exchange balance:` line comes free with profile
SPECS="$S" bash mechanism/run_placement.sh \
      $CODE_DIR/build_gpu/moby_solve $RUN_DIR/results_balance

# (b) skew/transfer split -- serialising, attribution only
BARRIER=1 SPECS="$S" bash mechanism/run_placement.sh \
      $CODE_DIR/build_gpu/moby_solve $RUN_DIR/results_barrier
```

Run these **inside one allocation** (an interactive `salloc -N4 ...
--gres=gpu:4 -t 1:00:00`, or add them to `submit_mechanism.sh`), so the two
passes and Phase 1 share one machine state.

### How to read it — decide this before looking

| `exchange balance` | `skew_barrier` vs `mpi_wait` (barrier pass) | conclusion |
|---|---|---|
| max/min ≈ 1 | wait ≫ skew | **transport**: the messages themselves are slow. Fix in the exchange (buffers, message structure, transport tuning). Partitioning will not help. |
| max/min ≫ 1, stable argmax | skew ≫ wait | **skew**: ranks arrive at wildly different times and the wait is idle. Fix upstream — load balance / partitioning. Find what the argmax rank has more of. |
| max/min ≫ 1 | wait ≫ skew | ranks differ in how much they transfer, not when they arrive → look at per-rank send volume (`exchange sizes` min/max). |

`rect` at 8 ranks / 4-per-node is the target: reproducibly ~708–753 us across
two jobs, two node sets and two mapping policies. `base` at the same placement
(101 us) is the control — run both and compare their balance lines.

## Rules that still bind

- **State the placement with every number.** It is worth up to 23x. `hosts.txt`
  and `--display-map` are captured per run; a run whose achieved node count
  differs from the request is marked `VOID` and excluded — check for those.
- **Same allocation for anything compared.** Machine state drifts; node-to-node
  variation is ~6 %.
- **Judge a run from `runtime.txt` drift**, not the final cumulative average
  (its last column is cumulative). All 23 runs so far sit at 0.992–0.997
  marginal/cumulative; anything outside that was disturbed.
- **`L2_div` must be identical** across placements of the same config and rank
  count. It has been in every run so far; if it is not, the run is broken.
- **Never quote a step time from a `BARRIER=1` run.**
- **Do not reopen** the halo copy kernel (at its ~51 %-of-peak coalescing floor)
  or redundant-computation schemes (blocked by the periodic-seam metric
  asymmetry).
- **Solver changes need the bit-exactness gate**, which lives on the workstation
  (`~/.moby_prof/gate_bitexact.sh`, reference tree in `~/.moby_ref/ref_src`) —
  it is NOT available on HoreKa. If this session changes solver code, say so
  explicitly and leave the gate to a workstation session. Diagnostics-only
  changes (timers, reductions, barriers) are bit-exact by construction but
  still get gated there before they are believed.

## Deliverable

A dated results file `overheadTest/results_horeka_<date>.md` in the repo
convention: the tables, the mechanism each number supports, and — explicitly —
what the numbers do **not** support. Then update the STATUS header of
`docs/next_session_2to1_penalty.md`. Commit the raw run directories too
(`git add -f`, they are under a `*.log` ignore) and push the branch.

If Phase 2 identifies the mechanism, **stop there and write it up**. Do not also
implement the fix in the same session: a partitioning or transport change needs
the workstation gate, and the analysis is the valuable half.
