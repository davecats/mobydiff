# 2:1 refinement — timing and profiling on HoreKa

A self-contained cluster package that measures what the 2:1 block-refinement
machinery actually costs, on real HPC hardware, across a rank sweep that crosses
the node boundary.

Everything measured so far was on 1–2 GPUs of a shared workstation
(RTX 3060 / A6000). That is enough to rank kernels against each other and it
produced the current optimisation track, but it cannot answer the questions that
decide whether the 2:1 machinery is production-ready:

- **What is the 2:1 overhead, like for like?** The small refined case has no
  matched single-level twin — its level-0 grid is half the production one — so
  every quoted figure so far is a per-cell comparison across two different
  grids. This package adds `refined_big_rect_jacobi`, which shares level-0 grid,
  block shape, refine box and solver with `rect_jacobi` **exactly**. Their
  per-cell ratio is the refinement machinery and nothing else.
- **What happens when the interface exchange leaves the node?** Every
  measurement to date was intra-node. The projection exchange is already the
  worst-scaling phase (71–75 % of all time lost at 2–4 ranks) and `mpi_wait`
  alone was 6–8 % of the step. The 4 → 8 rank step here is the first time it
  crosses an interconnect.
- **Does the Morton partitioning still balance at 8 and 16 ranks?** The xz-key
  fix (`4f04f6b`) cut peer traffic 20.5× and the per-rank send spread from 9.5×
  to 2.0× at 2–4 ranks. Whether the linear Z-order split stays balanced when
  ranks get thin is unmeasured — and *leaf-count balance is not exchange
  balance*: `MPI_Waitall` pays for the slowest peer.

## What it runs

Five configs, all cold-starting from their `.ini` (the `boundaryLayer` case is
analytic — **no restart field, no `moby_prepare` step**, so this whole package is
a few hundred kB), each profiled (`[output] profile = true`).

| config | grid (level 0) | blocks | leaves | cells | role |
|---|---|---|---|---|---|
| `base_jacobi` | 4096×176×192 | unblocked | — | 138.41 M | the unblocked rate; denominator of the block tax |
| `rect_jacobi` | 4096×176×192 | 64 44 48 | 1024 | 138.41 M | single-level baseline **and** the exact twin of the big refined case |
| `refined_big_rect_jacobi` | 4096×176×192 | 64 44 48 | 1792 | 242.22 M | **the headline case**: y-tile 1 refined in xz |
| `refined_yp82_rect_jacobi` | 2048×176×96 | 64 44 48 | 448 | 60.56 M | the existing small refined case; carries the local numbers over |
| `refined_yp82_rect_redblack` | 2048×176×96 | 64 44 48 | 448 | 60.56 M | red-black niter 3 vs Jacobi niter 6, at scale |

Rank sweep: **1, 2, 4** (intra-node, NVLink P2P) and **8, 16** (2 and 4 nodes).
`refined_big_rect_jacobi` holds ~39 GB of field state, so it starts at 4 ranks
(~9.7 GB/GPU on A100-40) and **will OOM at 1 rank** — that is why its rank list
differs.

The leaf counts above are **gates**, not documentation: the solver prints them
at init and `leaftable_test` was used to verify 448 and 1792 before shipping.
A different leaf or cell count means the refinement placement moved and the
comparison is void.

## Prerequisite: push the branch

`setup_and_run.sh` pins the commit so the cluster runs exactly the validated
code. As of writing, **the commit under test is not on the remote**:

```
COMMIT = a11e355a47e1535db4f3e9bae0dcc489eaec3567
         "Store the reciprocal diagonal and multiply. -4.3% of the step"
branch   optimiseBlockRefinement_parentBoundaryLayer   (313 commits ahead of origin/main, unpushed)
```

From the workstation, before running anything on the cluster:

```bash
git push -u origin optimiseBlockRefinement_parentBoundaryLayer
```

If pushing is not wanted, `rsync` the tree to `$CODE_DIR` instead and re-run —
the script accepts an existing `$CODE_DIR/.git` at the right commit and only
clones when the directory is absent. It stops with the exact fix rather than
silently testing a different commit.

## Running it

Transfer the package to the **workspace** (not `$HOME` — MPI-IO on the home
filesystem is slow and quota-limited), then run one command on a login node:

```bash
rsync -avP tutorials/turbulentBoundaryLayer/overheadTest/horeka/ \
      <user>@horeka.scc.kit.edu:/hkfs/work/workspace/scratch/<user>-<ws>/moby-2to1/

# on a HoreKa login node, from that directory:
WS=/hkfs/work/workspace/scratch/<user>-<ws> bash setup_and_run.sh
```

That builds parallel HDF5 (once), clones the code at the pinned commit into
`$WS/moby-2to1-code`, stages the run data and scripts into `$WS/moby-2to1-run`,
and submits the job. Monitor with `squeue --me` and
`tail -f $WS/moby-2to1-run/slurm-*.out`.

Edit before first use if needed: `--account` / `--partition` / `--mail-user` in
`submit.sh`, and `WS` (or `CODE_DIR`/`RUN_DIR`) above.

**Cost and wall clock.** 4 nodes × 6 h requested. 23 runs of 400 steps; at an
expected few tenths of a second per step the matrix is well under an hour, plus
~10 min HDF5 and ~10 min solver build on the first submission. The allocation is
sized for the 16-rank runs and the sweep deliberately runs the small rank counts
inside the same allocation, so the whole sweep sees one fixed machine state —
that is the only way the rank sweep is internally comparable.

**Resumable.** `run_matrix.sh` skips any run whose `run.log` exists. If the job
hits the wall clock, `sbatch submit.sh` again from `$RUN_DIR`. A run that failed
is renamed `run.FAILED.log` so it is retried rather than mistaken for complete.

## Output

`$RUN_DIR/results/` gets one directory per run — `config.ini` (exactly what ran),
`runtime.txt`, `run.log` — plus `provenance.txt` (commit, nodes, compiler, GPU)
and `summary.md` from `collect_profile.py`.

---

# Instructions for the Claude session on the HPC system

Start the session **on the cluster**, in `$RUN_DIR`, after the job has finished.
Everything below assumes the repo clone is at `$CODE_DIR`.

## First: establish that the numbers are trustworthy

Do this before reading a single timing. In order:

1. **Provenance.** `cat results/provenance.txt`. The commit must be
   `a11e355…`, `dirty` must be 0. If the clone is dirty, the results do not
   correspond to a commit — stop and say so.
2. **No failed runs.** `ls results/*/run.FAILED.log`. Any hit means the matrix is
   incomplete; report which, and why (read the log), before analysing the rest.
3. **The leaf/cell gates.** For each refined run,
   `grep "block refinement" results/<run>/run.log` must give **448** leaves for
   the small case and **1792** for the big one. A different count means the
   refinement placement moved; the comparison is void.
4. **Rank binding.** `grep -i "Cartesian dimensions" results/*/run.log` and check
   each run got the ranks its directory name claims. On a multi-node run confirm
   the GPUs were actually used — a silent fallback to one GPU per node reads as
   catastrophic scaling, not as an error. If in doubt, re-run one case with
   bindings printed:

   ```bash
   MPIRUN_EXTRA=--report-bindings CONFIGS=rect_jacobi RANKS_SMALL=8 \
       bash run_matrix.sh "$CODE_DIR/build_gpu/moby_solve" "$RUN_DIR/results_bind"
   ```
5. **Numerical sanity.** `tail -3 results/<run>/runtime.txt`: `L2_div` should be
   ~1e-5 and `Linf` bounded. A diverged run's timings are meaningless.
6. **Drift.** `runtime.txt`'s last column is a **cumulative** average. Judge each
   run from the slope over its last ~100 steps, not the final value; a rising or
   rise-then-fall average means contention or a transient, not a steady rate.
   This matters even on a batch node — the first steps carry allocation and
   first-touch costs.

## Then: the analysis, in this order

`python3 collect_profile.py results > summary.md` regenerates the report; read it
rather than re-deriving it, but verify its ratios against the raw logs for at
least one case.

1. **The 2:1 overhead, like for like.** `refined_big_rect_jacobi` vs
   `rect_jacobi` at equal rank count, **per cell**. This is the headline number
   the campaign exists to produce, and it has never been measured cleanly.
   Report it per rank count — if it grows with ranks, the overhead is in the
   exchange, not the kernels, and that is the finding.
2. **Where the 2:1 cost sits.** Difference the per-phase tables of those two runs
   at the same rank count. Expect the interface work in `proj/sweep`,
   `proj/apply` and the exchange buckets. Attribute the whole delta; an
   unexplained remainder means a phase is missing from the profiler's coverage
   (`step_timing: coverage … fraction` should be ~1.0).
3. **The node boundary.** Compare 4 ranks (intra-node) with 8 (2 nodes) for every
   config. The interesting quantity is not the speedup but **which bucket**
   absorbs the loss: `exch/mpi_wait` rising means the interconnect or a
   straggler; `exch/pack`/`unpack` rising means volume. Both were already known
   to grow *absolutely* with rank count on one node.
4. **Partition balance.** From the exchange-traffic table: `send pts/rank
   min..max`. A spread much above ~2× at 8 or 16 ranks reproduces, at higher rank
   counts, the defect the xz-key change fixed at 2–4 — and would be the next
   concrete optimisation target. Leaf-count balance is *not* exchange balance.
5. **Smoother, at scale.** `refined_yp82_rect_redblack` (niter 3) vs
   `refined_yp82_rect_jacobi` (niter 6). Red-black halves the exchange rounds, so
   its advantage should *grow* with rank count. Note the caveat the local work
   established: the two solvers need different `niter` to reach comparable
   quality, so s/step at fixed niter is not the comparison that matters — these
   configs already encode the calibrated pair.
6. **Block tax.** `rect_jacobi` vs `base_jacobi` per cell — the cost of blocking
   without refinement, which separates the block tax from the 2:1 tax.

## Writing it up

Follow the convention already in `overheadTest/`: a dated results file,
`results_horeka_<YYYY-MM-DD>.md`, in that directory, with the tables, the
mechanism for each number, and — explicitly — what the numbers do **not**
support. Then update the STATUS header of
`docs/next_session_2to1_penalty.md` and, if the ranking of targets changed,
`docs/next_session_jacobi_apply.md`.

## Rules that still apply on the cluster

- **Same-day, same-allocation A/B or nothing.** Machine state drifts; a number
  from another allocation is not a baseline. The sweep is inside one job for
  exactly this reason.
- **Never rebuild a binary while runs are in flight.**
- **Report buckets, not step times, if anything shared the node.** Check the job
  had exclusive nodes (it does, by SLURM allocation) before quoting absolute
  s/step.
- **An `ncu` ratio transfers to production only when the change moves bytes.**
  `ncu` locks the SM clock, so it overstates a pure compute reduction — this cost
  a 2× wrong projection once already.
- **Do not reopen** the halo copy kernel (at its ~51 %-of-peak coalescing floor)
  or redundant-computation schemes (blocked by the periodic-seam metric
  asymmetry).

## If you want to change the code and re-measure

The clone at `$CODE_DIR` is detached at the pinned commit. Work on a branch
there, rebuild with `./compile.sh gpu`, and run the matrix with a distinct
results directory:

```bash
bash run_matrix.sh "$CODE_DIR/build_gpu/moby_solve" "$RUN_DIR/results_<tag>"
```

Any change to the solver must still pass the bit-exactness gate before its
timing means anything — see `Verification` in `CLAUDE.md`. The gate harness lives
on the workstation (`~/.moby_prof/gate_bitexact.sh`), not here; a cluster-side
change should be gated back on the workstation before it is believed.
