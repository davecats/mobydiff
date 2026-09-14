# Handout — the between-iteration velocity exchange, and what is left after it

Written 2026-09-14, at `76f643f`, from the session that took **≈18 % off the
step** in five increments. **Read section 2 before planning anything: most of the
obvious targets are now closed by measurement, and one was closed by a
decision.**

## 0 — Start here

```bash
cd $WS/moby-2to1-code
git checkout optimiseBlockRefinement_parentBoundaryLayer
git pull --ff-only          # you need 76f643f or later
```

Read, in this order:

1. `tutorials/turbulentBoundaryLayer/overheadTest/results_stepwork_2026-09-14.md`
   — the last increment, and the two negative results that bound the next one
2. `results_apply_registers_2026-09-14.md` and
   `results_rdenom_registers_2026-09-14.md` — the register line, start to finish
3. `CLAUDE.md`, the `jacobi_apply` / register paragraphs
4. `horeka/exchange/PREREGISTERED_bodyforce.md` — the house style for writing the
   readings down first. It called the branch the trip change actually took.

## 1 — Where the step stands (8 ranks, measured, job 5145554)

| bucket | `rect_jacobi` 83.88 ms | `refined_yp82` 42.65 ms |
|---|---|---|
| `apply` | 31.83 ms **37.9 %** | 15.90 **37.3 %** |
| `momentum` | 19.27 **23.0 %** | 8.49 **19.9 %** |
| `sweep` | 12.47 **14.9 %** | 5.69 **13.3 %** |
| `vel_exchange` (in projection) | 7.41 **8.8 %** | 4.35 **10.2 %** |
| `phi_exchange` | 5.21 **6.2 %** | 3.81 **8.9 %** |
| `setup` (= `compute_rdenom`) | 3.67 **4.4 %** | 1.65 **3.9 %** |
| `apply_bc` (in projection) | 1.13 1.3 % | 0.84 2.0 % |
| `bodyforce` | 0.62 0.7 % | 0.45 1.1 % |
| `ibm_mu` | 0.02 — | 0.01 — |

Inside the exchange (13.21 / 8.50 ms per step total):
`local_copy` 7.47 / 3.18, `mpi_wait` 2.89 / 1.87, pack+unpack 2.61 / 2.31.

Register state, `cuobjdump -res-usage`, all `STACK:0 LOCAL:0` unless noted:

| kernel | regs | | kernel | regs |
|---|---|---|---|---|
| `jacobi_apply` k1 / k2 | 59 / **64** | | `step_momentum` predictor / copy | **128** / 66 |
| `jacobi_compute_phi` | **80** | | `add_eddy_viscosity_correction` | 136 |
| `compute_rdenom` | **80** | | `redblack_sweep` | 96 |
| `interface_correct` ×3 | 54 | | `airfoil_after_step` | 164 (`STACK:776`) |

## 2 — What is closed, and must not be reopened

| | why |
|---|---|
| **`jacobi_apply`** | Done. 64 registers, 46 % occupancy, 62 % of peak DRAM, **1.07x** its source-counted minimum. Fusing k1 into k2 pushes k2 back over the 64-register threshold — a measured bad trade. |
| **Eliding the `ibm%mu` multiply on body-free ranks** | **DECLINED BY THE USER, 2026-09-14, on readability grounds.** It is worth ~8 % of the step (mu is 3 of k2's 10.68 doubles/cell and `mu ≡ 1.0` exactly with no body, so skipping it is bit-exact), but it puts a uniform branch or a duplicated expression in the solver's hottest kernels. **Do not re-propose it.** |
| **The register lever on `step_momentum`** | Packing its nine Laplacian arrays into three moved the count **by nothing** (128 → 128). For that kernel the array bases are not the constraint; the ~30 live values of three fused component blocks are. |
| **`-gpu=maxregcount`** | It is a **target, not a cap**: globally it reached 102 on `step_momentum` only with `STACK:16`, and *raised* every other kernel toward the limit (`jacobi_apply` k2 **64 → 96**, `compute_rdenom` 80 → 100, `interface_correct` 54 → 80), undoing a day's work in one build. |
| **`local_copy`** | At its documented coalescing floor, and its volume is load-bearing: the redundant halo sweep is what makes results independent of `nb` and rank count. |
| **`mpi_wait`, overlap, partitioning** | Closed by the A0 probe and the 2026-09-10 §3 analysis. |
| **The exchange's per-launch cost** | Closed by `map(to: c)` (2026-09-11). 0–1 host-to-device copies per launch; per-launch cost at the floor. |
| **`ibm_mu`, the trip force** | Done, −99.5 % and −88.6 %. |

## 3 — Task 0 (do this first): re-run the campaign matrix

**Every ratio the campaign publishes now has a moved denominator.** ~18 % has
come off the step, unevenly across buckets, so the block tax, the strong-scaling
efficiencies and the 2:1 coarse-cell-equivalent in
`results_horeka_2026-09-14.md` are all stale. One 4-node job with the existing
`horeka/run_matrix.sh` + `submit.sh`; `accelerated` ran same-day on 2026-09-14
despite a three-day estimate, so submit early and do task 1 while it queues.

Do it **before** task 1, so the new baseline is published and task 1 is measured
against something current.

## 4 — Task 1: extend `sync_divergence_halos` to the multi-rank path

**This is the main task.** Between projection iterations only the divergence
stencil reads the velocity halo, so one plane per dimension of the normal
component suffices — not the full three-component 26-direction shell. The solver
**already does exactly this at one rank**; `pressure_projection` says why it
stops there:

> *Single-rank only: with peers the same planes would have to come over MPI, and
> the message is the whole point of the saving, so that needs the entry list
> partitioned first.*

**Size.** `proj vel_exchange` is 7.41 / 4.35 ms per step over 18 calls
(412 / 242 µs each) and **15 of the 18 are between-iteration**. A divergence-only
round moves roughly a tenth of the data but keeps the ~60 µs per-round launch
floor (pack + unpack + local copy at ~15–18 µs each after the `map(to: c)` fix),
so 412 → ~120 µs gives **≈4.5 / 2.5 ms per step, about 5 % of the step.**

**Why bit-exactness is nearly free here.** The single-rank path already uses the
reduced exchange, and *1 rank == 4 ranks is a standing gate that passes today* —
which is direct evidence that the reduced set delivers every halo value the next
divergence reads. The multi-rank version must deliver the same values, so the
existing rank-count gate is the right one to lean on. **Verify that claim against
the entry list before writing code**, do not assume it.

**Where the work is.** `comm.f90` — the most intricate module in the solver. The
entry list must be partitioned so a divergence round posts only the entries it
needs, in a canonical order both ends derive independently (the existing
invariant). Note the entries are already ordered *same-level-copies-first with
prefix counts* so the per-colour copy-only exchange is a prefix of the full one;
a second prefix, or a per-variable/direction mask, is the shape to look for.

**Watch for:** the launch count does not fall — pack, unpack and local copy are
one kernel each per round regardless of how many entries they cover — so the
saving is in bytes moved and in `mpi_wait`, not in launches. If the measured
saving is much larger than ~5 %, distrust it and find out what else changed.

## 5 — Task 2 (small, clean): `compute_rdenom`'s fp64 divide

`setup` is 4.4 / 3.9 % of the step and after the register cut the kernel is
**compute-limited**: SM 49.6 % against DRAM 39.8 %. The one arithmetic item left
is its per-cell `1.0d0/(...)`. The rdenom comment in `pressure_solver.f90` records
that removing exactly this divide from `jacobi_compute_phi` was worth **50.5 %**
of that kernel — but here the divide cannot simply move, because `rdenom` *is*
the reciprocal. **ncu roofline before any code**; the ceiling is the whole bucket.

## 6 — Task 3 (own session, own pre-registration): split `step_momentum`

23.0 / 19.9 % of the step, the largest item after `apply`, and genuinely
occupancy-limited: 128 registers, **23.96 % occupancy, 32.50 % of peak DRAM**
against its own sibling's 41.09 % / 74.22 % at 66 registers (job 5145549 — the
first ncu ever taken on it). Its traffic is **1.05x** its source-counted minimum,
so there is nothing to reclaim there either.

The only idea left is splitting the fused predictor into three per-component
kernels: far fewer live values each, at the cost of re-reading the shared
velocity components. **It is a real trade, not a free win** — the kernel is
already at 1.05x minimum traffic, so more traffic is the price of more occupancy,
and at 32.5 % of peak there is headroom but no guarantee. It is also the kernel
whose bit-exactness underpins every phase gate in this project.

**Do not start it inside another task.** It needs its own session, a
pre-registration that names the traffic increase it will accept, and the full
gate list.

## 7 — Rules and tooling (all of this was learned the hard way)

- **Which gate flags.** Production flags are the *tighter* test when the
  arithmetic source is byte-identical between the two binaries (a device index, a
  skipped kernel, a removed store). A change that **moves an expression** (a
  hoisted metric, a precomputed table, a reordered sum) needs
  `./compile.sh cpu_nofma` / `gpu_nofma` and
  `horeka/exchange/submit_nofma_gate.sh` on **both** sides — the compiler fuses
  `a*b+c` differently and the production gate fails at 1e-15..1e-12 with nothing
  wrong. Task 1 is a comm change: production flags.
- **The register loop is off-queue and free.** `cuobjdump -res-usage` on the
  linked binary, minutes per attempt. But `./compile.sh gpu` on a HoreKa login
  node poisons `build_gpu` (no GPU ⇒ multi-arch target ⇒ cannot link against the
  job's `cc80` objects); build into a separate dir with
  `-DOPENMP_OFFLOAD_FLAGS="-mp=gpu -gpu=cc80"`. See `horeka/README.md`.
- **Pin a queued job to a worktree**, never the live tree: a job that starts days
  later would otherwise build whatever the tree has become. Job 5145101 was
  cancelled and resubmitted as 5145120 for exactly this.
- **`dev_accelerated`** (2 nodes, 1 h) schedules same-day; the jobs in this
  campaign take 2–8 minutes.
- **Write the readings down first.** `PREREGISTERED_*.md`. It is not ceremony:
  on 2026-09-14 it named in advance the branch the trip change actually took,
  which is what pointed at the real fix.
- **Same allocation for anything compared.** Before and after inside one job.
- Job scripts to copy: `horeka/exchange/submit_bodyforce.sh` (build both,
  timeline A/B at 4 and 8 ranks, Pass G, 7-case suite) and `submit_trip.sh` (the
  same plus an ncu section).

**Housekeeping:** this session left six reference worktrees
(`moby-2to1-{headref,applynew,rdenomref,bfref,tripref,lapref}`). Keep one at
`76f643f` as the baseline for task 1 and `git worktree remove` the rest.

## 8 — Deliverable

A dated `overheadTest/results_*_2026-xx-xx.md` in the repo convention: the
before/after buckets at 4 and 8 ranks from one allocation, the gate outputs, the
graded pre-registration, and **explicitly what the numbers do not support**. Then
update `CLAUDE.md` and this file's STATUS header.

If task 1 turns out not to be worth ~5 %, **say so and stop there** — after it,
the remaining items are one small arithmetic question (task 2) and one large
uncertain refactor (task 3), and the honest position may be that this
optimisation track is finished.
