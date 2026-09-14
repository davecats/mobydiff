# Handout — `jacobi_apply` k2: can 88 registers become 64?

> **STATUS: CLOSED, 2026-09-14. Yes — by L2 alone, and the time followed.**
>
> Hoisting `face_grad_corr` out of the `collapse(4)` body (module arrays
> `cfLow`/`cfHigh`, `src/modules/pressure_solver.f90`) takes k2 from **88 to 64
> registers**, `STACK`/`LOCAL` still 0. ncu on the same node, both binaries:
> occupancy **29.6 → 46.0 %**, DRAM **42.5 → 60.5 % of peak**, k2 **7 565 →
> 5 309 us (−29.8 %)** with the traffic unchanged to the last digit printed
> (10.31 doubles/cell) — the mechanism §6 pre-registered, confirmed. `apply`
> falls **−21.8 % / −20.1 %** and the step **−8.6 % / −7.8 %** at 4 and 8 ranks,
> with `sweep` as an unmoved control. Every gate at `max_abs 0`.
>
> **L1 was not taken** and should not be: it costs a launch and the threshold it
> was for is already crossed. **L3/L4 are moot.**
>
> Full report, with what the numbers do not support and the next lever
> (`compute_rdenom`, 110 registers, the same hoist verbatim):
> `tutorials/turbulentBoundaryLayer/overheadTest/results_apply_registers_2026-09-14.md`.
> The 16-rank A/B (job 5145101) is still queued; §6 there is a prediction, not a
> measurement.
>
> Everything below is the handout as written, kept for the reasoning it records.


Written 2026-09-14 from the session that removed the exchange's per-launch
marshalling. **Read section 1 before planning anything: the obvious explanation
for this kernel's cost has already been measured and refuted.**

## Start here

```bash
cd $WS/moby-2to1-code
git fetch origin && git checkout optimiseBlockRefinement_parentBoundaryLayer
git pull --ff-only          # you need a438646 or later
```

Read, in this order:

1. `tutorials/turbulentBoundaryLayer/overheadTest/results_ncu_apply_2026-09-14.md`
   — the ncu measurement and what it refutes
2. `tutorials/turbulentBoundaryLayer/overheadTest/results_horeka_2026-09-14.md`
   — the post-fix phase table this sits in
3. `CLAUDE.md`

## 1 — What is closed, and must not be reopened

| | why |
|---|---|
| "`apply` wastes memory traffic" | **Refuted by ncu.** k2 moves 10.82 doubles/cell against a source-counted minimum of 8 (1.35x); `compute_phi` 1.10x, k1 1.07x. There is nothing to reclaim. |
| "`apply` has a coalescing defect" | **Refuted.** Load sectors/request 1.9–2.5, store 9.0 — the *same* as the control kernel. |
| The exchange | Done (`results_kernel_timeline_2026-09-11.md`). 0–1 host-to-device copies per launch, per-launch cost at the floor. Do not reopen it. |
| Shared memory or block size as the occupancy limit | 0 bytes shared, 128 threads/block, zero local memory. Registers are the only binding limit. |
| Any scaling table older than `results_horeka_2026-09-14.md` | Superseded by the `map(to: c)` fix. |

## 2 — The measurement this rests on

`jacobi_apply` is **33.1 %** of the refined 16-rank step (34.4 % single-level),
the largest single item, and it did not move when the exchange fix landed.
Job 5144931, A100, `refined_yp82` at 1 rank:

| kernel | us | **DRAM %peak** | **occupancy** | **registers** |
|---|---|---|---|---|
| `apply` k1 `p += phi·idt` (`F1L470`) | 1 543 | **64.8** | **45.3 %** | **59** |
| `compute_phi` (control) | 3 204 | 53.3 | 29.1 % | 94 |
| `apply` k2 face correction (`F1L502`) | **7 654** | **44.0** | **29.6 %** | **88** |

k1 and k2 share the grid, the block size, and most of their arrays. k1 has 29
fewer registers and gets 1.5x the occupancy and 1.5x the DRAM utilisation. **That
is the entire hypothesis: k2 is latency-bound because too few warps are
resident.**

### The occupancy arithmetic, which decides what is worth trying

A100: 65 536 registers/SM, 128 threads/block, so blocks/SM = ⌊65 536 / (128·R)⌋:

| registers | blocks/SM | threads | occupancy |
|---|---|---|---|
| 88 (now) | 5 | 640 | 31 % |
| 80 | 6 | 768 | 38 % |
| **72** | **7** | **896** | **44 %** |
| **64** | **8** | **1 024** | **50 %** |

**Going 88 → 80 buys one block; the thresholds that matter are 72 and 64.** A
change that removes four registers is worth nothing. Size any lever against this
table before implementing it.

## 3 — The fast loop: no job, no queue

Register counts are **not** emitted per object file (the kernels are compiled at
link time), but they are in the linked binary:

```bash
./compile.sh gpu
cuobjdump -res-usage build_gpu/moby_solve | grep -A1 -E 'jacobi_apply|compute_phi'
```

Verified against ncu: `REG:88` / `REG:59` / `REG:94`, exactly the profiler's
numbers. **Iterate here, minutes per attempt, and only submit a job once the
register count has actually crossed 72 or 64.**

`STACK` and `LOCAL` must stay **0** in that output. A register cut that spills is
a loss that looks like progress.

## 4 — Where the registers plausibly are

`jacobi_apply`'s second kernel (`src/modules/pressure_solver.f90`, the
`F1L502` region) is one `collapse(4)` loop containing **six** correction sites:

- three LOW faces, each `face_grad_corr(physLow(d,b), idx==1, d1?(idx,VAR_?,b), outLow(d), refd(d))`
  followed by an `is_interface` predicate and a masked update reading `phi` at two
  points and `ibm%mu`;
- three HIGH faces, each guarded by `idx == n .and. outHigh(d) .and. physHigh(d,b) == FACE_PHYS`,
  each with its own inlined `face_grad_corr` and its own `phi`/`mu` reads.

Ten array bases are live across the whole body: `phi`, `blk%q`, `ibm%mu`,
`blk%d1x/d1y/d1z`, `blk%physLow`, `blk%physHigh`, plus the mapped `outLow`,
`outHigh`, `refd`.

**The structural observation.** The three high-face branches execute on the
`i == nx`, `j == ny`, `k == nz` planes only — a ~1/nb fraction of threads, with
`nb = 64 44 48` — yet their operands are live for every thread in the kernel.

## 5 — The levers, in order

**L1 — split the high-face corrections into their own kernel.** This mirrors what
`interface_correct` already does (three planes per block, its own pass, for a
stated reason). It removes `physHigh`, `outHigh` and the `ip/jp/kp` reads from
k2's live set. Costs one extra launch: 18/step × ~15 us = **0.27 ms/step**, which
the register win must beat.
*Safety*: the low-face path writes `q(i)` for `i = 1..nx` and the high-face path
writes `q(nx+1)` — disjoint locations, no thread writes both. **Verify that
before splitting**, do not assume it.

**L2 — hoist `cf` out of the inner loops.** `face_grad_corr(physLow(1,b), i==1,
d1x(i,VAR_U,b), outLow(1), refd(1))` depends on `(b, i)` alone — not on `j` or
`k` — yet `collapse(4)` recomputes it in every thread. The same holds per
dimension. Precomputing three small per-`(b, index)` arrays removes the call and
its arguments from the loop body. Bit-exact by construction: identical
expression, identical inputs, computed once instead of nb² times.

**L3 — `-gpu=maxregcount`.** Global, and `step_momentum` is already at 128
registers; capping it at 64 would spill the largest kernel in the solver to buy
occupancy in a smaller one. Last resort, and only per-file if nvfortran allows it.

**L4 — fuse k1 into k2** (the separate ~3 % lever: one launch and one duplicated
`phi` read, ≈1.0 ms/step). **It pulls the opposite way on registers.** Do it
*after* L1/L2 and re-measure — if fusion pushes k2 back over a threshold in the
table above, it is a net loss.

## 6 — Pre-registered readings

Write the numbers down before running, as this campaign now does by habit.

| result | conclusion |
|---|---|
| registers cross 72, occupancy → ~44 %, DRAM → ~55–60 %, k2 time falls 20–30 % | the hypothesis holds; take L2 as well and re-measure at 4/8/16 ranks |
| registers fall but **time does not** | k2 was not occupancy-limited; the whole line is wrong and this handout should be closed, not extended |
| registers **do not fall** | the compiler is keeping them live for a reason the source reading missed — look at `face_grad_corr` inlining and at whether `outLow/outHigh/refd` are being kept in registers rather than constant memory |
| `STACK`/`LOCAL` become non-zero | spilling; revert regardless of what the clock says |

Expected size if L1 reaches 72 registers and DRAM utilisation follows k1's
example: k2 7 654 → ~5 600 us, i.e. **~2.5 ms/step at 16 ranks, 7–8 % of the
step.** That is the claim to falsify.

## 7 — Gates and rules

- **Bit-exact, both levers.** `max_abs 0` on the production Pass G (138 M and
  60 M points) and the 7-case suite, CPU and GPU.
  `horeka/exchange/submit_ab.sh` already runs exactly that shape against a
  worktree at the pre-change commit — point `REF_DIR` at it and reuse it.
- **Same allocation for anything compared**, and quote step times only from
  `horeka/exchange/submit_ab16.sh` (pass A, 8 and 16 ranks, before and after).
- **Check registers and `sm__warps_active`, not only the clock** — via
  `horeka/exchange/run_ncu.sh`, which needs `--bind-to none` (SLURM's 1-task
  cpuset defeats OpenMPI core binding and fails *before* the application starts,
  which looks exactly like a profiling-permissions denial and is not).
- `ncu --page raw` is a **wide** table; the solver's stdout lands in `ncu.csv`,
  not `ncu.log`. `collect_ncu.py` handles both — do not rewrite it from memory.
- `dev_accelerated` (2 nodes, 1 h) schedules same-day when `accelerated` is days
  out; chain jobs with `--dependency=afterany` because they share `build_gpu/`.

## 8 — Deliverable

A dated `overheadTest/results_*_2026-xx-xx.md` in the repo convention: the
register counts before and after from `cuobjdump`, the ncu occupancy/DRAM pair,
the step-time A/B at 8 and 16 ranks, the gate outputs, and **explicitly what the
numbers do not support.** Update `CLAUDE.md` and this file's STATUS header.

If L1 and L2 together do not cross 72 registers, **say so and stop** — the next
target is then `interface_correct`'s three kernels per call (54 launches/step,
already at the per-launch floor, so a kernel-count question worth ~0.5 ms/step),
or the 39-round exchange cadence, which is a numerics change and needs a
convergence argument rather than a bit-exactness gate.
