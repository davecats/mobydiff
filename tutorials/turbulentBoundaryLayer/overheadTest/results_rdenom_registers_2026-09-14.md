# The same lever on the other two projection kernels — and a gate that had to change

Jobs 5145507 (registers + ncu + 4/8-rank A/B + the production-flag gate) and
5145518 (the nofma gate). A100-SXM4-40GB, HoreKa. `ref` = `6708193` (the apply
register cut), `new` = `294cbb7`. Raw runs in `horeka/results_job5145507/` and
`results_job5145518/`. Readings pre-registered in
`horeka/exchange/PREREGISTERED_rdenom.md`, including — before the numbers landed
— why this was expected to be a weaker case than `jacobi_apply` k2.

Follows `results_apply_registers_2026-09-14.md`, which established the mechanism
by measurement. This applies it where the same argument holds:
`face_grad_denom` and the divergence metric `d1?(idx,VAR_P,b)` are functions of
the block and the face-normal index alone, exactly as `face_grad_corr` was, so
they join the static tables as `dnLow`/`dnHigh`/`d1P`. `jacobi_compute_phi` gains
nothing from the denominators but loses `d1x`/`d1y`/`d1z` to the single `d1P`.

## 1 — Registers (job 5145507, `cuobjdump` on the A100 build)

| kernel | ref | new | blocks/SM |
|---|---|---|---|
| `compute_rdenom` | **110** | **80** | 4 → 6 |
| `jacobi_compute_phi` | **94** | **80** | 5 → 6 |
| `jacobi_apply` k2 | 64 | 64 | 8 |
| `jacobi_apply` k1 | 59 | 59 | 8 |
| `redblack_sweep` | 96 | 96 | — |
| `step_momentum` | 128 | 128 | — |

`STACK:0 LOCAL:0` on all of them, both sides.

## 2 — ncu: and `compute_rdenom` turns out to have been the worst kernel in the projection

Both binaries, one node, one rank on `refined_yp82_rect_jacobi` (60.56 Mcell),
one full substage profiled (1 `compute_rdenom` + 6 × (`compute_phi`, apply k1,
apply k2)) — the probe now covers `compute_rdenom`, **which had never been
profiled at all**.

| kernel | us ref | **us new** | DRAM %peak | SM %peak | occupancy % | doubles/cell | vs min |
|---|---|---|---|---|---|---|---|
| `compute_rdenom` | 4 608 | **3 342 (−27.5 %)** | **28.9 → 39.8** | 41.7 → 49.6 | **23.7 → 35.0** | 4.27 | 1.07x |
| `jacobi_compute_phi` | 3 207 | **2 754 (−14.1 %)** | **53.8 → 62.3** | 47.1 → 52.4 | **29.1 → 35.0** | 5.51 | 1.10x |
| `jacobi_apply` k1 (control) | 1 541 | 1 542 (+0.1 %) | 64.9 → 64.8 | 66.2 | 45.5 | 3.21 | 1.07x |
| `jacobi_apply` k2 (control) | 5 352 | 5 355 (+0.1 %) | 62.1 → 62.1 | 40.6 | 46.1 | 10.68 | 1.33x |

Two things worth stating plainly.

**`compute_rdenom` was running at 28.9 % of peak DRAM** — worse than k2's 42.5 %
ever was, the lowest number this campaign has measured, on a kernel whose traffic
is 1.07x its source-counted minimum. It was called 3 times a step instead of 18,
so it never showed up in a bucket ranking, and nobody had pointed a profiler at
it. The register table found it, not the phase table.

**The two apply kernels are unmoved to 0.1 %** in the same run, on the same node,
with the same profiler settings. They are the control, and they did not move.

## 3 — Step times (job 5145507, before/after in one allocation)

| case | step | `sweep` | `setup` | `apply` (control) |
|---|---|---|---|---|
| `rect_jacobi` 4×1 | 177.57 → 172.87 ms **−2.65 %** | 28.26 → 24.71 **−12.6 %** | 9.04 → 7.28 **−19.5 %** | +0.2 % |
| `rect_jacobi` 8×2 | 93.74 → 91.28 **−2.63 %** | 14.21 → 12.50 **−12.0 %** | 4.53 → 3.68 **−18.9 %** | +0.5 % |
| `refined_yp82` 4×1 | 83.87 → 81.74 **−2.54 %** | 12.50 → 10.98 **−12.2 %** | 3.99 → 3.23 **−19.1 %** | +0.2 % |
| `refined_yp82` 8×2 | 47.66 → 46.63 **−2.15 %** | 6.48 → 5.69 **−12.1 %** | 2.05 → 1.65 **−19.6 %** | +0.1 % |

Graded against the pre-registration: `sweep` −5 to −12 % predicted, **−12.0 to
−12.6 %** measured (at the top of the band, marginally past it); `setup` −10 to
−25 %, **−18.9 to −19.6 %** (inside); step −1 to −2.5 %, **−2.15 to −2.65 %** (at
the top, marginally past it). That is the first row of the pre-registered table:
*the lever generalises at a smaller size; take it and stop.*

`apply` is the control here and moves +0.1 to +0.5 % — a hair the wrong way,
consistent across four runs and small enough that it is reported rather than
explained. The `before` column reproduces job 5145114's `after` column (the same
commit, a different allocation) to 0.3 %: a fourth allocation agreeing.

**16 ranks is not measured for this increment.** Scaling the per-bucket
reductions onto the measured post-apply 16-rank step gives ~1.26 ms of 55.02
(`rect`, −2.3 %) and ~0.59 of 30.78 (`refined`, −1.9 %) — arithmetic, flagged as
such.

## 4 — The gate had to change, and that is the other result here

**The production-flag gate FAILED**, at 1e-15 on velocities and 1e-12 on pressure
(job 5145507, every case except `lam30t`). It is not a defect and it is not
noise: it is FMA contraction. In the old code the compiler inlined
`face_grad_denom`, saw the literal `2.0d0*d1f`, and could fuse the following
multiply into it; with the value arriving as an opaque load from `dnLow` it fuses
differently. `CLAUDE.md`'s Verification section says exactly this and requires
`-Mnofma` / `-Mnofma -gpu=nofma` on **both** sides for a refactor comparison.

What made it a trap is that the cluster gate scripts were written for
`select_target_device`, a change whose arithmetic source was byte-identical
between the two binaries — there, production flags are a *strictly tighter* test
than nofma, and `run_mapgate.sh`'s header said so in capitals. The apply hoist
then passed at production flags too, by luck rather than by argument. This one
does not, and could not.

There was no nofma build on this machine at all: `compile.sh` only ever built
`build_cpu` and `build_gpu`, while `validation/` scripts have long expected
`build_cpu_nofma` / `build_gpu_nofma`. So `compile.sh` gained `cpu_nofma` and
`gpu_nofma` modes (extra flags empty for the production modes, which therefore
build byte-identically to before) and `submit_nofma_gate.sh` runs the gate with
them.

| gate, **both sides `-Mnofma -gpu=nofma`** (job 5145518) | result |
|---|---|
| Pass G, `rect_jacobi` (138 412 032 points) | **max_abs 0** |
| Pass G, `refined_yp82_rect_jacobi` (60 555 264 points) | **max_abs 0** |
| 7-case suite, 200 steps, GPU | **max_abs 0**, all 9 comparisons, incl. `nut`, `k`, `omega`, `gamma`, `rethetat` |
| CPU `build_cpu_nofma`, `min_channel` 200 steps | **max_abs 0** |

(The CPU **production** build was also max_abs 0 at 200 steps — the host compiler
happened to contract identically. That is luck, not evidence, and is recorded
only so the GPU-only failure is not mistaken for a GPU bug.)

## 5 — What these numbers do NOT support

- **Not that the kernels got cheaper.** Traffic is unchanged: 4.27 and 5.5
  doubles/cell, 1.07x and 1.10x their minima before and after. Only the rate
  moved.
- **Not a 16-rank measurement.** §3's last paragraph is arithmetic.
- **Not that `compute_rdenom` is now fixed.** At 80 registers it reaches the same
  35.0 % occupancy as `compute_phi` but only **39.8 % of peak DRAM against
  62.3 %**, and its SM throughput (49.6 %) now *exceeds* its DRAM throughput.
  It is no longer occupancy-limited; it is compute-limited — see §6.
- **Nothing about `redblack_sweep`** (96 registers, unchanged). The same hoist
  would apply, but its `atBnd` predicate is `idx == hi(d)`, the redundant-halo
  sweep window, not `idx == nb(d)`; that needs its own table, and it is not the
  production smoother.
- **Nothing about the exchange or `mpi_wait`.** Untouched, and their *share* of
  the step has grown again because the step is smaller.

## 6 — What is next

1. **`compute_rdenom`'s remaining cost is very likely the fp64 divide.** It is
   the one kernel in the projection that still does `1.0d0/(...)` per cell, and
   after the register cut its SM throughput (49.6 %) exceeds its DRAM throughput
   (39.8 %) — the signature of a compute limit. The rdenom comment in
   `pressure_solver.f90` records that removing exactly this divide from
   `jacobi_compute_phi` was worth **50.5 %** of that kernel. Here the divide
   cannot simply move — rdenom *is* the reciprocal — but 3 calls/step × ~3.3 ms
   is ~1.6 / 0.75 ms/step at 16 ranks and it is now the least efficient kernel in
   the projection. Worth an ncu roofline before any code.
2. `jacobi_compute_phi` and `compute_rdenom` are both at **80 registers**; 72 and
   64 remain unreached. Whether the next 8 exist in the source is unknown — the
   cheap structural bases are gone.
3. `interface_correct`'s three kernels per call (54 launches/step, at the
   per-launch floor) — a kernel-count question worth ~0.5 ms/step.
4. `mpi_wait`, ~10–12 % of the step and untouchable by partitioning or overlap.

## 7 — The register line, end to end

Across the two increments, on the 16-rank refined case (the apply part measured,
the second part scaled):

| kernel | registers | bucket |
|---|---|---|
| `jacobi_apply` k2 | 88 → 64 | `apply` −18.7 % (measured at 16 ranks) |
| `jacobi_compute_phi` | 94 → 80 | `sweep` −12 % |
| `compute_rdenom` | 110 → 80 | `setup` −19 % |

with every gate at `max_abs 0` and no kernel spilling. The step is **−6.5 %**
from the first increment (measured at 16 ranks) and about **−1.9 %** more from
the second.
