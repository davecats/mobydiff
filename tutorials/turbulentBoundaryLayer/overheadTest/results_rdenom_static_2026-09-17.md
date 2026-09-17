# `rdenom` was being recomputed three times a step to get the same answer

Job 5150041, A100-SXM4-40GB, HoreKa `dev_accelerated`, 100 steps, before/after in
one allocation at 1, 4 and 8 ranks with red-black as a control. `ref` =
`b9414bd` (the index list), `new` = `84e8265`. Raw runs in
`horeka/results_job5150041/`. Readings pre-registered in
`horeka/exchange/PREREGISTERED_rdenom_static.md`.

This is task 2 of `docs/next_session_divergence_halos.md` — but **not the change
that handout proposes**, and that is the finding.

## 1 — The handout's premise had expired

§5 frames task 2 as making `compute_rdenom`'s per-cell `1.0d0/(...)` cheaper,
"ncu roofline before any code". It inherits a comment in `pressure_solver.f90`
saying the metric tables are static

> *"... unlike rdenom: rdenom follows `ibm%mu`, which `update_ibm_mu` rewrites
> every substage"*

**That stopped being true on 2026-09-14**, in a different increment of the same
day (`results_stepwork_2026-09-14.md` §1): `update_ibm_mu` now returns
immediately on a rank holding no body, leaving `mu = 1.0` for the whole run. Two
changes landed a day apart and nobody joined them up. On every body-free case —
the production boundary layers, the channels, Beltrami — `rdenom` depends only
on `dnLow`, `dnHigh`, `d1P` and a `mu` that never changes again. It is as static
as the tables it is contrasted with.

So the lever is not the divide. It is `if (.not. rdenomStatic)`, with
`rdenomStatic` set from a new `ibm_mu_is_unit(ibm)` over the flag
`update_ibm_mu` already caches. A rank WITH a body keeps recomputing every
substage, because there `mu` really does follow `dt_gamma`.

## 2 — Measured

| case | step ref → new | | `setup` ref → new | |
|---|---|---|---|---|
| `rect_jacobi` 1 rank | 582.05 → 562.47 ms | **+3.36 %** | 28.238 → 8.670 | −69.3 % |
| `rect_jacobi` 4 ranks | 152.10 → 147.11 | **+3.28 %** | 7.329 → 2.418 | −67.0 % |
| `rect_jacobi` 8 ranks | 79.47 → 76.96 | **+3.16 %** | 3.686 → 1.218 | −67.0 % |
| `refined_yp82` 1 rank | 260.33 → 251.71 | **+3.31 %** | 12.376 → 3.799 | −69.3 % |
| `refined_yp82` 4 ranks | 72.10 → 69.96 | **+2.97 %** | 3.230 → 1.067 | −67.0 % |
| `refined_yp82` 8 ranks | 40.48 → 39.34 | **+2.80 %** | 1.644 → 0.537 | −67.4 % |
| **red-black 1 / 4 / 8** | **+0.07 / −0.16 / −0.01 %** | | 1.943 → 1.917, 0.534 → 0.526, 0.264 → 0.270 | **unmoved** |

**Red-black is the control and it behaves**: `redblack_projection` never calls
`compute_rdenom`, its `setup` moves by −1.5 to +2.3 % (noise on a 0.26–1.9 ms
bucket) and its step by ≤0.16 %.

## 3 — Why only −67 %, and why the real number is larger

The bucket did not fall by ~97 % as the pre-registration predicted, and the
reason is not that the skip half-works. It works exactly: on CPU, taking the run
from 10 to 40 steps drops `setup` per step by **3.90x while the TOTAL stays
constant** (10.5 vs 10.8 ms). `compute_rdenom` runs **once per run**.

What remains in the bracket is **one-time cost charged to every step by a
100-step benchmark**: the first-call allocation, zeroing and device mapping of
`phi`, `delta` and `rdenom` — for `rect_jacobi` at 1 rank, three arrays of
66×46×50×1024 doubles, ~3.7 GB to allocate, zero and map. That is 8.670 ms/step
× 100 steps = **867 ms of one-time work**, which is what the "after" column is
almost entirely made of.

Netting it out gives the kernel's own cost: **6.52 ms/call** on 138.4 M cells,
i.e. ~4 doubles/cell at ~0.68 TB/s — 44 % of the A100's peak, which matches the
39.8 % of peak `results_rdenom_registers_2026-09-14.md` measured for this kernel
independently. The arithmetic is consistent.

**So the measured +2.8 to +3.4 % is a lower bound.** Production runs are
thousands of steps, over which the one-time cost vanishes and the gain
approaches the whole bucket: **4.1 to 4.9 %** depending on case and rank count.
The 100-step benchmark understates this change specifically because the thing it
removes is per-step and the thing left behind is not.

## 4 — Graded against the pre-registration

| prediction | outcome |
|---|---|
| red-black `setup` unmoved | **PASS** — −1.5 to +2.3 % on a sub-2 ms bucket |
| `les_ibm` bit-exact (the body path still recomputes) | **PASS** |
| bit-exactness everywhere | **PASS** — Pass G on both production cases (138.4 M / 60.6 M points) and the 7-case suite, `max_abs 0`, production flags |
| `setup` −80 to −98 %, step −3.5 to −4.8 % | **MISSED** — −67 to −69 % and +2.8 to +3.4 % at 100 steps |

**The miss was avoidable and the pre-registration contained its own
refutation.** It quoted the CPU measurement — `3.166 → 1.073 ms/step`, which is
−66 % — and then predicted −80 to −98 % anyway, without asking why the CPU
number was not already −97 %. Had I asked, the one-time allocation would have
surfaced before the job ran instead of after. The predicted band was wrong; the
CPU number I already had was right.

## 5 — What this does not support

- **Nothing about cases WITH a body.** There `compute_rdenom` still runs three
  times a step and the only cost added is one cached boolean. `les_ibm` gates
  that it is still correct, not that it is still fast. **The handout's divide
  question survives intact for those cases** — it is simply now the minority.
- **Nothing about 16 ranks.** `setup` is a per-cell kernel so its share should
  be roughly rank-independent, and the 1/4/8 measurements are flat at −67 %;
  that is an expectation, not a measurement.
- **The 4.1–4.9 % asymptote is an inference**, from the CPU amortisation test
  plus arithmetic, not a long GPU run. A production-length run would settle it.
- **The campaign matrix (job 5150030) does not include this.** It was queued at
  `b9414bd` and answers the divergence question; its absolute ratios are ~3 %
  stale for Jacobi from the moment this landed, and red-black's are not.
