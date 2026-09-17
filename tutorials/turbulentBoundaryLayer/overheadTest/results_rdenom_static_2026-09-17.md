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

- **Nothing about cases WITH a body** — see the correction below, which is where
  that gap got closed.
- **Nothing about 16 ranks.** `setup` is a per-cell kernel so its share should
  be roughly rank-independent, and the 1/4/8 measurements are flat at −67 %;
  that is an expectation, not a measurement.
- **The 4.1–4.9 % asymptote is an inference**, from the CPU amortisation test
  plus arithmetic, not a long GPU run. A production-length run would settle it.
- **The campaign matrix (job 5150030) does not include this.** It was queued at
  `b9414bd` and answers the divergence question; its absolute ratios are ~3 %
  stale for Jacobi from the moment this landed, and red-black's are not.


## 6 — CORRECTION (2026-09-17, same day): body cases are not the minority

The section above called cases with a body "the minority". **That was my
assumption about the workload, not a fact, and the user corrected it:
production cases often have bodies.** The per-rank form helped only body-free
cases, so as written this change was aimed at the wrong half.

It generalises, and the generalisation is the same observation one level finer:
`mu = 1/(1 + dt*coef)` is **exactly** 1.0 wherever `coef` is zero, whatever `dt`
does, so the `dt`-dependence is confined to blocks that actually hold
coefficients. `ibm_body_blocks` (one device reduction per block, once per run)
replaces `ibm_mu_is_unit`; the first projection fills every block and then
narrows `rdenomBlocks` to the body ones. An empty list is the body-free case and
is exactly the behaviour described above.

Measured on CPU, `proj setup` ms/step:

| geometry | body blocks | `setup` |
|---|---|---|
| body-free (boundary layers, channels) | 0 / N | kernel never runs again |
| `les_ibm` — plane walls spanning the domain | **256 / 640** | 8.756 → 4.003 **−54 %** |
| `sailplane` at `nb = 10` — compact body, large domain | **48 / 4500** | ~−99 % expected |
| `sailplane` with `nb` unset | **1 / 1** | **no benefit** |

The last row is the real caveat and it is not a defect: with `nb` unset there is
one block per rank, the body touches it, and there is nothing to narrow. The
gain needs block granularity — which production cases set anyway (`CLAUDE.md`
recommends `nb = 32+`, and the airfoil cases run `refine_body` with thousands of
leaves). The fraction is **printed at init** for exactly this reason, like the
trip force's block list: a silent 100 % looks identical to a silent 0 %.

So the shape of the win is geometric, not binary: it is `1 − (body blocks / all
blocks)` of a 4–5 % bucket. A compact body in a large domain — the airfoil
case — keeps nearly all of it. Walls spanning the domain keep about half.

**The handout's divide question now survives only for the cells inside body
blocks**, which is a much smaller target than it was this morning, and smaller
the more finely the case is blocked.

### 6a — The GPU A/B, and what it does and does not show (jobs 5150243, 5150257)

**Job 5150243 found a regression, and it was mine.** `ibm_body_blocks` did one
device reduction PER BLOCK — 640 launches at init on `les_ibm`, **36.5 µs each**,
which over a 50-step run cost more than the narrowing saves: `setup` **+165 %**,
step **+7.2 %**. It is now a single kernel storing into a per-block flag (every
thread that writes writes the same 1, so the race is benign and needs no atomic
or array reduction).

It was invisible on CPU, where there are no launches, and CPU is where I
measured it. **A launch-cost question cannot be validated on the platform that
has no launches** — the same shape of error as extrapolating from two rank
counts.

**Job 5150243's sailplane row was not a regression either.** `h5maxdiff`
printed `max_abs=0 ** DIFFERS **`, which is contradictory on its face and is the
NaN signature: `d != 0.0` is true for NaN while `d > m` is false. The tutorial
ships `nsteps = 1` as a smoke test and diverges when driven 50 steps, in BOTH
binaries. Dropped from the A/B; its 48/4500 block fraction is an init-time
classification and stands.

**Job 5150257, corrected, 200 steps, `les_ibm`, `max_abs 0`:**

| | ref | new | |
|---|---|---|---|
| `setup` | 0.1539 | 0.1375 ms/step | **−10.6 %** |
| step | 6.209 | 6.198 ms/step | −0.18 % |

**That is a small win, and the honest reading is that `les_ibm` is close to the
worst case for this change.** Its walls span the domain (256/640 blocks, so the
ceiling is a 60 % cut), and `setup` is only **2.5 %** of its step against 4.6 %
on `rect_jacobi`, because WALE LES adds per-step work that `rdenom` does not
scale with. Ceiling here: 0.6 × 2.5 % ≈ **1.5 % of step**, and the rest of the
bucket is the one-time allocation. −0.18 % on the step is within noise of the
−10.6 % the bucket actually moved.

### 6b — What is still NOT measured

**The magnitude on a production-sized case with a COMPACT body.** That is the
airfoil shape — `sailplane` at `nb = 10` classifies 48/4500 blocks, so ~99 % of
a 4–5 % bucket should survive — but the mechanism is all that is demonstrated
for it. `les_ibm` shows the mechanism works and is positive; it is too small and
too wall-dominated to show the size.

Measuring it needs a stable, production-sized body case: `validation/naca0012`
or `sd7003`, which need `setup.sh` (STL generation + `moby_prepare`, minutes).
**Not attempted here.** Until it is, the body-case claim is: *correct, positive,
and bounded by `1 − (body blocks / all blocks)` times the `rdenom` share of the
step* — with only the small end of that range measured.
