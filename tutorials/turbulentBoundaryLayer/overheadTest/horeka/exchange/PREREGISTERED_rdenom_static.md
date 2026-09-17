# Pre-registered readings — rdenom is static on a body-free rank

Written **before** the job ran, at `fc816e1`. This is task 2 of
`docs/next_session_divergence_halos.md`, but **not the change that handout
proposes**, and the difference is the point.

## The handout's premise is stale

§5 frames task 2 as making `compute_rdenom`'s per-cell `1.0d0/(...)` cheaper,
with "ncu roofline before any code". That inherits a comment in
`pressure_solver.f90` which says the metric tables are static

> *"... unlike rdenom: rdenom follows `ibm%mu`, which `update_ibm_mu` rewrites
> every substage"*

**That stopped being true on 2026-09-14**, in a different increment of the same
day (`results_stepwork_2026-09-14.md` §1): `update_ibm_mu` now returns
immediately on a rank holding no body, leaving `mu = 1.0` for the entire run. So
on every body-free case — the production boundary layers, the channels,
Beltrami — `rdenom` depends only on `dnLow`, `dnHigh`, `d1P` and a `mu` that
never changes. **It is as static as the tables it is contrasted with, and it is
being recomputed three times a step to get the same answer.**

So the change is not to make the divide cheaper. It is to stop running the
kernel: compute `rdenom` once, then set `rdenomStatic` from a new
`ibm_mu_is_unit(ibm)` accessor over the flag `update_ibm_mu` already caches. A
rank WITH a body keeps recomputing every substage, because there `mu` really
does follow `dt_gamma`.

Bit-exact by construction and for the same reason as the `ibm_mu` increment:
same inputs, same expression, same result — the recomputation was writing the
values that were already there.

## The bucket being targeted (job 5149889, `list` column — the current code)

`proj setup`, which after the first call is `compute_rdenom` plus a handful of
early-returning host calls:

| case | `setup` ms/step | of step |
|---|---|---|
| `rect_jacobi` 1 rank | 28.267 | **4.9 %** |
| `rect_jacobi` 4 ranks | 7.317 | **4.8 %** |
| `rect_jacobi` 8 ranks | 3.678 | **4.6 %** |
| `refined_yp82` 4 ranks | 3.226 | **4.5 %** |
| `refined_yp82` 8 ranks | 1.642 | **4.1 %** |

On CPU (`min_channel`, 4 ranks, 10 steps) the bucket already reads
**3.166 → 1.073 ms/step with the call count unchanged at 30** — the bracket
still runs, the kernel does not.

## What the job must return

| result | conclusion |
|---|---|
| `setup` **−80 to −98 %** on the Jacobi cases, step **−3.5 to −4.8 %** | as designed; take it |
| **red-black `setup` unmoved** (0.265–1.928 ms/step) | the control: `redblack_projection` never calls `compute_rdenom`, so it must not move. If it does, something else is being measured |
| `les_ibm` bit-exact | the case WITH a body still recomputes every substage — the correctness half |
| `setup` falls much less than the CPU run suggests | the bucket is not mostly `compute_rdenom` on GPU; report the residual and what it is before quoting a win |
| step falls by materially more than `setup` did | distrust it; nothing else changed |
| any nonzero `max_abs` | a real difference. `les_ibm` covers the body path and the 7-case suite the rest |

**Gate flags: PRODUCTION.** No expression moves — the kernel is skipped, not
rewritten — so both sides contract identically and the production comparison is
the tighter test.

## What this cannot say

- **Nothing about cases WITH a body.** There `compute_rdenom` still runs in full,
  three times a step, and the only cost added is one cached boolean. `les_ibm`
  gates that it is still correct, not that it is still fast. **The handout's
  divide question survives untouched for those cases** — and is worth less,
  since they are the minority.
- **Nothing about 16 ranks** (`dev_accelerated` is 2 nodes). `setup` is a
  per-cell kernel, so its share should be roughly rank-independent; that is an
  expectation, not a measurement.
- **It does not close task 2 as written.** The per-cell divide is untouched. It
  makes it mostly moot for the cases the campaign runs.
