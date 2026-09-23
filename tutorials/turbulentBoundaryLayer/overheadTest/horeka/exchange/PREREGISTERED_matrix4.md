# Pre-registered readings — the rdenom work and the rough wall, full rank sweep

Written **before** the job ran, at `8fa0fc2`. Closes the two gaps
`results_horeka_2026-09-23.md` §4 names: that matrix predates the `rdenom` work
and the rough-wall benchmark.

## The three columns

| column | commit | runs |
|---|---|---|
| `ref` | `b9414bd` | the 5 body-free configs — **not** `rough_jacobi` |
| `mid` | `809759e` (`bench/rdenom-always`) | all 6 — HEAD with the narrowing alone disabled |
| `new` | `8fa0fc2` | all 6 |

**`ref` must not run `rough_jacobi`, and this is a trap rather than a
preference.** `config.f90` has no `case default`: an unknown key inside a known
section is **silently ignored**. A binary predating `7b2bc2a` therefore runs
`rough_jacobi` with `wall_shape = eggcarton` discarded, falling back to the
hardcoded 2D wavy wall — a different geometry, reported as a success. The column
is gated with `CONFIGS`.

## The four questions, and the pairing that answers each

| # | question | pairing |
|---|---|---|
| 1 | What is the `rdenom` narrowing worth on body-free cases, at every rank count? | `mid → new`, 5 configs |
| 2 | What is it worth on a case WITH a body, at every rank count? | `mid → new`, `rough_jacobi` |
| 3 | Is the roughness feature really dormant on body-free cases? | `ref → mid`, 5 configs — **a control** |
| 4 | What does the roughness cost, like for like? | `rough_jacobi` vs `rect_jacobi` within one column |

Question 3 is free and worth having: `ref → mid` differs by the (dormant)
roughness feature **plus** the `rdenomBlocks` indirection, which `mid` carries
but never benefits from. If that pairing is not ~0, the indirection costs
something and `mid → new` is flattered by exactly that amount.

## Predictions

| reading | prediction |
|---|---|
| Q1, `rdenom` on body-free, 1–8 ranks | **−3 to −4.5 %** of the step (measured −3.5 % at 4 and 8) |
| Q1 at 16 ranks | **smaller, −2 to −3.5 %**: `setup` is per-cell work removed, and at 8.65 Mcell/GPU the communication share is larger, so the same absolute saving is a smaller fraction |
| Q2, `rdenom` on `rough_jacobi` | **−2 to −2.8 %**, i.e. ~0.75× the body-free figure, because 25 % of its blocks hold the body (measured −2.46 / −2.32 % at 4 / 8) |
| Q3, `ref → mid` on body-free | **within ±1 %** of zero. A consistent negative across configs means the `rdenomBlocks` indirection is real and must be subtracted from Q1 |
| Q4, roughness cost | **+4 to +6 %** of the step against `rect_jacobi` (measured +5.5 / +5.0 % at 4 / 8) |
| any 16-rank difference below ~3 points | **not resolvable** — `results_horeka_2026-09-23.md` §0 measured that band by repeating one comparison across allocations. Anything inside it is reported as "no measurable effect", not as a number |

## What this campaign cannot say

- **Nothing physical about the rough wall.** `rough_jacobi` is a benchmark: the
  height and wavelengths were chosen to be resolved and to sit inside the first
  y-block, not to reproduce a published Re_τ or k_s⁺. Q4 is a *cost*, not a
  validation.
- **Nothing about compact bodies.** `rough_jacobi` holds its body in one y-block
  (25 % of blocks). An airfoil-shaped case (`sailplane` classifies 1.1 %) would
  sit near the body-free figure and is still unmeasured.
- **Nothing about bit-exactness.** That is carried by the gate jobs, not here;
  this matrix runs 200 steps and compares clocks only.
- **16 ranks remains one allocation.** The band above is the reason every
  16-rank claim here will be hedged, and no conclusion should rest on a single
  16-rank pair.
