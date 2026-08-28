# Red-black on the REFINED grid: iterations-to-residual and s/step

istmcetus GPU 1, 1 rank, 400 cold steps, `PROFILE=1`, all runs one session.
Case pair: `multiLevel_xz/refined_yp82_rect_{jacobi,redblack}.ini` — identical
grid, block shape, refinement and interface placement (448 leaves, 60.56 M
cells), only the smoother differs.

This is the §6 gate the single-level sizing could not give
(`results_smoother_2026-08-28.md`), because red-black gains a **per-colour
cross-level phi exchange** here that it does not have on one level.

`L2_div` after 400 steps is the residual proxy: it responds to `niter` on both
smoothers (Jacobi 2.93e-05 → 1.61e-05 for 6 → 12), so it measures projection
quality, not the flow.

## The measurements

| solver | niter | sor | `L2_div` | s/step |
|---|---|---|---|---|
| jacobi | 6 | 0.8 | 2.9336e-05 | 0.55388 |
| jacobi | 12 | 0.8 | 1.6053e-05 | 0.88812 |
| redblack | 3 | 1.5 | 5.1257e-05 | 0.45973 |
| redblack | 6 | 1.5 | 2.4878e-05 | 0.70033 |
| **redblack** | **6** | **1.7** | **1.9747e-05** | **0.70004** |

Fitted: Jacobi `L2_div ~ niter^-0.87` at 0.05571 s/iteration; red-black
`~ niter^-1.04` at 0.08020 s/iteration, both on a common 0.2196 s/step fixed
cost. **Red-black costs 1.440× per iteration here, against 1.308× on the
single-level grid** — that difference is the per-colour cross-level phi
exchange plus `interface_correct`.

## The omega sweep (§6's decisive number), niter = 6

| sor | 1.0 | 1.2 | 1.5 | **1.7** |
|---|---|---|---|---|
| `L2_div` | 2.5082e-05 | 2.7633e-05 | 2.4878e-05 | **1.9747e-05** |

Over-relaxation does help across an interface, and **`sor = 1.5` was not the
right setting** — 1.7 is 21 % better. Cost is flat in sor (0.6995–0.7003), as
expected.

**But 1.7 is at the stability edge**: at `niter = 3` it DIVERGED (`L2_div`
3.7e+50 by step 400) where `sor = 1.5, niter = 3` was fine. So the usable
operating point depends on `niter` as well as `sor`, and `sor = 1.7` needs
`niter ≥ 6` on this case. That is the §5 escalation-ladder trigger showing
itself; the ladder was not needed at `niter = 6`.

## Equal-residual comparison — the honest headline

Red-black at `niter = 6, sor = 1.7` reaches 1.9747e-05. Jacobi needs
`niter = 9.46` for the same residual — an **interpolation inside the measured
6–12 span**, not an extrapolation — costing 0.7465 s/step against red-black's
0.7000.

> **Red-black is 6.2 % faster at equal residual on the refined grid, against
> 26–51 % on the single-level grid. The win does not survive refinement.**

Under the *pressure-zero-mode* criterion instead (Jacobi needs `niter ≥ 12`,
red-black stays clean below 6) the comparison is 0.8881 → 0.7000, **21 %**. But
be careful quoting that: red-black at `niter = 6` leaves 1.23× Jacobi-12's
divergence residual, so the two are equal by the zero-mode criterion and *not*
equal by the residual one. Which number is right depends on which criterion the
production run actually needs.

## Is it the interface, or just the smaller grid?

Control: the same base resolution and block shape with refinement removed
(`refine_levels = 0`), `niter = 6`:

| | Jacobi | red-black (sor 1.7) | red-black advantage | cost premium |
|---|---|---|---|---|
| unrefined 2048×176×96 | 2.9395e-05, 0.31211 | 1.5354e-05, 0.34165 | **1.91×** residual | 1.095× |
| refined (60.56 M cells) | 2.9336e-05, 0.55388 | 1.9747e-05, 0.70004 | **1.49×** residual | 1.264× |

So refinement costs red-black on **both** axes, and the control separates it
from grid size: the per-iteration residual advantage falls 1.91× → 1.49×, and
the cost premium rises 1.095× → 1.264×. Together those turn a comfortable win
into a marginal one.

The mechanism is the one §5 of the red-black plan predicted: within a colour the
two sides of a level jump are relaxed *simultaneously* and patch each other, so
**cross-level coupling is additive (Jacobi-like) while everything inside a level
is Gauss-Seidel**. The interface is only ~0.3 % of the cells, but the mode it
leaves is not local — it shows up in the global residual.

## Where the time goes (s/step, niter = 12 Jacobi vs niter = 6 red-black sor 1.7)

| projection phase | Jacobi 12 | red-black 6 |
|---|---|---|
| sweep | 0.24577 | 0.38203 |
| apply (`interface_correct` on red-black) | 0.36074 | 0.00789 |
| phi_exchange | 0.02760 | 0.02550 |
| vel_exchange | 0.02094 | 0.05715 |
| apply_bc | 0.00983 | 0.00984 |
| **total** | **0.66555** | **0.48312** |

`interface_correct` is cheap (0.008 s/step — it is three plane kernels). The
red-black cost sits in the sweep and in `vel_exchange`, which is 2.7× Jacobi's
because red-black exchanges once per **colour**.

## What this changes

`docs/next_session_multirank_exchange.md` puts red-black first on the strength
of the single-level 26–51 %. On refined grids — which is what the 2:1 machinery
exists for — the honest figure is **6 % at equal residual, or 21 % under the
zero-mode criterion**, so it no longer obviously dominates the overlap and
partitioning work. Two things would move it back up, in order of cheapness:

1. **A better interface smoother.** §5's ladder — clamp `omega → min(omega,1)`
   on interface-adjacent cells, or level-ordered smoothing (relax coarse, patch,
   then fine, making the cross-level coupling multiplicative) — attacks exactly
   the 1.91× → 1.49× loss. Level ordering costs one extra phi exchange per
   colour.
2. **The partitioning fix** (`results_exchange_diag_2026-08-28.md`), which would
   make red-black's per-colour cross-level exchange a local copy instead of an
   MPI message at multi-rank.

Not measured here: multi-rank. Red-black's per-colour velocity exchange doubles
the round count against Jacobi's, and at 2 GPU ranks the exchange is already
17 % of the step, so the refined multi-rank figure is likely worse than this
single-rank 6 %.
