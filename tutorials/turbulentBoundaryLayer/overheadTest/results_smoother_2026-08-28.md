# Smoother cost at equal quality, A6000, 2026-08-28

Single-level production layout (`nb` unset, 138.4 M cells), istmcetus GPU 1,
400 cold steps, all four runs in one session with `PROFILE=1`. Each smoother at
its own relaxation factor (`sor = 0.8` damped Jacobi, `1.5` SOR), `NITER=` on
`run_overhead.sh`.

Red-black cannot cross a 2:1 interface yet — that is what
`docs/next_session_redblack_interface.md` R1 is for — so single-level is the
only grid where both smoothers run, and this is a *sizing* measurement for R1,
not a measurement of R1.

## The quality setting

User measurement (2026-08-28): **damped Jacobi needs `niter >= 12`** to keep the
pressure zero-mode away; **red-black stays clean below `niter = 6`** at similar
residuals. Everything in this directory has been running `niter = 6` Jacobi, so
**every recorded baseline in this repository is under-iterated** and understates
the true production cost.

The `L2_div` column below is an independent check on that ordering, and it
agrees:

| run | `L2_div` after 400 steps | s/step |
|---|---|---|
| Jacobi `niter = 6` | 1.845e-05 | 1.1758 |
| **Jacobi `niter = 12`** | 1.234e-05 | **1.8924** |
| **red-black `niter = 3`** | 1.653e-05 | **0.9243** |
| red-black `niter = 6` | 6.556e-06 | 1.3929 |

Red-black at `N` reaches a *better* residual than Jacobi at `2N`: `niter = 3`
beats Jacobi's 6 (1.65e-05 vs 1.85e-05) and `niter = 6` beats Jacobi's 12 by
1.9× (6.56e-06 vs 1.23e-05).

## The comparisons that matter

| at equal or better quality | Jacobi | red-black | gain |
|---|---|---|---|
| RB 3 vs Jacobi 6 | 1.1758 | 0.9243 | **−21.4 %** (and a better residual) |
| RB 6 vs Jacobi 12 | 1.8924 | 1.3929 | **−26.4 %** (and a 1.9× better residual) |

Interpolating red-black's projection (linear in `niter`: 0.1558 s/step per
iteration) against Jacobi at `niter = 12`:

| red-black `niter` | s/step | vs Jacobi 12 |
|---|---|---|
| 3 | 0.924 | −51.2 % |
| 4 | 1.088 | −42.5 % |
| 5 | 1.243 | −34.3 % |
| 6 | 1.393 | −26.4 % |

**So the prize is 26–51 % of the step**, depending where between 3 and 6
red-black is acceptable. That is larger than everything else on the multi-rank
list combined.

## The colour-sweep penalty is real, and the win survives it

At **equal** `niter = 6`: Jacobi 1.1758, red-black 1.3929 — red-black is
**18.5 % more expensive per iteration**, reproducing the 17 % measured on
2026-08-07 (1.4499 vs 1.2422). A colour sweep wastes bandwidth on the unused
colour. The gain above is `~2× fewer iterations` against that 1.185 penalty,
plus the fixed non-projection cost (~0.464 s/step, identical in all four runs)
which dilutes both.

## Where the time goes (s/step)

| | Jacobi 12 | RB 6 | RB 3 |
|---|---|---|---|
| momentum + ibm_mu + bodyforce + io | 0.464 | 0.456 | 0.455 |
| **projection** | **1.4281** | **0.9354** | **0.4675** |
| — `sweep` | 0.5719 | 0.9135 | 0.4564 |
| — `apply` | 0.8250 | — (in place) | — |
| — `phi_exchange` | 0.0086 | — (no phi) | — |
| — `vel_exchange` | 0.0065 | 0.0071 | 0.0036 |
| — `apply_bc` | 0.0147 | 0.0148 | 0.0074 |

Red-black has no separate `apply` (the correction is in place) and no `phi`
exchange at all on a single level, which is where the exchange-count halving
shows up: 12 rounds per substage at `niter = 6` against Jacobi's 24 at
`niter = 12`.

## Caveat — this does NOT size the refined case

Two things this measurement does not include:

1. **The 2:1 interface.** R1 gives each colour a cross-level phi exchange, which
   single-level red-black does not have. On the refined 2-rank GPU decomposition
   that is 7.39 MB per round over MPI (98.1 % of the phi entries are cross-level
   — `results_exchange_diag_2026-08-28.md`), so the refined win will be smaller
   than the numbers above and cannot be quoted from them.
2. **Multi-rank.** All four runs are 1 rank. Red-black's per-colour velocity
   exchange becomes an MPI round with peers.

The round accounting for the refined case is in
`docs/next_session_redblack_interface.md`: red-black at `N` costs the same
rounds as Jacobi at `2N` once the interface phi exchange is present, so at the
equal-quality settings the rounds are equal and the win is entirely compute.
