# The refined case at the production block shape, A6000, 2026-08-28

`multiLevel_xz/refined_yp82_rect_jacobi.ini` — the 2:1 xz-refined config at
`nb = 64 44 48` instead of the cubic `nb = 16` it shipped with. Phase 1 moved
the SINGLE-LEVEL production layout to that block shape and took 21.7 % off the
step; the refined case was never moved, so it was still paying the cubic tax.

Both refined configs were re-profiled **today with the same binary** (`0f91990`,
source identical to `2e74cc4`), because the yp100 profile on record predates
Phase 1 *and* Phase 3 and is not comparable — see "the stale profile" below.

## Gates on the config

Reported by the solver at init, and stated in the ini so a moved interface is
caught rather than silently timed:

| | leaves | cells |
|---|---|---|
| `refined_yp100_jacobi` (nb 16) | 15360 | 62.91 M |
| `refined_yp82_rect_jacobi` (nb 64 44 48) | **448** (192 coarse + 256 fine) | **60.56 M** |

Both hit their stated counts exactly.

## Result

| | s/step | ns/cell | overhead vs unblocked single level |
|---|---|---|---|
| `refined_yp100_jacobi` (nb 16) | 0.68544 | 10.895 | 1.280 |
| `refined_yp82_rect_jacobi` | **0.55189** | **9.114** | **1.071** |
| | **−19.5 %** | **−16.3 %** | |

The last column divides each case's ns/cell by `base_jacobi`'s 8.512 ns/cell
measured the same day (1.17811 s/step over 138.4 M cells) — i.e. what the run
costs above an ideal single-level, single-block rate at the same cell count.
**Moving the refined case to the production block shape removes about three
quarters of its excess cost, 28.0 % → 7.1 %.**

s/step and ns/cell differ because the two configs are not like-for-like:
`nb_y = 44` gives 4 y-tiles instead of 11, so the interface sits at y node 44
(y⁺ 82) instead of 48 (y⁺ 99) and the case carries 3.7 % fewer refined cells.
**Read the per-cell column.**

They are also not bit-comparable, and must not be: nb-independence does not hold
across a 2:1 interface (`validation/block_nb/README.md`), and these two place the
interface at different heights, so they solve different problems. The runtime
lines differ accordingly (L2_div 2.9336e-5 / 2.9209e-5) and that is expected,
not a regression.

## Where the 19.5 % came from (per step, 400 cold steps)

| step phase | nb 16 | rect | Δ s/step | Δ per cell |
|---|---|---|---|---|
| momentum | 0.16364 | 0.15380 | −6.0 % | −2.3 % |
| ibm_mu | 0.02303 | 0.01752 | −23.9 % | −21.0 % |
| bodyforce | 0.03895 | 0.03724 | −4.4 % | −0.7 % |
| vel_exchange | 0.01771 | 0.00553 | −68.8 % | −67.6 % |
| **projection** | 0.43453 | 0.33198 | **−23.6 %** | −20.6 % |
| io_stats | 0.00660 | 0.00497 | −24.7 % | −21.8 % |
| **total** | **0.68544** | **0.55189** | **−19.5 %** | **−16.3 %** |

| projection phase | nb 16 | rect | Δ s/step | share of the rect step |
|---|---|---|---|---|
| sweep | 0.14781 | 0.12138 | −17.9 % | 22.0 % |
| apply | 0.20309 | 0.17870 | −12.0 % | 32.4 % |
| phi_exchange | 0.03720 | 0.01327 | **−64.3 %** | 2.4 % |
| vel_exchange | 0.04004 | 0.01318 | **−67.1 %** | 2.4 % |
| apply_bc | 0.00560 | 0.00486 | −13.3 % | 0.9 % |
| `exch/local_copy` (all exchanges) | 0.07770 | 0.02518 | **−67.6 %** | 4.6 % |

Two mechanisms, both expected. The exchange collapses by ~2/3 because the halo
shell shrinks with the block shape: allocated-volume ratio `(18/16)³ = 1.4238`
against `(66/64)(46/44)(50/48) = 1.1230`, and halo cells are exactly what gets
copied. The volume kernels gain the smaller per-cell amounts (2–15 %) from the
same shrinkage.

## The stale profile, and what it got wrong

The yp100 profile previously on record (`runs/refined_yp100_jacobi_prof`,
0.86658 s/step) predates Phase 1 and Phase 3. Re-running it today gives 0.68544
— **20.9 % faster on the identical config**, almost all of it Phase 3's
`sync_divergence_halos` (`proj/vel_exchange` 0.13371 → 0.04004).

Anything derived from the old numbers was wrong in the same direction:

- projection halo exchange was read as **22 % of the refined step**. Fresh it is
  **11.3 %** at nb 16 and **4.8 %** at the production block shape.
- `proj/apply_bc` looked anomalous at 5× the single-level share (0.02731). Fresh
  it is 0.00560, in line with everything else. There is no anomaly to chase.

## What this means for the next target

At the production block shape the refined case's cost structure is now
essentially the single-level one:

| | refined (rect) | single-level `rect_jacobi` |
|---|---|---|
| `proj/apply` | 32.4 % | 33.4 % |
| `proj/sweep` | 22.0 % | 23.1 % |
| `step/momentum` | 27.9 % | 28.7 % |
| projection exchange | 4.8 % | 4.0 % |

**There is no longer a 2:1-specific cost centre at one rank.** The 2:1 machinery
itself costs ~3 % per cell (9.114 vs 8.809 ns/cell against single-level
`rect_jacobi`), consistent with the 3–8 % the README has always quoted for the
interface. Work on the refined case is now the same work as on the single-level
case — `jacobi_compute_phi` first (see `docs/next_session_jacobi_apply.md`).

The one place the exchange story could still be alive is **multi-rank**, which
none of this measures: `sync_divergence_halos` is single-rank only and with MPI
peers the projection falls back to the full velocity exchange 15 times per step.
At one rank that fallback is worth 0.027 s/step on the rect refined case; what
it costs over MPI is unmeasured, and measuring it needs a 4-rank profiled run of
this config.
