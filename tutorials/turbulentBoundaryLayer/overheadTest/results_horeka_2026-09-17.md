# The campaign re-measured, in three columns — and the divergence halo's real shape

Job 5147466, 4 nodes (hkn[0532,0627,0720,0733]), 52 min, the unchanged 23-run
matrix run **three times in one allocation**, 200 steps,
`--map-by numa --bind-to core`. Raw runs in `horeka/results_job5147466/`.

| column | commit | what it adds |
|---|---|---|
| `ref` | `55bee89` | the `new` column of `results_horeka_2026-09-14.md` — the **control** |
| `mid` | `3c2903a` | + the two register cuts and the three step-work increments (09-14) |
| `new` | `95312d7` | + the divergence-halo exchange (09-15) |

This supersedes `results_horeka_2026-09-14.md` §2–§4 and the ratios quoted from
it in `CLAUDE.md`. Collectors: `scaling_total.md` (`ref`→`new`),
`scaling_registers_stepwork.md` (`ref`→`mid`), `scaling_divhalo.md`
(`mid`→`new`).

## 0 — The control

`ref` is the same commit whose column `results_horeka_2026-09-14.md` published.
Its block tax here reads 1.049 / 1.091 / 1.069 / 1.042 / 1.088 at 1/2/4/8/16
ranks against 1.053 / 1.095 / 1.070 / 1.043 / 1.100 there — different nodes,
three days later. The matrix is sound and the other two columns can be read as
the changes and nothing else.

## 1 — The total: what has come off the step since `map(to: c)`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` | +18.0 % | +18.3 | +18.8 | +19.4 | **+21.3** |
| `rect_jacobi` | +19.0 | +21.2 | +21.2 | +21.8 | **+25.2** |
| `refined_yp82_rect_jacobi` | +18.3 | +20.5 | +20.3 | +21.7 | **+23.4** |
| `refined_yp82_rect_redblack` | +6.3 | +6.6 | +6.8 | +7.8 | **+5.4** |
| `refined_big_rect_jacobi` | – | – | +21.2 | +21.9 | **+24.0** |

**Red-black gets a fifth of what Jacobi gets**, and §4 says why.

## 2 — Split by increment

| config, 16 ranks | registers + step work | divergence halo | total |
|---|---|---|---|
| `base_jacobi` | +13.8 % | **+8.6** | +21.3 |
| `rect_jacobi` | +15.3 | **+11.8** | +25.2 |
| `refined_yp82_rect_jacobi` | +12.9 | **+12.0** | +23.4 |
| `refined_yp82_rect_redblack` | +6.1 | **−0.8** | +5.4 |
| `refined_big_rect_jacobi` | +16.2 | **+9.2** | +24.0 |

The two increments behave in **opposite directions with rank count**. The
register and step-work increments are per-cell work removed, so their fraction
*shrinks* as the per-GPU problem shrinks (`rect_jacobi` +19.8 % at 1 rank →
+15.3 % at 16). The divergence halo is message volume removed, so its fraction
*grows*: −1.0 % at 1 rank → **+11.8 % at 16**.

## 3 — REVISED ratios

Block tax (`rect_jacobi` against `base_jacobi`, inverse per-GPU throughput):

| ranks | published 09-14 | **now** |
|---|---|---|
| 1 | 1.053 | **1.036** |
| 2 | 1.095 | **1.053** |
| 4 | 1.070 | **1.038** |
| 8 | 1.043 | **1.010** |
| 16 | 1.100 | **1.033** |

Strong-scaling efficiency at 16 ranks: `base` 78 → **82 %**, `rect` 74 → **82**,
`refined` 60 → **65**, `redblack` 56 → **57**, `refined_big` 86 → **89**.

The 2:1 machinery, like for like (cost of the cells `refined_big` adds, against
the coarse cells beside them): **1.002x** at 4 ranks (was 1.004), **0.986x** at
8, **0.900x** at 16 (was 0.826). The 16-rank figure moves most because the
divergence saving is worth more to the case with fewer cells per message.

**Red-black against Jacobi — the largest revision.** Ratio of `redblack` to
`jacobi` step time on the same refined case:

| | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| published 09-14 | 0.76 | – | – | – | 0.808 |
| `ref` here | 0.759 | 0.748 | 0.761 | 0.782 | 0.804 |
| **`new` here** | **0.870** | **0.878** | **0.890** | **0.920** | **0.993** |

Red-black used to be ~20 % faster than Jacobi at 16 ranks. **It is now 0.7 %
faster — parity.** Nothing was done to red-black; Jacobi got faster and
red-black did not, because red-black does not use the projection's
between-iteration velocity round at all.

## 4 — What the divergence halo actually trades (the finding)

`results_divhalo_2026-09-15.md` measured this change at 4 and 8 ranks on two
nodes and called the `phi_exchange` side effect "observed, not explained". The
matrix explains it, and it is **larger and more consequential than that report
implies**.

The change does two things at once:

1. **It removes message volume from 15 of the 39 exchange rounds.** Per-round
   `mpi_wait` at 16 ranks: `rect` 170.6 → **80.0 µs**, `refined` 105.2 →
   **54.9**, `refined_big` 199.2 → **71.4** (−53 %, −48 %, −64 %).
2. **It reorders the entry list, and that costs every OTHER exchange.** The
   copies used to run all 26 directions of block 1, then all 26 of block 2;
   they now run the three `+axis` faces of every block, then the other 23 of
   every block. Same points, same values, worse destination locality.

**The clean measurement of (2) is at one rank**, where there is no MPI to
confound it — and a controlled pair falls out of the matrix:

| 1 rank, same solver, same change | blocks | `phi_exchange` | step |
|---|---|---|---|
| `base_jacobi` | **1** | −1.7 % | **−1.85 %** (a gain) |
| `rect_jacobi` | **1024** | **+11.9 %** | +0.96 % (a loss) |
| `refined_yp82_rect_jacobi` | 1024 | **+24.2 %** | +1.35 % (a loss) |

Same rank count, same binary pair; the only difference is how many entries
there are to reorder. With one block there is nothing to reorder, and the
divergence round alone is worth 1.85 %. With 1024 blocks the reordering costs
the scalar exchange 12–24 % and the change goes **net negative**.

So the trade is: **a fixed reordering cost on every exchange, against a large
message saving on 15 of 39 rounds.** Where messages dominate it wins by a lot
(+12 % at 16 ranks). Where they do not, it loses:

- **1 rank: −1.0 % (`rect`), −1.3 % (`refined`).**
- **Red-black at every rank count: −3.8 / −3.5 / −3.2 / −1.6 / −0.8 %.** Its
  `mpi_wait` per round is unchanged (96.2 → 97.4 µs at 16 ranks) because it
  never runs the reduced round; it pays the reordering and collects nothing.
  At 4 ranks its `proj vel_exchange` is **+23.1 %** and its `local_copy`
  **+30.7 %**, pure reordering cost.

### Correcting `results_divhalo_2026-09-15.md`

Two statements in that report do not survive this matrix:

- *"the largest single contributor is the device-local copy, not the message,
  which is why the change is worth something at one rank too."* **Wrong at one
  rank.** The −27 % `local_copy` it measured is real at 4 ranks, but at 1 rank
  the same bucket goes **+46 % / +64 %** and the step is a small net loss. Part
  of that is re-attribution — the old `dsSlot` kernels were deliberately not in
  `exch_timing` — but the step is the arbiter and it went the wrong way.
- *"`phi_exchange` … −1.7 % and −2.9 % at 8 ranks"*, read as the reordering
  sometimes helping. It does not: at 8 ranks the whole step got faster and
  ranks arrived more in sync, so `phi_exchange`'s share of `mpi_wait` fell. The
  reordering's own cost is the 1-rank number, and it is always a cost.

## 5 — What this does not support

- **Nothing about whether the reordering cost is necessary.** It is an artefact
  of making the divergence set a *prefix*. A design that leaves the entry order
  alone and drives the divergence round from an explicit index list would pay
  none of it — see §6 — but that is an argument, not a measurement.
- **Nothing about red-black's own scaling.** Its 16-rank efficiency (57 %) and
  its `mpi_wait` are unchanged by both increments; it was not optimised and is
  not claimed to be.
- **The 2:1 coarse-cell-equivalent below 1.0 at 8 and 16 ranks** still means
  what it meant: the refined case has fewer cells per launch and per message, so
  a per-launch or per-message saving helps it more than its coarse twin. It is
  not evidence that refinement is free.

## 6 — The concrete follow-up this measurement justifies

Keep the two-round entry enumeration exactly as it was, and drive the
divergence round from explicit index lists (`lDivEnt`/`sDivEnt`/`rDivEnt` plus
their point prefixes) rather than from a prefix of the reordered list. Roughly
six extra integer arrays; the divergence kernels are otherwise unchanged.

That was considered first and rejected on "more state" grounds. **The
measurement overturns those grounds**: the prefix costs 12–24 % of the scalar
exchange at 1 rank, turns red-black net negative at every rank count, and makes
the change a loss on single-rank runs. Expected value: recovers the red-black
regression entirely, turns the 1-rank result from −1.3 % into a gain, and adds
a little at 4–16 ranks where the reordering cost is currently swamped rather
than absent.
