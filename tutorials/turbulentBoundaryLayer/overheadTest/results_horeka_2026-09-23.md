# The divergence work across the full rank sweep — and one claim that does not survive it

Job 5150030, 4 nodes (hkn[0701,0709,0711,0733]), 50 min, the 23-run matrix at
**three commits in one allocation**, 200 steps, `--map-by numa --bind-to core`.
Raw runs in `horeka/results_job5150030/`.

| column | commit | |
|---|---|---|
| `ref` | `3c2903a` | before any divergence work |
| `mid` | `95312d7` | the prefix form |
| `new` | `b9414bd` | the index list |

This is the job that spent 2026-09-19 to 09-23 **held** at `Priority=0` after a
`user env retrieval failed` requeue — see `horeka/README.md`.

## 0 — The control

The `ref → mid` column here repeats a comparison job 5147466 already made in a
different allocation four days earlier. At 16 ranks:

| | this job | job 5147466 |
|---|---|---|
| `base_jacobi` | +8.3 % | +8.6 % |
| `rect_jacobi` | +13.5 | +11.8 |
| `refined_yp82` | +9.3 | +12.0 |
| `refined_yp82` red-black | +2.3 | −0.8 |
| `refined_big` | +8.8 | +9.2 |

Same sign, same order, but spread by up to **3 percentage points**. That is the
reproducibility of a 16-rank measurement across allocations, and §2 needs it.

## 1 — The whole divergence work, `3c2903a → b9414bd`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` | +0.3 % | +0.8 | +2.5 | +4.5 | **+6.5** |
| `rect_jacobi` | +0.0 | +3.9 | +4.2 | +5.2 | **+13.5** |
| `refined_yp82_rect_jacobi` | +0.4 | +3.6 | +3.8 | +7.0 | **+8.2** |
| `refined_yp82_rect_redblack` | **−0.0** | **−0.1** | **−0.0** | **−0.4** | **+0.5** |
| `refined_big_rect_jacobi` | – | – | +3.7 | +4.0 | **+9.3** |

**The red-black row is the point of the index list, and it is now measured
across the full sweep instead of chained from two A/Bs: neutral to within
0.5 % at every rank count.** The prefix form had cost it 0.8–3.8 %; that is
gone, and nothing was given up to remove it.

## 2 — CORRECTION: the index list does not measurably help at 16 ranks

`results_divlist_2026-09-17.md` §6 predicted the index list would "add a little
at 4–16 ranks", extrapolating a margin that was shrinking with rank count
(0.93 % at 4 ranks → 0.68 % at 8). Measured, `95312d7 → b9414bd`:

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` | +0.3 % | +0.1 | +0.1 | +0.3 | **−1.9** |
| `rect_jacobi` | +0.9 | +0.9 | +0.9 | +0.8 | **+0.1** |
| `refined_yp82_rect_jacobi` | +1.3 | +1.3 | +1.0 | +1.1 | **−1.2** |
| `refined_yp82_rect_redblack` | +3.4 | +3.1 | +3.0 | +1.7 | **−1.8** |
| `refined_big_rect_jacobi` | – | – | +1.1 | +1.0 | **+0.6** |

At 1–8 ranks it is a consistent **+0.8 to +1.3 %** (and +3 % for red-black,
which is the regression being handed back). **At 16 ranks every entry is inside
the ±2–3 point noise band §0 measured.** The honest statement is *no measurable
effect at 16 ranks*, not "a little". The extrapolation was labelled as one and
it was wrong in the direction the shrinking margin already hinted at.

This changes nothing about keeping the index list: it is neutral-to-positive
everywhere, it removes red-black's regression outright, and it is simpler than
the prefix (no `find_entry` in the divergence kernels).

## 3 — REVISED ratios

Block tax (`rect_jacobi` against `base_jacobi`, inverse per-GPU throughput):

| ranks | published 09-17 | **now** |
|---|---|---|
| 1 | 1.036 | **1.015** |
| 2 | 1.053 | **1.034** |
| 4 | 1.038 | **1.024** |
| 8 | 1.010 | **1.000** |
| 16 | 1.033 | **1.007** |

**The block tax at 16 ranks is now 1.007 — the blocked single-level case costs
essentially nothing against the unblocked reference, at every rank count.**

Strong-scaling efficiency at 16 ranks: `base` 75 → **80 %**, `rect` 70 → **80**,
`refined` 55 → **59**, `redblack` 52 → **52**, `refined_big` 83 → **88**.

The 2:1 machinery, like for like: **1.001x** coarse-cell-equivalents at 4 ranks,
0.984x at 8, **0.879x** at 16 (was 0.900x).

## 4 — What this does not support

- **It does not include the `rdenom` work** (`84e8265` onward) or the rough-wall
  benchmark. Those landed after the job was submitted; the `rdenom` increment is
  measured only at 1/4/8 ranks (jobs 5150041, 5151978) and adds a further
  ~3.5 % on body-free cases. **Another matrix pass is owed.**
- **16-rank differences below ~3 % are not resolvable** in a single allocation
  of this matrix. §0 is the evidence. Anything smaller needs repeated
  allocations, which no result here has.
- **Nothing about red-black's own scaling.** Its 16-rank efficiency is 52 %
  before and after; it was never optimised, only stopped from being taxed.
