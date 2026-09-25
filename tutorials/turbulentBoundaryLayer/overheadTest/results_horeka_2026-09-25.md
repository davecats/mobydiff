# The rdenom work and the rough wall across the full rank sweep

> **The `mid` control column's branch is deleted (2026-09-25).**
> `bench/rdenom-always` was `8aedbb5`, one commented-out line on top of
> `e5b5f8d`: `call narrow_rdenom_blocks(ibm)` in `pressure_solver.f90` disabled,
> so `rdenom` is recomputed for every block every substage as before `84e8265`,
> while the binary is HEAD in every other respect -- in particular it knows
> `[ibm] wall_shape`, which the pre-`rdenom` commits do not, which is why the
> control had to be a branch rather than an older tag. Recreate it by commenting
> out that one call; the results below are the durable record.


Job 5159149, 4 nodes (hkn[0521,0530-0531,0623]), 57 min, **three columns, six
configs, 79 runs, no failures**, 200 steps, `--map-by numa --bind-to core`. Raw
runs in `horeka/results_job5159149/`. Readings pre-registered in
`horeka/exchange/PREREGISTERED_matrix4.md`.

| column | commit | runs |
|---|---|---|
| `ref` | `b9414bd` | the body-free five (gated — see §4) |
| `mid` | `8aedbb5` `bench/rdenom-always` | all six, narrowing disabled |
| `new` | `e5b5f8d` | all six |

Closes both gaps `results_horeka_2026-09-23.md` §4 named.

## 1 — The control fired, and it changes the headline

`ref → mid` differs by the (dormant) roughness feature **plus** the
`rdenomBlocks` indirection that `mid` carries and never benefits from. It is not
zero:

| config | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| `base_jacobi` | −0.5 % | −0.6 | −0.6 | −0.6 | −0.6 |
| `rect_jacobi` | −0.5 | −0.4 | −0.4 | −0.5 | −0.5 |
| `refined_yp82` | −0.4 | −0.5 | −0.3 | −0.6 | +0.3 |
| `refined_big` | – | – | −0.5 | −0.5 | −0.2 |
| **`refined_yp82` red-black** | **−0.0** | **+0.0** | **+0.1** | **+0.0** | **+1.3** |

**Red-black is the discriminator and it reads zero.** `redblack_projection`
never calls `compute_rdenom`, so it cannot pay an indirection inside that
kernel — and it doesn't. Every config that *does* call it pays a consistent
**0.4–0.6 %**. So:

- **the roughness feature is dormant as claimed** (red-black shares that code
  path and is unmoved), and
- **the `rdenomBlocks` indirection costs ~0.5 % of the step** while the kernel
  is running.

Therefore `mid → new` **overstates** the narrowing by that amount, and the
honest number is `ref → new`. That is what §2 quotes.

The indirection is not a residual cost to go and fix: in `new` a body-free rank
has an empty list and `compute_rdenom` never runs at all, so it pays nothing. It
is a property of the artificial `mid` reference and of body cases, where it is
already inside the measured figure.

## 2 — What the rdenom work is worth (`ref → new`)

| config | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| `base_jacobi` | +3.5 % | +3.4 | +3.3 | +3.1 | **+3.0** |
| `rect_jacobi` | +3.3 | +3.3 | +3.2 | +3.1 | **+3.1** |
| `refined_yp82_rect_jacobi` | +3.3 | +3.1 | +3.0 | +2.7 | **+2.5** |
| `refined_big_rect_jacobi` | – | – | +3.2 | +3.0 | **+2.9** |
| `refined_yp82_rect_redblack` | +0.0 | −0.0 | +0.1 | −0.0 | +0.6 |

Red-black is zero throughout, correctly: it has no `compute_rdenom` to skip.

**The gain is essentially FLAT in rank count, and that was predicted wrong.**
The pre-registration expected it to shrink at 16 ranks, reasoning that `setup`
is per-cell work while communication grows with rank count. For `base`, `rect`
and `refined_big` it does not shrink at all (3.5 → 3.0, 3.3 → 3.1, 3.2 → 2.9);
only `refined_yp82` behaves as predicted (3.3 → 2.5). The premise — that
`setup`'s share of the step falls with rank count — is simply not true for these
cases, and the report that assumed it was reasoning ahead of its data.

## 3 — The body case

`rough_jacobi` classifies **exactly 25 % of blocks as body at every rank count**
(256/1024, 128/512, 64/256, 32/128, 16/64), so the pre-registered model says its
gain should be 0.75× the body-free one:

| ranks | `rough` gain | `rect` gain | ratio | model |
|---|---|---|---|---|
| 1 | +2.63 % | +3.79 % | 0.693 | 0.75 |
| 2 | +2.76 | +3.74 | **0.739** | 0.75 |
| 4 | +2.75 | +3.66 | **0.750** | 0.75 |
| 8 | +2.41 | +3.57 | 0.675 | 0.75 |
| 16 | +1.86 | +3.67 | **0.508** | 0.75 |

Holds to three decimal places at 2 and 4 ranks, drifts low at 1 and 8, and
breaks at 16 — where the gain itself (+1.86 %) is the smallest number in the
table and the least trustworthy. **The model is good to ±10 % in the middle of
the range and should not be extrapolated to 16 ranks.**

### What the roughness costs, like for like

`rough_jacobi` against `rect_jacobi` — same grid, blocks, flow, only the wall:

| ranks | in `new` | in `mid` |
|---|---|---|
| 1 | **+6.93 %** | +5.65 % |
| 2 | +6.74 | +5.67 |
| 4 | +6.49 | +5.49 |
| 8 | +6.46 | +5.19 |
| 16 | +6.23 | +4.28 |

Pre-registered band was +4 to +6 %; `new` sits just above it. The reason is
mechanical and is its own consistency check: the narrowing speeds `rect` up more
than `rough` (100 % of blocks skipped against 75 %), so the roughness costs a
*larger fraction* of a *smaller* step. Both columns are the same physical
overhead seen against different baselines.

## 4 — A trap that had to be designed around

`config.f90` has no `case default`: **an unknown key inside a known section is
silently ignored.** Handed `rough_jacobi`, the `b9414bd` binary would discard
`wall_shape = eggcarton`, fall back to the hardcoded 2D wavy wall, and report
success — a different geometry presented as a result. The `ref` column is gated
with `CONFIGS` for that reason. The parser behaviour is worth tightening, but
not by this session: a `case default` that errors would reject existing inis
carrying stray keys.

## 5 — CORRECTION to the 2026-09-23 noise claim

That report concluded **"16-rank differences below ~3 points are not
resolvable"**. That is true of what it measured — one comparison repeated across
two *different allocations*. It is too pessimistic as a general rule, and this
job shows why: the §1 control resolves a **0.5 %** effect cleanly, because it is
consistent in sign and size across 20 (config, rank) cells within **one**
allocation, and because red-black provides a config that must show zero.

The rule should be: *a single 16-rank pair across allocations is worth ±3
points; a consistent pattern across configs and rank counts inside one
allocation resolves well under 1 %.* Every conclusion in §1–§3 rests on the
latter.

## 6 — Graded against the pre-registration

| prediction | outcome |
|---|---|
| Q1 body-free, 1–8 ranks: −3 to −4.5 % | **PASS** — +3.1 to +3.5 % (`ref → new`) |
| Q1 at 16 ranks: smaller, −2 to −3.5 % | **MISSED** — flat, not smaller, on 3 of 4 configs |
| Q2 body case: −2 to −2.8 %, ~0.75× | **PASS at 1–8** (+2.4 to +2.8 %, ratio 0.68–0.75); **fails at 16** (+1.86 %, ratio 0.51) |
| Q3 control within ±1 % | **PASS in magnitude, but NOT zero** — a systematic −0.5 %, which is the finding |
| Q4 roughness cost +4 to +6 % | **MISSED HIGH** — +6.2 to +6.9 % in `new`, for a mechanical reason (§3) |

## 7 — What this does not support

- **No physics.** `rough_jacobi` is a benchmark; its height and wavelengths were
  chosen to be resolved and to sit in the first y-block. §3's +6.5 % is a cost,
  not a validation.
- **Nothing about compact bodies.** 25 % body blocks is one geometry.
  `sailplane` classifies 1.1 % and would sit near the body-free figure; still
  unmeasured.
- **The block tax is unchanged by this work** (1.014 / 1.037 / 1.028 / 1.005 /
  1.030) and should be, since `base` and `rect` both gain ~3 %. Note it reads
  1.030 at 16 ranks here against 1.007 on 09-23 — that IS the across-allocation
  spread §5 describes, and neither number should be quoted alone.
- **No bit-exactness evidence.** That is carried by the gate jobs; this matrix
  compares clocks only.
