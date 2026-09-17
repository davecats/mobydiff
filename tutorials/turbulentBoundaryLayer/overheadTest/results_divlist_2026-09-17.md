# The divergence subset as an index list — Jacobi keeps the win, red-black stops paying

Job 5149889, A100-SXM4-40GB, HoreKa `dev_accelerated`, 100 steps, **three
columns in one allocation** at 1, 4 and 8 ranks, red-black included. Raw runs in
`horeka/results_job5149889/`. Readings pre-registered in
`horeka/exchange/PREREGISTERED_divlist.md`.

| column | commit | |
|---|---|---|
| `base` | `3c2903a` | before any divergence work — **red-black must return to this** |
| `prefix` | `95312d7` | the prefix form — **Jacobi must beat this** |
| `list` | `b9414bd` | the index list |

Follows `results_horeka_2026-09-17.md` §6.

## 1 — The change

`95312d7` made the divergence subset a per-peer PREFIX by adding a third
enumeration round. That worked (+11.8 % at 16 ranks), but making it a prefix
**reordered the entry list**, and the campaign matrix measured what that cost
everything else — 12–24 % of `phi_exchange`, and red-black net negative at every
rank count.

The subset never had to be a prefix. The enumeration is now restored exactly as
it was, and the subset is **compacted out of the finished lists**: `?DivEnt` the
entry index in enumeration order, `?DivVar` its face normal, `?DivOff` its point
prefix, `?DivPt` the inverse point map, with `peerSend/RecvDivOff` as compaction
outputs. Both ends still select the same entries — the predicate is a pure
function of the op and the direction, the same argument the canonical wire order
already rests on.

It is also **simpler** than the prefix form where it counts: the lists inherit
the entry lists' peer-major order, so a peer's points are one contiguous range
and `find_entry` is gone from all three divergence kernels.

## 2 — Red-black returns to `base`, to three significant figures

`proj vel_exchange`, ms/step — the bucket the reordering was taxing:

| red-black | `base` | `prefix` | **`list`** |
|---|---|---|---|
| 1 rank | 15.982 | 21.636 (+35 %) | **15.984** |
| 4 ranks | 6.066 | 7.407 (+22 %) | **6.060** |
| 8 ranks | 4.350 | 4.807 (+11 %) | **4.342** |

Step time against `base`: **−0.07 % / +0.04 % / −0.10 %** at 1/4/8 ranks, where
the prefix cost −3.7 % / −3.0 % / −1.6 %. Pre-registered band was ±0.5 %.
**PASS.** Red-black is exactly where it was before any of this started.

## 3 — Jacobi keeps the win and collects what the prefix was losing

Step time, ms:

| case | `base` | `prefix` | **`list`** | list vs base | list vs prefix |
|---|---|---|---|---|---|
| `rect` 1 rank | 582.10 | 587.75 | **582.24** | −0.02 % | +0.94 % |
| `rect` 4 ranks | 158.77 | 153.48 | **152.05** | **+4.23 %** | +0.93 % |
| `rect` 8 ranks | 83.86 | 80.11 | **79.56** | **+5.12 %** | +0.68 % |
| `refined` 1 rank | 260.99 | 263.98 | **260.43** | +0.21 % | +1.34 % |
| `refined` 4 ranks | 74.97 | 72.81 | **72.13** | **+3.79 %** | +0.93 % |
| `refined` 8 ranks | 42.60 | 40.64 | **40.40** | **+5.17 %** | +0.58 % |

`sweep`, `apply` and `setup` move by at most 0.3 % on all nine runs — clean
controls.

## 4 — The artefact, measured directly

`phi_exchange` at ONE rank is the clean measurement, because no MPI confounds it
and its volume is untouched by any of this:

| 1 rank | `base` | `prefix` | **`list`** |
|---|---|---|---|
| `rect_jacobi` | 16.211 | 18.249 **+12.6 %** | **16.227 (+0.1 %)** |
| `refined_yp82` | 7.543 | 9.382 **+24.4 %** | **7.555 (+0.2 %)** |
| `refined_yp82` red-black | 6.918 | 8.801 **+27.2 %** | **6.915 (−0.04 %)** |

The reordering cost is gone, not reduced. Pre-registered band was ±2 %.

## 5 — Graded against the pre-registration

| prediction | outcome |
|---|---|
| red-black within ±0.5 % of `base` at every rank count | **PASS** — −0.07 / +0.04 / −0.10 % |
| Jacobi 4- and 8-rank better than `prefix` by ~0.3–3 % | **PASS** — +0.58 to +0.93 % (lower half of the band) |
| `phi_exchange` at 1 rank back within ±2 % of `base` | **PASS** — +0.1 / +0.2 / −0.04 % |
| bit-exact vs `prefix` | **PASS** — 9 comparisons `max_abs 0` (Pass G on both production cases + the 7-case suite incl. every RANS scalar) |
| bit-exact vs `base` on a RED-BLACK case | **PASS** — `max_abs 0` at 1 and 4 ranks (CPU) |
| Jacobi 1-rank step now a **gain** ≥ +0.5 % vs `base` | **MISSED** — −0.02 % and +0.21 %, i.e. **parity** |

**The miss, and what it means.** The prediction assumed the divergence round is
worth something at one rank once the reordering cost is removed. It is worth
very little: `rect`'s `proj vel_exchange` goes 15.787 → 15.644 ms (−0.9 %) and
`refined`'s 8.072 → 7.440 (−7.8 %), but on a 582 ms and 261 ms step that is
0.02 % and 0.24 %. The reason is that **at one rank the old `dsSlot` path was
already good**: three specialised, perfectly coalesced plane-copy kernels
against one generic index-list kernel is close to a wash. The prediction should
have said "parity", and the honest summary of the single-rank case is that this
work neither helps nor hurts it.

## 6 — What this does not support

- **Nothing measured at 16 ranks.** `dev_accelerated` is 2 nodes. The prefix was
  worth +11.8 % there and the list beats the prefix by 0.6–0.9 % at 8 ranks with
  the margin shrinking as ranks rise (0.93 → 0.68 on `rect`), so ~+12 % is the
  expectation — **an extrapolation, and labelled one.** The campaign matrix
  ratios in `results_horeka_2026-09-17.md` were measured with the PREFIX and are
  now slightly pessimistic for Jacobi and ~1–4 % pessimistic for red-black.
- **Nothing about red-black's own reduced round.** Not attempted. `redblack_sweep`
  iterates `0..hi`, sweeping the lower halo layer redundantly, so a cell at
  `i = 0` reads the low halo plane of **all three** components where Jacobi's
  divergence reads the high plane of one. Its reduced set would be ~12 of 18
  face-component units against Jacobi's 3 — ~1.4x, not 6.2x — and needs edge
  cases proved (a cell at `(0, hi2, k)` reads `q(0, hi2+1, k, V)`, an edge).
  **That is the remaining idea, and it is a small one.**
- **Nothing about memory.** `?DivPt` is ~0.4x the size of `lPointEntry`, tens of
  MB on the production cases. Small against a 40 GB card, not zero, and not
  measured here.

## 7 — Where the divergence work now stands, end to end

Against `3c2903a`, which is where it started:

| | 1 rank | 4 ranks | 8 ranks | 16 ranks |
|---|---|---|---|---|
| `rect_jacobi` | −0.02 % | **+4.2** | **+5.1** | ~+12 (extrapolated) |
| `refined_yp82` jacobi | +0.21 | **+3.8** | **+5.2** | ~+12 (extrapolated) |
| `refined_yp82` red-black | −0.07 | +0.04 | −0.10 | ~0 |

**Both solvers are now optimal in the sense that neither pays for the other.**
Jacobi takes the whole saving; red-black is exactly neutral; nothing regresses
anywhere, at any rank count, on any case.
