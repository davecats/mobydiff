# Pre-registered readings — the divergence subset as an index list, not a prefix

Written **before** the job ran, at `cb0a9c7`. Follows
`results_horeka_2026-09-17.md` §6, which is what justifies it.

## What is being changed, and why

`95312d7` made the divergence subset a per-peer PREFIX of the copy prefix by
adding a third enumeration round. That works, and at 16 ranks it is worth
+11.8 %. But making it a prefix **reorders the entry list**, and the matrix
measured what that costs everything else:

| | measured |
|---|---|
| `phi_exchange` at 1 rank, 1024 blocks | **+11.9 % / +24.2 %** |
| the same at 1 rank, `base_jacobi`, ONE block | −1.7 % (nothing to reorder) |
| red-black step, every rank count | **−3.8 / −3.5 / −3.2 / −1.6 / −0.8 %** |
| Jacobi step at 1 rank | **−1.0 / −1.3 %** |

The subset does not have to be a prefix. It can be an explicit index list —
`lDivEnt`/`sDivEnt`/`rDivEnt` with their own point prefixes and per-entry
component — leaving the enumeration EXACTLY as it was in `3c2903a`. The
divergence round then reads the same entries in the same order and delivers the
same values, while every other exchange goes back to byte-identical behaviour.

Three arrays per side (entry index, component, point prefix) plus a
point→index map, and `peerSendDivOff`/`peerRecvDivOff` become compaction
outputs rather than enumeration prefixes. **It also removes `find_entry` from
the divergence kernels**: the peer range is a property of the compact list, so
the binary search the prefix form needed is gone.

This design was considered first and rejected for carrying more state. The
measurement overturns that: the state is ~6 integer arrays, and the prefix
costs 12–24 % of a bucket that every solver configuration pays.

## Why red-black cannot simply share the reduced round

`redblack_sweep` iterates `i,j,k = 0..hi` — it sweeps the LOWER halo layer
redundantly with the owning neighbour, which is what makes its results
independent of the block and rank layout. A cell at `i = 0` reads
`q(0,j,k,U)`, `q(1,j,k,U)`, `q(0,j,k,V)`, `q(0,jp,k,V)`, `q(0,j,k,W)`,
`q(0,j,kp,W)`: the low halo plane of **all three** components. Jacobi's
divergence reads the high plane of **one**.

So red-black's own reduced set would be roughly the three low faces × 3
components + the three high faces × 1 = 12 of 18 face-component units, against
Jacobi's 3 of 18 — a ~1.4x cut, not 6.2x, and it needs edge cases proved (a
cell at `(0, hi2, k)` reads `q(0, hi2+1, k, V)`, an EDGE). **That is a separate
investigation and is not attempted here.** Red-black becomes optimal by paying
nothing, not by gaining.

## What the job must return

| result | conclusion |
|---|---|
| red-black step within **±0.5 %** of `3c2903a` at every rank count | as designed: the regression was entirely the reordering |
| Jacobi 1-rank step now a **gain** (`rect`/`refined` ≥ +0.5 % vs `3c2903a`) | the reordering cost is gone and the divergence saving shows through |
| Jacobi 4- and 8-rank step **better than `95312d7`** by ~0.3–3 % | recovers what the prefix was costing `phi_exchange` |
| `phi_exchange` at 1 rank back to within ±2 % of `3c2903a` | the direct measurement of the artefact being removed |
| bit-exact vs **`95312d7`** (same values, different bookkeeping) | the divergence round is unchanged |
| bit-exact vs **`3c2903a`** on a RED-BLACK case | the enumeration really is restored, not merely similar |
| red-black still slower than `3c2903a` | the reordering was not the whole story — report it and look again |
| any Jacobi case slower than `95312d7` | the index list costs more than the prefix saved; report both and keep the prefix |

**Gate flags: PRODUCTION.** No expression moves; values are copied, never
recomputed. Both comparisons above are `max_abs 0` or the change is wrong.

## What this cannot say

- **Nothing about red-black's own reduced round.** Not attempted, and the
  ~1.4x ceiling above is an estimate from the read set, not a measurement.
- **Nothing new about 16 ranks for Jacobi.** The prefix already wins there and
  the reordering cost is swamped by the message saving; the index list should
  add a little, and "a little" is the prediction.
- **Nothing about the memory cost being free.** The point→index map is ~0.4x
  the size of `lPointEntry`, tens of MB on the production cases. Small against
  a 40 GB card, not zero.
