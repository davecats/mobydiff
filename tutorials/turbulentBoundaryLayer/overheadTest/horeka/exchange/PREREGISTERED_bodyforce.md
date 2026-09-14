# Pre-registered readings — the trip force and the body-free `mu` refresh

Written **before** job `moby_bodyforce` ran. Unlike the two register increments,
this one is not a rate improvement: it is **work that was being done and did not
need to be**, so the prediction is sharper and a miss is more damaging.

## What is being removed

**1. `fill_trip_kernel` ran over the whole domain every substage.** It wrote
`f_u = f_w = 0` (always zero for a trip, and already zero from `bf%f = 0.0d0` at
allocation) and evaluated `f_v` everywhere, while its own `ex < -50` cutoff
confines the envelope to |x−x₀| < lx√50 = 28.3 and |y| < ly√50 = 7.07 in a
750×100 domain. Now: the v component only, over a block list built once at init
from the same cutoff evaluated at each block's closest point. On a 256-cell
scaled-down twin, 8 of 64 blocks are selected.

**2. `update_ibm_mu` ran an fp64 divide per ghost-inclusive cell × 3 components ×
3 substages to compute `mu ≡ 1`** on cases with no immersed body. Now: one
device reduction on the first call decides whether any coefficient is non-zero on
this rank, and if none is, the kernel is skipped for the rest of the run —
`init_ibm` already wrote `mu = 1.0`, and `1/(1+dt·0)` is exactly 1.0.

## The buckets being targeted (job 5145507, `after` column — the current code)

| case | `bodyforce` | `ibm_mu` | sum |
|---|---|---|---|
| `rect_jacobi` 4×1 | 5.512 ms (3.2 %) | 8.347 ms (4.8 %) | **8.0 %** |
| `rect_jacobi` 8×2 | 3.584 (3.9 %) | 4.171 (4.6 %) | **8.5 %** |
| `refined_yp82` 4×1 | 3.079 (3.8 %) | 3.646 (4.5 %) | **8.2 %** |
| `refined_yp82` 8×2 | 2.237 (4.8 %) | 1.850 (4.0 %) | **8.8 %** |

## What the job must return

| result | conclusion |
|---|---|
| `ibm_mu` → ~0 (a few µs of bracket), `bodyforce` → the selected-block fraction, step **−7 to −9 %** | as designed; take it |
| `ibm_mu` → ~0 but `bodyforce` barely moves | the selected-block fraction is much larger than the scaled-down twin's 8/64, i.e. the envelope test is too loose for the STRETCHED full grid — report the fraction and reconsider, do not quote a win |
| either bucket unchanged | the guard is not firing; the case is not the one being measured, or the flag is wrong. Investigate before anything else. |
| the step falls by materially MORE than the two buckets | something else changed; distrust it and look for it |
| any nonzero `max_abs` | a real difference. Pass G runs the actual tripped boundary-layer cases and `les_ibm` the actual IBM path, so either would catch its own half. |

**Gate flags: PRODUCTION, not nofma** — and that is a deliberate call under the
rule written in `run_mapgate.sh` two increments ago. Neither change MOVES an
expression: `fill_trip_kernel`'s surviving arithmetic is textually identical and
`update_ibm_mu`'s kernel is untouched (skipped, not rewritten). So both sides
contract identically and the production-flag comparison is the tighter test. If
it fails at 1e-16 the rule was wrong here and the nofma gate is the fallback —
that outcome is worth recording either way.

## What this job cannot say

- **Nothing about 16 ranks** (`accelerated` queue); the two buckets' shares at 8
  ranks are the basis for any extrapolation, flagged as such.
- **Nothing about cases WITH a body.** There, `update_ibm_mu` still runs in full
  and the only cost added is one device reduction, once. `les_ibm` gates that it
  is still correct, not that it is still fast.
- **These are not 2:1 costs.** Both twins shrink, so the block-tax,
  strong-scaling and coarse-cell-equivalent ratios all move and must be
  re-measured rather than rescaled.
