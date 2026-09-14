# `jacobi_apply` is occupancy-limited, not wasteful — and my own inference was wrong

Job 5144931, hkn0403, A100-SXM4-40GB, `refined_yp82_rect_jacobi` at 1 rank
(60.56 Mcell), 12 profiled launches after two full steps. `jacobi_compute_phi` is
the **control**: the one kernel already known to sit near its limit, so the others
are read against it rather than against a spec sheet.

Context: `results_horeka_2026-09-14.md` measured `apply` at **33.1 %** of the
post-fix refined 16-rank step against `sweep`'s 10.3 %, both called 18 times.

| kernel | us | doubles/cell | min | vs min | ld sect/req | st sect/req | L2 hit % | **DRAM %peak** | SM %peak | **occupancy %** | **registers** |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `jacobi_compute_phi` | 3 204 | 5.48 | 5 | 1.10x | 2.49 | 9.00 | 51.5 | 53.3 | 47.0 | 29.1 | 94 |
| `jacobi_apply` k1 (`p += phi·idt`) | 1 543 | 3.21 | 3 | 1.07x | 1.90 | 9.00 | 49.9 | **64.8** | 66.1 | **45.3** | **59** |
| `jacobi_apply` k2 (face correction) | 7 654 | 10.82 | 8 | **1.35x** | 2.13 | 8.97 | 51.9 | **44.0** | 32.1 | **29.6** | **88** |

## What this refutes — mine

From the phase table alone I inferred that k2 "moves ~17.6 doubles/cell where
~10 suffice", i.e. wasted roughly half its traffic, and wrote that an
uncoalesced access pattern was the likely cause and possibly a several-fold
prize. **That is wrong.** k2 moves **10.82** doubles/cell against a minimum of 8
— 1.35x, and the two other kernels are at 1.07–1.10x. Load sectors/request is
1.9–2.5 and store 9.0 across all three, the *same* for the control, so there is
no coalescing defect to find. **There is no wasted traffic to reclaim.**

The inference failed because it assumed the kernels run at a fixed ~1.3 TB/s
effective bandwidth and converted time into bytes. They do not: DRAM utilisation
is 44–65 % of peak and varies by kernel, which is the entire finding.

## What it is instead

`apply` costs 2.87x `compute_phi` (9 197 us against 3 204) for two compounding
reasons, and neither is waste:

- **it genuinely moves 1.98x the bytes** — 3.21 + 10.82 = 14.0 doubles/cell
  against 5.48 — because correcting three velocity faces from `phi` and `mu` is
  simply more work than forming a divergence;
- **it runs at 0.83x the memory efficiency** (44.0 % of peak DRAM against 53.3 %).

1.98 × 1.21 = 2.40, which is k2/compute_phi = 7 654/3 204 = 2.39 exactly.

And the efficiency gap is **register-limited occupancy**. At 128 threads/block,
88 registers/thread allows 5 blocks/SM = 640 of 2 048 threads:

| kernel | registers | occupancy | DRAM %peak |
|---|---|---|---|
| `apply` k1 | **59** | **45.3 %** | **64.8** |
| `compute_phi` | 94 | 29.1 % | 53.3 |
| `apply` k2 | 88 | 29.6 % | 44.0 |

k1 — the same arrays, the same grid, the same block size, 29 fewer registers —
gets 1.5x the occupancy and 1.5x the DRAM utilisation. Local memory is zero in
all three, so there is headroom to trade registers for occupancy rather than for
spills.

## Two levers, measured rather than assumed

| lever | what it removes | estimated |
|---|---|---|
| **Cut k2's register count** 88 → ≤64 (8 blocks/SM, ~50 % occupancy) | lifts DRAM utilisation toward k1's 65 % | if k2 reaches 60 %: 7 654 → ~5 600 us, i.e. **~2.5 ms/step at 16 ranks, 7–8 % of the step** |
| **Fuse k1 into k2** | one launch and one duplicated `phi` read (k1 reads `phi`, k2 reads it again) | ~0.5 ms/step of traffic + 0.27 ms of launches ≈ **1.0 ms/step, 3 %** |

The register lever is the larger and the less obvious one. How to take it is an
open question — `-gpu=maxregcount` is global and risks spilling every other
kernel, so the candidates are restructuring k2 (its interface/outlet branch
variables are live across the whole body) or splitting it. **Whatever is tried
must be checked against `launch__registers_per_thread` and
`sm__warps_active`, not against a step time alone**: a change that cuts
registers by spilling to local memory will look like progress on paper and lose
on the clock.

## What this does not say

- **Not that k2 is badly written.** 1.35x minimum traffic with no coalescing
  defect is a competent memory-bound kernel; it is limited by how many warps the
  SM can keep in flight, not by how it accesses memory.
- **Not that 50 % occupancy is achievable.** 88 → 64 registers is a target, not a
  plan, and nothing here says the arithmetic fits in 64.
- **One rank, one configuration, 12 launches, one GPU.** The per-kernel
  efficiency is a per-rank property so it should carry, but it was measured at
  one point.
