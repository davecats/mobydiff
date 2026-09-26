# The campaign on the consolidated head

Job 5163909, 2026-09-26, **dev_accelerated, 2 nodes** (hkn[0401,0403]),
39 minutes, 44 runs, no failures, 200 steps, `--map-by numa --bind-to core`.
Raw in `horeka/results_job5163909/`.

Two columns in one allocation:

| column | commit | what it is |
|---|---|---|
| `ref` | `8fa0fc2` | the head `results_horeka_2026-09-25.md` measured |
| `new` | `2a0ed47` | the consolidated head: the `boundaryLayer` and `scalar` merges, and the jacobi-interface feature port |

**Ranks 1/2/4/8 only.** `accelerated` was backlogged to 09-28/09-29, so the
sweep ran on `dev_accelerated` (2 nodes = 8 GPUs); **16 ranks is submitted
separately as job 5163915 into the same results directory** and is NOT in the
tables below. Every 16-rank figure in the 09-25 report therefore still stands
un-rechecked.

## THE REF COLUMN RUNS `convection = skew`, and that is the whole design

The shipped configs used to carry `[flow] convection = skew`. The port's S3
lockdown made that key an error, so it was stripped from them — which means the
two columns would otherwise have compared skew-symmetric convection against
DIVERGENCE-form convection and called the difference "the consolidation".
`config.f90` has no `case default`, so neither binary would have complained: the
new one rejects the key it no longer knows, the old one silently ignores its
absence. The ref column is therefore run from a staged config copy with the key
appended, which is what `8fa0fc2` was published with.

**So the numbers below do NOT include the cost of hardwiring skew.** Both
columns compute the same convection. What they measure is everything else the
consolidation brought.

## Finding 1 — the consolidation costs 1.3–2.0 % of the step, everywhere

| config | n=1 | n=2 | n=4 | n=8 |
|---|---|---|---|---|
| `base_jacobi` | −2.0 % | −1.9 % | −1.7 % | −1.7 % |
| `rect_jacobi` | −1.8 % | −1.5 % | −1.5 % | −1.4 % |
| `refined_yp82_rect_jacobi` | −1.6 % | −1.3 % | −1.4 % | −1.5 % |
| `refined_yp82_rect_redblack` | −1.8 % | −1.9 % | −1.7 % | −1.6 % |
| `refined_big_rect_jacobi` | — | — | −1.4 % | −1.4 % |

A loss, on every config and every rank count, in a band narrow enough
(1.3–2.0 %) to be one cause rather than several, and **flat in rank count**, so
per-cell or per-step work rather than communication.

**LOCALISED to the momentum predictor** by the phase tables, which both columns
carry (`[output] profile` is on for campaign runs). `rect_jacobi` at 4 ranks:

| bucket | ref | new | delta |
|---|---|---|---|
| `momentum` | 38.721 ms | 41.366 ms | **+2.644 ms** |
| `projection` | 102.911 ms | 102.568 ms | −0.343 ms |
| `bodyforce` | 0.626 ms | 0.548 ms | −0.078 ms |
| everything else | | | under ±0.01 ms |
| **total_measured** | 145.884 ms | 148.122 ms | **+2.238 ms** |

So it is the fused predictor kernel and nothing else — the projection actually
got marginally faster and partly offsets it. **My first explanation was wrong
and is recorded here so it is not repeated**: I attributed this to the scalar
merge putting dormant code on the step path. It is not that. The `scalar` bucket
reads 0.0007 ms, and the scalar merge does not touch the predictor kernel at
all.

Inside that kernel the ONLY change between the two commits is the S3 lockdown
removing the three `if (skew)` guards. **Both columns execute the skew
corrections** (the ref config carries the key — verified in the staged ini and
in the run's own `config.ini`), so the arithmetic is identical and the guard was
taken on both sides.

**MECHANISM, from ncu (job 5163917, `run_ncu.sh KERNEL_RE=step_momentum`, both
binaries on one node).** The predictor kernel:

| | ref (`8fa0fc2`) | new (`b624d16`) |
|---|---|---|
| time | 15,226 us | **16,938 us (+11.2 %)** |
| doubles/cell | 16.87 | 16.87 |
| vs source-counted minimum | 1.05x | 1.05x |
| ld / st sectors per request | 3.62 / 8.67 | 3.50 / 8.67 |
| L2 hit | 68.14 % | 68.15 % |
| **occupancy** | 23.76 % | **23.69 %** |
| DRAM % of peak | 34.51 | **31.02** |
| SM % of peak | 40.33 | **34.09** |
| registers | 128 | **100** |

Read it in this order:

1. **Traffic is byte-for-byte identical** — same doubles/cell, same ratio to the
   source-counted minimum, same sectors per request, same L2 hit rate. It is not
   memory, and it is not extra work.
2. **Occupancy did not move** (23.76 → 23.69) even though registers fell 128 →
   100. So the 28 freed registers bought NOTHING: this kernel's occupancy is
   capped by something else (grid/block shape under `collapse(4)`), not by
   registers. **That retracts the inference in the first version of this
   section**, which read "registers down, therefore occupancy up" — occupancy
   was never measured there, only assumed, and it was wrong.
3. **Both utilisation axes FELL** — DRAM 34.51 → 31.02 % of peak, SM 40.33 →
   34.09 % — while the time rose. Same bytes, same occupancy, less issued per
   unit time: the kernel got MORE LATENCY-BOUND.

So the compiler spent the freed registers badly. Removing the three `if (skew)`
guards let it cut register pressure by 28, which gained no occupancy because
occupancy is capped elsewhere, and the shorter live ranges lengthened dependency
chains — which a kernel running 15 of 48 warps per SM cannot hide. The
predictor was ALREADY the solver's occupancy-limited kernel; it had no slack to
absorb this.

**The cheap fix, not yet taken:** put the corrections back inside a condition
the compiler cannot fold (a mapped logical that is always true), which restores
the ref binary's 128-register schedule at the cost of one never-taken branch.
That is a one-line change worth ~1.5 % of the step, and it needs its own
before/after ncu plus the nofma suite — an optimisation increment, not a
footnote to this campaign, so it was left for a session that can gate it.
Recommended AGAINST: `-gpu=maxregcount`, which CLAUDE.md records as a TARGET
rather than a cap — it raises every other kernel toward it.

## Finding 2 — the cross-config ratios barely moved

These are the numbers the campaign exists for, and they are ratios BETWEEN
configs at one head, so they are what goes stale when the step time moves.
`new` column, against the 09-25 publication:

| quantity | 09-25 (to 16 ranks) | 09-26 `new` (to 8 ranks) |
|---|---|---|
| block tax, 1/2/4/8 | 1.015 / 1.034 / 1.024 / 1.000 | **1.010 / 1.032 / 1.025 / 0.996** |
| 2:1 machinery, coarse-cell-equivalents | 1.001x at 4 ranks | **0.997x at 4, 0.979x at 8** |
| red-black / Jacobi, 1/2/4/8 | — (0.808 at 16) | **0.883 / 0.894 / 0.905 / 0.937** |
| strong scaling at 8 ranks | — | base 90 %, rect 92 %, refined 80 %, red-black 76 %, big 96 % |

So the consolidation's 1.5 % is nearly common-mode: it moves every config by
about the same fraction, which is why the ratios survive it. The block tax is
**within 0.005 of its published value at every rank count**, and the 2:1
machinery still costs what the coarse cells beside it cost.

Red-black is still the faster smoother at these rank counts (0.88–0.94 of
Jacobi's step) and closing as ranks rise, which is the same trend that reached
parity at 16 ranks on 09-17. Whether it still does is job 5163915's question.

## What this does not answer

* **16 ranks.** Job 5163915. Until it lands, every inter-node figure in the
  09-25 report is the most recent one there is.
* **The cost of skew itself.** Deliberately excluded (see above). Measuring it
  needs a pre-lockdown binary run BOTH ways, which is a different job and a
  different question — and the branch's own figure (~+3.4 % s/step) was measured
  on another case and machine, so it should not be carried over.
* Nothing on the predictor regression — Finding 1 now carries the ncu
  measurement and the mechanism. What remains is to TAKE the fix and gate it.
