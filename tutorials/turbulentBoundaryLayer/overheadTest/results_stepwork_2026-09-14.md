# Work that did not need doing — and where the register lever stops

Jobs 5145541, 5145549 and 5145554, A100-SXM4-40GB, HoreKa, `dev_accelerated`,
before/after inside each allocation, 100 steps. Raw runs in
`horeka/results_job5145541/`, `results_job5145549/`, `results_job5145554/`.
Readings pre-registered in `horeka/exchange/PREREGISTERED_bodyforce.md`.

Follows the two register increments of the same day
(`results_apply_registers_2026-09-14.md`, `results_rdenom_registers_2026-09-14.md`).
Those made kernels faster. **These three mostly make kernels not run**, which is
why the total is larger than either.

**Cumulative, from `cd0f279` to `03f0e1b`:**

| case | step |
|---|---|
| `rect_jacobi` 4×1 | 172.63 → 158.75 ms **−8.04 %** |
| `rect_jacobi` 8×2 | 91.19 → 83.88 **−8.02 %** |
| `refined_yp82` 4×1 | 81.67 → 75.00 **−8.16 %** |
| `refined_yp82` 8×2 | 46.48 → 42.65 **−8.23 %** |

Every gate below is `max_abs 0` at **production flags** — all three changes leave
the surviving arithmetic textually identical, so by the rule set in
`run_mapgate.sh` the production comparison is the tighter one, and it passed.
Nine comparisons per job: Pass G on both production cases (138 M and 60 M points)
and the 7-case suite including every RANS scalar.

## 1 — `update_ibm_mu` on a case with no body (job 5145541)

It ran an fp64 **divide per ghost-inclusive cell × 3 components × 3 substages** to
compute `mu = 1/(1+dt·0) = 1` — which `init_ibm` had already written. One device
reduction on the first call now decides whether any coefficient is non-zero **on
this rank**, and if none is, the kernel never runs again.

| case | `ibm_mu` |
|---|---|
| `rect_jacobi` 4×1 | 8.329 → 0.036 ms **−99.6 %** |
| `rect_jacobi` 8×2 | 4.149 → 0.019 **−99.6 %** |
| `refined_yp82` 4×1 | 3.635 → 0.017 **−99.5 %** |
| `refined_yp82` 8×2 | 1.857 → 0.009 **−99.5 %** |

The residual is the profiler bracket. The answer is local and correctly so: `mu`
is pointwise, so a rank with no body has `mu = 1` whatever other ranks hold, and
the fields are bit-identical either way. `les_ibm` gates that a case WITH a body
still computes it.

## 2 — The trip force, in two steps, the first of which the pre-registration called

`fill_trip_kernel` refreshed the WHOLE domain every substage: it wrote
`f_u = f_w = 0` (always zero for a trip) and evaluated `f_v` everywhere, although
its own `ex < -50` cutoff confines the envelope to |x−x₀| < 28.3, |y| < 7.07 in a
750×100 domain.

**First half** (job 5145541): write only `f_v`, and only over the blocks the
envelope can reach — a list built once at init from the same cutoff evaluated at
each block's closest point, so skipping a block writes back the zeros it already
holds.

| case | `bodyforce` |
|---|---|
| `rect_jacobi` 4×1 | 5.497 → 2.279 ms **−58.5 %** |
| `refined_yp82` 8×2 | 2.242 → 1.791 **−20.1 %** |

**This is the second branch of the pre-registration, not the first**, which said:
*"`bodyforce` barely moves ⇒ the selected-block fraction is much larger than the
scaled-down twin's 8/64 — report the fraction and reconsider, do not quote a
win."* It is reported: the solver prints **40 of 256 blocks** at 4 ranks and 40 of
128 at 8 (26/112 and 26/56 refined). The trip sits at low x, so the rank that
owns it owns a sixth to a half of the work. And those blocks are exactly the ones
paying the **24-mode Fourier sum — 48 transcendentals per cell** — which removing
stores does nothing about.

**Second half** (job 5145549): `g_k(z)` and `g_{k+1}(z)` depend on z alone, so
they are tabulated per (k, listed block); and they change only when the random
walk redraws, once per `trip_ts` ≈ hundreds of steps, so the table is refreshed
then, not per substage. Only the blend `(1−b)g_k + b g_{k+1}` stays per cell,
because b(t) does change.

| case | `bodyforce` | total from the start |
|---|---|---|
| `rect_jacobi` 4×1 | 2.272 → 0.626 ms **−72.4 %** | 5.497 → 0.626 **−88.6 %** |
| `rect_jacobi` 8×2 | 2.246 → 0.621 **−72.3 %** | 3.562 → 0.621 **−82.6 %** |
| `refined_yp82` 4×1 | 1.800 → 0.459 **−74.5 %** | 3.063 → 0.459 **−85.0 %** |
| `refined_yp82` 8×2 | 1.791 → 0.461 **−74.3 %** | 2.242 → 0.461 **−79.4 %** |

Gated additionally with `trip_ts = 0.005`, so the walk advances four times in 20
steps — the caching path a normal run does not reach for 4000 steps.

## 3 — `step_momentum`: the measurement is the result

**First ever ncu on the solver's largest kernel** (job 5145549, `refined_yp82` at
1 rank, 60.56 Mcell):

| kernel | us | doubles/cell | min | vs min | DRAM %peak | occupancy % | registers |
|---|---|---|---|---|---|---|---|
| `step_momentum` predictor | 16 168 | 16.87 | 16 | **1.05x** | **32.50** | **23.96** | **128** |
| `step_momentum` qs→q copy | 2 800 | 6.67 | 6 | 1.11x | **74.22** | **41.09** | **66** |

It is occupancy-limited exactly as `jacobi_apply` k2 and `compute_rdenom` were,
and more starkly: its own sibling — same grid, same block size, 62 fewer
registers — gets 1.7x the occupancy and **2.3x the DRAM utilisation**. Traffic is
1.05x its source-counted minimum, so again there is nothing to reclaim.

**But the lever that worked twice today does not work here.**

- **Packing the nine Laplacian coefficient arrays into three changed the register
  count by nothing: 128 before, 128 after.** For this kernel the array bases are
  not the binding constraint — the ~30 live values of three fused component
  blocks are. That is a real limit on the pattern, and it is the reason this
  section exists.
- **A register cap reaches 102 only by spilling.** `-gpu=maxregcount:102` puts
  `step_momentum` at 102 with `STACK:16`, and applied globally it *raises* every
  kernel toward the cap — `jacobi_apply` k2 **64 → 96**, `compute_rdenom`
  80 → 100, `compute_phi` 80 → 96, `interface_correct` 54 → 80 — undoing the
  day's projection work, while `add_eddy_viscosity_correction` spills 32 bytes
  and `redblack_sweep` 8. The handout's L3 warning is confirmed and sharpened:
  the flag behaves as a *target*, not only a cap.

So the register line stops at `step_momentum`. What the pack did buy, measured
rather than assumed (job 5145554) because a change to the most
bit-exactness-critical kernel in the solver should not stay on a hypothesis:

| case | `momentum` | step |
|---|---|---|
| `rect_jacobi` 4×1 | 39.75 → 38.62 ms **−2.82 %** | −0.65 % |
| `rect_jacobi` 8×2 | 19.88 → 19.27 **−3.05 %** | −0.81 % |
| `refined_yp82` 4×1 | 17.36 → 16.84 **−3.00 %** | −0.50 % |
| `refined_yp82` 8×2 | 8.773 → 8.490 **−3.23 %** | −0.65 % |

and ncu says where it came from: 16 168 → 15 221 us (−5.9 %), DRAM
32.50 → 34.52 %, **registers and occupancy unchanged** (128, 23.8 %). It is the
contiguity alone — the three stencil coefficients for one index now share a
sector — which is exactly the hypothesis the job was run to test. Kept.

## 4 — What these numbers do NOT support

- **These are not 2:1-refinement costs.** Both twins shrink, so the block tax,
  strong-scaling and coarse-cell-equivalent ratios all move. They must be
  **re-measured, not rescaled**, before any of them is quoted again.
- **The trip result is case-specific in scope but not in kind.** It applies to the
  `boundaryLayer` case, which is what this whole campaign benchmarks and how the
  production boundary layer is tripped. Channels and airfoils have no trip.
- **`ibm_mu`'s win does not exist on cases with a body.** There the kernel still
  runs in full; the only thing added is one device reduction, once.
- **Nothing at 16 ranks.** All three jobs are 4 and 8 ranks on `dev_accelerated`.
  The two buckets' shares at 8 ranks are the basis for any extrapolation.
- **`step_momentum` is not fixed.** It is still 128 registers at 24 % occupancy
  and 34.5 % of peak DRAM — the least efficient large kernel in the solver — and
  this report says only that the two cheap levers do not move it.

## 5 — What is next

1. **`step_momentum` needs a different idea, not another register trick.** The
   one not yet tried is splitting the fused predictor into three per-component
   kernels: far fewer live values each, at the cost of re-reading the shared
   velocity components. At 34.5 % of peak DRAM there is headroom for more
   traffic, but it is a real refactor of the kernel whose bit-exactness underpins
   every phase gate in this project, and it should not be started without an
   allocation to measure it in.
2. `compute_rdenom`'s fp64 divide — now compute-limited (SM 49.6 % > DRAM
   39.8 %), the remaining ~3 % bucket.
3. The exchanges, ~21–26 % of the step, with `mpi_wait` inside — closed for
   overlap and partitioning.
4. **Re-run the campaign matrix.** Roughly 15 % has come off the step today and
   every published ratio in `results_horeka_2026-09-14.md` has a moved
   denominator.
