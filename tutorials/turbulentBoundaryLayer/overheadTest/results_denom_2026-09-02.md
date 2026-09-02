# Hoisting the projection diagonal out of the Jacobi iteration

istmcetus, refined production case (`refined_yp82_rect_jacobi.ini`, 448 leaves,
60.56 M cells), 400 steps a side, A/B against `2c2d6bd`, both production flags,
one session, **zero foreign compute apps before and after every run**. Runtime
diagnostic lines identical between the binaries in every pairing.

## The observation

`jacobi_compute_phi` recomputed the projection diagonal every iteration. It is a
function of the face kinds, the grid metrics and `ibm%mu` only — and
`update_ibm_mu` runs once per RK substage, *before* `pressure_projection` — so
it is **identical across all `niter` iterations**. It was being formed 6 (or 12)
times per substage to produce the same numbers.

`compute_denom` now forms it once per substage into `pdenom`; the iteration
kernel reads it. That removes six branchy `face_grad_denom` calls, six `mu`
loads and ~11 fp64 operations per cell per iteration, which is the right lever
for a kernel `ncu` measured at **82 % SM throughput against 46 % DRAM**, 128
registers, 33 % occupancy — on a GPU whose fp64 rate is a small fraction of its
fp32 rate.

The division is untouched and the stored value is the old expression moved
verbatim, so the change is **bit-exact**, not a reciprocal-multiply trade.

## Results

| config | step ref | step cand | | `proj/sweep` |
|---|---|---|---|---|
| 1 GPU rank, niter 6 | 0.558443 | 0.539894 | **−3.32 %** | −26.4 % |
| 2 GPU ranks, niter 6 | 0.310020 | 0.301745 | **−2.67 %** | −25.8 % |
| **2 GPU ranks, niter 12** | 0.500011 | 0.474426 | **−5.12 %** | −26.3 % |

`proj/apply` is unchanged in every case (0.0910 → 0.0911, 0.1797 → 0.1797,
0.1823 → 0.1824), which is the control: only the sweep was touched.

**The sweep falls a consistent ~26 %.** What varies is how much of that survives
to the step, because `compute_denom` gives some back:

| | sweep saved | `compute_denom` cost | net |
|---|---|---|---|
| 1 rank, niter 6 | 0.032769 | 0.013943 | 0.018826 (3.37 %) |
| 2 ranks, niter 6 | 0.015517 | 0.007015 | 0.008502 (2.74 %) |
| 2 ranks, niter 12 | 0.032196 | 0.006881 | 0.025315 (5.06 %) |

## Why it is worth more at the setting production should use

The saving is **per iteration** (18 sweeps/step at `niter = 6`, 36 at 12); the
cost is **per substage** (3 `compute_denom` calls/step either way). So the gain
grows with `niter`: **−2.67 % at niter 6, −5.12 % at niter 12**, and 12 is the
setting the convergence measurement says damped Jacobi actually needs to keep
the pressure zero-mode away.

`compute_denom` is not cheap in absolute terms — it writes `pdenom` (484 MB per
call on this case) and reads `mu` (1.45 GB), traffic the old code never paid.
That is why ~45 % of the sweep saving is given back at `niter = 6` and only
~21 % at 12.

## Cost

One extra field-sized array (`pdenom`, same bounds as `phi`; ~544 MB on this
case). If memory ever becomes the binding constraint on a larger case, this is a
knob to reconsider — the alternative is recomputing, which is what it replaced.

## Gates

7-case suite bit-exact (`max_abs 0`, `-Mnofma` / `-gpu=nofma`) CPU **and** GPU;
`validation/block_nb` CPU uniform-flow gates still EXACT at two refinement
levels, plus 1 rank == 4 ranks.

## Not done

The red-black sweep recomputes the same diagonal per colour, and would benefit
from the same hoist — but it also needs the six *per-face* metrics individually
for its in-place corrections, not just their sum, so it is six arrays rather
than one. Not attempted; it is a different trade.
