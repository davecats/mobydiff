# HoreKa, 2026-10-06: the fresh profile of the current head (second cleaning pass)

Job **5184720** (`dev_accelerated`, 2 nodes hkn[0401,0403], 8 x A100-40GB,
36 min, nvhpc 25.3, `horeka/exchange/submit_profile2.sh`, head `b1f41be`).
ONE column, no comparison: six configs at 1/2/4/8 ranks + `min_channel`
(200 steps, phase buckets), the nsys per-kernel timeline at 1 rank on five
configs (`collect_kernels.py`, niter-aware since this job), ncu on every
kernel of one `refined_yp82_rect_jacobi` step (the `rough_jacobi` ncu leg
skipped past its 3 steps: `LAUNCH_SKIP=260` next time).

Raw tables: `horeka/results_job5184720/{phases,kernels,ncu}.md`. **The
write-up is section 7 of `docs/cleanup_review_2026-10-06.md`**; its headline:
the trip's `add_bodyforce_correction` streams the whole domain for all three
components every substage and is 8.5 % of the step on every boundary-layer
config, while the force is `v` only on the trip's block list -- the fill was
optimised on 2026-09-14, the correction was not. Step times reproduce the
2026-10-01 matrix (`base_jacobi` 75.83 ms at 8 ranks vs 75.93); block tax
1.006/1.028/1.023/0.994 at 1/2/4/8; red-black `zero_phi` is 9.3 % of that
config's device time.
