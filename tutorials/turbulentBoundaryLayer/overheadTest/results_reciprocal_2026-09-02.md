# Reciprocal diagonal: −4.3 % of the step, and an ncu lesson

istmcetus, refined production case, 2 GPU ranks, 400 steps a side, A/B against
`b90563b` (the diagonal hoist), both production flags, zero foreign compute apps
before and after every run.

## Result

| | ref | cand | | `proj/sweep` |
|---|---|---|---|---|
| niter 12 | 0.472928 | 0.452720 | **−4.27 %** | 0.088758 → 0.066232 (−25.4 %) |
| niter 6 | 0.301875 | 0.292633 | **−3.06 %** | 0.044298 → 0.033075 (−25.3 %) |

`compute_rdenom` costs +0.0016 s/step for the added `1/x` (once per substage).
`proj/apply` is unchanged, as expected.

Cumulative over the two projection changes, from the pre-hoist baseline:

| | before hoist | now | |
|---|---|---|---|
| niter 12 | 0.500011 | 0.452720 | **−9.46 %** |
| niter 6 | 0.310020 | 0.292633 | **−5.61 %** |

## The ncu projection was too optimistic, and the reason is methodological

`ncu` measured the division at **50.5 %** of `compute_phi` and I projected
≈ 9.6 % of the step. Production gives **−25.4 % on the sweep and −4.27 % on the
step** — a factor of two less.

**`ncu` locks the SM clock to 1.40 GHz.** That makes compute artificially
expensive relative to memory, so it *overstates* a pure compute reduction. The
earlier hoist's −26.2 % transferred to production exactly because that change
removed **traffic and compute together**; this one removes only compute, and it
does not transfer.

Rule for the next time: an ncu ratio transfers when the change moves bytes;
when it only moves instructions, treat the ncu figure as an upper bound and
measure before quoting.

## Numerical effect, measured

Not bit-exact by construction. Against the pre-change binary, 20 steps:

| case | `un` | `pn` | relative |
|---|---|---|---|
| min_channel | 9.2e-14 | 3.0e-12 | u 5.1e-15, p 1.6e-11 |
| les_ibm | 7.1e-14 | 1.4e-12 | |
| beltrami_y | 1.1e-14 | 4.4e-14 | |
| turb180 | 1.1e-14 | 4.4e-14 | (k 4.4e-16, nut 1.2e-15) |

Round-off throughout — velocities a few ulps, pressure larger because it
accumulates the increment over iterations and steps. No case differs
structurally. The runtime `L2_div` is **identical to printed precision** in both
A/B pairings (1.60527692E-05 at niter 12, 2.93355685E-05 at niter 6), so the
change is invisible at diagnostic resolution.

## Gates

Bit-exactness against a pre-change binary is deliberately given up for the
projection. Everything else holds, and exactly:

- uniform oblique flow through a 3-level patch: **EXACT** (0.0, `pn` spread 0.0),
  CPU and GPU, both refinement modes — uniform flow gives `div = 0`, so
  `phi = -omega*0*rdenom` is exactly zero either way;
- **1 rank == 4 ranks EXACT**, `validation/block_nb` and
  `validation/redblack_interface`;
- **CPU == GPU EXACT** (`max_abs 0`, nofma) on the refined min_channel;
- red-black + 2:1 gates pass.

Going forward the 7-case suite works again as a refactor gate: it builds both
sides from the current tree, so the reference binaries were regenerated at this
commit and future changes baseline against it.
