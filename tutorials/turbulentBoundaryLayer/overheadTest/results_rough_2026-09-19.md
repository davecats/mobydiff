# A rough boundary layer, and the body-case number that was missing

Job 5151978, A100-SXM4-40GB, HoreKa `dev_accelerated`, **200 steps**, before/
after in one allocation at 4 and 8 ranks. `new` = `7b2bc2a`, `ref` =
`03305e4` (branch `bench/rdenom-always`). Raw runs in
`horeka/results_job5151978/`.

`results_rdenom_static_2026-09-17.md` §6b ended by admitting that the body-case
magnitude was not measured: `les_ibm` is too small and too wall-dominated, and
the sailplane diverges past its smoke-test length. **This closes that**, on a
case shaped like the production target.

## 1 — The case

`configs/rough_jacobi.ini` is `rect_jacobi`'s grid, blocks and flow exactly
(4096×176×192, `nb = 64 44 48`, `lx = 750`), plus a 3D sinusoidal wall:
`h·sin(kx x + φx)·sin(kz z + φz)`, the egg-carton roughness of MacDonald, Chung,
Hutchins, Ooi & Sandberg (JFM 2017), at `[ibm] wall_shape = eggcarton`.

It exists because **every other config here is body-free**, so none of them can
measure anything that depends on the immersed boundary — and the `rdenom`
narrowing is exactly that.

**It is a BENCHMARK, not a validated physics case.** The height and wavelengths
are chosen to be resolved by this grid and to sit inside the first y-block, not
to reproduce a published Re_τ or k_s⁺. Nothing physical should be quoted from it.

## 2 — The measurement

| case | step ref → new | | `setup` ref → new | |
|---|---|---|---|---|
| `rect_jacobi` 4 ranks (body-free control) | 151.617 → 146.265 ms | **−3.53 %** | 6.808 → 1.225 | −82.0 % |
| `rect_jacobi` 8 ranks | 79.279 → 76.497 | **−3.51 %** | 3.417 → 0.615 | −82.0 % |
| **`rough_jacobi` 4 ranks** | 159.927 → 155.990 | **−2.46 %** | 6.801 → 2.689 | **−60.5 %** |
| **`rough_jacobi` 8 ranks** | 83.278 → 81.345 | **−2.32 %** | 3.402 → 1.382 | **−59.4 %** |

**Body blocks: 64 of 256 per rank — 25 %**, exactly the fraction the scaled twin
predicted, because the roughness lies entirely within the bottom y-block.

## 3 — The prediction was quantitative and it held

The claim was that the saving scales as `1 − (body blocks / all blocks)`:

| | predicted | measured |
|---|---|---|
| `setup` drop, body case | 0.75 × 82.0 % = **61.5 %** | **60.5 % / 59.4 %** |
| absolute saving, rough ÷ rect | **0.75** | **0.736** |

Both within about a percentage point, from the block fraction alone. The model
is not a story fitted afterwards — it was written down in
`PREREGISTERED_rdenom_static.md` and the block count was measured before this
job ran.

## 4 — The one-time residual, now measured rather than inferred

`results_rdenom_static_2026-09-17.md` §3 argued that what remains in `setup`
after the change is **first-call allocation charged to every step by a short
benchmark**, and inferred it from a CPU amortisation test. Two GPU runs of
different length now settle it directly (`rect_jacobi`, 8 ranks):

| run length | `setup` after | implied one-time total |
|---|---|---|
| 100 steps (job 5150041) | 1.218 ms/step | **121.8 ms** |
| 200 steps (this job) | 0.615 ms/step | **123.0 ms** |

Constant to 1 %. The residual is fixed cost, not per-step work, so the fraction
keeps falling with run length: the same change reads **−67 %** of the bucket at
100 steps and **−82 %** at 200. Production runs are thousands of steps.

## 5 — Where the `rdenom` work now stands

| case shape | body blocks | step, 200 steps |
|---|---|---|
| body-free (boundary layers, channels) | 0 % | **−3.5 %** |
| rough wall, roughness in one y-block | 25 % | **−2.4 %** |
| `les_ibm`, walls spanning the domain, small case | 40 % | −0.18 % (its `setup` is only 2.5 % of a much heavier step) |
| compact body (`sailplane`, `nb = 10`) | 1.1 % | not measured — would be near the body-free figure |

All of it bit-exact: the 7-case suite at `max_abs 0` including every RANS scalar,
and `wavychannel`/`les_ibm`/`min_channel` against the pre-roughness commit.

## 6 — What this does not support

- **No physics.** See §1. The roughness is not calibrated and the case has no
  validated statistics.
- **Nothing about passive scalars**, which do not exist in the solver. A
  rough-wall scalar-transport benchmark — the MacDonald et al. forced-convection
  configuration — needs generic scalar transport first: a module, a fused RK3
  kernel, BCs, IBM Dirichlet coupling, io/restart and validation. Its own track.
- **Nothing about 16 ranks**; `dev_accelerated` is 2 nodes. `setup` is a per-cell
  kernel and its share is flat across 4 and 8 ranks here, so rank-independence
  is an expectation, not a measurement.
- **The compact-body end of the range is still unmeasured.** `sailplane`
  classifies 48/4500 but cannot be driven far enough to time; an airfoil case
  would settle it.
