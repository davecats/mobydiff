# Block-overhead timing baseline

The measurement case that the block/2:1-interface optimisation work hangs on
(`docs/next_session_block_overhead.md`). It answers two questions on the
production boundary-layer grid:

1. **what does the equal-size block lattice cost** when it buys nothing
   (single level, `nb` set vs `nb` unset), and
2. **what does 2:1 refinement give back** (the same grid with the freestream at
   half resolution in x and z).

It is a *timing* case, not a physics case: 400 cold-start steps on a 138 M-cell
grid produce no meaningful flow. Nothing here should be read as a result.

## Layout

```
singleLevel/base_redblack.ini      production layout (nb unset), red-black SOR
singleLevel/nb16_redblack.ini      same grid, nb = 16 (33792 blocks)
singleLevel/base_jacobi.ini        production layout, damped Jacobi
singleLevel/nb16_jacobi.ini        same grid, nb = 16, damped Jacobi
multiLevel_xz/refined_yp100_jacobi.ini
                                   half-resolution freestream + wall band
                                   refined 2:1 in x,z (15360 leaves, 62.9 M cells)
run_overhead.sh                    runs them into runs/<name>/ (gitignored)
summarise.py                       the two rates, the ratios and the gates
```

Both solvers are timed on purpose. Red-black sweeps the open halo layer
**redundantly** (this is what makes results independent of `nb`); damped Jacobi
does not. Equal overheads under the two ⇒ the tax is not redundant arithmetic.

## Running

```bash
module load toolkits/nvhpc/25.9
./run_overhead.sh                    # all five, ~10-15 min each on an A6000
./summarise.py                       # re-print the table without re-running
```

`CUDA_VISIBLE_DEVICES` defaults to **1** (istmcetus GPU 1). corax and cetus GPU 0
belong to the production campaign — check `nvidia-smi` before launching, the
machines are shared. Absolute times on cetus are ~3.1× corax's; **only ratios
transfer between machines.**

`summarise.py` prints two independent rates per run — the solver's own `chron`
loop timer and the rate differenced out of `runtime.txt` (whose last column is a
*cumulative* average, so consecutive lines have to be differenced to get the
marginal cost and to drop the start-up transient). They should agree to a few
percent; a larger gap means the run was disturbed and the number is unusable.

## Reference measurements (A6000, 2026-08-06)

| | s/step | ratio |
|---|---|---|
| `nb` unset, redblack (production layout) | 1.456 | 1.000 |
| `nb = 16`, redblack | 2.052 | **1.409** |
| `nb` unset, jacobi | 1.240 | 1.000 |
| `nb = 16`, jacobi | 1.753 | **1.413** |
| 2:1 xz refined (y⁺ ~100), `nb = 16`, jacobi | 0.859 | 0.490 (cells 0.455) |

**The cost model.** 1.409 measured against `(18/16)³ = 1.4238` predicted: cost
tracks the *allocated* volume, so the overhead of any block shape is
`(1+2/nb_x)(1+2/nb_y)(1+2/nb_z)` and does not need re-measuring per case. The
lattice costs 41 %; the 2:1 interface costs 3–8 %.

**The free gate.** The block decomposition is result-invariant by design
(Phase 1: "results EXACTLY independent of nb and rank count"), and the runs
confirm it — `nb = 16` and `nb` unset produced identical runtime lines (jacobi
L2_div 1.73814896E-05 / Linf 1.00322464E+00; redblack 6.37658051E-06 /
1.00344109E+00). `summarise.py` checks this automatically. Any change to
blocking or storage must reproduce the same fields bit-for-bit.

## Interface placement

`refine_dims = xz` never subdivides y, so a 2:1 interface can only sit on a
**y-tile boundary** — every `nb_y` cells. With `nb = 16` on the production y line
(blayer, ny = 176, `resolved_height` 36, Δy⁺_wall 0.223) the choices are:

| y-tiles refined | y node | y | y⁺ | cells / production |
|---|---|---|---|---|
| 1 | 16 | 0.357 | 7.7 | 0.318 |
| 2 | 32 | 1.813 | 39.3 | 0.386 |
| **3** | **48** | **4.569** | **99.0** | **0.455** |
| 4 | 64 | 8.221 | 178.2 | 0.523 |
| 5 | 80 | 12.432 | 269.5 | 0.591 |
| 6 | 96 | 17.015 | 368.8 | 0.659 |
| 8 | 128 | 26.968 | 584.5 | 0.795 |
| 11 | 176 | 100.000 | 2167.5 | 1.000 |

(y⁺ at the developed c_f = 0.00464, i.e. u_τ = 0.0482, ν = 1/450.)

The shipped configuration is the **3-tile / y⁺ 99** row. This is the trade-off
per-direction `nb` has to weigh: a larger `nb_y` is cheaper but leaves fewer
legal interface heights.

## Provenance

Rebuilt 2026-08-07. The original directory was written on 2026-08-06 but never
committed and did not survive; the measurements above are the ones it recorded
(`docs/next_session_block_overhead.md`). The configs were reconstructed from
`production.ini` and pinned by three invariants the original reported, all of
which the rebuild reproduces exactly:

- the refined case's **15360 leaves** (verified with `build_cpu/leaftable_test`:
  6144 coarse + 9216 fine),
- its **0.455 cell ratio** (`34.60 M × (8/11 + 4·3/11) / 138.41 M = 0.4545`),
- **y⁺ ~100** at the interface (y node 48 → y⁺ 99.0 on the production y line).

The one thing the rebuild cannot reproduce is the original's absolute seconds,
which is why every comparison here is a ratio.
