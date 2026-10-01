# HoreKa, 2026-09-30 / 10-01: the drivers under the case-file contract, the predictor register guard, and a node-boundary finding

Three jobs on `dev_accelerated` (2 nodes, 8 x A100-40GB, nvhpc 25.3), all
200 steps, `--map-by numa --bind-to core` unless stated:

| job | what | binaries |
|---|---|---|
| 5172140 | `submit_matrix5.sh`, `CONFIGS="base_jacobi rect_jacobi"`, ranks 1/2/4/8 | new `0811811` (step 8 head), ref `8fa0fc2` (pre-step-7, builds inline) |
| 5172166 | `submit_predictor.sh` (PREREGISTERED_predictor.md) | new `1cf8d58` (guard), ref `0811811` |
| 5173620 | rank-placement test, `base_jacobi` at 8 ranks | `0811811` |

## 1. The drivers prepare the case file (handout item 2) — PASS

Every `new` run directory holds `overhead.case.h5` and `prepare.log` beside
`run.log`; no `ref` directory does, which is correct (8fa0fc2 predates the
contract and builds inline). No `PREPARE FAILED`, no `run.FAILED.log`.

| config | n | leaves (prepare.log) | expected |
|---|---|---|---|
| `base_jacobi` (nb unset) | 1 / 2 / 4 / 8 | 1 / 2 / 4 / 8 | one block per rank: layouts 1 1 1, 2 1 1, 2 2 1, 2 2 2 |
| `rect_jacobi` (nb 64 44 48) | 1 / 2 / 4 / 8 | 1024 at every count | 64 x 4 x 4 |

Step time against the 2026-09-26 matrix (`new` = `ae37ebf` there), s/step:

| config | n | 09-26 | today | delta |
|---|---|---|---|---|
| `base_jacobi` | 1 | 0.56280 | 0.56224 | −0.1 % |
| | 2 | 0.27975 | 0.28021 | +0.2 % |
| | 4 | 0.14447 | 0.14422 | −0.2 % |
| | **8** | **0.07775** | **0.13200** | **+69.8 %** (section 3) |
| `rect_jacobi` | 1 | 0.56844 | 0.56803 | −0.1 % |
| | 2 | 0.28866 | 0.28830 | −0.1 % |
| | 4 | 0.14812 | 0.14789 | −0.2 % |
| | 8 | 0.07746 | 0.07717 | −0.4 % |

The `ref` column of the same job reproduces its 09-26 values to 0.1–0.5 %
(`base_jacobi` n=8: 0.07662 vs 0.07642), so the allocation is comparable.
Seven of eight runs sit inside the matrix's noise. Everything that landed
between the two heads (step 6 boundary rows, step 7 prepare, F5 `apply_bc`
once per projection, F1) is therefore invisible on these body-free configs
at 1–8 ranks, with one exception.

## 2. The predictor register guard — every pre-registered number hit

ncu, `refined_yp82_rect_jacobi`, both binaries on one node, the fused
predictor kernel (and its sibling as the control):

| | ref (no guard) | new (guard) | predicted |
|---|---|---|---|
| registers | 100 | **128** | 100 → 128 |
| occupancy | 23.70 % | 23.77 % | unchanged within 0.1 |
| doubles/cell | 16.84 | 16.84 | identical |
| ld / st sectors per request | 3.50 / 8.67 | 3.62 / 8.67 | identical |
| DRAM % of peak | 30.98 | **34.47** | 31 → ~34.5 |
| SM % of peak | 34.08 | **40.33** | 34 → ~40 |
| kernel time | 16,937 us | **15,221 us (−10.1 %)** | −9 to −11 % |
| control kernel (66 regs) | 2,757 us | 2,757 us | unmoved |

The new column IS the 8fa0fc2 row of `results_horeka_2026-09-26.md`
(15,226 us, 34.51, 40.33, 128): the guard restores that schedule exactly.

The step, both columns in one allocation:

| run | step ref → new | | momentum bucket |
|---|---|---|---|
| `base_jacobi` n4 | 144.341 → 141.522 ms | −1.95 % | 43.71 → 40.60 (−7.1 %) |
| `base_jacobi` n8 | 132.012 → 130.168 | −1.40 % | 21.54 → 19.98 (−7.2 %) |
| `rect_jacobi` n4 | 148.142 → 145.695 | −1.65 % | 41.44 → 38.75 (−6.5 %) |
| `rect_jacobi` n8 | 77.190 → 76.059 | −1.47 % | 20.71 → 19.37 (−6.5 %) |
| `refined_yp82` n4 | 70.012 → 68.980 | −1.47 % | 18.06 → 16.90 (−6.4 %) |
| `refined_yp82` n8 | 39.023 → 38.432 | −1.51 % | 9.07 → 8.48 (−6.5 %) |

Predicted: step −1.3 to −2.0 %, momentum −6 to −7 %, flat in rank count.
Controls on `rect` n4: projection 101.928 → 102.174 ms (+0.24 %), `apply`
63.336 → 63.495, `sweep` 24.762 → 24.837, body force 0.550 → 0.550. The
1.3–2.0 % the consolidation cost is recovered, by one never-taken branch.

Equivalence was gated on the workstation before submission (nofma, against
`~/step8_ref_binaries` = `0811811`): 7-case suite and 9-case scalar suite,
CPU and GPU, every dataset `max_abs 0`.

Compiler note: nvhpc 25.9 on the workstation (cc86) folds harder than the
cluster's 25.3 (cc80) — the unguarded kernel is 78 registers there, 100
here; both go to 128 with the guard, `STACK 0`.

## 3. FINDING: an nb-less case puts a different face on the node boundary since step 7

`base_jacobi` at 8 ranks: 132.0 ms/step, reproduced three times in two
jobs, against 77.7 on 09-26 and 76.6 for 8fa0fc2 in the same allocation.
Momentum is unchanged; the whole difference is exchange wait:

| | 8fa0fc2 (inline rank box) | `0811811` (case file) |
|---|---|---|
| MPI layout | 2 2 2 | 2 2 2 |
| `mpi_wait` | 3.26 ms/step | **58.8 ms/step** |
| wait balance (max/min over ranks) | 1.70 | 1.04 (every rank waits equally: transfer, not skew) |
| projection | 52.6 ms | 95.0 ms |

**Mechanism.** The old rank box gave rank r the box at its Cartesian
coordinates (rank = 4x + 2y + z, x slowest). Since step 7 an nb-less case
is one block per rank numbered along the Morton curve and split linearly,
and `morton_key` puts x in the LOWEST bit (id = x + 2y + 4z; block table
dumped from the case file: ids 0–3 at z = 0, 4–7 at z = 96). With ranks
filled node by node, the node boundary used to cut x and now cuts z. On
this domain (4096 x 176 x 192, blocks 2048 x 88 x 96) an x-face is
88 x 96 = 8,448 cells and a z-face 2048 x 88 = 180,224 — and z is
periodic, so each rank reaches its off-node z neighbour through both
faces: ~43x the cross-node cells per exchange.

**Decisive test (job 5173620):** same binary, same case file, only the
rank placement changed.

| placement | node boundary cuts | step | `mpi_wait` |
|---|---|---|---|
| `--map-by numa` (ranks 0–3 / 4–7) | z | 131.99 ms | 58.85 ms |
| `--map-by node` (even / odd) | x | **77.58 ms** | **3.43 ms** |
| `--map-by numa`, repeated | z | 132.02 ms | 58.88 ms |

**Consequences.** The fields are untouched (results are independent of
ownership by construction). What changed is a performance property of
nb-less multi-node runs only: every config with an explicit `nb` already
used the Morton split and is unaffected (`rect_jacobi` n=8: −0.4 %). In the
campaign, `base_jacobi` is the unblocked reference of the block tax, so
**the block tax at 8 ranks and above may not be formed from a post-step-7
`base_jacobi` run with the default placement** (0.0772/0.1320 would read
0.58). Until the ordering question below is decided, run `base_jacobi`
with `MPIRUN_EXTRA`-level placement that cuts x, or quote 1–4 ranks only.

**Not fixed here — it is a design decision.** Options: (a) number the
nb-auto lattice in Cartesian rank order instead of Morton (the file
carries `block_nb_auto`, but the leaf-table readers check Morton order);
(b) choose the Morton bit significance by block-face area so the top-level
cut takes the smallest face (changes the canonical ids of every multi-block
case with equal bit counts per direction, and its mobygeom mirror);
(c) leave the ordering and give the benchmark an explicit placement. The
xz-mode key reorder of 2026-08-28 (peer traffic 20.5x lower) is the
precedent for (b).

## 4. Housekeeping found on the way

- `run_ncu.sh` judged success on a literal `jacobi` in the CSV, so a
  `KERNEL_RE=step_momentum` profile printed FAILED over complete data
  (job 5163917 did). It now greps its own regex.
- The three bit-exactness suite drivers deleted `.${pfx}.case.h5` before
  re-preparing; the solver's default name has no leading dot, so a case
  file from an earlier run with the same MODE label was reused and then
  refused as stale. Fixed (`0d7f4ab`).
- The solver banner's `commit:` is a hard-coded string in `init.f90`
  (`7aa1c7b`), identical in every build: it identifies nothing. Not changed.
