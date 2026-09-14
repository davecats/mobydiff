# The between-iteration velocity exchange, reduced to what the divergence reads

Job 5145805 (A/B at 4 and 8 ranks on both production cases, Pass G, the 7-case
suite — one allocation) and 5145806 (Pass G re-run). A100-SXM4-40GB, HoreKa,
`dev_accelerated`, 100 steps. `ref` = `3c2903a`, `new` = `95312d7`. Raw runs in
`horeka/results_job5145805/`. Readings written down first, in
`horeka/exchange/PREREGISTERED_divhalo.md`.

Task 1 of `docs/next_session_divergence_halos.md`.

## 1 — The change

15 of the 18 velocity halo rounds per step sit between projection iterations,
and the only velocity halo anything reads before the next divergence is
`q(nb+1)` of the component NORMAL to that face. The solver already knew this and
already exploited it — **at one rank.** With peers it fell back to the full
26-direction same-level shell of three components, because those planes would
have to come over MPI and the entry list was not partitioned for it.

It is now. `entry_round` gives the enumeration a third round: pure `+axis`
same-level face COPY entries first, then the remaining same-level copies, then
the cross-level entries. A divergence round is therefore a per-peer **prefix of
the copy prefix** — the trick `copyOnly` already used, one level down — derived
identically on both ends of every message without negotiation, and it carries
**one variable per entry** (`lDivVar`/`sDivVar`/`rDivVar`, the face normal)
instead of three. `dsSlot` and its three little kernels went with it: the
same-rank half is that same prefix, run by one kernel, peers or not.

## 2 — Why the reduced set is complete: the entry list, not a gate

The handout asked for this to be checked against the entry list rather than
assumed. It holds by construction.

1. **What is read.** `jacobi_compute_phi` forms the divergence over
   `i,j,k = 1..nb`, so the only halo cells it touches are
   `q(nb(d)+1, ·, ·, VAR_d)` at *interior* tangential indices. `jacobi_apply`
   reads `phi`, not the velocity halo; `interface_correct` is documented
   interior-tangential-only; `apply_bc` works on physical ghosts the block owns.
2. **What writes it.** From `entry_boxes`: `off(d) = +1` gives
   `dstLo(d) = nb(d)+1, ext(d) = 1`; `off(d) = -1` writes the plane 0;
   `off(d) = 0` reaches `nb(d)+1` only through the tangential extension, i.e.
   only at *halo* indices of the other dims. So inside the range the divergence
   reads, **the pure `+axis` face entries are the only writers** — edges and
   corners land at tangential index 0 or nb+1, outside it.
3. **Cross-level ones write nothing there anyway.** At a 2:1 `+axis` face
   `interface_normal_dim` returns `d` and unpack skips `var == rNrm`, because
   the low-side block owns that face and reconstructs it from its own pressure.
   Excluding them is not an approximation.

So: pure `+axis` same-level face COPY entries **are** the writers of the
divergence halo. Everything else the full copy round wrote is either not read
until the next substage — after the last iteration's full exchange restores it —
or is written with an identical value.

## 3 — Volume: the partition selected what it was meant to

Measured, both production cases, at both rank counts:

| case | wire, `div` vs `copy-only` | device-local, `div` vs `copy-only` |
|---|---|---|
| `rect_jacobi` | 1.914 vs 12.365 MB/round → **6.46x** | 7.19 M vs 44.62 M pts → **6.20x** |
| `refined_yp82` | 1.199 vs 7.627 → **6.36x** | 2.47 M vs 15.15 M → **6.14x** |

Pre-registered band: 1/6.0 to 1/6.5. **PASS**, and it lands on the 6.24x the
entry arithmetic predicts for `nb = 64 44 48` (8000 values against 49896).

## 4 — Time (job 5145805, before/after in one allocation)

| case | `proj vel_exchange` | step |
|---|---|---|
| `rect_jacobi` 4×1 | 12.295 → 5.799 ms **−52.8 %** | 158.19 → 153.11 **−3.21 %** |
| `rect_jacobi` 8×2 | 7.471 → 3.473 **−53.5 %** | 83.83 → 80.02 **−4.55 %** |
| `refined_yp82` 4×1 | 6.201 → 3.351 **−46.0 %** | 74.91 → 72.71 **−2.93 %** |
| `refined_yp82` 8×2 | 4.352 → 2.321 **−46.7 %** | 42.68 → 40.58 **−4.93 %** |

**Controls.** `apply`, `sweep` and `setup` move by at most 0.4 % on every case,
and the `before` column reproduces job 5145554's `after` column (refined 8×2
`proj vel_exchange` 4.352 here against 4.353 there, 0.01 %), so the allocation
and the binaries are comparable.

**The launch count did not fall, as predicted**: 39 pack, 39 unpack, 39 local
copy and 24 cross-level copy per step on both sides. The saving is bytes and
wait, which is what the exchange buckets say it is:

| bucket, 8 ranks | `rect` before → after | `refined` before → after |
|---|---|---|
| `local_copy` | 7.472 → 5.437 **−27.2 %** | 3.195 → 2.324 **−27.3 %** |
| `mpi_wait` | 2.906 → 1.962 **−32.5 %** | 1.756 → 1.336 **−23.9 %** |
| `pack` | 1.280 → 0.713 **−44.3 %** | 1.235 → 0.700 **−43.3 %** |
| `unpack` | 1.305 → 0.983 **−24.7 %** | 1.084 → 0.869 **−19.8 %** |

Those four sum to −3.86 ms on `rect` 8×2 against a `proj vel_exchange` fall of
−4.00 ms, and to −2.04 against −2.03 on `refined` 8×2 — the bucket arithmetic
closes. **The largest single contributor is the device-local copy, not the
message**, which is why the change is worth something at one rank too.

Per between-iteration round (15 of the 18): **267 µs saved** on `rect` 8×2,
**135 µs** on `refined` 8×2, **433 / 190 µs** at 4 ranks. The handout forecast
412 → ~120 µs on `rect` 8×2, i.e. ~292 µs; measured 267.

## 5 — Graded against the pre-registration

| prediction | outcome |
|---|---|
| `proj vel_exchange` −40 to −60 % | **PASS** — −52.8 / −53.5 / −46.0 / −46.7 |
| `div nv=1` at 1/6.0 to 1/6.5 of `copy-only nv=3` | **PASS** — 6.46x and 6.36x |
| step −3.5 to −5.5 % | **PASS at 8 ranks** (−4.55, −4.93); **MISSED LOW at 4 ranks** (−3.21, −2.93) |
| any nonzero `max_abs` | none — see §7 |

**The 4-rank miss has a named cause, and it is not this change's own bucket.**
`phi_exchange` — whose volume the change does not touch at all — went
**+7.3 % (`rect`) and +9.6 % (`refined`) at 4 ranks**, while going −1.7 % and
−2.9 % at 8. That is +0.49 and +0.42 ms — **7 % of the `rect` 4-rank gain and
15 % of the `refined` one**, and it is most of why those two fell short.

The entry REORDERING is the only thing that can explain it. The scalar exchange
walks points in entry order, and where the copies used to run all 26 directions
of block 1, then all 26 of block 2, they now run the three `+axis` faces of
every block and then the other 23 of every block. Same points, same values,
different scatter locality. Why that costs ~8 % at 4 ranks and *pays* ~2 % at 8
is not established here; the per-rank point count is 2x larger at 4 ranks, so a
cache-residency explanation is available but untested. **It is recorded as an
observed redistribution, not explained.**

## 6 — What these numbers do not support

- **Nothing about 16 ranks.** `dev_accelerated` is 2 nodes. The message is
  smaller, so `mpi_wait` plausibly benefits more where the network matters more,
  but the 8-rank `mpi_wait` gain (−32 % / −24 %) is the basis for any
  extrapolation and must be labelled one.
- **Nothing about a launch saving.** Round counts are identical on both sides;
  anyone reading this as fewer kernels is reading it wrong.
- **Nothing about red-black.** It uses per-colour copy-only exchanges, not this
  path. Its volume is unchanged — but its wire layout is not, because the
  reordering touches every exchange, which is why it has its own gate (§7).
- **The campaign ratios move again.** Both twins shrink and they shrink
  unevenly (−4.55 % against −4.93 % at 8 ranks), so the block tax, the
  strong-scaling efficiencies and the 2:1 coarse-cell-equivalent all have a
  moved denominator on top of the one job 5145798 is already re-measuring.

## 7 — Gates

**Production flags on both sides, which is the tighter test here**: the values
are copied, never recomputed, and every surviving kernel's arithmetic is
textually what it was, so no expression moves and both sides contract
identically (the rule in `run_mapgate.sh`).

| gate | result |
|---|---|
| 7-case suite, 200 steps, GPU (job 5145805) | **max_abs 0**, all seven |
| `min_channel` 4 ranks / 1 rank (blocks + 2:1 + Chebyshev) | 0 / 0 |
| `beltrami_slaby` | 0 (4 datasets) |
| `les_ibm` (file IBM + WALE) | 0 (5 datasets, incl. `nut`) |
| `turb180`, `wf180_y30` | 0 (7 datasets, incl. `k`, `omega`) |
| `lam30t` | 0 (9 datasets, incl. `gamma`, `rethetat`) |
| Pass G, production cases (job 5145806) | **max_abs 0** — `rect` 138 412 032 pts, `refined` 60 555 264 pts |
| `min_channel` 1 == 2 == 3 == 4 ranks, CPU | **exactly equal** — four peer topologies |
| `validation/redblack_interface` 1 and 4 ranks, CPU | **max_abs 0** |

The rank-count identity is the load-bearing one: it is the standing gate the
single-rank reduced path was validated against, and it now compares the reduced
path against itself across four different partitions of the same blocks.

**A tooling failure worth recording.** Job 5145805's Pass G produced no
comparison: `h5maxdiff` is a build product, not a tracked file, and the pinned
worktree did not have it — while `run_exchange.sh` deletes its snapshots whether
or not the comparison ran, so the gate was *lost* rather than merely unreported.
`submit_divhalo.sh` now builds `h5maxdiff` before the passes, and 5145806 re-runs
the pass. Pinning a job to a fresh worktree costs the build products the live
tree happens to carry.

## 8 — What is left

`docs/next_session_divergence_halos.md` §5 (`compute_rdenom`'s fp64 divide,
~4 % of the step, needs ncu first) and §6 (splitting `step_momentum`, its own
session). Plus, new from this job: the `phi_exchange` locality effect above — if
the 4-rank sign is real and reproducible, the scalar exchange's entry order is
worth one measurement, because it is a free 0.5 ms if the ordering can be
decoupled from the prefix requirement.
