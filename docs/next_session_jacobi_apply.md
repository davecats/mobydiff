# Next session — `jacobi_apply`, the largest single phase in the step

Handout. **Measure first.** The previous four phases of this track all turned on
a measurement that contradicted the obvious guess; assume this one will too.

## Why this target

Per-phase on the boundary-layer production grid (A6000, 400 cold steps, 1 rank,
`nb = 64 44 48`, 1.2596 s/step — `overheadTest/results_phases01_breakdown_2026-08-27.md`):

| phase | s/step | share |
|---|---|---|
| **`proj/apply`** | **0.4072** | **32.3 %** |
| `step/momentum` | 0.3505 | 27.8 % |
| `proj/sweep` | 0.2806 | 22.3 % |
| `proj/vel_exchange` | 0.0664 | 5.3 % |
| `proj/phi_exchange` | 0.0229 | 1.8 % |
| `exch/local_copy` (all exchanges) | 0.0976 | 7.8 % |

`jacobi_apply` is the biggest single item in the step — larger than the fused
momentum predictor, and 6x the exchange that Phase 3 just optimised. **The
exchange story is over**: everything the block lattice costs in traffic has been
attacked (block tax 1.4297 -> ~1.04). What is left is kernel efficiency.

Note the phase shares moved a lot during 2026-08-27, so **re-profile before
committing to anything** (`PROFILE=1 ./run_overhead.sh`, then `phase_table.py`).
The numbers above are post-Phase-1 but PRE the Phase-3 divergence sync.

## What the kernel does

`pressure_solver.f90 jacobi_apply` — TWO `target teams distribute parallel do
collapse(4)` launches over the same interior index space, called `nIter` times
per RK substage (18x per step at niter 6):

1. **Pressure update**: `q(i,j,k,VAR_P) += phi(i,j,k)*idt`. Pure streaming:
   read p, read phi, write p.
2. **Velocity face corrections**: per cell, three LOW-face corrections
   `q_face += (phi_below - phi_self)*face_grad_corr(...)*mu`, each guarded by
   `if (cf /= 0.0d0)`, plus three conditional HIGH-face corrections for owned
   2:1-interface and outlet faces.

## Hypotheses, in the order they are worth testing

None of these is established. Rank them by what `ncu` says, not by this list.

1. **`phi` is read twice** — once per launch, over the whole interior. Fusing
   the two passes into one kernel removes an `nb^3` read of `phi` and one launch.
   Cheap to try, and bit-exact *if* the fused order keeps each cell's arithmetic
   identical (the pressure update and the face corrections touch disjoint
   outputs, so fusing is a scheduling change, not a numerics one). Rough ceiling:
   `phi` is 1 of ~16 arrays touched, so ~6 % of `apply` = ~2 % of the step.
2. **Register pressure / occupancy**, the copy-kernel failure mode. `apply`
   carries `face_grad_corr` calls, six branches and the interface/outlet
   predicates. If `Block Limit Registers` caps occupancy the way it did for the
   halo gather (4 blocks/SM, 33 %), splitting the rare interface/outlet
   high-face work into its own small kernel would leave a lean common path.
3. **Branch divergence**: `if (cf /= 0.0d0)` and the `i == nx .and.
   is_interface(...)` tests are uniform across almost every warp (interior
   cells) but not at block faces.
4. **Store efficiency**: the three face writes are to `VAR_U/V/W` planes
   `(nb+2)^3` apart, the same structural issue that put a 1.84x floor under the
   halo copy. Unlike that case the accesses here are contiguous in `i`, so this
   is probably NOT the limiter — verify rather than assume.

## Method

```
mpirun -n 1 ncu --kernel-name regex:jacobi_apply --launch-skip 12 --launch-count 4 \
    --section SpeedOfLight --section MemoryWorkloadAnalysis --section Occupancy \
    --section WarpStateStats <bin> <ini>
```

- **Performance counters are permission-blocked on the local workstation**
  (`ERR_NVGPUCTRPERM`). `ncu` works on **istmcetus**; GPU 1 is the free one.
- A small case profiles fine and iterates in seconds: the production ini with
  `nx` and `nz` cut 8x (`~/.moby_prof/rect_small.ini` pattern; keep `nb`
  dividing the reduced grid).
- Read `Block Limit Registers` and `Achieved Occupancy` FIRST. For the halo
  gather those two, not bandwidth, were the whole story.

## Gates

The standard discipline: **bit-exact, `max_abs 0`, CPU AND GPU** on the 7-case
suite (`scratchpad gate_bitexact.sh` pattern — reference binary from a detached
worktree at HEAD, both sides `-Mnofma` / `-gpu=nofma`), plus
`validation/block_nb/run_gates.sh` and 1 rank == 4 ranks.

Fusing kernels or splitting off a rare path is a scheduling change and MUST be
bit-exact. If it is not, something reordered the arithmetic — find it rather
than accepting a tolerance.

## Hard-won process notes from the 2026-08-27 session

- **Same-day A/B or nothing.** Machine state drifted ~5 % over three weeks and
  faked a gain. `run_overhead.sh` takes `SUFFIX=<tag>`; build the reference from
  a worktree at HEAD and run both binaries the same afternoon.
- **Read `runtime.txt`'s drift, not the final average.** A contended run reads
  as a fixed cost; a rising or rise-then-fall cumulative average means
  contention, not your change.
- **Never rebuild a binary while a timing matrix is running.** Done twice; the
  second time it also left an orphaned solver on the GPU (kill on the SOLVER
  cmdline, then verify with `pgrep` + `nvidia-smi`).
- **Read the complete pass/fail list before theorising.** A partial gate list
  produced a confident and wrong "chebyshev correlates perfectly" diagnosis.
- **`1 rank == 4 ranks` does not catch everything.** It cancels anything both
  rank counts share — it missed the periodic-seam metric asymmetry entirely.

## NEXT-SESSION PROMPT

> Read `docs/next_session_jacobi_apply.md` and CLAUDE.md. Branch
> `optimiseBlockRefinement_parentBoundaryLayer`, HEAD `730d0ca`; Phases 0, 1 and
> 3 (vel_exchange) are DONE and committed — block tax 1.4297 -> 1.0369,
> production step 1.6923 -> 1.2192 s/step, all bit-exact. **Build fresh nofma
> reference binaries at HEAD before touching any code** (detached worktree +
> `build_nofma.sh` pattern; recreate it from the doc if the scratchpad is gone).
>
> Target: **`jacobi_apply`, 32.3 % of the step** — the largest single phase,
> bigger than the fused momentum predictor, and never examined. The exchange
> work is finished; what is left is kernel efficiency.
>
> **MEASURE BEFORE CHANGING ANYTHING, in this order.** (1) Re-profile:
> `PROFILE=1 ./run_overhead.sh` in `tutorials/turbulentBoundaryLayer/overheadTest`
> then `phase_table.py` — the shares moved a lot on 2026-08-27 and the numbers in
> this doc predate the divergence sync. (2) `ncu` the kernel on **istmcetus**
> (performance counters are permission-blocked on the local workstation,
> `ERR_NVGPUCTRPERM`; GPU 1 is the free one), recipe in the doc. Read
> `Block Limit Registers` and `Achieved Occupancy` FIRST. (3) Only then pick a
> fix, from the ranked hypotheses in the doc or a better one the profile
> suggests. This ordering is not ceremony: this track's obvious guess was wrong
> twice — the halo copy was register-limited, not bandwidth-limited, and R0's
> predicted bit-exactness failed on a periodic-seam metric asymmetry.
>
> Gates for any change, non-negotiable: the 7-case suite bit-exact (`max_abs 0`,
> `-Mnofma` / `-gpu=nofma`, **CPU AND GPU**), `validation/block_nb/run_gates.sh`,
> and 1 rank == 4 ranks. Fusing kernels or splitting off a rare path is a
> scheduling change and MUST be bit-exact; if it is not, something reordered the
> arithmetic — find it rather than accepting a tolerance.
>
> Report the gain as a **same-day A/B**: `SUFFIX=<tag>` + `BIN=<ref>` on
> `rect_jacobi.ini` (the production layout) with `base_jacobi.ini` as the
> control, both binaries the same afternoon. Machine state drifts ~5 % over
> weeks and has already faked a gain once. Judge a timing delta from
> `runtime.txt`'s drift within the run, not the final cumulative average, and
> **never rebuild a binary while a timing matrix is running**.
>
> Do NOT reopen: the halo copy kernel (at its ~51 %-of-peak coalescing floor,
> `docs/next_session_block_overhead.md`), or R0 / redundant-computation schemes
> (blocked by the periodic-seam metric asymmetry — fixing that is its own
> non-bit-exact numerics change and buys < 1.8 %).
>
> Stop cleanly after any increment: one commit, its own gate, STATUS header in
> the handout, and a results file in `overheadTest/`.

## Related

- `docs/next_session_block_overhead.md` — STATUS headers for Phases 0, 1, 3 and
  the reverted 3a, with the cost model and its caveats.
- `validation/block_nb/README.md` — the 2:1 interface is NOT nb-independent;
  refined layouts must be gated on uniform-flow preservation.
- Phase 3a's failure: block metrics are not bitwise equal across a PERIODIC
  seam. Any future scheme that recomputes instead of exchanging hits this.
