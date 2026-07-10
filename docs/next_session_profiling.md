# Next session — profile & optimise (the refined channel)

Branch `claude/jacobi-interface`. Goal: **measure where the step time goes now**,
on the GPU, for the 2:1-refined channel, then attack the dominant cost. The
suspected bottleneck is halo exchange, but the last hard profile is **stale** and
must be redone before touching anything (see "Why the old numbers are stale").

This is an optimisation session: every change must be **bit-exact** vs the current
solver (it changes *how* work is scheduled, not *what* is computed). Build
`-Mnofma` / `-gpu=nofma`, compare with `tools/compare_fields.py` (max_abs 0), CPU
AND GPU, on the standard suite (min_channel, les_ibm channel + refine_body,
Beltrami y-slab). An overlap/batching refactor that is not bit-exact is a bug.

## Why the old numbers are stale (read first)

`[[refinement-perf-profile]]` (memory) recorded, at commit `a432af7` via the
`MOBY_PHASETIME` hook, a ~560 ms/step split: projection 63%, **reflux 23%**,
momentum 9%, vel-exchange 4% — "halo exchange ~62% of the step". Two things have
changed since, so this is a starting hypothesis, NOT a measurement to optimise
against:

1. **`momentum_reflux` is GONE** (production-config lockdown, commit `02c8927`
   area). The whole reflux (23% of the step, all halo exchange) no longer exists.
   The step is now dominated by the projection.
2. **`MOBY_PHASETIME` is GONE** (the MOBY_* hook cleanup). There is currently NO
   per-phase timer in the code. `chron` times the whole loop; `les_prof` times the
   three LES sub-phases only.

So step one is to **re-instate a lightweight phase timer and re-profile the
current code**, then decide the target from fresh data.

## Current step structure (what to time)

Per time step: 3 RK substages, each (main.f90, ~L205-234):

```
update_ibm_mu
update_bodyforce            (no-op unless [force] custom; cheap)
[turbulence, per model — timed by the existing turb_timing profiler:
  les:   update_sgs_viscosity -> turb%nut
  rans:  rans_substage (pin + ghosts + k/omega exchange_scalar_halos +
         fused transport + nut assembly; +gamma/Re_thetat when transition)
  iddes: update_sgs_viscosity(nut_sgs) + compute_iddes_fd (a second
         gradient-tensor sweep) + rans_substage + blend_iddes_nut
 then exchange_scalar_halos(nut)]
momentum(...)              predictor kernel (+ eddy-viscosity/body-force corrections)
apply_bc
exchange_halos([u,v,w], syncface=.true.)         ! 1 full velocity exchange
pressure_projection(...)    ! the loop below
```

`pressure_projection` (pressure_solver.f90, ~L115-135) loops `nIter` times, and
**each iteration does two halo exchanges**:

```
do iIter = 1, nIter
    jacobi_compute_phi(...)
    [cheb: cheb_combine(...)]
    exchange_scalar_halos(phi, ifaceRow=.true.)   ! (1) phi halo
    jacobi_apply(...)
    apply_bc(...)
    if (last) exchange_halos([u,v,w,p])            ! (2) full 4-var
    else      exchange_halos([u,v,w], interp=.false.) !     3-var velocity
end do
```

So the projection alone is `nIter` × (1 scalar + 1 velocity) exchanges per
substage. With `nIter=6`, that is 12 exchanges/substage × 3 = **36 exchanges/step
in the projection**, plus 3 syncface velocity exchanges — the reflux's 27/step are
gone. On a single rank `nPeers=0`, so these are device-side `copy_local_entries`
gather kernels over every leaf (12k+ leaves at nb=8), not MPI. On multi-rank they
are real messages.

The exchange machinery is already refactored (comm.f90): one weighted gather per
entry (`entry_gather_map`), entries ordered same-level-copy-first with prefix
counts, same-rank entries fused into one flat device copy kernel, off-rank into one
message per peer. The nonblocking API `start_halo_exchange` / `finish_halo_exchange`
exists but `exchange_halos` currently calls them back-to-back (no overlap).

## Step 1 — re-instate a phase timer and profile

Two options; do the first, keep it minimal and REMOVABLE (do not resurrect the
MOBY_* env-hook sprawl):

- **In-code phase timer.** Reuse `les_wall_seconds()` (system_clock wrapper in
  les.f90) — or a tiny `chron`-style accumulator — to bracket, per substage:
  `momentum`, `syncface exchange`, and inside the projection `jacobi kernels`
  (compute_phi + cheb_combine + jacobi_apply) vs `projection exchange`
  (phi-scalar + velocity). Print ms/step per phase at loop end (like
  `write_les_profile`). Keep it behind one clean flag or a compile guard, or just
  a `[output] profile = true` config key — NOT a scatter of env vars. Bracket on
  the host; on GPU insert `!$omp target update` only if you need a sync point
  (prefer wrapping whole phases so kernel-launch async cost shows honestly, then
  add one sync at the phase boundary).
- **GPU profiler.** `nsys`/`ncu` are NOT on PATH in this env (check
  `/opt/nvidia/hpc_sdk/*/*/*/profilers`); if available, `nsys profile --stats=true`
  on a ~20-step run gives a kernel/mem breakdown. The memory notes a prior nsys run
  showed `cuMemcpyDtoHAsync` 43% but with an incomplete timeline — trust the
  in-code phase timer over nsys here.

**Profile case:** the refined channel. Use `tutorials/min_channel/input.ini`
(128×64×8, nb=8, two wall bands refined to level 1, niter=6 chebyshev) — small and
representative — or scale to `validation/channel_interface/refined_y110.ini`
(128×64×128) for a heavier, more realistic split. GPU, ~50-100 steps after a short
warm-up. Sweep **nb = 8, 16, 32** (fewer/larger blocks) and note how the
exchange fraction moves — the Phase-1 redundant-halo-sweep + per-block entry
overhead scales with block count, so larger nb should shrink exchange cost (see the
Phase-1 note in CLAUDE.md: nb=32 +19%, nb=16 +49% halo overhead vs nb=default).

Report the fresh split (projection kernels vs projection exchange vs momentum vs
syncface). THAT decides the target.

## Step 2 — optimise the dominant cost

First, DON'T repeat two things already settled in `[[refinement-perf-profile]]`:
- **Entry-map precompute already landed** (commit `0b06bec`, ~12% faster,
  bit-exact): per-point→entry maps replaced the per-point binary search in the
  copy/pack/unpack kernels. The exchange kernels are already this-optimised.
- **Merging the per-iteration phi + velocity exchange was found NOT feasible**:
  `jacobi_apply` sits between them — the phi halo must be exchanged BEFORE
  `jacobi_apply` (it reads phi-gradient halos), the velocity exchange comes AFTER
  (feeds the next iteration's divergence). They cannot fuse. Do not re-try this.

Remaining hypotheses, in likely-payoff order (re-rank against fresh Step-1 data):

1. **Nonblocking overlap (core/shell split) — the real structural win.** Overlap
   the projection's halo exchange with the interior Jacobi compute: split each
   sweep into core points (stencil never touches a halo) and a one-cell shell;
   launch `start_halo_exchange`, compute the core while it is in flight, then the
   shell after `finish_halo_exchange`. Keep the update formula in ONE shared point
   kernel — do not duplicate the stencil. `docs/nonblocking_overlap_strategy.md`
   sketches this but predates the Chebyshev-Jacobi solver (it still describes the
   red-black solver — UPDATE that doc). This is where the exchange cost actually
   hides now that the reflux is gone.
2. **Shrink / remove the redundant open-halo sweep.** Per `[[refinement-perf-profile]]`
   #3: nb=8 (the wall-band granularity: a 24-cell band = 3 blocks of 8, which nb=16
   can't align) gives ~(10/8)³−1 ≈ 95% redundant-halo overhead on the volume
   kernels (jacobi/momentum). Larger blocks are NOT a general fix (the band forces
   small nb); the real lever is Phase 4 — overlap or eliminate the redundant halo
   sweep. Quantify the sweep's share in Step 1 before investing.
3. **Do the full velocity exchange less often.** Every Jacobi iteration re-exchanges
   all 3 velocity components (`interp=.false.`); only the last does the 4-var.
   Check whether the intermediate velocity halos are genuinely consumed before the
   next `jacobi_compute_phi` or whether the mid-loop velocity exchange can be
   thinned. Subtle (`jacobi_apply` reads velocity halos) — verify against the
   divergence gate, not just non-blowup. Lower-confidence; measure first.
4. **copy_local_entries kernel efficiency over many leaves** (single-rank path):
   the same-rank gather is one flat device kernel over all entries; check occupancy
   / launch overhead at 12k+ leaves.

## Gates (build `-Mnofma` / `-gpu=nofma`)

- **Bit-exact** vs the current binary (this is a scheduling change) on: min_channel
  (blocks+2:1+cheb), les_ibm channel (file IBM+WALE) + refine_body, Beltrami
  y-slab — CPU AND GPU, `tools/compare_fields.py` max_abs 0. Reuse the harness from
  the body-force session (scratchpad `gate/`: short inis, reference runs from the
  pre-change binary, then compare). On multi-rank, also 1-rank == N-rank bit-exact
  (the exchange is where a bug would hide).
- **Speed:** report ms/step before/after on the profile case (GPU), and confirm the
  optimisation actually moved the phase it targeted.
- Both `compile.sh cpu` and `compile.sh gpu` green after each step.

## Watch for

- The projection is a damped-Jacobi / Chebyshev-Jacobi smoother now (NOT red-black
  SOR); the old overlap doc and some comments still say "red-black". Fix as you go.
- The 2:1 interface transfer lives inside the exchange (`entry_gather_map`,
  RESTRICT/PROLONG, const-1/2 injection, ghost blend). Any exchange refactor must
  preserve it exactly — the interface is validated and locked (CLAUDE.md
  "Production-config lockdown"); do not change the numerics, only the scheduling.
- Body force adds one cheap correction kernel per substage only when `[force]` is
  enabled; irrelevant to the profile unless a `custom` force is on.
- Keep any new timer minimal and easy to remove — do not recreate the MOBY_* hook
  sprawl that was just deleted.
- RANS (T2-T4) adds per-substage `exchange_scalar_halos` calls on top of nut:
  k + omega always, gamma + Re_thetat with `[rans] transition = true` — 2-4 extra
  scalar exchanges per substage in the halo-bound budget. The anticipated remedy
  (user-approved 2026-07-10, sequenced AFTER T5 so the turbulence feature set
  re-gates once): an AUGMENTED q — the transported RANS scalars become extra
  cell-centred variable slots of `blk%q`, allocated only under RANS, so one
  batched `exchange_halos([...])` call carries them all (they ride the validated
  cell-centred 2:1 class like p/nut; nut itself stays put — its consumer chain is
  untouchable; the RK scratch/oldrhs pairs stay separate arrays). Pure bit-exact
  refactor, justified by the fresh profile of a RANS/IDDES case, not assumed.

## NEXT-SESSION PROMPT (updated 2026-07-10, post-IDDES-T5)

> Read `docs/next_session_profiling.md` and CLAUDE.md. Branch
> `claude/jacobi-interface`. RE-PROFILE the GPU step, then optimise the dominant
> cost. The last hard profile is stale (reflux removed; `MOBY_PHASETIME` deleted);
> the turbulence portion already has a timer (`turb_timing`, turbulence.f90) —
> extend the same pattern, do NOT add env hooks. STEP 1: add a minimal, removable
> phase timer behind one clean `[output]` config key, bracketing per substage:
> momentum / syncface exchange / projection-kernels / projection-exchange (and
> keep the existing turb_timing subdivision for the turbulence block). Profile on
> GPU, ~100 steps after warm-up: (a) the 2:1-refined channel
> `tutorials/min_channel/input.ini`, sweeping nb = 8/16/32; (b) a heavier refined
> case (`validation/channel_interface/refined_y110.ini` if present); (c) ONE
> RANS/IDDES case (`validation/iddes/iddes180.ini` 20-step or
> `validation/rans_sst/turb180.ini`) to measure the per-substage
> exchange_scalar_halos burden (k+omega+nut, +2 with transition, +fd's second
> gradient sweep under iddes) — this decides whether the user-approved AUGMENTED-q
> batching increment (RANS scalars as extra cell-centred `blk%q` slots, one
> batched exchange; nut and its consumer chain stay put) is worth its re-gating
> cost. STEP 2, from the fresh split only: the likely wins are, in order, the
> core/shell nonblocking overlap of the projection's per-iteration exchanges
> (start_halo_exchange / compute core / finish / shell — update
> `docs/nonblocking_overlap_strategy.md` off the red-black solver first), the
> redundant open-halo sweep share, thinning the mid-loop velocity exchange
> (subtle — verify against the divergence gate), and copy_local_entries occupancy
> at 12k+ leaves. Do NOT re-try the phi+velocity exchange merge (settled
> infeasible: jacobi_apply sits between them) and do NOT re-do the entry-map
> precompute (landed). Every change is a SCHEDULING change and must be bit-exact
> vs the pre-change binary (`-Mnofma`/`-gpu=nofma`, compare_fields max_abs 0,
> CPU AND GPU, 1==N ranks) on the standard list — min_channel, les_ibm
> ± refine_body, Beltrami y-slab (5 steps), and since RANS scalars ride the
> exchange add turb180 + lam30t (all fields incl. k/omega/nut/gamma/rethetat);
> an exchange-scheduling bug hides in the multi-rank path, so gate 4-rank
> explicitly. Report ms/step before/after per phase; confirm the optimisation
> moved the phase it targeted. Do NOT touch the locked 2:1-interface numerics
> inside the exchange. Make a plan first; execute after. Deferred elsewhere:
> the flat-plate inlet increment and transition/wall_function under iddes
> (docs/next_session_iddes.md).
