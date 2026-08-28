# Next session — 2:1 interfaces on the red-black SOR projection

## NOTE 2026-08-28 — the communication argument does NOT hold; the compute one does

Measured on the refined boundary-layer config at 2 GPU ranks
(`overheadTest/results_exchange_diag_2026-08-28.md`). Rounds per substage with
a 2:1 interface present:

| | per iteration | iters | phi | velocity | total |
|---|---|---|---|---|---|
| Jacobi `niter = 6` | 1 phi + 1 velocity | 6 | 6 | 6 | 12 |
| red-black `niter = 3` | 2 colours × (1 phi + 1 velocity) | 3 | 6 | 6 | **12** |

Red-black needs a communication **per colour**, and §4 gives each colour both a
cross-level phi exchange and a velocity exchange, so halving the iterations
cancels exactly against doubling the exchanges per iteration. On a single-level
grid red-black needs no phi and the count does halve (6 against 12) — but that
is where red-black already works.

Bytes do not save either, on the current decomposition: red-black's phi is
cross-level-only, and cross-level is **98.1 %** of the full phi exchange because
the 2:1 interface always lands on a rank boundary (the `refine_dims = xz` Morton
order puts the y tile in the high bits, so every fine block precedes every
coarse one; the cross-level peer count is identical at 2 and 4 ranks).

**The case to make is compute, and it is bigger.** `niter = 3` is 6 half-sweeps
against Jacobi's 6 `compute_phi` + 6 `apply` passes — roughly half the
projection arithmetic. `sweep` + `apply` are 46 % of the 2-rank step, so the
prize is of order 20 % of the step, several times anything the overlap work can
reach, and it helps single-rank runs equally.

**Do this before writing the 150–250 lines**: measure iterations-to-residual for
red-black vs damped Jacobi on the existing **single-level** path, where both
smoothers already run. It costs no new code and it is the number the whole
increment rests on ("3 SOR ≈ 6 Jacobi" is an expectation, not a measurement).
§6 already lists it as a gate; promote it to a precondition.

R0 (§3) remains CLOSED — attempted and reverted in `e77c75e`, not bit-exact
because block metrics differ across a periodic seam.

Branch `boundaryLayer` (or `claude/jacobi-interface`). Goal: **remove the
`[pressure] solver = redblack` × 2:1-refinement mutual exclusion**
(`pressure_solver.f90:98-108`) by giving the red-black sweep the one mechanism it
is missing at a level jump — cross-level transmission of the pressure increment
inside the projection loop.

Motivation is twofold. (i) SOR (`omega` up to ~1.5, Gauss-Seidel) converges
substantially faster per iteration than damped Jacobi (`omega <= 0.8`), and the
projection dominates the step. (ii) The boundary-layer campaign found
`chebyshev + niter=6 + Dirichlet-p outlet` to be **unstable on long steady runs**
(2-dx pressure mode, e-fold ~36 t.u.; see memory `boundary-layer-case` and
`docs/next_session_boundary_layer.md`), leaving plain Jacobi — stable but slow —
as the only option for that class of case. Red-black is the third option, and it
should not be restricted to single-level grids.

This is **not** a re-derivation of the interface numerics. It is a change of
*smoother* on an operator that is already validated.

## 0. Why this is tractable now (read before touching `interface_projection_derivation.md`)

`docs/interface_projection_derivation.md` §6c-§6g reads like a graveyard for
red-black + 2:1: mean-preserving prolong, pressure agglomeration, additive
Schwarz — all tried, all failed. **Those were failures of the OPERATOR, not of
the smoother.** They were fought on `claude/blocks`, where `D` and `G` were not
an adjoint pair at the interface (velocity injected across the jump, above-owns
ownership, one-sided relaxation). §1 of that same document states the criterion
that matters: *a fix is correct iff it keeps `G = -D^T`*.

The operator that criterion demands has since been built and validated (Phases
3b/3c/3d-file, the const-1/2 lockdown, the turbulent-channel and LES/IBM
campaigns):

- **composite gradient metric** `face_grad` (`pressure_solver.f90:453`):
  `(2/3)d1` on the fine side, `(4/3)d1` on the coarse side — the same physical
  `1/d` over the coarse+fine half-cell gap, each side formed from its own `d1`;
- **low-side-owns-face** ownership plus **symmetric relaxation on own copies**:
  both sides carry the shared face in the denominator and correct their own copy
  (the one-sided owner-overwrite variant is the unconditionally unstable one —
  `docs/block_refinement_strategy.md` §6);
- **mean-preserving cross-level transfer of the increment**: `ifaceRow` restrict
  of the fine interface ROW into the coarse ghost, injection prolong into the
  fine ghost (`comm.f90 copy_local_scalar_entries`; note the pressure BLEND
  weights `lWp/lWpDst` are **not** applied to scalars — the raw coarse value is
  what the composite metric expects);
- **conservative reconciliation** by the final full exchange.

So `L = DG` is already SPD and consistent at the interface. Red-black SOR is a
different *splitting* of that same `L`, and converges for any `0 < omega < 2` on
an SPD system. What follows is plumbing, not numerics.

**Do not** revisit the MP-prolong / agglomeration work of §6c-§6g. It targets a
defect the current operator does not have.

## 1. Two things that already work in `redblack_sweep`

1. **The interface metric is already correct.** `redblack_sweep` builds
   `gLo1..gHi3` from `face_grad_denom`, which defers to `face_grad` for every
   non-outlet face kind — including `FACE_COARSE` / `FACE_FINE`. Both the
   diagonal and the in-place correction therefore already carry the composite
   `1/d`.
2. **The redundant halo sweep is already off at interfaces.**
   `iLo/jLo/kLo = merge(1, 0, physLow(d,b) /= 0)` (`pressure_solver.f90:596-598`)
   suppresses the halo-layer sweep on any non-`FACE_OPEN` face, and
   `FACE_COARSE = 3` / `FACE_FINE = 4` are non-zero. Sweeping a halo cell that
   belongs to another level would be meaningless; the existing predicate already
   prevents it. Because the rule is keyed on **face kinds** — a property of the
   leaf table, not of the rank split — the `1 rank == 4 ranks` identity survives
   untouched.

## 2. The one missing mechanism

In-place SOR never materialises the pressure increment: it folds it straight into
`p` and the six faces. Same-level coupling still crosses block boundaries by two
routes, neither involving an increment array:

- the **redundant halo sweep** applies the neighbour cell's half of the
  correction to this block's own copy of the shared face, locally;
- the per-colour `exchange_halos(..., interp=.false.)` refreshes the halo copies.

**Across a level jump neither route exists**: the per-colour exchange is
same-level-copies-only by construction, and there is no `phi` to
restrict/prolong. The two sides would relax blind to each other, the shared face
would carry only half its gradient, and the projection would converge to the
wrong fixed point (or ring at the interface).

Jacobi solves exactly this with `exchange_scalar_halos(c, phi, blk,
ifaceRow=.true.)` (`pressure_solver.f90:186`) followed by the two-sided face
correction in `jacobi_apply`.

**Design rule this session establishes, for both smoothers:**

> `phi` crosses a block boundary **only to cross a level jump**.

## 3. Increment R0 (recommended first) — Jacobi's `phi` halo by redundant computation

Independent of red-black, cheaper, and it proves the redundant-halo consistency
argument on a path that is already gated by a large bit-exact suite. It also
turns §2's design rule into a fact on the Jacobi side.

**What `jacobi_apply` actually reads from the `phi` halo:**

- low faces → `phi(0,j,k)`, `phi(i,0,k)`, `phi(i,j,0)` with the *other two
  indices interior*: the three low-halo **face planes only**, never edges or
  corners;
- the high halo plane `phi(nb+1,...)` → **only** at outlet faces (mirrored
  locally by `apply_scalar_bc`, no comm) and at 2:1 interface faces.

The three low-halo face planes are exactly the cells red-black's redundant sweep
already computes (`nLowerHaloDirections == 1`, `pressure_solver.f90:605-609`). So
extend `jacobi_compute_phi` over that halo layer and **delete the same-level part
of the scalar exchange**.

**Why it is consistent.** Block A's halo cell at index 0 in direction `d` is
block B's cell at index `nb`. Its far face in `-d` is *interior to B* (never
wall/closed/interface), and its face kinds in the other two directions are shared
with A because the blocks are lattice-aligned. The `atBnd` tests `i == 1` /
`i == nb` therefore reproduce the owner's stencil exactly. (This is why `nb >= 4`
matters, and it is the same argument that makes the red-black velocity sweep
rank-independent.) Where the halo cell is *not* same-level, `iLo` is already 1 and
the case falls through to the cross-level exchange — precisely where it belongs.

**Payoff.** The Jacobi loop drops from two exchanges per iteration to one on
single-level grids. With `niter = 6-12` that is a large share of the projection,
and `docs/next_session_profiling.md` already names the per-iteration exchanges as
the suspected dominant cost. Extra compute: divergence + denominator on `3*nb^2`
cells against `nb^3` interior — **+9% of one kernel at nb=32**, +19% at nb=16.

**Catches.**

- **Chebyshev**: `cheb_combine` must run over the same halo layer so `delta` and
  `z` exist there; the recursion at a halo cell is identical to the owner's, so
  it stays consistent. Verify on a chebyshev case explicitly.
- The halo `phi` needs the halo **velocity** to be current: it is, from the
  previous iteration's velocity exchange. Do not reorder.
- Expected **bit-exact** (same expression, same inputs, same operand order) — but
  that is a claim to measure, not to assume. Arrange the halo-cell arithmetic to
  be literally the same source expression.

**Gate R0**: the standard 7-case suite bit-exact (`-Mnofma` / `-gpu=nofma`,
`max_abs 0`), CPU AND GPU, *including* a chebyshev case and a 2:1-refined case
(min_channel), plus `1 rank == 4 ranks` EXACT. Report the measured s/step change
on the refined channel and on the boundary-layer production grid.

## 4. Increment R1 — red-black + 2:1: "sweep + interface patch"

Give red-black a `phi` array, allocated **only when the case has 2:1
interfaces**, and apply interface faces in a separate post-colour kernel that
reuses the Jacobi formula verbatim.

Per colour, `redblack_projection` becomes:

```
if (hasInterface) zero phi
call redblack_sweep(color)          ! stores phi(cell) for this colour;
                                    ! is_interface faces stay in denom but are
                                    ! NOT corrected in place
if (hasInterface) then
    call exchange_scalar_halos(c, phi, blk, ifaceRow=.true.)   ! cross-level suffix only
    call interface_correct(blk, ibm, ...)                       ! the shared routine
end if
call apply_bc(blk, bc)
call exchange_halos(c, blk, [U,V,W(,P on the last colour)], interp=.false.)
```

**Why it is correct:**

- `interface_correct` is **byte-identical in formula to the interface branches of
  `jacobi_apply`**: `q_face += (phi_nbr - phi_self) * face_grad * mu`, i.e. the
  low-face branch at index 1 and the `is_interface` high-face branch at index
  `nb`. Factor it out of `jacobi_apply` and call it from both solvers — one
  implementation of the interface correction, per the no-duplication convention.
- **Zeroing `phi` per colour makes the per-colour halves sum to the full
  two-sided correction.** In the red pass only red cells hold a non-zero
  increment, so each interface face receives the red half; the black pass
  supplies the other half. No double counting.
- **Conservation holds half-sweep by half-sweep.** Each half is separately
  mean-preserving: the coarse ghost is the mean of the fine interface row
  (`ifaceRow`), the fine ghost is the injected coarse value, and
  `(4/3)d1_c == (2/3)d1_f`. Whichever subset of fine cells is active in a colour,
  the coarse face motion equals the mean of the fine sub-face motions.
- **Only the cross-level entries are needed** (§2). The entry lists are already
  ordered same-level-copies-first with prefix counts (`nLocalCopyPts`,
  `peerSendCopyOff` / `peerRecvCopyOff`), so a *suffix-only* phase for
  `exchange_scalar_halos` mirrors the existing `copyOnly` path. Running the full
  scalar exchange is also correct, just more traffic — do that first if it
  shortens the first working version.
- **Patch interior cells only** (`1..nb` in all three indices). Tangential-halo
  copies of an interface face are repaired by the same-level velocity exchange
  that follows — which is why the patch must run **before** it. Edge/corner `phi`
  ghosts (plain injection, not `ifaceRow`) are then never read.
- Masking `is_interface` faces out of the in-place sweep correction (for interior
  AND halo cells) means no stray writes into tangential halo copies of interface
  faces, and keeps the interface correction in exactly one place.

**Alternative considered**: leave the self-half in the in-place sweep and let the
patch add only the neighbour half. Mathematically identical (the correction is
additively decomposable: `(phi_below - phi_self) = phi_below + (-phi_self)`), but
it splits the interface formula across two kernels. Prefer the full-patch form.

### Concrete edits

| File | Change |
|---|---|
| `pressure_solver.f90` | Narrow the init guard (98-108): drop the refinement rejection, keep the Chebyshev rejection and the even-global-size-in-periodic check. Chebyshev stays exclusive — it needs a stationary linear operator. |
| `pressure_solver.f90` | Allocate `phi` for the red-black path only when interfaces exist (`block_refine_levels > 1 .or. refine_body .or. refine_nboxes > 0`, or better: a `blk` flag set from the actual leaf levels). |
| `pressure_solver.f90` | `redblack_sweep`: store `phi(i,j,k,b)`; mask `is_interface` faces out of the six in-place corrections (keep them in `denom`). |
| `pressure_solver.f90` | Factor the interface-face correction out of `jacobi_apply` into `interface_correct`; call from both paths. Jacobi must stay bit-exact through this refactor. |
| `pressure_solver.f90` | `redblack_projection`: add zero → phi exchange → patch inside the colour loop, all `hasInterface`-gated. |
| `comm.f90` | (optional, perf) `phase`/`crossOnly` argument on `exchange_scalar_halos` for the cross-level suffix, mirroring `copy_local_entries(phase=2)` and the `copyOnly` MPI prefix logic. |
| `config.f90` / `docs/configuration.md` | Update the `[pressure] solver` documentation once the exclusion is gone. |

Rough size: 150-250 lines in `pressure_solver.f90`, no new numerics.

## 5. The real risk, and the escalation ladder

Within one colour the two sides of a jump are relaxed **simultaneously** and patch
each other, so the cross-level coupling is **additive (Jacobi-like) while
everything within a level is Gauss-Seidel**. That is the BCM regime
(`block_refinement_strategy.md` §6: symmetric corrections, own copies,
stage-frozen ghosts) which is the *stable* one — the unstable variant was the
one-sided owner-overwrite. But it is untested at `omega = 1.5`: over-relaxation is
safe for the GS part and unproven for the additive interface rows.

If the interface rows ring, in increasing order of cost:

1. Clamp `omega -> min(omega, 1)` on interface-adjacent cells. A local,
   face-kind-keyed rule, so rank/`nb` independence is untouched.
2. **Level-ordered smoothing**: relax coarse-level cells, patch, then fine-level
   cells — the cross-level coupling becomes multiplicative/GS. One extra phi
   exchange per colour; this is the textbook AMR smoother ordering.
3. Extra interface-only smoothing passes (cheap: interface rows are a small
   fraction of the cells).
4. Hybrid: damped-Jacobi rows at the interface, SOR everywhere else.

Diagnose with the *residual divergence after a converged projection*, split
interior / coarse band / fine band — **not** the velocity change
(`jacobi_interface_handout.md`, "Method note"). A growing per-colour corrector
change at fixed `niter` is the signature of a non-contractive splitting.

## 6. Gates (all mandatory before declaring R1 done)

Dormancy / regression:

- **Single-level red-black bit-exact** vs `7e1e4b3` (`-Mnofma` / `-gpu=nofma`,
  `max_abs 0`), CPU AND GPU — guaranteed by the `hasInterface` gating.
- **Jacobi path bit-exact** through the `interface_correct` refactor on the
  standard 7-case suite, CPU AND GPU.

Interface correctness (the same ladder every interface phase has used):

- **Uniform / oblique flow through a 3D refined patch: EXACT (0.0, incl. `pn`)**
  with `solver = redblack` — the decisive transfer gate
  (`validation/refine2d/`, `validation/multilevel_body/`).
- **Global mass residual** with a refined patch → round-off.
- **Beltrami y-slab** interface regression; laminar channel-patch vs
  uniform-fine convergence order ≈ 2 (`validation/channel_interface/`).
- **`1 rank == 4 ranks` EXACT** and **CPU == GPU** with redblack + refinement.

Physics:

- **Developed Re_tau 180 wall-band channel**: no interface band in `u'`/`v'` (the
  validated reflux-off signature), compared against the *Jacobi* solution of the
  same case, not only against uniform-128.
- `refine_body` stability over ~2000 steps on a body case (naca0012 or sd7003 at
  reduced level).

Performance / motivation:

- Iterations-to-residual and s/step vs Jacobi and Chebyshev at equal work, on the
  refined channel.
- The `omega` stability limit **with** an interface present (sweep 0.8 → 1.7),
  reported explicitly — it is the number that decides whether red-black is worth
  using there.
- The boundary-layer motivation: does `redblack + niter=6 + Dirichlet-p outlet`
  stay clean where `chebyshev` grew the 2-dx mode? (`lowtrip_srdgridy_finex_redblack`
  is the single-level version of this test; repeat with a refined patch.)

## 7. Deferred / out of scope: unequal block sizes across an interface

Recorded here because it comes up whenever the block lattice is discussed.

- **Equal *number* of blocks per side is already not required.** A coarse block
  face is fed by up to 4 fine sub-entries (2 per edge, 1 per corner) in a fixed
  child order; counts differ by 4 (8 by volume) at every interface today.
- **Equal block *size* is required — but not by the interface.**
  `blk%nb(1:3)` is a single global value and `blk%q(0:nx+1,0:ny+1,0:nz+1,NVAR,
  nBlocks)` is one rectangular array over all blocks (`blocks.f90:62,308`). That
  invariant buys the `collapse(4)` volume kernels, the closed-form Morton
  `zorder_owner/start/count` and load balance, the leaf-table arithmetic mirrored
  in `moby_prepare`, and the `(nBlocksGlobal, nb^3)` io datasets.
- The interface **numerics** are indifferent to it: `face_grad`, low-side
  ownership, the `ifaceRow` restrict and the extent-based entries with per-dim
  affine gather maps are all per-face / per-cell. What they genuinely require is
  (a) **dyadic cell nesting** (midpoint subdivision, so fine sub-faces tile the
  coarse face exactly and restrict is a mean) and (b) the level change lying **on
  block faces** (a block has one level; face kinds are per block face).

Relaxation options, by cost: (1) tune `nb` — halo overhead scales ~`3/nb`
(measured: nb=32 +19%, nb=16 +49% vs the default layout) — and use
`refine_dims = xz`; (2) **per-level block size** (`nb_l`), the smallest
structural change that lets a coarse zone use big blocks — implement as one array
(or offset range) per level with volume kernels launched per level, NOT by
padding to `nb_max` (fine blocks dominate, 8x memory waste); (3) fully
variable-size patches (Berger-Oliger) — ragged arrays, neighbour search instead
of `neighbor_origin`, new load balance and file format: a different code.

Note that the usual felt inefficiency is **refinement granularity** (a whole
block refines, amplified by 2:1 smoothing and the one-block 26-neighbour buffer),
which is attacked by smaller `nb` or a tighter buffer rule — trading against halo
overhead in the opposite direction.

## 8. Suggested order of work

1. **R0** (§3): Jacobi `phi` halo by redundant computation. Bit-exact gate on the
   existing suite + chebyshev + refined case. Commit separately; report the
   s/step win.
2. **R1a** (§4): factor `interface_correct` out of `jacobi_apply`. Pure refactor,
   Jacobi bit-exact. Commit separately.
3. **R1b** (§4): red-black `phi` + zero + cross-level exchange + patch; drop the
   init guard. Run the §6 ladder, starting with uniform-flow-through-a-patch
   (0.0) — if that is not exact, stop and fix the transfer before looking at
   anything else.
4. **R1c**: the `omega` sweep with an interface, and the escalation ladder (§5)
   only if needed.
5. Update `CLAUDE.md` (the `[pressure] solver` line and the "mutually EXCLUSIVE
   with 2:1 refinement" claim), `docs/configuration.md`, and add a STATUS header
   here.

## References

- `docs/interface_projection_derivation.md` — §1 (the `G = -D^T` criterion) is
  the lens; §6c-§6g are **obsolete** (old operator).
- `docs/block_refinement_strategy.md` §6 — ownership, symmetric relaxation,
  stage-frozen ghosts, conservative reconciliation, red-black parity across a
  jump.
- `docs/jacobi_interface_handout.md` — the diagnosis method (residual divergence,
  banded).
- `docs/next_session_profiling.md` — why the exchange count matters (R0).
- `docs/next_session_boundary_layer.md` + memory `boundary-layer-case` — the
  chebyshev/outlet instability that motivates a third smoother.
