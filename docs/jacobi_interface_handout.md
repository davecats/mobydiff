# 2:1 interface on a damped-Jacobi projection — handout

Branch **`claude/jacobi-interface`** (forked from `cbd07d0`). This is a clean
restart of the 2:1 coarse–fine interface treatment on a deliberately simple
pressure solver, after every red-black scheme failed the local mass-conservation
gate. Read this first; then `docs/interface_projection_derivation.md` and
`docs/interface_review.md` for the literature analysis carried over from the
prior work.

## Why we restarted

On `claude/blocks` the production composite (above-owns) and the experimental
agglomeration / mean-preserving (MP) schemes were all driven with **red-black
SOR**. None achieved *local* mass conservation at the 2:1 interface; the
converged interface divergence floored at O(0.7) and the velocity trajectory was
worse than baseline. The strategy now: **simplify the smoother to damped Jacobi**
(no colouring) so the interface formulation is tractable, get it consistent on
the simplest cases, then re-introduce a faster scheme.

## What is in place (commits, oldest→newest)

- `56f121c` **Restart base** — config-driven `[flow] initial = beltrami|tgv`
  exact initial condition (replaces the `MOBY_BELTRAMI` env hack), set at each
  component's staggered coordinate. Carried over the Beltrami/TGV validation
  cases and the diagnostics (`check_beltrami`, `check_tgv`,
  `interface_diagnostics`, `plot_beltrami_fields`/`_slices`, `plot_tgv_error`)
  and the literature docs.
- `c15f416` **Damped-Jacobi projection** (replaces red-black SOR). Per iteration:
  `phi = -omega*div/denom` from the frozen velocity → `exchange_scalar_halos(phi)`
  (cross-level RESTRICT/PROLONG carries the coupling) → apply `p += phi/dt_gamma`
  and the per-face velocity corrections (race-free, no colouring). **Must be
  damped: `omega < 1`** (config key still `sor`, default 0.8). Verified
  empirically: `omega ≥ 1.1` diverges (checkerboard mode); Jacobi cannot be
  over-relaxed (SOR needs Gauss-Seidel/updated values, which red-black provides).
  Single-block Beltrami converges to the discretisation floor.
- `f1637e1` **Composite 2:1 interface stencil** — at a 2:1 face the gradient
  metric is the coarse-fine `1/d` (`face_grad`: `(2/3)d1` fine side, `(4/3)d1`
  coarse side; same physical `1/d`). Stops the interface from contaminating the
  interior (refined-patch interior error 54× better than the interface-agnostic
  sweep). This is a genuine improvement but **not yet consistent** (see below).
- `b43e604` **`MOBY_PROJONLY` / `MOBY_PREDONLY`** — skip momentum / skip
  projection, to partition the interface error between the correction and the
  predictor.
- `447b6db` **`MOBY_DIVDUMP`** — dump the discrete divergence the solver sees,
  before (`div_pre`) and after (`div_post`) the first projection, as companion
  field files. Captured at `rkStage==1` (projonly runs the projection at every
  substage, so a later capture would measure an already-multiply-projected field).

All of the above are **off by default and bit-exact otherwise**. CPU build is the
reference here (`build_cpu`, `-Mnofma`); GPU not yet re-validated on this branch.

## The verified diagnosis (the crux)

Using the exact-relaxation test (`MOBY_PROJONLY` projects the analytically
divergence-free Beltrami; `MOBY_DIVDUMP` reports `D` of the field):

slab_y, 32³, niter=1000:

| | interior | coarse band | **fine band** |
|---|---|---|---|
| `D·u_exact` (exact div-free field) | **0.0 exact** | 1.2e-2 | **0.997 — O(1)** |
| `D·u_after` (1 converged projection) | 1.7e-4 | 1.5e-3 | **0.124** |

**The operators are inconsistent at the fine interface.** The discrete divergence
of a divergence-free field is O(1) on the fine interface cells (it is exactly 0
in the interior — the Beltrami has zero discrete divergence on a uniform staggered
grid), and a *fully converged* projection cannot drive it to zero: `G` does not
undo what `D` measures, so `D` and `G` are not a consistent (adjoint) pair there.

**Root cause:** the fine cell's interface *normal velocity* is a **coarse-
resolution** value — it is filled by the velocity PROLONG / exchange instead of
being fine-owned. The fine divergence then differences a coarse-resolution flux
against a fine-resolution one: error `O(h_coarse)·(1/h_fine) = O(1)`. (Earlier
guesses — a low-k tangential-averaging O(h) error, and a fine-side checkerboard —
were a minor coarse-side effect and a symptom respectively, not the cause. A
coarse-side tangential-curvature correction `tang_curv` was tried and reverted:
it did nothing, because the defect is this O(1) fine-side one.)

Method note: **measure the residual divergence `D·u_after`, not the velocity
change** — the change conflates "the projection legitimately removing the exact
field's nonzero discrete interface divergence" with "the operators are broken".

## The fix (next step, scoped)

Make `D` and `G` consistent at the interface, per the literature (Almgren–Bell–
Colella JCP 1998; Martin–Colella JCP 2000 — quadratic coarse-fine interpolation;
Guittet–Theillard–Gibou JCP 2015 — stable non-graded MAC) and the old branch's
`mp_fine_owns`:

- **Only the pressure crosses the interface by PROLONG.** The interface normal
  velocity is **fine-owned at fine resolution**: the fine reconstructs its
  interface face from its own predictor plus the prolonged coarse *pressure*
  gradient (`face_grad`/`ifGrad`), and that face is **RESTRICTED to the coarse**
  (average of the 4 fine sub-faces) for conservation.
- The velocity exchange must **stop overwriting the fine-owned interface face**
  (prolong fills only the deep halo needed by the momentum stencil). In this
  Jacobi codebase the initial full exchange and the final full exchange currently
  prolong the coarse velocity onto the fine interface face — that is the corruption.
- Re-measure `D·u_after`: success = fine-band residual → discretisation floor
  (O(h²)), not O(0.1).

## How to run the gate

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
cmake --build build_cpu -j        # build_cpu is -Mnofma; CPU is the reference
# single-block sanity (no interface): projection change must be ~0
# refined / slab consistency test (the decisive one):
#   MOBY_PROJONLY=1 MOBY_DIVDUMP=1 mpirun -x MOBY_PROJONLY -x MOBY_DIVDUMP -n 1 \
#     build_cpu/main input.ini       # input from validation/beltrami/slab_y.ini,
#                                     # nx=32 nb=4 nsteps=1 niter=1000 field_interval=0
# then compare r_divpre_<step>.h5 / r_divpost_<step>.h5 (pn = divergence),
# split coarse/fine interface band vs interior.
```

Cases: `validation/beltrami/{uniform,slab_y,refined_fast}.ini` (set `initial =
beltrami`, `sor = 0.8` for Jacobi). `slab_y` isolates a single flat interface
direction; `refined_fast` is the 3D patch (edges + corners). Diagnostics:
`tools/interface_diagnostics.py` (band/interior error + roughness),
`tools/plot_beltrami_fields.py` (solver/exact/error slices). Beware: run cases
**one at a time** — `mpirun -n 1` still spawns a small process tree; do not
`pkill build_cpu/main` while a run is live (it matches the running mpirun and you
will kill your own job — this bit us repeatedly).

## What we are NOT chasing

- The broad low-k **pressure** error is a Jacobi under-convergence artifact
  (present single-block too) — ignore it; it clears with more iterations / a
  faster solver later.
- **Long-time Beltrami stability** is a chaotic, uninformative metric (Davide).
  Use the 1-step exact-relaxation `D·u_after` consistency gate instead.
