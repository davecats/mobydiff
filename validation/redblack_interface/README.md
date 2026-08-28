# Red-black SOR across a 2:1 interface (R1)

Gates for `[pressure] solver = redblack` on a refined grid — the mutual
exclusion with 2:1 refinement was removed in R1
(`docs/next_session_redblack_interface.md` §4).

```bash
./run_gates.sh cpu     # includes the 1 rank == 4 ranks gate
./run_gates.sh gpu
```

## What the mechanism is

In-place SOR never materialises the pressure increment, so across a level jump
the two sides would relax blind to each other. R1 gives red-black the one thing
it was missing: with an interface present the sweep also stores its increment in
`phi`, a per-colour cross-level exchange transmits it (`ifaceRow` restrict of
the fine interface row into the coarse ghost, injection into the fine ghost),
and `interface_correct` — the same routine the Jacobi path uses — applies

    q_face += (phi_neighbour - phi_self) * face_grad * mu

Interface faces stay in the sweep's DENOMINATOR but are not corrected in place.
`phi` is zeroed per **colour**, not per iteration, so each interface face
receives the red half in one pass and the black half in the other; the two sum
to the full two-sided correction with no double counting. Single-level runs
allocate no `phi` and take exactly the old code path.

## The blind gate — read before adding tests here

The uniform-flow-through-a-patch test, which the R1 plan called the decisive
transfer gate, **cannot see the interface patch at all**. Uniform flow makes the
divergence, and therefore `phi`, exactly zero, so `interface_correct` adds
nothing. Measured with the patch deliberately disabled: max deviation still
`0.000e+00`.

It remains a real gate on the halo TRANSFER operators — a per-direction indexing
mistake in the entry generation, the gather maps or the ghost blend still shows
up — and is kept for that. But it must not be quoted as evidence that the patch
works.

The gate with power over the patch is the **refined channel**. Negative control,
same case, 20 steps:

| | max‖u‖ | max‖p‖ |
|---|---|---|
| `interface_correct` enabled | 1.82e+01 | 2.05e-01 |
| `interface_correct` disabled | 4.42e+04 | 2.33e+11 |

## What these gates cover

- red-black + 2:1 runs and stays finite (and 200 steps by hand, stable);
- **1 rank == 4 ranks EXACT** on red-black + 2:1 — this exercises the
  cross-level phi exchange over MPI against the same-rank local-copy path;
- uniform oblique flow through a 3-level patch EXACT (transfer only).

Verified by hand alongside these, not automated here:

- **CPU == GPU EXACT** (`max_abs 0`, nofma binaries) on red-black + 2:1;
- **single-level red-black bit-exact** against the pre-R1 binary — the
  `hasIface` gating means an unrefined run is untouched;
- the **Jacobi 7-case suite bit-exact**, CPU and GPU, through the R1a extraction
  of `interface_correct` out of `jacobi_apply`.

## NOT yet done — R1 is not complete

`docs/next_session_redblack_interface.md` §6 lists these as mandatory before
declaring R1 done, and they have not been run:

- global mass residual with a refined patch → round-off;
- Beltrami y-slab interface regression, and the laminar channel-patch
  convergence order ≈ 2 against a uniform-fine reference (note
  `validation/beltrami/run_beltrami.sh` is stale: it drives a `MOBY_BELTRAMI`
  env hook removed in the 2026-06-30 cleanup);
- developed Re_τ 180 wall-band channel: no interface band in u′/v′, compared
  against the **Jacobi** solution of the same case;
- `refine_body` stability over ~2000 steps on a body case;
- the **`omega` stability limit with an interface present** (sweep 0.8 → 1.7).
  Within one colour the two sides of a jump are relaxed simultaneously and patch
  each other, so cross-level coupling is additive while everything within a
  level is Gauss-Seidel. That is the stable BCM regime, but it is untested at
  `omega = 1.5`, and §5 has the escalation ladder if the interface rows ring.
  These gates run at `sor = 1.5` and are stable at 20–200 steps, which is
  evidence but not the sweep.
- performance: iterations-to-residual and s/step against Jacobi **on a refined
  grid**. The 26–51 % sizing in
  `tutorials/turbulentBoundaryLayer/overheadTest/results_smoother_2026-08-28.md`
  is single-level and does not include the per-colour cross-level phi exchange
  R1 adds.
