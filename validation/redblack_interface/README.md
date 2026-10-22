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

## Re-validated 2026-10-21, and the R1 list closed

RE-VALIDATED 2026-10-21, commit 532ef54 (banner), production flags, workstation
CPU (52 s, 4 ranks) and RTX 3060 (6 s): `run_gates.sh cpu` and `gpu` ALL PASS
(finite: max|un| 18.23, 1 rank == 4 ranks EXACT, the 3-level transfer gate
0.0 / 0.0). The red-black suite case `../block_nb/beltrami_rect_redblack.ini`
is max_abs 0 against the pass-2 reference set (CPU and GPU, nofma).
Settings changed for the rerun: none.

The items of the list below were run the same day (`run_r1_items.sh`,
decision 6 of the re-validation handout):

- **the developed Re_tau 180 wall-band channel against Chebyshev** --
  `../channel_interface/developed/README.md`, "The red-black twin": the
  same statistics to the sampling scatter (<= 0.5 %), no interface band in
  either, red-black 6 at 0.143 s/step against Chebyshev 12 at 0.120 on an
  A6000. **Red-black SOR across a 2:1 interface works** and is admissible
  on refined cases (decision 1); it is not the faster projection there.
- **the omega stability limit with an interface present** (`run_r1_items.sh
  cpu`, `refined_channel.ini` for 200 steps on 4 ranks): stable at sor 0.8,
  1.0, 1.2, 1.4, 1.5 and 1.6 (max|un| 18.28 throughout), DIVERGES at 1.7
  (max|un| 3e44). The production 1.5 has margin; the limit lies in (1.6, 1.7).
- **the Beltrami y-slab regression** (`run_r1_items.sh gpu`, istmcetus A6000):
  the exact decaying flow through a y-normal 2:1 slab at t = 2, L2 velocity
  error 3.095e-3 with red-black 6 / sor 1.5 against 3.088e-3 with Chebyshev
  12 -- ratio 1.002 (tolerance 1.5): the per-colour transfer adds no
  interface error over the Jacobi path's (which is itself the const-1/2
  transfer's, `../beltrami/README.md`).
- **`refine_body` over 2000 steps** (the les_ibm `channel_ibm_refine.ini`,
  2:1 around an immersed wall, red-black 6 / sor 1.5): finite, max|un| 22.3.
- **the global mass residual with a refined patch**: not run -- the `divsum`
  tool left with the `MOBY_*` hooks; the surviving mass statement is the
  body_outlet checker's interior rms divergence (1.0e-4 on the pocket case).
- **performance**: the measured 19 % per-step cost on the band channel above
  is the number; the exchange-round reduction stays with the performance
  handout.

## The R1 list as it stood before 2026-10-21

`docs/next_session_redblack_interface.md` §6 lists these as mandatory before
declaring R1 done, and they have not been run:

- global mass residual with a refined patch → round-off;
- Beltrami y-slab interface regression, and the laminar channel-patch
  convergence order ≈ 2 against a uniform-fine reference (the Beltrami inis
  select their IC with `[flow] initial = beltrami`; `run_beltrami.sh` no longer
  relies on the removed `MOBY_BELTRAMI` hook, updated 2026-09-26);
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
