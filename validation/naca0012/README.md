# A3 INCREMENT 2 — NACA 0012 full-turbulent SST sanity (Re_c = 1e5)

Quasi-2D airfoil case through the full A0-A2 + A3 stack: `[case] name =
airfoil` freestream composition (inlets x_min/y_min/y_max at the aoa
angle, outlet x_max), file-based IBM from an extruded NACA 0012 STL,
`refine_body` at 5 levels, `[turbulence] model = rans` (SST, transition
off, resolved walls), penalization C_L/C_D as the runtime statistic, and
the INCREMENT-1 scalar inlet values active at the three inlet faces.

Grid: 12c x 12c x 0.1875c, base 512x512x8 (nb = 8, Delta0 = 0.0234c),
refine_levels = 4 -> Delta = 1.465e-3c at the surface (y+_1 ~ 2-4 at the
turbulent midchord). LE at (4.5, 6.0); the freestream carries the angle,
the geometry never rotates. dtmax = 4e-4: the eddy-viscosity correction
is explicit and the Peclet limiter molecular-only — see the ini comment.

Gate (sanity band, not tight — first-order penalization D_eff ~ D + h,
Dirichlet blockage, y+ resolution): |C_L(0)| < 0.02, lift slope over
aoa = 0/4/8 in [0.06, 0.13]/deg (2pi rad = 0.1097/deg), C_D(0) in
[0.008, 0.035] (XFOIL fully-turbulent NACA 0012 at Re 1e5: C_D0 ~
0.012-0.015, C_L(4) ~ 0.40-0.44).

## FOUND while gating (2026-07-14): penalization forces need --keep-buried

The first aoa = 4 run converged STEADY to C_L = 0.018 / C_D = 0.0166 while
the FLOW carried the full lift: box-contour circulation Gamma = -0.185
(contour-independent) => Kutta-Joukowski C_L ~ 0.37, classic suction-side
speedup, zero reversed flow. The A2 reduction is faithful (an independent
numpy sum over coef*u*V reproduces the solver to all digits): the SINK
itself is incomplete. refine_body removes leaves buried inside the body;
their FACE_CLOSED faces are exact zero-flux walls that absorb the core
pressure loading with NO coef bookkeeping, so sum(coef u dV) keeps only
the thin-shell contribution — the friction-dominated drag survives, the
pressure-dominated lift vanishes. The validated cylinder was immune by
accident: its legacy (non-tiled) coefficient file has no block_active
table, so all blocks were kept. Fix: `mobygeom block-table --keep-buried`
(zeroed buried masks; the solver's own builder then keeps the core too).

## Workflow

```bash
./setup.sh              # STL + grid.h5 + 5-level block-table file (~1 h)
./run_sweep.sh          # aoa 0/4/8 sequential GPU runs + check_naca.py
./run_sweep.sh 4        # a single angle
```
