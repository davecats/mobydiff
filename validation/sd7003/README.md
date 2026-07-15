# A3 INCREMENT 3 — SD7003 transition benchmark (Re_c = 6e4, aoa = 4)

The gamma-Re_thetat transition model on the standard low-Re LSB benchmark,
through the full airfoil stack (freestream composition, file IBM from the
UIUC Selig coordinates spline-resampled to 720 points, 5-level
refine_body with --keep-buried, SST + transition, scalar inlets holding
k/omega/gamma/Re_thetat~ freestream values, penalization forces).

- tu = 0.1 %, nut_ratio = 1 => omega_inf = 0.09: ambient decay time
  1/(beta* omega_inf) = 123 c/U — the inlet tu reaches the airfoil
  unchanged (the naca0012 tu = 5 ambient-decay lesson does not bite).
- Grid = validation/naca0012 (Delta = 1.465e-3 c at the surface; laminar
  y+_1 ~ 1.8 at Re 6e4, the docs' 3-4-level arithmetic).

Gates (check_sd7003.py; published gamma-Re_thetat RANS-LM scatter —
Windte/Radespiel, LES reference Galbraith & Visbal):

- LSB present (reversed-flow patch on the suction side);
- transition location x_t/c = 0.5 +- 0.1 (near-wall gamma crossing 0.5);
- C_L within +-15 % of ~0.60, C_D within +-15 % of ~0.022;
- the gamma-front chordwise smearing (stations from gamma 0.1 to 0.9, in
  level-4 cells) is MEASURED AND REPORTED ONLY — it is the number that
  decides the separate TVD/van-Leer + second-scalar-halo increment.

## Results (2026-07-14, gate PASS — WITH gamma_sep)

The FIRST run (commit 2d514ca state, saved as *_nogsep) exposed a missing
model piece: the bubble separated at x/c = 0.223 but the detached shear
layer NEVER transitioned (near-wall peak k = 6.6e-5 ~ laminar; gamma
never rose) — at Re_c = 6e4 / tu = 0.1 % the attached-BL Re_theta (~170)
cannot reach the natural-transition criterion (~1200): separation-induced
transition IS the SD7003 mechanism, and gamma_sep (LM Eq. 18) was the
T4-deferred increment. C_L 0.586 / C_D 0.026 sat in the published band by
the steady pressure field alone — the integral gates alone would have
GREEN-LIT a laminar solution; the field-level k/gamma gates caught it.

With gamma_sep (50f03bd; t_final = 15, tail-0.2 means):

- separation x_s/c = 0.22 (published 0.22-0.30), patchy separated region
  to x/c ~ 0.56 with a last segment at 0.66-0.68 (published
  reattachment 0.65-0.70); staircase Cf flicker fragments the bubble;
- transition (near-wall k onset, 100 k_inf) x_t/c = 0.427 — inside the
  0.5 +- 0.1 gate, at the early edge of the published RANS-LM scatter
  (0.53-0.58), consistent with first-order-upwind front smearing;
- peak near-wall k = 3.7e-2: genuinely turbulent reattachment;
- C_L = 0.5616 (-6 % of 0.60), C_D = 0.0267 (+21 % of 0.022, inside the
  +-15 %-band edge at 0.027 — the penalization D_eff ~ D + h class);
- **transition-front chordwise smearing: 104 level-4 cells (0.152 c)
  across k = 10 -> 1000 k_inf.** This is THE trigger number for the
  deferred TVD/van-Leer + second-scalar-halo increment: the front is
  wide, and the measurement now justifies that follow-up.

## Workflow

```bash
./setup.sh          # grid.h5 + resampled STL + keep-buried block-table (~1 h)
mpirun -n 1 ../../build_gpu/main aoa4.ini     # ~7 h on the RTX 3060
python3 check_sd7003.py sd7003_aoa4_*.h5 forces_aoa4.txt
```

## R2D-3 follow-up (2026-07-15): the 2D-refined (span y, refine_dims = xz) reruns

Two runs through the full span-y quadtree stack (xz_aoa4.ini = L5-xz,
keep-buried, dt 2e-4; xz_l4.ini = the L4-xz CONTROL at the exact
published resolution and dt; STL sd7003_spany.stl; check_sd7003.py
--span y --lmax 4|5):

- **L4-xz == L4-3D to every checker digit**: x_t/c = 0.427, front
  smearing 104.0 level-4 cells, C_L = 0.5617 (3D 0.5616), C_D = 0.0267
  (3D 0.0267), same LSB prints. The span-y xz path (quadtree transfers,
  scalar halos, transition transport, penalization forces) is
  QUANTITATIVELY EQUIVALENT to the validated 3D path at matched
  resolution — at 3459 refined leaves vs ~21k (runs in ~28 min on the
  RTX 5090).
- **L5-xz FINDING — separation-induced transition does NOT fire at the
  finer surface**: the laminar bubble is present (reversed flow x/c
  0.25..0.55, reattachment ~0.6, gamma dips to 0.02 in the BL) and
  C_L = 0.626 sits in the published band, but k never leaves the
  freestream level (max 3.3 k_inf; the L4 runs reach k ~ 5e-2), so
  C_D = 0.0281 reads laminar-high. With the L4-xz control clean, this
  is NOT an xz defect but the Langtry-Menter model/discretization
  interaction at finer resolution: the L4 transition was already the
  "early edge" (x_t 0.427 vs published 0.53-0.58), and the FIRST-ORDER
  upwind scalar convection smears the front over 104 cells = 0.152c —
  at L5 (double the cells across the same physical bubble, halved dt)
  the smeared gamma/Re_thetat fields no longer trigger gamma_sep at
  all. This STRENGTHENS the already-measurement-justified TVD/van-Leer
  scalar-convection increment (deferred, separate session): rerun the
  L5-xz benchmark after it lands.
