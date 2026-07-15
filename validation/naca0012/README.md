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

## LE ripple follow-up (2026-07-14): the runs ARE second-order — verified

Diagnosis of the USE_IBM_SECONDORDER / coefficient-file question (user
prompt), correcting the earlier "binary mask" statement:

- `USE_IBM_SECONDORDER` guards only the ANALYTIC `set_ibm_coeff` path
  (ibm.f90): fluid cells adjacent to the body get the graded coefficient
  sum((d0-d)/d)/d0^2 / re per solid neighbour, d from bisection to the
  surface — the sharp-interface second-order Laplacian correction (cf.
  arXiv:2506.14328): the wall-distance-weighted coefficient corrects the
  viscous stencil so the effective no-slip plane sits AT the surface.
- mobygeom's `stl_coeff_tile_from_mesh` — shared by `stl-ibm-coeff`
  (legacy cylinder file) and `block-table` (the A3 files) — writes the
  SAME graded values unconditionally (`stl_segment_distances` is the STL
  bisection analog). Verified in the files: ibm_coeff_n0012.h5 carries
  graded cells in the whole first fluid ring (median ~4 = the half-cell
  crossing value at Delta4, range 1.5e-2..5.4e4; 190 graded cells per LE
  block), ibm_coeff_re40.h5 carries 920 per component. The solver reads
  file coefficients verbatim.
- CONCLUSION: every A3 run (and the validated cylinder) already used the
  second-order body description; the flag is moot for file-based IBM.
  The LE ripple fan (~1 % rms u upstream of the nose, 40 % omega
  cell-to-cell at 2 cells decaying to 0 at 40) is the RESIDUAL of the
  second-order scheme at a strongly curved staircase surface — no
  regeneration or rerun warranted; the SD7003 transition gates remain
  the sentinel, and a calibrated smoothed-mask/Brinkman treatment stays
  the (post-A3) escalation increment if they fail.

## Results (2026-07-14, gate PASS)

t_final = 10 (tail-0.2 means), keep-buried coefficient file, dt = 4e-4:

| alpha | C_L     | C_D    |
|-------|---------|--------|
| 0     | -0.0013 | 0.0186 |
| 4     | +0.3838 | 0.0209 |
| 8     | +0.7447 | 0.0290 |

Lift slope 0.0932/deg = 85 % of 2pi (the expected Dirichlet-far-field
blockage reduction at 12c); |C_L(0)| = 0.0013 (grid/staircase asymmetry
negligible); C_D(0) = 0.0186 vs XFOIL fully-turbulent ~0.013-0.015 (the
first-order penalization D_eff ~ D + h class + y+_1 ~ 2-4 resolution).
aoa4 rerun history: the first (removed-core) run gave C_L = 0.018 — see
the keep-buried finding above. Runs executed across three hosts (local
RTX 3060 / istmcetus A6000 / istmcorax RTX 5090 — cc120 build dir
build_gpu_corax).

## LE fan root cause (2026-07-14, user OpenFOAM counter-evidence): the
## CENTRAL-SCHEME cell-Reynolds parasite, not the staircase

The user implemented the same graded IBM in OpenFOAM and saw NO fan —
refuting "inherent staircase residual". Full diagnosis (evidence in
slice_le_new.npz / slice_fanre1000.npz / the cylinder snapshot):

- the fan's dominant ripple wavelength is 2.44 cells with 68 % of its
  energy near the grid Nyquist — the stationary parasite of CENTRED
  convection at high cell-Reynolds number (u Delta4/nu = 146 here);
- upstream decay: cylinder Re 40 (cell-Re 1.25) ripple dies within 8
  cells (0.0038 -> 0.0004 -> 1e-5); airfoil Re 1e5 fan persists 64+
  cells (0.023 -> 0.016 -> 0.010 -> 0.005 -> 0.003);
- the full 2x2 matrix (user-requested deconfounding of Re vs turbulence
  model; upstream strip rms at 4/8/16/32/64 cells):
  Re 1e5 + rans   0.0228 0.0160 0.0101 0.0052 0.0026  (fan)
  Re 1e5 + none   0.0271 0.0208 0.0148 0.0100 0.0042  (fan, slightly
                  STRONGER — the RANS ambient nut ~ 10 nu was damping it)
  Re 1e3 + none   0.0025 0.0003 ...                    (no fan)
  Re 1e3 + rans   0.0025 0.0003 ...                    (no fan, identical)
  The fan correlates purely with the cell-Reynolds number and is
  INDEPENDENT of the turbulence model (the 1e3 rows' 32-64-cell
  residual ~0.002 is large-scale starting transient, ~41 % short-wave);
- niter is NOT the mechanism (6 vs 12 clean-p: fan rms identical to 4
  digits; more iterations only damp the TEMPORAL ringing).

The IBM forcing is the (rough) SEED; the centred scheme at cell-Re >> 2
is the AMPLIFIER that radiates it upstream; OpenFOAM's upwind-biased/
limited div schemes damp the parasite. Mitigation plan (user decision
2026-07-14, above TVD): (R1) refinement ground truth, then a band-filter
or band-dissipation correction pass. A one-sided convective-stencil
variant was implemented and REVERTED (user decision: too invasive in the
hottest kernel; also its first form — substituting the ADVECTING pair
members — blew up in 1500 steps because it breaks discrete continuity of
the advecting flux; only transported-value substitution is admissible).

## R1 refinement ground truth (2026-07-15): one level collapses the fan

The ring cells' local cell-Re scales QUADRATICALLY with resolution
(cell-Re_ring ~ S Delta^2/nu ~ 60 / 7.6 / 1.9 at L4/L5/L6). Fan strip
rms at FIXED physical distances upstream of the nose (y band
[5.95,6.05], high-pass 5 native cells; R1 files are remove-buried:
forces meaningless by construction, fan metric only; L5 run t_final 1.5
on the RTX 5090, 65094 leaves):

  case            0.006c   0.012c   0.023c   0.047c   0.094c  shortwave
  L4 RANS (base)  0.0228   0.0160   0.0101   0.0052   0.0026     68 %
  L5 RANS (R1a)   0.0015   0.0016   0.0021   0.0027   0.0025     33 %

One refinement level cuts the near-nose fan 15x; the flat ~0.002
residual is the shared large-scale transient floor, not a fan (both
cases meet at the far strip). The parasite is effectively dead already
at ring cell-Re ~ 7.6.

R1 CONTROLS (user-requested; all at t = 1.5 for fair comparison):

  case                      0.006c   0.012c   0.023c   0.047c   0.094c
  L4 dt4e-4 (t-twin)        0.0172   0.0113   0.0057   0.0022   0.0020
  L4 dt6e-5 (dt control)    0.0179   0.0114   0.0055   0.0024   0.0021
  L4 + fan-box (fine MEDIUM,
    same L4 ring; box
    4.35-4.55 x 5.92-6.08)  0.0170   0.0113   0.0053   0.0008   0.0000
  L5 (fine ring + medium)   0.0015   0.0016   0.0021   0.0027   0.0025

Verdicts: (1) dt INNOCENT — a 6.7x smaller step changes every strip by
< 4 % (the gap to the t = 10 baseline is fan growth in time, not dt);
(2) the near fan (<= 0.023c) is controlled by the RING SEED, not the
medium: refining the medium alone changes nothing, refining the ring
collapses it 11x — the cell-Peclet-of-the-medium hypothesis is refuted
for the core fan; (3) the FAR strips were largely LEVEL-INTERFACE
artifacts, not parasite: with the whole analysis window at one level
they drop to exactly 0.0 (the twin's and L5's far-strip floor ~0.002
coincides with L4->L3->L2 interfaces inside the window — the earlier
"fan reaches 64+ cells" was partly this). (4) L6 (ring cell-Re ~ 1.9)
BLEW UP under model = none (under-resolved DNS at Re 1e5: the fine band
resolves real shear-layer instabilities the coarse wake cannot carry;
RANS OOM'd at 121M cells on 49 GB) — dropped as redundant: L5 already
sits at the transient floor. mobygeom block-table now combines
--refine-box WITH body classification (the solver ini must carry the
same [blocks] refine key for the builder cross-check).

## Band filter ([ibm] band_filter, 2026-07-15): the production option

Per the user decision (refinement = the reference answer; an optional
near-zero-cost-when-off filter for everything else): a 3-point low-pass
qs += (theta/4)(q_{i-1} - 2 q_i + q_{i+1}) per direction, applied ONLY on
a compressed list of near-body fluid DOFs (band_width = 3 cells; solid
marker from the coefficients, dilated with a halo exchange per pass so
the band crosses block boundaries exactly; per-direction bits exclude
any solid read — the v1 lesson; physically pinned faces skipped).
qs-only operator splitting (NO oldrhs term), increment x mu (the A2
force bookkeeping stays exact). OFF = the pass is never called and
nothing is allocated or mapped: bit-exact and zero cost by construction
(7-case suite ALL BIT-EXACT).

Fan bench (L4 grid, dt 4e-4, t = 1.5, RANS; 1587456 band DOFs = 2.1 %):

  case                   0.006c   0.012c   0.023c   0.047c   0.094c
  L4 no filter           0.0172   0.0113   0.0057   0.0022   0.0020
  L4 + filter th=0.5     0.0061   0.0025   0.0012   0.0019   0.0020
  L5 ground truth        0.0015   0.0016   0.0021   0.0027   0.0025

At/below the refinement ground truth from 0.012c outward; the 0.006c
residual (4x the L5 floor) sits at the band edge (strip at 4 cells,
band 3). STABILITY: cells filtered in all three directions amplify by
1 - 3 theta => theta >= 2/3 is UNSTABLE (theta = 1.0 blew up within 40
steps, C_L -> 1e147); config hard-errors above 0.6. COST/BIAS: les_ibm
channel 200 steps stable (73728 band DOFs); the Re 40 cylinder
(D/h = 32, 6912 band DOFs) pays +3.9 % C_D (1.7587 vs 1.6925 same-t) —
the filter thickens the coarse near-wall profile; on fine-ring grids
the bias should shrink with the band/BL ratio. RECOMMENDATION: default
OFF; enable for LE-sensitive production runs at moderate resolution;
prefer one more refinement level for benchmark-grade results. Follow-up
candidates: filtered-NACA force check, width/theta tuning, per-DOF
theta scaling by active direction count.

PROCESS LESSON: remote build dirs must be rebuilt too — the first
"filter" bench ran corax's stale binary, which silently ignored the
unknown ini keys and reproduced the twin bitwise.

## Workflow

```bash
./setup.sh              # STL + grid.h5 + 5-level block-table file (~1 h)
./run_sweep.sh          # aoa 0/4/8 sequential GPU runs + check_naca.py
./run_sweep.sh 4        # a single angle
```
