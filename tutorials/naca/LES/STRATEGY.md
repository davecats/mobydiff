# NACA 0012 AoA 5 — wall-resolved LES: strategy

STATUS (2026-07-29): strategy agreed with the user — **Re_c = 4e5**
(matching the RANS seed), WALE instead of ADM-RT and the zero-gradient
outlet accepted, the production run goes to a CLUSTER (this machine only
prepares: features, gates, case file, launch scripts). The Schlatter &
Orlu trip is PORTED verbatim from the validated boundary-layer
implementation in `~/Codes/mobydiff.bl` (see F2 below).
Seed case: the validated RANS tutorial `../rans/` (NACA 0012, AoA 5,
Re_c = 4e5, OpenFOAM-matched, C_L 0.514 / C_D 0.0130).
Discretisation template: Atzori et al., Flow Turb. Combust. 105 (2020),
Sect. 2.2 (`../atzori.pdf`) — the well-resolved-LES setup of Vinuesa et
al. (2018) / Negi et al. (2018) for a NACA 4412 at 5 deg.

## 1. What the paper prescribes (Sect. 2.2)

| Item | Paper value |
|---|---|
| Domain (streamwise x vertical x span) | 6c x 4c x 0.2c |
| Airfoil placement | LE 2c from the front boundary, TE 3c from the rear |
| Span | periodic |
| Far-field BCs (front/top/bottom) | Dirichlet velocity sampled from an auxiliary k-omega SST RANS (200c domain); < 1% off the LES mean |
| Outlet | Dong et al. (2014) outflow |
| Near-wall resolution (wall units) | tangential Delta_xt+ = 18, wall-normal Delta_yn+ = 0.64 (first point) to 11 (in the BL), span Delta_z+ = 9 |
| SGS | ADM-RT relaxation-term filter (Schlatter 2004) |
| Tripping | localized volume force (Schlatter & Orlu 2012), both sides at x/c = 0.1 |
| Transient protocol | init from the RANS field; polynomial order ramped P5 (4 flow-overs) -> P7 (2) -> P11, then 2 more before sampling |
| Statistics | 10-15 flow-over times (c/U_inf) for converged stresses |
| Their cost | ~220M grid points, ~1M CPU-h per 10 flow-overs (Cray XC40) |

The resolution criteria are wall-unit based (they come from Negi et al.'s
Re 4e5 LES-vs-DNS validation), so they transfer to any Re — the physical
spacings below are what changes.

## 2. Viscous scales, measured from our RANS case

From `../rans/cpcf_c11_aoa5_final_cf.dat` (Re_c = 4e5, turbulent after the
x/c = 0.1 trip): l* = nu/u_tau =

| station x/c | 0.05-0.15 | 0.2-0.4 | 0.4-0.6 | 0.6-0.8 | 0.8-1.0 |
|---|---|---|---|---|---|
| l*/c at Re 4e5 | 4.5e-5 | 4.7e-5 | 5.6e-5 | 6.6e-5 | 9.1e-5 |
| l*/c at Re 2e5 (cf ~ scaled) | ~8.5e-5 | ~8.8e-5 | ~1.1e-4 | ~1.3e-4 | ~1.7e-4 |

Design point: the strict fore/mid stations (l* = 4.7e-5 c at 4e5,
~8.8e-5 c at 2e5).

## 3. The two structural translations to mobydiff

**(a) Isotropic in-plane cells.** Our xz-quadtree AMR gives isotropic
cells per level in the chord-lift plane; tangential and wall-normal
directions share one spacing that rotates around the profile. The paper's
0.64(+)-normal x 18(+)-tangential anisotropy is not representable. The
binding constraint for a penalized IBM without a wall model is the
wall-normal one: the first fluid cell centre must sit in the viscous
sublayer (the RANS campaign already established the y+ ~ 1.5-3 practice;
C10 at y+ ~ 5 paid +18% drag). So the finest band is chosen by
Delta+ ~ 3 (first centre y+ ~ 1.5), which over-resolves the tangential
target 18 by ~6x — the price of isotropy, and a step toward the later DNS.
The BL bulk (above the wall band) only needs Delta+ <= 11, i.e. two
levels coarser.

**(b) Span multiplies everything.** In refine_dims = xz the span keeps a
single global uniform line: Delta_y+ = 9 applies to the entire domain,
and every in-plane cell is multiplied by ny. This makes the paper's SMALL
6c x 4c box structurally necessary (on the RANS 128c x 96c box the far
field alone would cost ~44M cells at base level). Adopting the paper's
domain is therefore both faithful and required.

## 4. Reynolds number: Re_c = 4e5 (DECIDED, user 2026-07-29)

The LES matches the validated RANS seed exactly (direct overlay of
C_p/C_f/C_L/C_D and BL profiles on the `../rans/` OpenFOAM-anchored
case; the RANS field supplies far-field BCs and the initial condition
with no re-run). The production run goes to a cluster; the size below
is the resource request, not a local-feasibility constraint.

Sizing at Re 4e5 (from the Sect. 2 l* table):
- finest wall band Delta = c/8192 = 1.22e-4 c (L11 off Delta0 = c/4):
  Delta+ = 2.7 fore, first cell centre y+ ~ 1.3;
- BL bulk L9 (Delta+ 10.8 fore — the paper's <= 11 in-BL target);
- span ny = 480, Delta_y = 4.17e-4 c, Delta_y+ = 9.3/9.0 at the
  fore/mid stations (the paper's 9);
- estimated total ~6e8 cells (in-plane ~1.25e6 across levels x 480
  span planes) — the paper had 2.2e8 at Re 2e5 with anisotropic
  spectral elements; isotropy in-plane + doubled Re explain the ratio;
- memory ~100 GB of device-resident state -> e.g. 2 x 80 GB or
  4 x 40 GB cluster GPUs (3D MPI decomposition + Z-order block split
  are already the production path);
- dt ~ 5e-5 (CFL at Delta11, suction-peak |u| ~ 1.7); a 15-flow-over
  statistics campaign is ~3e5 steps, plus the staged transient.

## 5. Proposed case definition (Re 4e5)

- Geometry: NACA 0012 (closed TE, `make_airfoil_stl.py`, span y,
  Ly = 0.2c), AoA 5 via inlet-flow angle as in the RANS case.
- Domain: lx = 6, lz = 4, ly = 0.2; nose at (x, z) = (2.0, 2.0)
  (LE 2c from inlet, TE 3c from outlet, mid-height).
- Grid: base Delta0 = c/4 -> nx = 24, nz = 16, ny = 480; nb = 8
  (base lattice 3 x 2 blocks x 60 y-tiles); refine_dims = xz,
  refine_body = true, refine_levels = 11, **keep_buried** (the runtime
  penalization C_L(t)/C_D(t) is the convergence monitor — the A2
  landmine: removal invalidates it; CV forces stay as the cross-check
  on the mean).
- Protection boxes (concentric, wake-skewed, RANS-case style): L9 over
  the full BL envelope + near wake, L8 wake to ~1c behind the TE,
  grading outward; L10/L11 from refine_body + cascade.
- BCs: x_min inlet + z_min/z_max Dirichlet velocity from the sampled
  RANS far field (feature F1 below); x_max the existing validated
  outlet (BC_OUTFLOW + Dirichlet p; the paper's Dong outflow is an
  accepted difference — monitor wake exit in the pilot); y periodic.
- SGS: `[turbulence] model = les`, `[les] model = wale, ibm_aware =
  true` (the combination validated in validation/channel_interface/les
  and les_ibm, incl. across 2:1 interfaces). The paper's ADM-RT
  relaxation filter is a documented modelling difference.
- Trip: the ported Schlatter & Orlu forcing at x/c = 0.1, both sides
  (feature F2/F2b below), same station as the RANS kpin/ktrip and the
  paper.
- Time: dt = dtmax = 5e-5 (CFL 0.8 at Delta11 with |u|max ~ 1.7);
  niter = 18, accel = chebyshev (the RANS-case projection settings —
  the pn-cleanliness argument for CV borders carries over).

## 6. Solver/tooling gaps (each its own gated increment, dormant paths bit-exact on the 7-case suite CPU+GPU)

- **F1 — far-field Dirichlet from file.** Extend the A0 per-point BC
  machinery (`_profile = parabola` already stores per-point values) with
  `_profile = file` reading (u,v,w) per boundary point for
  inlet/patch faces. Gate: a uniform-freestream file reproduces the
  existing freestream case bit-exact; the RANS-sampled ring recovers the
  imposed circulation-consistent angle.
- **F2 — trip forcing: PORTED (2026-07-29).** `[force] type = trip`
  brought over verbatim from the validated boundary-layer implementation
  in `~/Codes/mobydiff.bl` (bodyforce.f90 module wholesale + the
  dns_type trip_* fields + the [force] trip_* config keys; call sites
  were already identical). Schlatter & Orlu (2012): f_v = amp *
  exp(-((x-x0)/lx)^2 - (y/ly)^2) * g(z,t), g a unit-rms nmodes-Fourier
  random spanwise function regenerated every ts with a C^1 smooth step
  in time; deterministic seed, host random walk + device fill kernel.
  Gate: dormant bit-exact 7-case suite (CPU done locally; GPU leg
  pending a free GPU — all hosts busy with production runs / the local
  3060 driver awaits reboot).
- **F2b — trip orientation for the airfoil.** The ported trip is in
  BL axes: it forces the y-component with the Gaussian in (x,y) about
  y = 0 and the random modes along z. The quasi-2D airfoil runs span=y
  (the only quadtree is refine_dims = xz), so the trip must force w
  with modes along y, centred on each airfoil surface (two strips).
  Increment: config-selected trip component/span axis + repeatable
  strip centre (x0, z0), DEFAULTING to the .bl convention so existing
  inis and the ported arithmetic stay untouched. Gate: axis-permuted
  twin == the .bl orientation under coordinate swap, + dormant
  bit-exact.
- **F3 — online statistics.** Device-resident running accumulators
  (mean u,v,w,p,nut + second moments u_iu_j, at least) written in the
  block-table snapshot layout; avoids the ~5 GB/snapshot flood that
  offline averaging would need. Gate: equals offline snapshot averaging
  to round-off on a short run.
- **T-BC — boundary sampler.** Python tool: read a converged RANS
  snapshot (chunked, compare_snapshots.py machinery — NOT
  compare_fields.py, the L11 OOM landmine), interpolate (u,w) onto the
  LES boundary rings, write the F1 file.
- **T-IC — cross-grid restart.** Generalize interp_restart.py to
  resample a block-table snapshot from the RANS grid/domain onto the
  LES layout (different domain, offsets, levels).
- **T-POST — surface_cp_cf LES mode** (mean-field input, no k -> pure
  anchored fit); cv_forces on the mean field + Reynolds-stress border
  flux from F3; spanwise two-point correlations (span-sufficiency
  check); BL profile extraction at x/c stations.

## 7. Execution phases (simulations run on the CLUSTER; this machine
## prepares — features, gates, case file, inis, launch scripts)

1. **P0 — features**: F1, F2b, F3 landed and gated individually
   (F2 already ported). Local CPU gates + GPU gates when a GPU frees.
2. **P1 — boundary + IC extraction**: T-BC samples the converged
   `../rans/` Re 4e5 field on the LES boundary rings; T-IC resamples
   it onto the LES layout. No new RANS run needed.
3. **P2 — prepare**: STL (span 0.2c) + moby_prepare on the 6x4x0.2
   domain with keep_buried; gates: 1==4 ranks identical file, mask/
   dwall sanity, leaf count vs the ~6e8-cell estimate (fixes the
   cluster resource request).
4. **P3 — pilot** (ny = 64, finest L10): full pipeline shakedown on
   the cluster — trip amplitude calibration, transition location,
   outlet behaviour with a turbulent wake, force monitor, stats
   module, s/step measurement. Prepared here as a ready-to-submit
   case; small enough for a single cluster GPU.
5. **P4 — production, staged like the paper**: IC from P1 at an
   L9-capped grid (~4 flow-overs) -> deepen L10 (~2) -> L11 (~2, via
   T-IC/interp_restart) -> statistics over 10-15 flow-overs (~3e5
   steps at dt 5e-5). Multi-GPU MPI; rank count fixed by the P2 leaf
   table and the cluster's per-GPU memory.
6. **P5 — validation report**: C_p/C_f and C_L/C_D overlaid on the
   `../rans/` case, its OpenFOAM reference and XFOIL; u+/
   Reynolds-stress BL profiles at x/c stations; spanwise two-point
   correlations (span sufficiency); convergence appendix (paper
   style).

## 8. Known risks / accepted differences

- Staircase IBM roughness: laminar Cf ran 1.5-1.7x Thwaites in the RANS
  campaign; the trip at x/c = 0.1 keeps the laminar run short, but the
  LES transition may be staircase-assisted — the pilot measures it.
- Isotropic wall band over-resolves tangentially (6x vs the paper) —
  cost, not accuracy, and it shortens the later DNS step.
- Const-1/2 2:1 interface transfer is mildly dissipative for the
  smallest scales crossing to coarse levels (~5% v'/w' energy, the
  validated residual); BL-bulk levels are chosen so the interfaces sit
  outside the near-wall production region.
- Outlet: zero-gradient + p-Dirichlet instead of Dong; the wake crosses
  it 3c downstream. Pilot gate. (User-accepted 2026-07-29.)
- WALE instead of ADM-RT relaxation filtering. (User-accepted
  2026-07-29.)
