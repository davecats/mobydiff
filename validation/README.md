# Test cases — index

> **PLAN (2026-09-26; the handout is `docs/next_session_revalidation.md`,
> 2026-10-21):** after a round of code simplifications and a thorough
> code analysis (done 2026-10-21), every gate, validation case and tutorial
> listed here is re-gone-through ONE BY ONE (re-run, re-measure, update or
> retire). Until a case's row in that handout's ledger says DONE, the numbers
> recorded in its README are what the ledger says they are: niter-6 values
> for the long Chebyshev cases (all Chebyshev inis moved to `niter = 12` on
> 2026-09-26), pre-F1/F7/outlet values for most, and the inis that declare
> plain Jacobi still carry `niter = 6` (moved when their case is rerun).

The repository has exactly two homes for runnable cases:

- `tutorials/` — user-facing cases (how to set up and run a real flow);
- `validation/` — every regression / verification gate (this directory).

There is no separate `smoke/` or `profile_*/` directory any more (cleaned up
2026-09-26; everything removed is recoverable from git history, see the end of
this file). Each subdirectory has its own README with commands and recorded
results; this file only says WHAT each case exercises.

## The standard regression suite

Since numerics review step 7 every run reads a CASE FILE (`[case] file`,
default `<field_prefix>.case.h5`) that `moby_prepare` writes from the same ini
-- the solver never builds one and refuses a missing or stale file. Every
driver here therefore runs `tools/prepare_if_missing.sh <ranks> <solver> <ini>`
before a solve: it resolves the file the solver will read and prepares it
only when it does not exist (the inis that name a committed or generated
reference file are never recomputed), with the moby_prepare next to the
solver binary and the SAME rank count as the solve (an unset `[blocks] nb`
is one block per prepare rank). A suite run leaves one `<prefix>.case.h5`
per output prefix beside the snapshots; the `les_ibm` case reads the
committed `ibm_coeff_case.h5` (the exact conversion of the legacy
`ibm_coeff.h5`, which left the tree on 2026-09-30).
A driver that derives a solve ini from a case ini points it at a case file
with `tools/set_case_file.sh <ini> <case.h5>` (writes `file =` into `[case]`,
adding the section when absent); `[ibm] coeff_file` is no longer a key.

The cases below are the "7-case suite" (+ the scalar legs) used for every
bit-exactness gate (`-Mnofma` / `-gpu=nofma` on both sides, `tools/h5maxdiff`
or `tools/compare_fields.py --tolerance 0`). Driver:
`tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange/run_portgate.sh`
(or `run_mapgate.sh`); `validation/scalar/run_bitexact*.sh` on a workstation.

| case | what it covers |
|---|---|
| `tutorials/min_channel/input_gpu.ini` (1 and 4 ranks) | blocks + 2:1 interface + Chebyshev-Jacobi projection, multi-rank exchange |
| `validation/beltrami/slab_y.ini` | exact NS solution through a y-normal 2:1 slab |
| `validation/channel_interface/les_ibm/channel_ibm.ini` (+ `_refine`) | file-based IBM + WALE LES (± refine_body / 2:1) |
| `validation/rans_sst/turb180.ini` | k-omega SST, resolved walls |
| `validation/rans_sst/wf180_y30.ini` | SST wall functions |
| `validation/rans_sst/lam30t.ini` | gamma-Re_theta transition |
| `validation/scalar/{conduction,prsweep,wave}.ini` | passive scalar transport (name `s1` explicitly: h5maxdiff's default list skips scalars) |
| `validation/scalar/{ibmwavy,ibmwavyr}.ini` | scalar + analytic IBM (dirichlet/adiabatic), ± refine_body |
| `validation/scalar/{mixed,mixedles}.ini` | scalar + analytic IBM, one scalar in EACH wall mode (dirichlet/adiabatic/conjugate), ± WALE: both transport kernels in one run (2026-10-20) |

## validation/

| directory | exercises | driver | status |
|---|---|---|---|
| `blasius/` | laminar ZPG Blasius boundary layer on a one-sided natural (stretched) y line, Blasius u+v inlet, Dirichlet-p outlets at x_max and the top: the analytic similarity field is minted into a restart and must be HELD to t = 2000 (theta, H, u/v profiles vs the independent ODE); the laminar stretched-line gate the F1 flux-form stencil was validated on | `run_blasius.sh`, `compare_blasius.py` | re-validated 2026-10-21 |
| `beltrami/` | 3D Beltrami/ABC exact decaying NS solution (`[flow] initial = beltrami`): 2nd-order convergence (`uniform.ini`, nx 32/64/128) and the 2:1 interface artifact for x/y/z slabs (`slab_{x,y,z}.ini`) and a 3D patch with edges+corners (`refined.ini`) | `run_beltrami.sh`, `tools/check_beltrami.py` | re-validated 2026-10-22 |
| `interface_decay/` | stability of the 2:1 interface: white noise on a 3D refined patch must decay everywhere (catches interface-localized growth that smooth-flow gates miss). ~35 s on 4 CPU ranks | README, `tools/check_interface_decay.py` | re-validated 2026-10-21 |
| `block_nb/` | per-direction `[blocks] nb`: single-level nb-independence (bit-exact across tilings and rank counts) and uniform oblique flow through a 3-level patch with non-cubic nb; `beltrami_rect.ini` and its red-black twin `beltrami_rect_redblack.ini` (a real flow through the patch, 2:1 faces of unequal extents; the suite's only red-black case) are cases 8 and 9 of `scalar/run_bitexact.sh` | `run_gates.sh` | re-validated 2026-10-21 |
| `refine2d/` | `refine_dims = xz` quadtree: all-refined == doubled-resolution twin (`allref_xz`/`twin_xz`), uniform flow through xz patches (`uniform_xz`, `unibody_xz`), Beltrami order pair (`bp_xz_32/64`), RANS scalars across xz interfaces (`rans_xz`) | `run_gates.sh`, `setup.sh` | re-validated 2026-10-21/22 |
| `redblack_interface/` | `[pressure] solver = redblack` across a 2:1 interface (R1); 1 rank == 4 ranks | `run_gates.sh` | re-validated 2026-10-21 |
| `multilevel_body/` | 3-level `refine_body` on a cylinder: uniform oblique flow EXACT across all interface levels (zero-coef twin), per-level dwall vs an exact prism reference | `setup.sh`, `run_gates.sh` | re-validated 2026-10-21 |
| `channel_interface/` | turbulent Re_tau 180 channel across the 2:1 interface: wall bands (`reference`/`uniform128`/`refined_y55`/`refined_y110`), developed statistics (`developed/`), edge/corner core patch (`core_patch/`), LES across refinement (`les/`; its `run_wall_jump.sh` is the one-step CPU gate of a level jump that touches a wall, seen by the SGS model), LES<->IBM coupling (`les_ibm/`, suite case). Long GPU runs | `run_validation.sh`, per-subdir READMEs | re-validated 2026-10-21/22 |
| `freestream/` | inlet/outlet patch types + Dirichlet-pressure outlet, the outlet face PREDICTED (2026-10-01): oblique freestream exact, inflow/outflow Poiseuille vs periodic (last-cell p = G dx/2), Lamb-Oseen vortex exit, 1 == 4 ranks, config contradictions; the outlet on a LOW face (a mirrored pair, Poiseuille and vortex run in -x) equal to the x_max image to round-off; uniform flow recovered from rest through a refined patch touching an outlet; restart == continuous at max_abs 0; bit-exactness on the outlet cases | `run_gates.sh`, `run_restart_outlet.sh`, `run_bitexact_outlet.sh` | re-validated 2026-10-21 |
| `cylinder/` | `[case] name = airfoil` on a cylinder: control-volume C_L/C_D (Re 40 steady drag, Re 100 shedding/Strouhal), empty-domain zero force; the stored pressure and the forces at two step sizes and three projections (`run_outlet_mode.sh`), St against the outlet distance and the blockage (`shedding_fit.py`) | `setup.sh`, `check_cylinder.py` | re-validated 2026-10-22 |
| `rans_geometry/` | SST wall distance / y_eff / wall-cell classification: file IBM (flat, ± refine_body) vs closed form, analytic wavy wall (± refine_body) vs scipy | `setup.sh`, `check_rans_geometry.py` | re-validated 2026-10-21 |
| `rans_sst/` | k-omega SST (laminar decay, Re_tau 180/395, IBM channel, refined bands), wall functions (y+ sweep, IBM), transition (lam30t, laminart, turb180t) | `run_gates.sh`, `check_gates.sh` | re-validated 2026-10-21 |
| `body_outlet/` | an immersed body crossing an inlet AND an outlet plane (B1/B2, 2026-10-05): Poiseuille between immersed slabs in an inflow/outflow channel vs its periodic twin (outlet-face profile, p pinned at G dx/2, solid faces on the plane at round-off, mass balance), the outlet on a LOW face (mirror), stream along y and z (permutation), the body ending one cell before the plane / starting at it, 1 == N ranks, restart == continuous, a zeroed-coefficient twin == body-free at max_abs 0; since 2026-10-06 the `pocket` group: a rough wall's trough ON the outlet plane (the laminar form of the production rough case's blow-up, stable only with the total-pressure datum on backflow) | `run_gates.sh` | re-validated 2026-10-21 |
| `bodyforce/` | the volumetric body force (2026-10-08): the correction kernel visits the force's SUPPORT only (the trip: v on its block list; every other source: the whole domain). `trip_small.ini` (a 512 x 176 x 48 boundary layer, trip ON on 7 of 32 blocks) REF == NEW at max_abs 0 and its trip-off twin must differ (a vacuous gate once ran the GENERIC case, which is why a foreign `[case.<x>]` section is a config error now); `force_sine.ini` (generic `[force]` sine profile on three components) REF == NEW | `run_gates.sh` | re-validated 2026-10-21 |
| `rans_inlet/` | RANS scalar freestream values at declared inlet faces, 1 == 4 ranks | `run_gates.sh` | re-validated 2026-10-21 |
| `iddes/` | IDDES (SST + WALE blend): developed channel log layer, fd shielding profile, fd_force limits, IBM stability | `run_gates.sh`, `check_gates.py` | re-validated 2026-10-22 |
| `courant_limit/` | the convective limit of the RK3 step, bracketed: a uniform flow along one / two / three directions in a periodic box, seeded with a 1e-6 perturbation, a fixed step just below and just above `dt x sum_d |u_d|/h_d = sqrt(3)`; `cflmax` bounds that sum, so the same value is stable or not for the three flows (1.60 / 1.85) | `run_gate.sh` | re-validated 2026-10-21 |
| `diffusion_limit/` | the explicit-diffusion limit of the RK3 step, bracketed: decaying Beltrami flow at Re = 1 in a periodic box on grids with three / two / one fine direction, a fixed step just below and just above `dt x sum_d nu/h_d^2 = 0.628` (decay below, blow-up and the solver's one-time warning above); `pecletmax` bounds that sum, so the same value is stable or not on all three grids (0.60 / 0.66) | `run_gate.sh` | re-validated 2026-10-21 |
| `penalization/` | the penalization factor in isolation (numerics review F7): uniform velocity + uniform body force + uniform moderate coefficient in a periodic box integrates `du/dt = f - lambda u`, and a Dirichlet-penalized uniform scalar its twin `ds/dt = S - lambda (s - s_body)`; the rational third-order factors must converge at order >= 2.7 and reach the steady fixed point to round-off (pass an older binary as `REF`: implicit Euler reads order 1, the step-9 exponential round-off for the velocity) | `run_decay.sh` | re-validated 2026-10-21 |
| `prepare/` | `moby_prepare`: analytic geometry bit-exact vs inline (P0: wavy, wavy_refine, wavysolid), STL geometry vs mobygeom references + periodic shift invariance (P1; the flat-slab references come from `rans_geometry/setup.sh`), big geometries (P1b: now the sailplane leg only), the configured analytic wall drives the classification (gate 4 of `run_gates.sh`, `check_eggclass.py`, 2026-10-06); and the block order -- a legacy-order case file and a re-prepared (minimum-surface order) one give identical fields, tiles and restarts across the two orders, the solver's partition line equals `tools/partition_analysis.py` | `run_gates.sh`, `run_gates_stl.sh`, `run_gates_big.sh`, `run_gates_order.sh` | re-validated 2026-10-22 (P1b retired) |
| `scalar/` | passive scalars S0-S5a: uniform preservation, conservation, conduction, MMS, Pr sweep, determinism, LES/SST closures, IBM dirichlet/adiabatic, heated cylinder, in-solver statistics, thermal wall function (domain + IBM walls), `count = 0` bit-exactness | `run_gates*.sh`, `run_bitexact*.sh` | re-validated 2026-10-22 |
| `conjugate/` | conjugate heat transfer at the IBM interface (C1-C3, F1-F5): multi-material slabs, contact resistance, oblique/curved interfaces, transient capacity, Nusselt diagnostic, channel with conducting walls, annular pipe geometry | `run_gates_c{1,2,3}.sh`, `run_gates_pipe.sh` | re-validated 2026-10-21 |

## tutorials/

| tutorial | exercises | status |
|---|---|---|
| `channel_kmm180/` | Re_tau 180 channel (KMM), natural-stretched y, `[case] name = channel` cold start | re-validated 2026-10-21 (statistics block) |
| `min_channel/` | minimal-span refined channel; the suite's blocks + 2:1 + Chebyshev case (`input_gpu*.ini`) | re-validated 2026-10-21 |
| `cht/` | conjugate heat transfer: Flageul turbulent channel (`channel/`) and Neuhauser pipe (`pipe/`) vs published references | re-validated 2026-10-22 (channel a991a1c, pipe 532ef54) |
| `naca/rans/` | NACA 0012, alpha 5, Re_c 4e5 SST vs OpenFOAM (C_L, Cp, Cf), control-volume forces | re-validated 2026-10-24 (532ef54, from scratch on the nose grid on HoreKa) |
| `naca/LES/` | strategy note only (`STRATEGY.md`, the planned NACA LES); no runnable case yet | pending |
| `sailplane/` | file-based IBM from a CAD STL through `moby_prepare` | re-validated 2026-10-21 |
| `turbulentBoundaryLayer/` | spatially developing TBL with the CaNS/SIMSON trip, Blasius inlet, outlets: `smooth/` (the code comparison) and `rough_macdonald/` (egg-carton roughness through the outlet, Hama Delta U+ 3.41 vs 3.72); `overheadTest/` holds the performance campaign (configs, drivers, results) | stands (2026-10-05/06, after every change that touches it) |

## Removed on 2026-09-26 (recover with `git show <commit>^:<path>`)

- `profile_200_20_gpu_{ibm,noibm}/` — nsys/ncu outputs of a pre-block-era binary.
- `smoke/` — old run outputs (h5/xdmf/vtk) and inis of the pre-block code;
  superseded by `min_channel`, `freestream/pois_io` and `prepare/wavy*`.
- `tutorials/wavychannel/` — despite the name, a 5-step plain channel with no IBM.
- `tutorials/naca/*.py` — post-processing of the removed polar sweep
  (`naca/rans/postProcess/` carries the live versions).
- `validation/momentum_interface/`, `channel_interface/interface_benchmark/`,
  `beltrami/{uniform32,slab_y_diag,slab_y_diag500}.ini` — needed the removed
  `MOBY_*` diagnostic hooks (RETIRED records).
- `validation/channel_interface_mfu/` — unvalidated salvage from `claude/blocks`,
  never run on this code; `channel_interface/` is the validated case.
- `validation/channel_interface/VALIDATION_CASES.md`, `developed/run_reflux_study.sh`
  — index/driver for the removed reflux and interface toggles.
- `validation/beltrami/{PENDING_TESTS.md,refined_fast.ini}` — the `claude/blocks`
  agglomeration study (`MOBY_AGGLOM_*`), code gone.
- `validation/taylor_green/` (+ `tools/check_tgv.py`, `tools/plot_tgv_error.py`)
  — 2D TGV, subsumed by the 3D Beltrami gates.
- `validation/naca0012/`, `validation/sd7003/` — the NACA 0012 SST sanity /
  L5-xz fan bench and the SD7003 transition benchmark. Their force gates used
  the removed penalization integral and they carried no `cv_box`; recorded
  results stay quoted in `refine2d/` and `prepare/`. The live airfoil case is
  `tutorials/naca/rans/`.
- `validation/refine2d/{gate_leaftable.sh,compare_leaftable.py,leaf_grid.ini,box_xz.ini}`
  — needed the deleted `mobygrid`; `box_xz` was a historical error-stop check.

Moved: `tutorials/interface_decay/` -> `validation/interface_decay/`.
