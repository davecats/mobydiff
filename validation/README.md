# Test cases — index

> **PLAN (2026-09-26):** after a round of code simplifications and a thorough
> code analysis, every gate, validation case and tutorial listed here is to be
> re-gone-through ONE BY ONE (re-run, re-measure, update or retire). Until then,
> numbers recorded in the per-case READMEs for cases with long runs are niter-6
> values: all Chebyshev inis moved to `niter = 12` on 2026-09-26 and only the
> cheap gates have been re-measured (see CLAUDE.md "Test-case cleanup").

The repository has exactly two homes for runnable cases:

- `tutorials/` — user-facing cases (how to set up and run a real flow);
- `validation/` — every regression / verification gate (this directory).

There is no separate `smoke/` or `profile_*/` directory any more (cleaned up
2026-09-26; everything removed is recoverable from git history, see the end of
this file). Each subdirectory has its own README with commands and recorded
results; this file only says WHAT each case exercises.

## The standard regression suite

Since numerics review step 7 every run reads a CASE FILE (`[case] file`,
default `<field_prefix>.case.h5`, built in-process when absent), so a suite
run leaves one `<prefix>.case.h5` per output prefix beside the snapshots; the
`les_ibm` case reads the committed `ibm_coeff_case.h5` (the exact conversion
of the legacy `ibm_coeff.h5`, which the solver no longer reads).

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

## validation/

| directory | exercises | driver |
|---|---|---|
| `beltrami/` | 3D Beltrami/ABC exact decaying NS solution (`[flow] initial = beltrami`): 2nd-order convergence (`uniform.ini`, nx 32/64/128) and the 2:1 interface artifact for x/y/z slabs (`slab_{x,y,z}.ini`) and a 3D patch with edges+corners (`refined.ini`) | `run_beltrami.sh`, `tools/check_beltrami.py` |
| `interface_decay/` | stability of the 2:1 interface: white noise on a 3D refined patch must decay everywhere (catches interface-localized growth that smooth-flow gates miss). ~35 s on 4 CPU ranks | README, `tools/check_interface_decay.py` |
| `block_nb/` | per-direction `[blocks] nb`: single-level nb-independence (bit-exact across tilings and rank counts) and uniform oblique flow through a 3-level patch with non-cubic nb | `run_gates.sh` |
| `refine2d/` | `refine_dims = xz` quadtree: all-refined == doubled-resolution twin (`allref_xz`/`twin_xz`), uniform flow through xz patches (`uniform_xz`, `unibody_xz`), Beltrami order pair (`bp_xz_32/64`), RANS scalars across xz interfaces (`rans_xz`) | README |
| `redblack_interface/` | `[pressure] solver = redblack` across a 2:1 interface (R1); 1 rank == 4 ranks | `run_gates.sh` |
| `multilevel_body/` | 3-level `refine_body` on a cylinder: uniform oblique flow EXACT across all interface levels (zero-coef twin), per-level dwall vs an exact prism reference | `setup.sh`, `run_gates.sh` |
| `channel_interface/` | turbulent Re_tau 180 channel across the 2:1 interface: wall bands (`reference`/`uniform128`/`refined_y55`/`refined_y110`), developed statistics (`developed/`), edge/corner core patch (`core_patch/`), LES across refinement (`les/`), LES<->IBM coupling (`les_ibm/`, suite case). Long GPU runs | `run_validation.sh`, per-subdir READMEs |
| `freestream/` | inlet/outlet patch types + Dirichlet-pressure outlet: oblique freestream exact, inflow/outflow Poiseuille vs periodic, Lamb-Oseen vortex exit, 1 == 4 ranks, config contradictions | `run_gates.sh` |
| `cylinder/` | `[case] name = airfoil` on a cylinder: control-volume C_L/C_D (Re 40 steady drag, Re 100 shedding/Strouhal), empty-domain zero force | `setup.sh`, `check_cylinder.py` |
| `rans_geometry/` | SST wall distance / y_eff / wall-cell classification: file IBM (flat, ± refine_body) vs closed form, analytic wavy wall (± refine_body) vs scipy | `setup.sh`, `check_rans_geometry.py` |
| `rans_sst/` | k-omega SST (laminar decay, Re_tau 180/395, IBM channel, refined bands), wall functions (y+ sweep, IBM), transition (lam30t, laminart, turb180t) | `run_gates.sh`, `check_gates.sh` |
| `rans_inlet/` | RANS scalar freestream values at declared inlet faces, 1 == 4 ranks | `run_gates.sh` |
| `iddes/` | IDDES (SST + WALE blend): developed channel log layer, fd shielding profile, fd_force limits, IBM stability | `run_gates.sh`, `check_gates.py` |
| `prepare/` | `moby_prepare`: analytic geometry bit-exact vs inline (P0: wavy, wavy_refine, wavysolid), STL geometry vs mobygeom references + periodic shift invariance (P1; the flat-slab references come from `rans_geometry/setup.sh`), big geometries (P1b: now the sailplane leg only), and step 7 -- prepared body-free / nb-less / analytic cases solved from the case file == the reference inline solve, on 1 and 4 ranks and on the GPU | `run_gates.sh`, `run_gates_stl.sh`, `run_gates_big.sh`, `run_gates_step7.sh` |
| `scalar/` | passive scalars S0-S5a: uniform preservation, conservation, conduction, MMS, Pr sweep, determinism, LES/SST closures, IBM dirichlet/adiabatic, heated cylinder, in-solver statistics, thermal wall function (domain + IBM walls), `count = 0` bit-exactness | `run_gates*.sh`, `run_bitexact*.sh` |
| `conjugate/` | conjugate heat transfer at the IBM interface (C1-C3, F1-F5): multi-material slabs, contact resistance, oblique/curved interfaces, transient capacity, Nusselt diagnostic, channel with conducting walls, annular pipe geometry | `run_gates_c{1,2,3}.sh`, `run_gates_pipe.sh` |

## tutorials/

| tutorial | exercises |
|---|---|
| `channel_kmm180/` | Re_tau 180 channel (KMM), natural-stretched y, `[case] name = channel` cold start |
| `min_channel/` | minimal-span refined channel; the suite's blocks + 2:1 + Chebyshev case (`input_gpu*.ini`) |
| `cht/` | conjugate heat transfer: Flageul turbulent channel (`channel/`) and Neuhauser pipe (`pipe/`) vs published references |
| `naca/rans/` | NACA 0012, alpha 5, Re_c 4e5 SST vs OpenFOAM (C_L, Cp, Cf), control-volume forces |
| `naca/LES/` | strategy note only (`STRATEGY.md`, the planned NACA LES); no runnable case yet |
| `sailplane/` | file-based IBM from a CAD STL through `moby_prepare` |
| `turbulentBoundaryLayer/` | spatially developing TBL with the CaNS/SIMSON trip, Blasius inlet, outlets; `overheadTest/` holds the performance campaign (configs, drivers, results) |

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
