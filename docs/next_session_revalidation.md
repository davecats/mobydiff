# Next session: the one-by-one numerics re-validation of `validation/README.md`

STATUS: handout written 2026-10-21 on head `a6f50aa` (the four decisions of
that day committed: F13, F14, the nut read-back, the CLAUDE.md log move).
**NOTHING RE-RUN YET; the seven decisions of section 4 were TAKEN by the
user on 2026-10-21 (the same day) and are folded into sections 2.1, 3.4
and 7.** This is the second half of the user's standing order
of 2026-09-26 ("first the code simplifications and a thorough analysis of the
code; THEN re-go through every gate, validation case and tutorial ONE BY ONE
-- re-run, re-measure, update or retire"); the first half closed with
`docs/next_session_code_simplification.md`. Section 7 is the LEDGER: one row
per case, to be updated as the work proceeds. A case is DONE when its row
says so AND its README carries the dated block of section 2.4; the pass is
done when no row says PENDING.

Read first: `validation/README.md` (the index and what each case exercises),
CLAUDE.md "State of the solver" (gate rules, landmines, hosts),
`docs/cleanup_review_2026-10-06.md` section 6 (the driver-staleness list,
each a one-line fix to take when the case is touched), and the README of the
case at hand. The survey behind this handout (every README read in full on
2026-10-21) is condensed in sections 3 and 7; line numbers there are of
`a6f50aa`.

## 1. Why every number on file is suspect, and which change touches what

The numbers in the per-case READMEs were recorded between 2026-06-29 and
2026-10-06. Since the oldest of them the following changed the arithmetic or
the settings. Per case, section 7 names which of these it has NOT been
measured after.

| change | date | moves |
|---|---|---|
| Chebyshev `niter` 6 -> 12 (every Chebyshev ini, 51 files) | 2026-09-26 | every long Chebyshev case. NOT the inis that declare plain Jacobi (no `accel` line) -- those still say `niter = 6` (section 2.1) -- nor the red-black cases (KMM180 niter 3, boundary layer 6, sailplane 20), the NACA (18), the exactness gates with their own counts (Beltrami 200, conjugate/scalar 2/40/50, prepare 3) |
| F6 KMM180 `natural_dyw_plus = 0.5` | 2026-09-27 | the KMM180 tutorial (no statistics recorded since) |
| F2 advective conjugate default | 2026-09-27/28 | conjugate `conserve` (re-baselined), the CHT pipe (re-run at niter 6), the CHT channel (NOT re-run) |
| F4 sailplane BCs | 2026-09-27 | the sailplane (one-step leg re-run, no number recorded) |
| F1 flux-form viscous stencil | 2026-09-30 | every STRETCHED line: rans_sst natural y (turb180 re-measured: u_tau 1.0008 -> 1.0000), KMM180 (statistics within scatter), Blasius (unchanged to the digit), the CHT channel's stretched y, the NACA, the boundary layer (geometric y; re-run after). Uniform lines are an identity |
| F7 rational penalization factor | 2026-10-01 | every body: les_ibm (`a_wale` re-measured at 20 steps only), cylinder (re-measured), scalar S3/S4 body gates (re-measured), NACA (restart check only), CHT pipe and channel (NOT re-run), sailplane, rans_sst `ibm180*`, iddes `ibm`, multilevel_body, refine2d `unibody_xz`, body_outlet and bodyforce (measured after) |
| the predicted outlet face, the total-pressure datum on backflow | 2026-10-01 / 10-06 | every outlet: freestream, Blasius, cylinder, body_outlet, the boundary layer (all measured after); the NACA (restart check only); the sailplane (not re-run) |
| `pecletmax` and `cflmax` bound the SUM over directions | 2026-10-02 | the step of every adaptive case; Courant- or diffusion-bound ones take a smaller step by design (turb180, lam30t, conduction, prsweep, turbsst, pois_io moved then) |
| `nut` read back on a RANS/IDDES restart (`b50ea9f`) | 2026-10-21 | every restart FROM a RANS or IDDES snapshot: the NACA restart leg, iddes `fd1`/`toggles`, scalar `wfsst` (from `turbsst_186711.h5`), rans_sst's own restart checks. A restart from a non-RANS IC (`ibm180*` from the les_ibm IC) reinitialises nut as before |
| the Courant reduction over the interior (`431a4f5`) | 2026-10-21 | a Courant-bound cold start FROM REST with an inlet only; no suite case, measured |
| F13 restart precedence (`7ee2776`) | 2026-10-21 | behaviour, not numbers: a restart whose ini leaves a run-control key at its default while the file holds another value STOPS (section 2.2) |

## 2. Rules for every re-run

### 2.1 Settings (the ini the case is re-run with)

- **The projection of a rerun case (decision 1): Chebyshev `niter = 12` if
  it runs stably on the case, else plain Jacobi `niter = 12`; and red-black
  SOR at `niter = 3` or `6` becomes an admissible production choice for
  refined cases IF the R1 band channel of decision 6 shows it works across
  the 2:1 interface (until then SOR stays where it is declared).** Record
  which one the case ended up with and why in its README block. The
  2026-09-26 sweep moved the Chebyshev inis only; the inis that declare
  plain Jacobi still carry `niter = 6` and move to `accel = chebyshev`,
  `niter = 12` WHEN THEIR CASE IS RERUN (falling back to Jacobi 12 on an
  instability, which is itself a finding to record): all twelve
  `rans_sst` channels (`turb180 turb395 turb180t laminar laminart lam30
  lam30t wf180_y05/15/22/30/45`), their scalar twins in `validation/scalar`
  (`turbsst wfsst wferr wfs180_y05/15/30/45 ibmwf180 ibmwf1000 conduction`),
  `refine2d/rans_xz.ini`, `cylinder/empty.ini`, `freestream/oblique.ini`.
  Keep what the other inis declare: red-black where set (KMM180, the
  boundary layer, the sailplane), NACA 18, the exactness/stability gates
  with their own counts (Beltrami 200, `interface_decay` 3, conjugate and
  scalar 2/40/50, prepare 3, `bodyforce/trip_small` jacobi 6 -- it is a
  bit-exactness gate of the force kernel).
- **`cflmax = 1.0` on every rerun** (the user's 2026-10-05 decision; 0.8
  is 46 % of the RK3 limit since `cflmax` bounds the sum); delete the `sor`
  line where it is ignored (Chebyshev / Jacobi); do NOT sweep the inis for
  either -- change the keys of the case being rerun. The
  turbulentBoundaryLayer tutorial keeps `convection = divergence` (the one
  CaNS-comparison exception).
- A case whose ini is in one of the bit-exactness suites (`min_channel`,
  `les_ibm[_refine]`, `beltrami_yslab`, `turb180`, `wf180_y30`, `lam30t`,
  the eleven scalar cases, the outlet seven, `beltrami_rect[_redblack]`)
  changes the suite's case for every future gate when its ini changes: run
  that suite once after the change, NEW = the head at nofma, REF =
  `~/pass2_ref_binaries` at nofma, and expect `max_abs 0` still (both
  binaries read the new ini; the keys are old).

### 2.2 Restarts

- **F13 stops a restart whose ini leaves a run-control key at its default
  while the file holds another value** (`re`, `forcing_*`, `cflmax`,
  `pecletmax`, `dtmax`, `niter`, `sor`, `[ibm] enabled`, the BC rows of
  non-periodic faces). The message names the key and says to set it. The
  grid keys must agree outright (a field cannot be restarted on another
  grid). Set every run-control key explicitly in a restart ini; an explicit
  key wins silently, which is how `niter = 12` or `cflmax = 1.0` is applied
  to an archived snapshot. The archived snapshots in use all read: every
  metadata attribute but `natural_dyw_plus` is REQUIRED by the reader
  (`field_hdf5.c fdm_h5_read_metadata`), so a file old enough to lack one
  would stop with "could not read restart metadata", and none of the files
  of section 2.5 does.
- **Zero `pn` when restarting a niter-6 snapshot** (the cylinder README's
  clean-p recipe, with h5py on a COPY: the snapshot's pressure was converged
  at 6, the first steps at 12 otherwise see a pressure jump the forces and
  the stored-p diagnostics remember). A solver-minted template carries the
  OUTLET PLANE of its own run (`un_xmax` planes): an IC tool must write or
  delete it; a restart without the plane falls back to the zero-gradient
  value (freestream README L148-151, 1.9e-4).
- A RANS/IDDES snapshot now restarts with its own `nut` (restart ==
  continuous at `max_abs 0` on turb180 / wf180_y30 / lam30t and the IDDES
  pair, 2026-10-21). A RANS snapshot cold-started BEFORE 2026-08-05 is not
  bit-reproducible with the current binary (the k halo IC fix); its physics
  is.
- The final snapshot of a `t_final`-terminated run is a BAD restart (its
  trimmed last step; scalar README L1336): restart from a cadence snapshot.
- Mint a restart template from the SAME ini whose faces the run will use
  (the Blasius lesson: a sibling variant's template transplants its boundary
  rows; F13 now stops that case instead of running it wrong).

### 2.3 Binaries, bookkeeping, hosts

- Physics numbers at PRODUCTION flags; nofma only for bit-exactness
  statements. ONE binary for every leg of a case; the banner's commit hash
  goes into the README block. Build the GPU hosts from the SAME commit
  (istmcorax: `~/compile_corax.sh gpu` -> `build_gpu_corax`; istmcetus runs
  the local `build_gpu`, cc86). Name the device of every timing (RTX 3060
  here, A6000 istmcetus, RTX 5090 istmcorax, A100 HoreKa).
- Prepare with the CPU build on the rank count the solve uses (an unset
  `nb` is one block per prepare rank); drivers call
  `tools/prepare_if_missing.sh`. Five drivers default to the `build_*/main`
  symlink (freestream, iddes, multilevel_body, rans_sst, rans_inlet
  `run_gates.sh`): it still resolves; change it to `moby_solve` when the
  driver is touched. `les_ibm/setup.sh`, `cylinder/setup.sh`,
  `rans_geometry/setup.sh` default to `/home/davide/ibmc/bin/python` (not
  here: `~/ibmc`).
- GPU statistics/heat files are not bit-reproducible run to run (the
  atomic accumulation order; 4e-16); the fields are. A statistics identity
  is a CPU statement.
- `module load toolkits/nvhpc/25.9` in the shell of every `mpirun`;
  `OMPI_MCA_hwloc_base_binding_policy=none`; absolute `NEW`/`REF`/`BIN`
  paths; always `mpirun`; never edit a driver while it runs; never run
  `freestream/run_gates.sh` and `run_bitexact_outlet.sh` at once (the first
  regenerates the ICs the second restarts from); stage binaries in a
  scratch `bin/<label>/` so the tree can be rebuilt meanwhile; kill batches
  by PID from `pgrep -x moby_solve`, never `pkill -f`.
- The conjugate gate drivers REWRITE tracked record files
  (`validation/conjugate/*.dat`, `*.heat.txt`): `git status` after a batch.
- Hosts on 2026-10-21: the workstation (RTX 3060, 12 GB: tiers A and B, two
  streams); istmcetus (2x A6000 48 GB, both IDLE; `nohup` + status file,
  `numactl --cpunodebind` per GPU, GPU0 = cpus 0-63, GPU1 = 64-127; no h5py
  there); istmcorax (RTX 5090 32 GB, IDLE; 2.6x an A6000 on KMM180); HoreKa
  (4x A100 nodes; `dev_accelerated` one job per user in ~40-min resumable
  chunks, the user's own jobs on `accelerated`; the boundary layer only).
  Launch the tier-C runs on the remote hosts FIRST (they take days) and read
  them out last.

### 2.4 The record

Per case directory, ONE dated block at the top of the results in its README:

```
RE-VALIDATED 2026-MM-DD, commit <hash> (banner), production flags,
<device>, <ranks> ranks, <wall time>: <the gate numbers, one line each,
old value -> new value where it moved, PASS/FAIL against the recorded
tolerance>. Settings changed for the rerun: <keys>.
```

Older results stay as a superseded line ("2026-07-08 at niter 6: ...") only
where the mechanism they document is still cited; otherwise they go. Stale
descriptions of removed behaviour (the "niter = 6" boxes, "reflux",
"penalization C_D", "READY TO RUN", "HISTORICAL" numbers that no longer
reproduce) are deleted, not annotated. Then: the row in section 7 (status,
date), the status column of `validation/README.md` (add one: `re-validated
2026-MM-DD` / `RETIRED` / `stands`), and one worklog bullet per directory
in `docs/history/worklog_2026.md` (the numbers that moved and why, the
decisions, the files added or removed). A case RETIRED gets a two-line
header in its driver or README: what it covered, what covers it now, and
`git show <commit>^:<path>` for the record.

### 2.5 Inputs NOT in git

Symlink them into the tree (untracked; the symlinks of 2026-09-27 for the
suites exist: `validation/scalar/IC_turbles.h5`, `IC_turbslab.h5`,
`channel_interface/les_ibm/IC_refine.h5`, `ibm_coeff_blocks.h5`,
`rans_geometry/ibm_coeff_blocks_l{1,2}.h5`). Before declaring a case
failed, grep its log for `MPI_ERR_NO_SUCH_FILE`.

| file | lives in | needed by |
|---|---|---|
| `tutorials/channel_kmm180/channel_kmm180_restart.h5`, `_dywplus0p5.h5` (361 MB each) | `~/Codes/mobydiff.scalar/tutorials/channel_kmm180/` | KMM180 (`_dywplus0p5` is the F6 line, 2026-09-30 legs); `tools/make_channel_restart.py` for every channel_interface IC (`developed/IC.h5`, `core_patch/IC.h5`, `les/*_IC.h5`, `refine2d/IC_xz.h5`) |
| `validation/channel_interface/developed/IC.h5` (210 MB), `IC_xz.h5` | `mobydiff.scalar/validation/channel_interface/developed/` | developed, refine2d xz channel |
| `validation/channel_interface/core_patch/base_IC.h5` (34 MB), `IC.h5` (41 MB) | `mobydiff.scalar/.../core_patch/` (`IC.h5` in `mobydiff.bl`) | core_patch |
| `validation/channel_interface/les/uniform_IC.h5` (6.3 MB), `patch_IC.h5`, `reference_IC.h5`, `runs/{uniform,slab,patch,reference}` | `mobydiff.scalar/validation/channel_interface/les/` | les; iddes `iddes180` (its IC) and the WALE statistics reference |
| `validation/rans_sst/turb180_132565.h5` (the converged T2 field) | `mobydiff.scalar/validation/rans_sst/` | iddes `fd1`/`toggles` (`run_gates.sh` hard-codes the step: regenerate from the new turb180 run and read `latest`) |
| `validation/iddes/iddes180_dev_10000.h5` | `mobydiff.scalar/validation/iddes/` | iddes `toggles` (or regenerate) |
| `validation/scalar/turbles_67600.h5`, `turbsst_186711.h5` | `mobydiff.scalar/validation/scalar/` | S4 `s4stats`, `wfsst` |
| `validation/cylinder/cvpz_20301.h5`, `cyl_re40_20001.h5` (67 MB) | `mobydiff.scalar/validation/cylinder/` | cylinder (clean-p Re 40 state), scalar S3/S4 `cyl` (`cylheat` from the Re 40 steady field with `pn` zeroed) |
| `validation/multilevel_body/ibm_coeff_ml3_zero.h5`, `refine2d/ibm_ml3_xz_zero.h5` | `mobydiff.scalar/validation/{multilevel_body,refine2d}/` | regenerable: `multilevel_body/setup.sh` + `zero_coef.py` |
| `validation/conjugate/{pipe_prod.h5 (2.8 GB), IC_prod.h5, p_pr_settle2_30000.h5 (5.6 GB each), chan_conj*.h5}` | `mobydiff.scalar/validation/conjugate/` (the 2026-09-20 pipe campaign; the F2 re-run in `~/pipe_rerun_f2/`) | CHT pipe stats leg; `chan_conj` |
| `validation/conjugate/cht/cht149_flageul.h5` (539 MB), `cht149_nb32.h5` | `mobydiff.scalar/validation/conjugate/cht/` | the CHT channel (a developed seed candidate; the README says "no small one to ship") |
| `tutorials/naca/rans/c11_aoa5_460000.h5` (460 MB, t = 29) | `~/Codes/mobydiff/tutorials/naca/rans/` | NACA baseline restart; `c11_nose_660000.h5` (production, t = 34) NOT FOUND on this filesystem (2026-10-21) -- decision 2 |
| `tutorials/turbulentBoundaryLayer/restart_field.h5` (5.6 GB, step 1587500), `production_stats.h5` | in the tree, untracked; the HoreKa run dirs `rough_tbl_run`, `rough_io_run` | the boundary layer (stands, section 3.4) |
| Neuhauser tars (7.2 GB), `~/ibmc` venv (zarr, trimesh, shapely, h5py, matplotlib) | `~/ibmc`; the tars' location is not in the README | conjugate pipe pin-down, `run_gates_pipe.sh PY=`, cylinder STL regeneration, pipe post-processing |
| the retired mobygeom + venv | `tools/mobygeom*` (kept as the cross-implementation reference) | prepare P1b, refine2d gate (b), `rans_geometry/setup.sh`'s original path -- RETIRE those legs (section 3.5) |

## 3. The order of work

Four tiers. Tier C is launched on the remote hosts at the START of the
first session (two A6000s + the 5090 are idle) and read out at the end;
tiers A and B fill the workstation meanwhile. `rans_sst` goes first in any
case (the user: "the RANS rans_sst cases first, the nut read-back changed
their restart").

### 3.1 FIRST: `validation/rans_sst` (+ the scalar RANS twins)

What it is: the T2/T3/T4 physics gates of k-omega SST (resolved walls Re_tau
180/395, laminar Re_tau 10, the file-IBM channel `ibm180`, the refined wall
bands `refine180`), wall functions (`wf180_y05..y45`, `ibm180wf`),
gamma-Re_theta transition (`laminart lam30 lam30t turb180t`). Numbers on
file: 2026-07-08/09 (niter 6, pre-F1, pre-F7, pre-sum-limiter); turb180
alone re-measured 2026-09-30 at plain Jacobi 12 (u_tau 1.0000, log-law
0.049, U+ 18.16; recorded in `docs/history/next_session_after_step7.md`
only). Last changes that touch it: F1 (natural y of turb180/turb395/
turb180t), F7 and the les_ibm IC (`ibm180*`), the nut read-back (its
restart), the ini moves of 2.1.

Do:
1. Inis: the twelve plain-Jacobi channels -> `accel = chebyshev`, `niter =
   12`; all sixteen `cflmax = 0.5` -> `1.0`; drop `sor`. Check `dtmax`
   stays the binding limit on the channels (5e-2; the sum limiter at cflmax
   1.0 on 8x96x8 is far above it) so the step does not change with cflmax.
2. `RANKS=4 BIN=<abs>/build_cpu/moby_solve ./run_gates.sh` for the channels
   (seconds to minutes each: 8x96x8 to 16x64x8, 4000-7500 steps at dtmax);
   `ibm180` and `ibm180wf` on the GPU (`BIN=<abs>/build_gpu/moby_solve`,
   200 000 steps each, 18.6 ms/step on the RTX 3060 in 2026-07 = ~1 h each;
   or istmcetus). Then `./check_gates.sh` (python3 + h5py). Expected: T2
   PASS with turb180 u_tau 1.0000 (F1), turb395 and the wf series within
   their tolerances (uniform-y wf cases: F1 is an identity there, the
   Chebyshev switch moves the steady state only in the last digits), the
   transition cases unchanged in kind (lam30t laminarises, turb180t stays
   turbulent).
3. Restart == continuous at the DEVELOPED state, not only at step 10: restart
   turb180 from its t = 100 cadence snapshot to t = 200 and compare with the
   continuous run's final field (`tools/h5maxdiff`, expect `max_abs 0` on
   `un vn wn pn k omega nut`); the same for `ibm180wf` from a mid-run
   snapshot (20 steps suffice). This is the gate of the nut read-back on a
   converged RANS field; the 2026-10-21 gate was ten steps from a cold start.
4. Regenerate `turb180_<last>.h5` for iddes (section 2.5) and point
   `iddes/run_gates.sh` at `latest`.
5. The 9-case suite once (the three suite inis changed), both binaries at
   nofma: `max_abs 0`.
6. README: the dated block replaces the 2026-07-08/09 results; the 2026-08-05
   cold-start IC story stays (mechanism, still cited); the "RE-BASELINE
   REQUIRED" paragraph (about `~/s5a_ref_binaries`) goes.
7. The scalar twins right after, same ini moves: `validation/scalar`
   `run_gates_s2.sh sst kays wferr` (turbsst 186 711 steps on 4 CPU ranks),
   `run_gates_s5.sh unit ref sweep` (the thermal wall function), the
   `ibmwf180`/`ibmwf1000` budgets (ibmwf1000 on a GPU), `wfsst` (restarts
   from `turbsst_186711.h5`, a RANS snapshot: the nut read-back is live
   here). The S4 `heat` checker's stale penalization term (its graded
   column reads 3.8e+25 against the solver's 0.25 while the 1 % `--balance`
   closure passes; found 2026-10-20, identical on the reference binary) is
   fixed in this step: `check_scalar_stats.py heat` must use the RATIONAL
   factor and the case file's coefficient column.

### 3.2 Tier A: the cheap gates (one workstation session, CPU + GPU streams)

Each is minutes; run them in two batch scripts (CPU 4 ranks / GPU) with
absolute paths; record the dated block even where the numbers do not move.

- `courant_limit/run_gate.sh`, `diffusion_limit/run_gate.sh`: 2026-10-02
  numbers; the Courant reduction changed on 10-21 (gate 11/11 that day, not
  recorded in the README). Record.
- `penalization/run_decay.sh` (2026-10-01, post-F7; `REF=` an older binary
  shows order 1): record.
- `freestream/run_gates.sh` all groups + `run_restart_outlet.sh` +
  `run_bitexact_outlet.sh` (1 and 4 ranks): 2026-10-01 numbers; the
  `refined` group starts from REST through a patch touching an outlet -- the
  one gate shape the interior-only Courant reduction can touch if it is
  Courant-bound (it runs at `dtmax 2.5e-3`: the reduction moves nothing
  while dtmax binds; check the printed dt is the cap). The vortex gate is
  informational (no threshold): give it one (the 2026-10-01 reflected
  fraction 9.80e-5, monotone) or say so in the index.
- `body_outlet/run_gates.sh all` + `pocket` + `gpu` (2026-10-05/06): record;
  the rough/smooth boundary-layer cut-downs stay laminar at t = 300 and are
  not a gate (README L144-147): say so, the production rough case on HoreKa
  (2026-10-08) is the measured one.
- `bodyforce/run_gates.sh` (no README!): write a ten-line README from the
  driver's header and the index row; record.
- `prepare/run_gates.sh` (P0 29/29 on 2026-10-21), `run_gates_stl.sh` (P1,
  2026-07-16; the sphere leg SKIPS without the venv: run it with `~/ibmc`),
  `run_gates_order.sh` (2026-10-01, `OLD=~/step9b_ref_binaries`, the last
  legacy-order writer): record. P1b RETIRE (section 3.5).
- `multilevel_body/setup.sh` + `run_gates.sh` (2026-09-26, pre-F7: the body
  twin is zero-coefficient, so F7 cannot move it; record), `rans_geometry`
  (flat 2026-09-26 via moby_prepare; wavy 2026-07-07: re-run both modes with
  `check_rans_geometry.py`; the `setup.sh` mobygeom path is documented as
  retired, moby_prepare is the generator), `rans_inlet/run_gates.sh`
  (2026-09-26), `block_nb/run_gates.sh cpu|gpu` (2026-09-26; its FINDING
  table of 2026-08-27 -- the 2:1 interface is not nb-independent at 1e-4 --
  re-measure once; `beltrami_rect[_redblack]` have no recorded result:
  record the suite's), `redblack_interface/run_gates.sh cpu|gpu` (no dated
  result; the R1 "NOT yet done" list: DECIDED (6) -- run the cheap items,
  mass residual, Beltrami y-slab order, `refine_body` 2000 steps, the omega
  sweep 0.8-1.7, here, AND the Re_tau 180 band channel vs Jacobi as a
  tier-C run on istmcetus: it is the test that admits SOR as a production
  projection across a 2:1 interface, decision 1; the performance item is
  deferred to the performance handout), `interface_decay`
  (no recorded result, `niter = 3`: run, record the contraction),
  `refine2d`: `uniform_xz`/`unibody_xz` (2026-09-26) + the Beltrami pair
  `bp_xz_32/64` + `rans_xz` (2026-07-15, niter 6/50: re-run, expect the
  order 2.83 within 0.05), gate (c) `allref_xz == twin_xz` at tolerance 0.
- `beltrami/run_beltrami.sh cpu|gpu`: NO README exists (the index row is
  the only description); run the uniform convergence (nx 32/64/128, order
  2) and the three slabs + `refined`, and WRITE the README with the numbers.
- `blasius/run_blasius.sh` (2026-10-01, GPU 1 rank, 41 s; ONE-RANK gate):
  record; the stale "NO accel" comment in the ini goes.
- `tutorials/min_channel`: a suite case; the README's HISTORICAL blow-up
  table and its "niter 6" sentence go, the README says what the case is
  for today (the suite's blocks + 2:1 + Chebyshev case, 1 == 4 ranks).
- `tutorials/sailplane`: prepare (2 ranks, minutes) + the one-step solve on
  the GPU + the P1b sailplane leg (the only P1b leg kept): record the
  bit-exactness; no physics number exists or is needed.
- scalar S1 (`run_gates.sh` all groups; `cond_64` 19 min GPU), S3 all groups
  except `cyl` (2026-10-01), S4 all groups except `cyl` (2026-08-07, need the
  turbles snapshots: section 2.5), S5a (with 3.1): record. The `uniform`,
  `det`, `restart`, `noeffect` identities must read exactly 0.
- conjugate C1 `run_gates_c1.sh all` (61 legs, 2026-10-21 PASS in the
  simplification gates; numbers of C3 era), C2 (checker-only legs have no
  binary dependence; the `dt`/`cylinder` legs do), C3 `fraction transient
  nusselt budget channel` (the `channel` leg is `chan_conj`: 50 min on one
  CPU rank, under a minute on the GPU), `run_gates_pipe.sh all` (2026-09-15,
  `~/f1_ref_binaries`): record. The 2026-09-27 `conserve` drift (1.286e-12)
  and the C3-era drifts (6.9e-18) are two code states not reconciled in the
  README: the new block states ONE number.

### 3.3 Tier B: the GPU-hour class (workstation or one A6000)

- `channel_interface/les_ibm`: `run_ibm_les.py --case all` (`a_wale`,
  `b_none`, `c_refine`; 25 600 steps each, 24 ms/step on the RTX 3060 =
  ~10 min each), `measure_nut.py`, `ibm_les_stats.py`. On file: 2026-06-30
  (niter 6) + the 2026-10-01 `a_wale` 20-step F7 check; the rational factor
  never ran the statistics. Expect bulk U near 15.07 (the niter-12 value
  seen on 2026-10-01), nut step ~2x, the solid-nut and band ratios as
  recorded.
- `iddes/run_gates.sh` all groups (2026-07-10, niter 6): `iddes180` t 0..25
  on 64x48x64 (50 000 steps, 0.029 s/step on the RTX 3060 = ~25 min; needs
  `uniform_IC.h5` and the WALE `runs/uniform/stats` reference, section 2.5),
  `fd0`/`fd1` (from the NEW turb180 field), `ibm`, `ranks`, `toggles`;
  `check_gates.py`. The IDDES restart now reads nut (and `fd` is restart
  state): `fd1` restarts from a RANS field. Expect the fd profile (1.000
  through y+ 15; band means 0.87 / 0.03) and the log-layer deviation
  (0.5 % vs WALE) to hold.
- `cylinder`: `setup.sh` (moby_prepare), Re 40 steady C_D (12 000 steps, 88
  ms/step on an A6000 = 18 min; the gate entry still says 1.6924 ± 7e-4 from
  the pre-outlet series while the outlet section measured 1.69682 ± 1e-5:
  ONE number after this run), Re 100 Strouhal (20 000 steps, 30 min; the
  gate entry quotes the 2026-07-13 PENALIZATION St while the control-volume
  one is 0.1744 at 16 x 16: the gate becomes the control-volume number and
  `check_cylinder.py`'s "penalization C_D" header goes), `empty`, force
  determinism 1 == 4 ranks. `run_f7_gates.sh` and `run_outlet_mode.sh` are
  one-off investigation drivers needing `~/step8_ref_binaries` and a sibling
  checkout: RETIRED header, kept under "how the numbers were produced".
  Then scalar S3/S4 `cyl` (`cylheat` from the NEW Re 40 state with `pn`
  zeroed; Nu 3.37 in the literature band [2.8, 4.0]).
- scalar S2 `det` and the `s4stats` group need a `turbles` cadence snapshot
  (section 2.5); S2 `les`/`band` are tier C.

### 3.4 Tier C: the long tier (remote hosts; launch first)

Estimates scale the recorded per-step times (section 7 has each); total
~150-250 GPU-h, the two CHT tutorials most of it, as the 2026-09-26 plan
said.

| case | run | where | cost | notes |
|---|---|---|---|---|
| CHT pipe `tutorials/cht/pipe` | FIRST a timing on istmcorax (the 5090 holds the 22.7 GB): ~200 steps from the settled state with Chebyshev `niter = 12` and with red-black `solver = redblack, sor = 1.5, niter = 6`, same binary, s/step of each (and that both are stable over the window); THEN `run_pipe.sh stats` with the FASTER one on HoreKa, from `p_pr_settle2_30000.h5` (5.6 GB, scp'd to the workspace with the case file and the ini; `cflmax 1.0` set explicitly; `pn` zeroed in a copy), chunked `dev_accelerated` or `accelerated`; post-process locally (no h5py on HoreKa) | istmcorax (timing) + HoreKa (run) | 168 000 steps: 24-26 h on the 5090 at 0.52-0.56 s/step; on one A100 of the order of a day, chunked | on file: the 2026-09-28 F2 re-run AT NITER 6 ("the binary is the ONLY change"), pre-F1/F7; the ini says 12 now. DECIDED (4): benchmark Chebyshev 12 vs SOR 6, run the fastest on HoreKa |
| CHT channel `tutorials/cht/channel` | `run_flageul.sh` from a developed Re_tau 149 seed (`cht149_flageul.h5` 539 MB in the sibling checkout is the candidate; `SEED=`), develop 80 000 + settle 65 000 + clean 650 000 steps; `check_cht.py`, `asset/compare_flageul.py` | HoreKa (`dev_accelerated` chunks or `accelerated`, the user allows both; the seed and the case file go up by scp; the projection per decision 1) | 715 000 steps at 0.562 s/step on the RTX 3060 = 112 h; a few days chunked on A100s, 4 ranks | undated, niter unknown (6 by its date), pre-F1 (stretched y), pre-F2, pre-F7. The one CHT comparison against published data at a flat interface. DECIDED (3): rerun on HoreKa. `check_cht.py` prints FAIL on two understood checks (README L134-151): re-tolerance after the run |
| `channel_interface/developed` | `run_developed.py --arch gpu --ranks 2`: transient 16 000 + statistics 64 000 steps; `tools/plot_channel_stats.py` | istmcetus (2 ranks = the two A6000s; `dims = 2 1 1`, never split y) | ~0.6 s/step on an RTX-class GPU = ~14 h, half on two | NO numbers recorded at all (its reflux study is settled history): this run produces the first niter-12 statistics of the band-refined channel. Needs `IC.h5` (section 2.5) and the uniform-fine reference (`run_reference.sh`, 0.29 s/step, 6 GB) |
| `channel_interface` main (`reference`, `refined_y110`, `refined_y55`) | `run_validation.sh gpu 1`: transient 5 + statistics 20 turnovers each; `tools/channel_interface_validation.py` | istmcetus, one case per GPU, the reference last (256x128x256: 80 000 steps at ~0.3 s/step on an A6000 = ~7 h; the refined ones ~3 h each) | ~15 GPU-h | NO numbers recorded; the README's solver box says niter 6. DECIDED (7): run and record |
| `channel_interface/core_patch` | gates 1, 2, 4 (`core_patch.ini` 250 steps, `base.ini`), gate 5 `run_gate5.py --arch gpu --ranks 2` (t 5..25, 201 snapshots, ~8 GB) | istmcetus | gate 5 ~5 h | 2026-06-29, niter 6; gate 3 needed `MOBY_STEPDIV` (removed): RETIRE gate 3. The base-128 pressure null mode (pn std 2e5) was blamed on niter 6: the run says whether it is gone |
| `channel_interface/les` | `run_les.py --arch gpu --case all` (`uniform`, `slab`, `patch` 64^3 + WALE, `reference` 128x64x128 DNS), `les_stats.py`; `run_wall_jump.sh` (2026-10-02, seconds) | istmcetus / corax | the LES cases ~1-2 h each, the reference several hours | 2026-06-30, niter 6 ("the under-converged niter=6 projection" named in the p column); needs the four ICs (section 2.5); the "READY TO RUN" header goes. Also the iddes WALE reference |
| `refine2d` xz channel | `developed/run_developed.py --refine-dims xz` from `IC_xz.h5`, `analyze_xz_channel.py` | istmcorax | 80 000 steps, ~4 h on the 5090 (2026-07) | 2026-07-15, "niter = 6 throughout"; the last xz validation with turbulence |
| KMM180 `tutorials/channel_kmm180` | `input.ini` cold start (`[case] name = channel`, red-black 3, `t_final 200`, dt 3.125e-4 = 640 000 steps) OR a two-leg restart from `channel_kmm180_restart_dywplus0p5.h5` as on 2026-09-30 (30 t.u. discarded + 20 averaged, 64 000 steps) | istmcorax (0.0208 s/step = 3.7 h full; 22 min for the two-leg form) | 0.4-4 h | NO statistics recorded in the README ("the case is on the one-by-one re-validation list"); the F1 legs' statistics are in `asset/f1_2026-09-30/` (U+ 18.51, u'+ 2.65, both ~1 % from KMM). Record a statistics block in the README from the two-leg form; the full cold start only if the user wants the tutorial's own developed field regenerated |
| NACA `tutorials/naca/rans` | `run_case.sh scratch` on the NOSE grid (stages L10 -> L11 -> nose per the README protocol, `interp_restart.py` between them; prepare ~40 min on 20 CPU ranks, done on the workstation and scp'd), the stopping table's tolerances; then `postprocess.sh`, `cv_forces.py`, `compare_openfoam.py` | HoreKa (A100s, chunked; the user: "use HoreKa if necessary to speed up") -- the production file `c11_nose_660000.h5` is NOT on this filesystem, so the scratch run also regenerates the tutorial's own restart | 2026-08: ~33 h stage 1 on an RTX 3060 + ~46 h stages 2-3 on an A6000; on A100s a day or two chunked | 2026-08-08 at niter 18 (not 6 -> 12), pre-F1/F7/outlet/nut; the 2026-10-01 restart check read C_D 0.01275 on either binary against the table's 0.01294 ("one scatter below", unresolved): the scratch run settles it. DECIDED (2): re-converge from scratch on the nose grid |
| scalar S2 `les`, `band` | `run_gates_s2.sh les band` (`turbles` then `turbslab`: relax + re-target + statistics legs, 1 rank, pinned `dims 1 1 1`) | istmcetus (one GPU each) | 64x48x64 at dtmax 5e-4 to t 40 + relax: ~2 h each on an A6000 (the README says "hours to days" for the CPU form) | 2026-08-07, niter 6; the band excess sign question (-0.0265 vs +0.0012, a normalisation effect, README L1300, L1384) is settled by stating ONE normalisation in the block. DECIDED (5): rerun both on istmcetus |
| boundary layer `tutorials/turbulentBoundaryLayer` | -- | HoreKa | 500 000 steps at 0.18 s/step on 4 A100 = 25 h | **STANDS**: measured 2026-10-05/06 on main with the predicted outlet, post-F1, red-black (no niter change), `cflmax 1.0` and `convection = divergence` already by the user's decision; the only later change (the interior Courant reduction) cannot touch a restart from a developed field. Record "stands" in the ledger and the index; no run |

### 3.5 Retire (with a header) -- no decision needed

`prepare/run_gates_big.sh` (P1b; mobygeom references, the airfoil rows
cannot be re-run; the sailplane leg moves into `tutorials/sailplane`'s
record), `refine2d` gate (b) (mobygeom + mobygrid), `core_patch` gate 3
(`MOBY_STEPDIV`), `cylinder/run_f7_gates.sh` + `run_outlet_mode.sh`
(investigation drivers; their findings stay in the README), the conjugate
C3 transient control rows (already retired 2026-10-06, keep the note),
`scalar/compile_nofma.sh` (duplicates `./compile.sh cpu_nofma`), the
`rans_geometry/setup.sh` mobygeom path (moby_prepare is the generator since
2026-09-26), `min_channel`'s historical blow-up table, the NACA polar sweep
(removed 2026-09-26; needs a `cv_box` and a re-measured reference if
revived -- not in this pass).

## 4. DECIDED by the user on 2026-10-21 (asked one by one)

1. **Projection of the plain-Jacobi `niter = 6` inis when rerun:** Chebyshev
   `niter = 12` if it works stably on the case, otherwise plain Jacobi
   `niter = 12`; and if red-black SOR is observed to work across the 2:1
   interface too (to be tested: decision 6), SOR at `niter = 3` or `6` is
   also acceptable. Rule in section 2.1.
2. **NACA: re-converge from scratch on the nose grid**, on HoreKa if that
   speeds it up (the production restart file is lost; the scratch run
   regenerates it). Row in 3.4.
3. **CHT channel: rerun on HoreKa**, `accelerated` allowed as well as
   `dev_accelerated`. Row in 3.4.
4. **CHT pipe: first time Chebyshev-Jacobi `niter = 12` against red-black SOR
   `niter = 6` on the case; run the stats leg with the FASTER one, on
   HoreKa.** Row in 3.4 (the timing on istmcorax, the only local device that
   holds the 22.7 GB).
5. **Scalar S2 LES legs (`turbles`, `turbslab`): rerun both on istmcetus.**
6. **Red-black R1 list: run the cheap items AND the Re_tau 180 band channel
   vs Jacobi**; the performance item waits for the performance handout.
   Decision 1 makes this the admission test of SOR on refined cases.
7. **The three main `channel_interface` cases: run and record** (istmcetus,
   one per A6000, the reference last).

Host plan after the decisions: HoreKa takes the NACA scratch run, the CHT
channel and the CHT pipe (three chains; `dev_accelerated` is one job per
user, `accelerated` is allowed for these -- the user's own jobs sit there
too; two chains that build the same worktree must not overlap). istmcetus
takes `developed/`, the three main channel cases, the S2 LES pair,
`channel_interface/les`, `core_patch` gate 5 and the R1 band channel (two
A6000s, so queue them in a status-file chain per GPU). istmcorax takes the
pipe timing, the refine2d xz channel and KMM180. The workstation takes
tiers A and B.

## 5. Text to fix on the way (found in the survey; fix when the case is touched)

- `validation/README.md` header PLAN note -> points here; add the status
  column. `tutorials/cht/README.md` says the pipe conductivity pair agrees
  to 0.1 %, the pipe README says 0.2 %. `tutorials/min_channel/README.md`
  says the inis are at niter 6 (they are at 12). `channel_interface/
  README.md` solver box: niter 6. `channel_interface/les/README.md` L54
  "READY TO RUN". `cylinder/README.md` gate entry L68 (1.6924) vs the outlet
  section (1.69682); the clean-p protocol text says "production niter = 6".
  `iddes/README.md` example uses `build_cpu/main`. `scalar/README.md`
  L992-1000 lists the 7-case suite under `run_bitexact_s3.sh` with a typo;
  `run_bitexact_s3.sh:21` cites `~/s2_ref` (gone). `conjugate/README.md`
  two conservation drifts (2.4). `rans_sst/README.md` "RE-BASELINE
  REQUIRED" (`~/s5a_ref_binaries`). `prepare/run_gates_order.sh` writes
  `cflmax = 0.8` into a generated ini. `channel_interface/run_validation.sh`
  header, `developed/run_reference.sh`, `run_uniform128.sh` describe removed
  behaviour. The five `main`-symlink driver defaults and the three
  `/home/davide/ibmc` venv defaults (2.3).
- Gate leftovers: ~600 untracked `bx_*`, `i12*`, `s7_*` outputs in
  `validation/beltrami`, `tutorials/min_channel`, `rans_sst`, `les_ibm`
  (gitignored; delete per directory when its case is done).

## 6. Not in this handout

The performance ideas measured but not taken (the predictor split by
component, the red-black reduced exchange round, 7-5, the `apply_bc`
face-row skip) and the feature tail (S5b TVD + second halo, C4 corners, the
IDDES flat-plate inlet, transition under IDDES, the MacDonald rough-wall
forced convection) -- CLAUDE.md "Open, in order" (2) and (3), their own
handouts after this one. A gate that FAILS in this pass is a finding: fix
the root (never a compensating tweak of the gate), gate the fix at nofma
against `~/pass2_ref_binaries`, cut a new reference set from the commit if
a suite moves -- the rules of CLAUDE.md "Gate rules".

## 7. The ledger

Status: PENDING / DONE yyyy-mm-dd / STANDS (measured after every change
that can touch it; recorded, no run) / RETIRED. "Last" is the newest
recorded result; "since" lists the changes of section 1 it has NOT been
measured after (`12` = niter 6 -> 12, `F1`, `F2`, `F7`, `out` = the outlet,
`sum` = the sum limiter, `nut`, `cou` = the interior Courant reduction,
`ini` = the plain-Jacobi ini still at 6).

| case | last | since | tier | where | status |
|---|---|---|---|---|---|
| `rans_sst` T2/T3/T4 | 2026-10-21 (ALL PASS; wf180_y30 near-centre 3.02 % -> tolerance 0.035; ibm180wf == ibm180 and both restart gates max_abs 0) | -- | first | WS CPU + GPU | DONE 2026-10-21 |
| `scalar` RANS twins (S2 sst/kays/wferr, S5a, ibmwf, wfsst) | 2026-10-22 (ALL PASS at Chebyshev 12; the ibmwf budgets are the S4 heat gate's instrument, 2026-10-01 stands) | -- | first | WS | DONE 2026-10-22 |
| `courant_limit` | 2026-10-21 (12/12, unchanged) | -- | A | WS CPU | DONE 2026-10-21 |
| `diffusion_limit` | 2026-10-21 (12/12, unchanged) | -- | A | WS CPU | DONE 2026-10-21 |
| `penalization` | 2026-10-21 (to the digit; the step-9 REF leg no longer runs) | -- | A | WS | DONE 2026-10-21 |
| `freestream` (all groups, restart, bit-exact outlet) | 2026-10-21 (ALL PASS, numbers unchanged) | -- | A | WS | DONE 2026-10-21 |
| `body_outlet` (all, pocket, gpu) | 2026-10-21 (ALL PASS to the digit; gpu/twin at nofma) | -- | A | WS | DONE 2026-10-21 |
| `bodyforce` (README written) | 2026-10-21 (ALL PASS) | -- | A | WS | DONE 2026-10-21 |
| `prepare` P0 / P1 / order | 2026-10-22 (29/29, 15/15, 43/43 at nofma) | -- | A | WS | DONE 2026-10-22 |
| `prepare` P1b | 2026-07-17 (mobygeom) | -- | -- | -- | RETIRED 2026-10-21 (header; the sailplane leg is the tutorial's record) |
| `multilevel_body` | 2026-10-21 (ALL PASS, numbers to the digit) | -- | A | WS | DONE 2026-10-21 |
| `rans_geometry` | 2026-10-21 (flat 0.0, wavy 2.3e-11, sphere test PASS) | -- | A | WS CPU | DONE 2026-10-21 |
| `rans_inlet` | 2026-10-21 (0.14 % / 2.15 %, ranks EXACT) | -- | A | WS | DONE 2026-10-21 (CPU; the GPU ulp line of 2026-07-13 stands) |
| `block_nb` (+ the beltrami_rect suite cases) | 2026-10-21 (cpu + gpu ALL PASS; the suite pair max_abs 0; the FINDING table not re-measured) | -- | A | WS | DONE 2026-10-21 |
| `redblack_interface` (+ the R1 cheap items) | 2026-10-21/22 (cpu + gpu ALL PASS; omega limit in (1.6, 1.7); y-slab rb6/cheb12 error ratio 1.002; refine_body 2000 steps finite) | -- | A | WS + cetus | DONE 2026-10-22 |
| `redblack_interface` R1 band channel (Re_tau 180 vs Jacobi) | 2026-10-21 (same statistics to <= 0.5 %, no band; 19 % slower per step than Chebyshev 12) | -- | C | cetus | DONE 2026-10-21 (SOR admissible on 2:1) |
| `interface_decay` | 2026-10-21 (contraction 7.7x velocity / 10.9x pressure) | -- | A | WS GPU | DONE 2026-10-21 |
| `refine2d` uniform/unibody / Beltrami pair + rans_xz / xz channel | 2026-10-21 (driver `run_gates.sh`: uniform EXACT, twin EXACT, order 2.73, rans_xz EXACT) / xz channel on istmcetus (corax reboots) | -- | A / A / C | WS / WS / cetus | DONE 2026-10-21 (xz channel pending) |
| `beltrami` (README written) | 2026-10-22 (order 2.00 with Chebyshev 200; plain Jacobi 200 floors at 128^3; slabs 3.08e-3 symmetric) | -- | A | WS | DONE 2026-10-22 |
| `blasius` | 2026-10-21 (theta 1.66 %, H 0.82 %; ini moved to Chebyshev 12) | -- | A | WS GPU 1 rank | DONE 2026-10-21 |
| `tutorials/min_channel` | README rewritten 2026-10-21 (the suite case) | suite case | A | -- | DONE 2026-10-21 |
| `tutorials/sailplane` | 2026-10-21 (prepare 2 ranks; 1 == 2 GPU ranks max_abs 0) | -- | A | WS GPU | DONE 2026-10-21 |
| `scalar` S1 / S3 / S4 (not cyl) | 2026-10-22 (ALL PASS; the tile-order checker fix; heat == solver to 1e-15) | -- | A | WS | DONE 2026-10-22 |
| `scalar` S3/S4 `cyl` | 2026-10-22 (Nu 3.3673, CV cross-check PASS; S4 cyl 4/4) | -- | B | WS GPU | DONE 2026-10-22 |
| `scalar` S2 `les`/`band` | 2026-10-21 (both PASS on an A6000 in 3 + 12 min; band excess +0.0086 in the checker's normalisation) | -- | C | cetus | DONE 2026-10-21 |
| `conjugate` C1 / C2 / C3 / pipe gates | 2026-10-21 (61 / 75 / all / 36 PASS; drift 1.288e-12 stated once) | -- | A | WS | DONE 2026-10-21 |
| `channel_interface/les_ibm` | 2026-10-22 (ALL PASS; bulk U 15.01, nut step 2.00, band/core 0.05) | -- | B | WS GPU | DONE 2026-10-22 |
| `iddes` | 2026-10-22 (ALL PASS; fd 0.872/0.034, log layer 0.3 %, fd1 drift 1e-12 with the nut read-back) | -- | B | WS GPU | DONE 2026-10-22 |
| `cylinder` | 2026-10-22 (C_D 1.6968 +- 3e-5; St 0.1734, C_D 1.4385; empty exact; forces 1 == 4 ranks) | -- | B | WS GPU | DONE 2026-10-22 |
| `channel_interface/developed` | 2026-10-21 (first niter-12 statistics, within ~1 % of KMM, no band; red-black twin equal) | -- | C | cetus | DONE 2026-10-21 |
| `channel_interface` main (3 cases) | none | -- | C | cetus | PENDING (decided 7: run) |
| `channel_interface/core_patch` | 2026-10-22 (gate 5 PASS within 2 % incl. p; the base-128 null mode gone) | -- | C | cetus | DONE 2026-10-22 (gate 3 retired) |
| `channel_interface/les` (+ wall jump 2026-10-02) | 2026-10-22 (four cases PASS; the niter-6 pressure null mode gone) | -- | C | cetus | DONE 2026-10-22 (the wall-jump gate stands at 2026-10-02) |
| `tutorials/channel_kmm180` | 2026-10-21 two-leg statistics block (5090; within ~1 % of KMM, scatter vs the F1 legs) | -- | C | corax | DONE 2026-10-21 |
| `tutorials/naca/rans` | 2026-08-08; restart check 2026-10-01 | F1, F7, out (check only), nut | C | HoreKa (scratch, nose grid) | PENDING (decided 2) |
| `tutorials/cht/pipe` | 2026-09-28 at niter 6 | 12, F1, F7 | C | corax (timing) + HoreKa | PENDING (decided 4: fastest of Cheb 12 / SOR 6) |
| `tutorials/cht/channel` | undated | 12, F1, F2, F7 | C | HoreKa | PENDING (decided 3) |
| `tutorials/turbulentBoundaryLayer` (+ rough) | 2026-10-05/06, 10-08 | (cou: cannot touch it) | -- | HoreKa | STANDS |
