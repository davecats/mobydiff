# Next session: after numerics review step 7 — loose ends, then step 8

STATUS (2026-09-30): **loose ends 1, 2, 4, 5 DONE and gated; 3, 6, 7
are notes; step 8 (F1) IN PROGRESS — see the STEP 8 block at the end.**
The step-7 commits were rebased onto origin before the push (the pre-rebase
ids cited below: `d2249b1 -> 59ea085` (7-0), `5d95937 -> 91f01cc`,
`18b95c1 -> fc7f05e`, `3461610 -> 0d2d18f`, `da64b37 -> 13fb213`).

- **1 (HoreKa drivers) DONE.** Every launch site (`run_matrix.sh`,
  `exchange/run_{exchange,portgate,mapgate,timeline,ncu}.sh`,
  `submit_rdblk.sh`, `submit_port_gate.sh` pass 2) prepares the case file in
  the run directory before its `mpirun`, on the SAME rank count. The helper is
  looked up BESIDE THE BINARY'S OWN CHECKOUT
  (`<dir of moby_solve>/../tools/prepare_if_missing.sh`), not through
  `MOBY_ROOT` or `$HERE` (the submit scripts stage the drivers into a run
  directory and half of them do not pass `MOBY_ROOT`): the `ref` column of a
  two-binary job resolves to the ref checkout's helper, or to nothing for a
  pre-step-7 checkout, which builds inline -- the helper arrived in the same
  commit that made the solver refuse a missing file, so "helper present" is
  exactly "case file needed". A local reference set must therefore carry
  `tools/prepare_if_missing.sh` (the map-gate test failed once because it did
  not). Exercised LOCALLY on the CPU build: `run_matrix.sh` (1 and 2 ranks,
  `overhead.case.h5` + `prepare.log` per run directory, 1600 leaves split
  800/rank) and `run_mapgate.sh` (min_channel + Beltrami, new vs the step7b
  reference, `max_abs 0`). **NOT verified on the cluster** -- run one short
  `submit.sh` before quoting any campaign number (README says so too).
- **2 (alias retirement) DONE.** Reference set `~/step7b_ref_binaries`
  (commit `784ee0b`, PROVENANCE inside; one directory per build + the
  helper, `REF=~/step7b_ref_binaries/build_cpu/moby_solve`). Then: the 24
  committed inis moved to `[case] file` (the `results_job*/` copies are
  records and were left alone), `[ibm] coeff_file` deleted from config.f90
  (`seen%case_file`/`seen%ibm_coeff_file` and the both-set error with it),
  `dns%ibm_coeff_file` RENAMED `dns%case_file` (20 sites), the helper's
  alias branch dropped, and every driver that wrote the alias into a
  generated ini rewritten: the conjugate/pipe drivers pass the FULL ini to
  prepare (its `[case] file` names the output; the explicit second argument
  still wins -- multilevel's zero twin relies on that), the STL-swap drivers
  insert `stl_file` under `[ibm]` instead of replacing the alias line, and
  the four that inject a case-file key use the new `tools/set_case_file.sh
  <ini> <case.h5>` (writes `file =` into `[case]`, adds the section when
  absent, replaces an existing line, never touches `[restart] file`).
  `run_ibm_les.py` got a section-aware key setter for the same reason
  (`[case] file` and `[restart] file` share a key name; its bare `^file =`
  regex would have rewritten both). Docs: configuration.md (`file` row in
  `[case]`, `coeff_file` row gone), running.md, tutorials.md, CLAUDE.md.
  GATES, all at PRODUCTION flags (nothing recomputed), ref = step7b:
  moved-ini gate (iddes_ibm, multilevel dwall/uniform, rans_geometry
  flat_l1/flat_refine incl. the ransgeom dump, ibm180, ibm180wf) `max_abs 0`
  on every dataset, CPU AND GPU; 7-case suite `max_abs 0` CPU 1 rank, CPU 4
  ranks, GPU; 9-case scalar suite `max_abs 0` CPU; P0 prepare gates 26/26
  and P1 STL gates 15/15 (both use `set_case_file.sh` through their twins);
  step-7 prepare gates vs the 7-0 INLINE reference (`~/step7_ref_binaries`,
  the only set that can still build inline): 44/44 PASS (the driver's `short_ini` now goes through `set_case_file.sh`).
- **3, 6 (campaign restarts, sibling checkouts)**: unchanged, notes only.
  The turb180 RANS channel was re-run in full at niter 12 on the consolidated
  head as the step-8 baseline (see below) and PASSES its T2 gate: u_tau
  1.0008, log-law 0.049 (tol 0.06), U+ centreline 18.16 vs DNS 18.20 -- one
  long-tier re-measurement discharged.
- **4 DONE**: `ibmwavy` at `cpu_nofma` vs `gpu_nofma`, 200 steps from ONE
  prepared case file: `max_abs 0` on un/vn/wn/pn/theta/phi. The 1.1e-16
  production-flag `theta` difference is FMA contraction, not a branch.
- **5 DONE**: `measure_nut.py` reads the case file's `coef_blocks` tiles
  (axes (x, y, z, var) per leaf, as `fdm_h5_case_append_coef` writes them --
  NOT the (z, y, x) of the snapshot blocks; this geometry varies only in y so
  the two orders are indistinguishable on it, the writer settles it). Gate:
  the ported solid mask is IDENTICAL to the legacy loader's (69632 of 327680
  cells) and the script's output on the 51 archived `a_wale` stats snapshots
  (mobydiff.scalar) is unchanged to the printed digit. `ibm_coeff.h5` (25 MB)
  left the tree; README/ini/setup comments follow.
- **7**: plan entry only, untouched.

FOUND ALONG THE WAY — first misdiagnosed, then resolved the same day:
**the B0 Blasius precursor case reproduces on `main`, and its IC generator
is fine.** Recovered from `4f819fc^` (`topbc_outlet/blasius2d.ini` == the
original `d82a97c` ini; niter 6 -> 12, `[case] file` added). My first
attempt minted the one-step TEMPLATE from the Dirichlet-top VARIANT's ini,
and the analytic-IC run then diverged by t ~ 14 (dt collapsing to 2e-6,
|p| 1e10), identically to every digit on `6db60ef`, `59ea085` and HEAD,
while the uniform start was stable. The snapshot metadata explains it: a
restart file carries its `bc_type`/`bc_value` rows and the ini overrides
only the rows it sets explicitly, so the outlet-ini run restarted with the
top face still pinned to the variant's Dirichlet u = 1, v = 1 -- a uniform
inflow through the top, fatal on any binary. LESSON, general: **mint a
restart template from the SAME ini whose faces the run will use**; a
template from a sibling variant silently transplants its boundary rows.
Re-minted from the outlet ini, the documented flow works on HEAD (F1
binary, GPU): plain Jacobi niter 12 -> theta 1.33 %, H 0.46 %, du/Ue
2.4e-3, dv 0.125 (gate PASS; B0 recorded 1.13 / 0.34 / 2.4e-3 / 0.066 on
the pre-review code at niter 6); niter 6 -> 1.72 % / 0.70 %; dt at its 0.5
cap throughout, |p| < 1e-2. The 2026-09-27 F3 note's "main is stable"
therefore stands. Scratch: the session's `scratchpad/blasius/outlet_ic/`.
PROJECTION SENSITIVITY on this case (same IC, HEAD binary, t = 2000):
red-black SOR 1.5 niter 6 -> theta 1.55 %, H 0.58 %; Chebyshev-Jacobi
niter 12 -> 1.55 % / 0.58 % (identical to four digits, i.e. the converged
projection); Jacobi 12 -> 1.33 / 0.46; Jacobi 6 -> 1.72 / 0.70. All PASS
the 2 % / 2 % / 1 % / 15 % bands; the residual differences between the
Jacobi runs are the projection's own truncation on the outlet faces. The
Chebyshev run holds dt = 0.5 to t = 2000 with |p| < 5e-3 -- the B0-era
instability (e-fold ~36 t.u.) is absent, so the ini's "NO accel" comment
is a stale warning on main.
**F1 ON THIS CASE (the laminar stretched re-validation, DONE):** the same
IC on the pre-F1 `~/step7b_ref_binaries` GPU binary gives theta 1.33 % /
H 0.45 % (niter 12) and 1.71 % / 0.70 % (niter 6) -- unchanged to the
printed digit -- and the t = 2000 fields differ by 3.4e-5 in u, 1.4e-6 in
v, 2.9e-6 in p (u scale 1). The flux-form change is invisible at the gate
level on a laminar one-sided natural line.

## The contract that now holds (one paragraph)

Every run reads a CASE FILE: `[case] file`, else `<field_prefix>.case.h5`
(a dot: gate drivers glob snapshots as `<prefix>_*.h5`). `moby_prepare
input.ini` writes it (`case_file_name`, config.f90, is the one resolver both
executables use); `moby_solve` only reads it — a missing file stops with the
`mpirun -n N moby_prepare input.ini` line printed, a stale one with the
offending `key = value` pair printed (the `case_inputs` echo attribute;
geometry rule: an ini naming no STL accepts the file's geometry lines).
Preparing and solving are separate executables with nothing repeated
between them (the user's decision of 2026-09-30; the 7-3 in-process
auto-prepare is gone). `[blocks] nb` stays optional: unset = one block per
rank of the PREPARE run's layout, so prepare on the rank count you solve
with. Gate drivers run `tools/prepare_if_missing.sh <ranks> <solver> <ini>`
before a solve: it prepares only when the resolved file is absent (never
recompute a committed or generated reference file), with the `moby_prepare`
next to the solver (none beside a pre-step-7 reference binary, which builds
inline). The rank-box layout, the legacy global-layout coefficient reader,
the inline classify dispatch and the device coefficient kernel are gone.

## Loose ends, in order

**1. The HoreKa campaign drivers were NOT updated (no cluster access from
this workstation).** Every solver launch there runs a config with no case
file and will stop at "case file not found":
`overheadTest/horeka/run_matrix.sh:103`, `exchange/run_exchange.sh:63`,
`exchange/run_portgate.sh:55`, `exchange/run_mapgate.sh:74`, and whatever
`run_timeline.sh` / `run_ncu.sh` / `setup_and_run.sh` launch. Add
`"$ROOT/tools/prepare_if_missing.sh" "$ranks" <solver> config.ini` in the
run directory before each `mpirun` — SAME rank count (`base_jacobi` and the
other nb-less configs are one block per rank by the nb rule; a file
prepared on 4 ranks and solved on 16 stops with "rank owns no blocks"), and
a fresh run directory per rank count, which the matrix already makes. The
`ref` column of a two-binary matrix runs a pre-step-7 binary: the helper
finds no `moby_prepare` beside it and does nothing, which is right. The
campaign gate scripts (`submit_*_gate.sh`) compare a `ref` and a `new`
binary on the 7-case suite through `run_portgate.sh`: same edit. Timing
note for the reports: the prepare is outside the timed step loop, but the
first-call allocation residual the 2026-09-19 note describes is unchanged.
Verify with one short `submit.sh` before any campaign number is quoted.

**2. The `[ibm] coeff_file` alias.** 22 live inis still use it (cht/channel,
cht/pipe, naca/rans, cylinder, iddes, multilevel_body, rans_geometry,
rans_sst ibm180*, refine2d, scalar/uniform3, conjugate/*, les_ibm). It
prints a note per run and is kept ONLY because the 7-0 reference binary
does not know `[case] file` and the suites run both binaries on the same
ini. Once a reference set is cut from `da64b37` or later (do that first:
`~/step7b_ref_binaries`, PROVENANCE, all four builds — the gate helper then
prepares on the ref side too, since a `moby_prepare` sits beside it), move
the inis to `[case] file` and delete the alias from config.f90
(`seen%ibm_coeff_file`, the note, the both-set error).

**3. Gates that need campaign restarts absent on this machine** (all
pre-existing, none re-gated today): scalar S3 leg (s) and S4 entirely
(`turbles_67600.h5`, from `run_gates_s2.sh les`), the cylinder gates
(`cyl_re40_20001.h5`), the les_ibm campaign, the CHT tutorials. They live in
the sibling checkouts (memory: `generated-inputs-sibling-checkouts`). When
one is run again, it runs on the separated executables for the first time:
its driver has the helper (S4, cylinder README, cht run scripts prepare
explicitly), but nobody has watched it work.

**4. The one non-zero GPU comparison.** `ibmwavy` / `ibmwavyr` on the GPU vs
the CPU reference differ in `theta` alone at 1.1e-16, and the same
difference exists between the OLD GPU binary and the CPU reference: the
scalar transport's CPU-vs-GPU production-flag arithmetic, not step 7. It
was never gated at nofma. A `cpu_nofma` vs `gpu_nofma` run of `ibmwavy`
would say whether it is FMA contraction (expected) or a real branch.

**5. les_ibm's legacy `ibm_coeff.h5` (25 MB, committed)** is kept only
because `measure_nut.py` reads its `coef` dataset. Port that script to the
block layout (`ibm_coeff_case.h5`, `coef_blocks`) and drop the legacy file
from the tree (git history keeps it; the conversion is exact, gated
`max_abs 0` at 7-3).

**6. Sibling checkouts and staging trees** (`mobydiff.scalar`,
`mobydiff.bl`, the HoreKa `moby-2to1-run` tree) predate the contract: their
drivers launch the solver without a case file. Anything re-run from there
needs the helper or an explicit prepare.

**7. Increment 7-5 (plan entry only).** A per-leaf weight column written by
prepare (body block, cut cells × the conjugate rate share, trip blocks,
refined blocks) and a weighted prefix sum replacing the closed-form
`zorder_start/count`; results stay rank-count independent because the
split moves work, never arithmetic. Not started; not needed for step 8.

## STEP 8 — F1, the flux-form viscous stencil: DONE 2026-09-30 (two re-validations owed)

Implemented as the plan entry says (`slice_grid_direction`: cell-centred
directions `lapM = 1/(hm·width)`, `lapP = 1/(hp·width)`; the staggered
direction keeps the Taylor spelling, which is the same number). Gated by the
new `src/test_stencil.f90` (`stencil_test`, 48/48 identities to ≤ 6.5e-14),
the nofma suites CPU+GPU (uniform lines round-off, stretched truncation-level
-- every number in the review's section 10 step 8 MEASURED block) and the
converged turb180 RANS channel on its natural line (T2 gate PASS, u_tau
1.0008 -> 1.0000, baseline vs F1 in the session scratch `turb180_base/` and
`turb180_f1/`). `docs/numerical-methods.md` states the flux form.
DONE the same day: (a) the Blasius precursor gate (finding above: unchanged
to the printed digit, fields 3.4e-5 apart); (b) KMM180 developed statistics,
pre-F1 vs F1 on istmcetus (two-leg runs from the archived dyw+ 0.5 restart,
30 t.u. discarded + 20 t.u. averaged, 64 000 steps each): every profile
within the 20-t.u. sampling scatter (U+ centreline 18.515 vs 18.510, u'+
peak 2.652 vs 2.623, -<u'v'>+ peak 0.7192 vs 0.7193; max |dU+| 0.041 at
y+ 21), u_tau 1.0000 both, both ~1 % from KMM DNS --
`tutorials/channel_kmm180/asset/f1_2026-09-30/`. F1 is therefore validated
on three stretched cases (turb180 RANS, Blasius laminar, KMM180 DNS) and the
operator identities; NOTHING is owed for step 8. (The archived reference
statistics were NOT used: they are the t = 0..5 cold-start transient.)
**The next reference set must be cut AFTER this commit** (the suites now
read ulps against step7b on uniform lines, so step7b is no longer a
max_abs-0 reference for anything downstream of the momentum predictor).

## (plan text, as written before the work)

Review section 9 (F1) and section 10 step 8. `slice_grid_direction` uses
the Taylor three-point `2/(hm(hm+hp))` weights for every variable, which is
the conservative flux form only where the point is the midpoint of its
control volume; on stretched lines the molecular viscous term is
non-conservative and inconsistent with the (flux-form) SGS and scalar
diffusion (KMM180 natural line: first-cell weight 0.913× the flux weight).
The fix is one formula, `lapM = d1/hm, lapP = d1/hp`; it is a NUMERICS
change (truncation-level moves on stretched cases) and the review schedules
it with the stretched-case re-validation so the RANS channel and Blasius
numbers are re-measured once. Uniform-line cases must stay `max_abs 0`
(the two forms coincide there) — that is the bit-exactness half of its gate.
Needs the new reference set of item 2 and remote GPU runs for the long
cases.

## Landmines, inherited

- Build environment: `source /etc/profile.d/lmod.sh && module load
  toolkits/nvhpc/25.9` in the SAME shell as `./compile.sh` AND as every
  `mpirun` — the system `/usr/bin/mpirun` dies in MPI_Init.
- Gate batches (memory `gate-batch-workflow`): absolute `NEW=`/`REF=`
  paths (the suites `cd` into case directories); never rebuild a build
  directory while a batch uses its binary; the 9-case suite is 1-rank only
  (`turbles`/`turbslab` pin `[mpi] dims`); `pkill -f` patterns must be
  anchored (`^python3 …`) or they kill the calling shell.
- `[mpi] dims` in a prepare input describes the SOLVE (the P0 wavy inis
  pin `1 1 1`); prepare does not check it against its own rank count.
- A case file is a build product: `*.h5` is gitignored, and the
  `tutorials/sailplane/*.h5` exception is gone. Keep it that way.
- Do not run a bare `moby_prepare` on an ini whose `[case] file` /
  `coeff_file` names a committed or generated reference file: it
  overwrites it. That is why the drivers go through the helper.
- `validation/redblack_interface/.gate_*.ini` and `.cmp_rank.txt` are that
  driver's untracked scratch (pre-existing, not gitignored).
- `tools/compare_fields.py` reassembles onto the finest lattice (OOM on
  deep-refinement snapshots): use `validation/prepare/compare_snapshots.py`.
