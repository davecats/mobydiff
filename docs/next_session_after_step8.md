# Next session: after numerics review step 8 — reference set, cluster check, then step 9

STATUS: **NOT STARTED (handout written 2026-09-30, HEAD `35864a3`).** The
previous handout (`docs/next_session_after_step7.md`) is closed: its loose
ends 1, 2, 4, 5 and step 8 (F1) are DONE and gated, 3, 6, 7 are notes.
Read `CLAUDE.md` (the step-7 bullet states the case-file CONTRACT; the
"F1 IS DONE" paragraph states what step 8 changed), then work the items
below IN ORDER. Item 1 is a precondition for every later gate; item 2 is
the one thing the previous session could not do from this workstation.
Do not start steps 10-12.

## What holds now (one paragraph each)

**The case-file contract.** Every run reads `[case] file` (default
`<field_prefix>.case.h5`), written by `mpirun -n N moby_prepare input.ini`
on the SAME rank count as the solve (an unset `[blocks] nb` is one block per
prepare rank). `moby_solve` stops on a missing or stale file. The
`[ibm] coeff_file` alias is GONE (2026-09-30): every ini and driver uses
`[case] file`, drivers that derive a solve ini call
`tools/set_case_file.sh <ini> <case.h5>`, and every gate driver runs
`tools/prepare_if_missing.sh <ranks> <solver> <ini>` first. A local
reference set is a directory per build PLUS `tools/prepare_if_missing.sh`
(`~/step7b_ref_binaries` is the pattern; the flat `moby_solve_cpu` sets are
only good for pre-step-7 binaries that build inline).

**F1, the flux-form viscous stencil (step 8, `0737b88`).** In a
cell-centred direction `slice_grid_direction` builds `lapM = 1/(hm*width)`,
`lapP = 1/(hp*width)`; the face-staggered direction keeps the Taylor
spelling, which is the same number there. Gated by `stencil_test` (48/48
telescoping + symmetry identities), the nofma suites (uniform lines
ROUND-OFF, u <= 8e-14 / p <= 1.5e-12, four cases exactly 0; stretched cases
truncation-level), the converged turb180 RANS channel (T2 PASS, u_tau
1.0008 -> 1.0000), the B0 Blasius laminar case (unchanged to the printed
digit) and KMM180 developed statistics (within the 20-t.u. sampling scatter,
`tutorials/channel_kmm180/asset/f1_2026-09-30/`). CONSEQUENCE for gating:
**`~/step7b_ref_binaries` is no longer a max_abs-0 reference for anything
downstream of the momentum predictor** -- it reads ulps on uniform lines and
truncation-level on stretched ones. Hence item 1.

## Items, in order

**1. Cut the post-F1 reference set.** `~/step8_ref_binaries` from HEAD
(`35864a3`; source last changed at `0737b88`), all four builds
(`./compile.sh {cpu,gpu,cpu_nofma,gpu_nofma}` in a shell with the module
loaded, from a clean tree -- `git status --porcelain -- src tools
CMakeLists.txt compile.sh` empty), laid out as `build_<mode>/{moby_solve,
moby_prepare}` + `tools/prepare_if_missing.sh` + a PROVENANCE.txt naming the
commit and the flags (copy `~/step7b_ref_binaries/PROVENANCE.txt` and
edit). Do NOT rebuild any build directory while a gate batch uses it
(memory `gate-batch-workflow`). Everything after this gates against
`~/step8_ref_binaries`; retire the older sets from the READMEs' examples as
you touch them.

**2. Verify the HoreKa drivers on the cluster.** `run_matrix.sh`,
`exchange/run_{exchange,portgate,mapgate,timeline,ncu}.sh`,
`submit_rdblk.sh` and `submit_port_gate.sh` pass 2 now prepare the case file
in the run directory before each `mpirun` (helper looked up beside the
binary's own checkout). Exercised only LOCALLY on the CPU build. Run ONE
short `submit.sh` (or `submit_matrix5.sh` with `NSTEPS=20`, `CONFIGS=base_jacobi
rect_jacobi`) and confirm `overhead.case.h5` + `prepare.log` appear beside
`run.log`, the leaf count printed by prepare matches the rank count's
expectation (nb-less configs: one block per rank), and the step time matches
the 2026-09-26 matrix to within its noise. Only then may a campaign number be
quoted again. If the `ref` column is a pre-step-7 checkout it builds inline
and its run directory has no case file -- that is correct, not a failure.
ALSO on the cluster: the momentum-predictor regression of the consolidated
head (1.3-2.0 % of the step, `results_horeka_2026-09-26.md`; the compiler
spent the freed registers on longer dependency chains) has a one-line fix
proposed and NOT taken -- a mapped always-true logical around the three
corrections -- which needs its own before/after ncu + nofma suite. Same
allocation, if the queue allows.

**3. Step 9 — F7, the exact penalization factor** (review section 10; its
text is authoritative). `update_ibm_mu` produces two factors, `muA = B/(lam
dt + B)` on `q` and `muB = 1/(lam dt + B)` on the substage increment, with
`B = lam dt / expm1(lam dt)` and `B = 0` for `lam dt > 700`; the predictor
applies them; the projection keeps `muB`. Body-free: `lam = 0` gives
`B = 1`, both factors exactly 1.0, so body-free cases are bit-exact BY
CONSTRUCTION -- gate that at max_abs 0 on the 7-case and 9-case suites
minus the body cases, against `~/step8_ref_binaries`, CPU and GPU, at
production flags if the body-free arithmetic is textually untouched, nofma
otherwise. Body cases move: re-measure the cylinder Re 40 C_D (its case
file must be re-prepared through `validation/cylinder/setup.sh`; the
Re 40 restart `cyl_re40_20001.h5` lives in a sibling checkout -- memory
`generated-inputs-sibling-checkouts`) and the les_ibm law of the wall
(`validation/channel_interface/les_ibm/`, 4000 steps on a GPU), and run a
`dt`-halving study on the cylinder drag: the cut-cell time error should
drop faster than before. One session plus a remote GPU.

**4. (When time allows) Put the B0 Blasius case back in the tree.** It was
deleted in the 2026-09-26 tutorial consolidation and recovered from
`4f819fc^` (`topbc_outlet/blasius2d.ini` == the original `d82a97c` ini);
it runs and passes on main (theta 1.33 %, H 0.46 % at Jacobi niter 12;
red-black niter 6 and Chebyshev niter 12 both 1.55 % / 0.58 %). Re-add
`blasius2d.ini`, `template.ini` (minted from the SAME ini -- see the
landmine below), `make_blasius_ic.py`, `compare_blasius.py` (fix its
`sys.path` insert to `tools/`) under `tutorials/turbulentBoundaryLayer/`
or `validation/blasius/`, with its README stating the four-projection
table; it is the laminar stretched-line gate the F1 work used. Update the
stale "NO accel = chebyshev" comment in the ini: Chebyshev is stable here
to t = 2000 at dt 0.5 on main (the B0-era mode is gone, F3 retraction).

**5. Notes carried forward, unchanged.** (a) Gates that need campaign
restarts absent on this machine (scalar S3 leg (s) and S4, the cylinder
gates, the les_ibm campaign, the CHT tutorials) live in the sibling
checkouts; each runs on the separated executables for the first time when
it is next run -- its driver has the helper, nobody has watched it work.
(b) Sibling checkouts and the HoreKa `moby-2to1-run` tree predate the
contract; anything re-run from there needs the helper or an explicit
prepare. (c) Increment 7-5 (per-leaf weight column + weighted Morton split)
is a plan entry, not started. (d) The one-by-one re-measurement pass
(validation/README.md is the checklist; the long tier is still at niter-6
values) resumes after step 9; this session discharged turb180 (T2 PASS at
niter 12) and the Blasius case from that list.

## Landmines, inherited and new

- **A restart file carries its boundary rows** (`bc_type`/`bc_value` in the
  snapshot metadata) and the ini overrides only the rows it sets explicitly.
  Mint a restart template from the SAME ini whose faces the run will use: a
  template from a sibling variant silently transplants its boundary rows
  (the outlet-top Blasius case restarted with a uniform v = 1 inflow through
  the top and diverged by t = 14 on every binary; one afternoon lost).
- **Two separate `mpirun -n 1` on istmcetus both bind to core 0** and one
  rank then sits on the wrong NUMA node (GPU0 = cpus 0-63, GPU1 = 64-127):
  67-71 % GPU utilisation instead of 98 %. Launch with `--bind-to none
  numactl --cpunodebind=<node>` (or `taskset -c`), or re-pin a running rank
  with `taskset -pc` (safe: affinity never touches the arithmetic). The
  memory `remote-hosts` has the measured per-step costs: KMM180 0.021 s/step
  on the RTX 5090, 0.054 (pinned: ~0.04) on one A6000.
- **The archived KMM180 statistics file is the t = 0..5 cold-start
  transient**, not a developed reference; a comparison against it measures
  the transient. Use the 20-t.u. developed windows in
  `tutorials/channel_kmm180/asset/f1_2026-09-30/` (pre-F1 and F1) as the
  reference pair from now on, and restart from the archived dyw+ 0.5 twin
  with 30 t.u. discarded.
- `*.h5` and `*.png` are gitignored tree-wide: an asset you mean to commit
  needs `git add -f`, and `git add <dir>` skips them SILENTLY.
- Build environment: `source /etc/profile.d/lmod.sh && module load
  toolkits/nvhpc/25.9` in the SAME shell as `./compile.sh` and every `mpirun`.
- Gate batches (memory `gate-batch-workflow`): absolute `NEW=`/`REF=` paths;
  never rebuild while a batch runs; the 9-case suite is 1-rank only; anchor
  `pkill -f` patterns with `^` (an unanchored one kills the calling shell,
  exit 144 -- it did, twice today); a `VAR=... && cmd &` list backgrounds the
  assignment too.
- `run_gates_step7.sh` compares against the reference's INLINE solve, so its
  REF must be the 7-0 set `~/step7_ref_binaries/moby_solve_cpu`, the only
  one that still builds inline; every other gate takes the step8 set.
- `[mpi] dims` in a prepare input describes the SOLVE; prepare does not
  check it against its own rank count.
- Do not run a bare `moby_prepare` on an ini whose `[case] file` names a
  committed or generated reference file: it overwrites it.
- `tools/compare_fields.py` reassembles onto the finest lattice (OOM on
  deep-refinement snapshots): use `validation/prepare/compare_snapshots.py`.
- `validation/redblack_interface/.gate_*.ini` and `.cmp_rank.txt` are that
  driver's untracked scratch (pre-existing).
