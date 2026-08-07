# Verification debt + the host/device staleness audit (branch `scalar`)

STATUS: **DONE 2026-08-07.** Written 2026-08-06 at the end of the session
that created the debt (commits `2c90aea`, `7d54556`, `a69a615`); both halves
were executed on 2026-08-07. Results are in §A (the audit) and §B (the gate
re-runs) below; §1 and §2 are kept verbatim as the specification they were
measured against.

**The audit found NO unprotected site in the solver.** Every host-side
consumer of a device-mapped array is protected by an explicit `target
update` or host-authoritative by construction, and each verdict was
MEASURED, not argued: the escape hatches were dropped or the host copy
perturbed in an instrumented GPU build, and the GPU output was checked to
move (or not). One documentation gap was found and fixed — the `[force]
type = custom` user hook invited a host-side `blk%q` read with no warning
that the host copy is stale inside the time loop (`bodyforce.f90`,
comment-only). The class itself is now a CLAUDE.md coding convention.

**The gate re-runs found two things, and neither was the two fixes.** Every
group reproduces its recorded numbers except the S2 `band` excess, which
turned out to be a normalisation in the METRIC rather than a change in the
solver: `cmd_band` compares raw `theta'_rms`, so the two campaigns' 3.6 %
wall-flux difference (the control had not converged) lands in the ratio;
normalised, the core matches the control to 0.05 % and the residue is a
localized 2.6 % DEFICIT at the interface — the const-1/2 restriction's known
dissipation, the opposite sign from a band. Run-down and the suggested
one-line metric fix in `validation/scalar/README.md`.

And: `run_gates_s2.sh les` NaN'd 54
steps into its statistics leg, and chasing that turned up a SOLVER-level
trap: a run that stops on `t_final` takes one extra step whose `dt` is the
accumulated round-off in `t_current` (2.46e-11 against a 1e-12 stopping
tolerance), and the projection's pressure on that step is amplified by
`1/dt` — |pn| = **1.5e6** at step 60001 against **9.1** at step 60000, with
the velocity identical to 8.5e-6. The FINAL snapshot of any
`t_final`-terminated run is therefore a bad restart. **FIXED on both
levels**: the gate driver retargets from the last PERIODIC snapshot
(`les_legs`/`last_periodic`), and `trim_dt_for_final_time` now zeroes a
`remaining` that is a negligible FRACTION of the step (`< 1e-6 dt`), so the
loop exits instead of taking it. Gated three ways in §B1 — inert on all 32
bit-exactness case-runs, effective on the exact leg that produced the
finding (10400 steps not 10401, `max|pn|` 9.09 not 1.5e6, and the last real
step bit-identical), and not over-triggering on a genuine final partial
step.

Prerequisite reading: the STATUS header of `docs/next_session_scalar.md`
(what the two 2026-08-05 fixes were), `validation/rans_sst/README.md` (the
cold-start IC write-up, which is where the landmine is documented in full),
and the "Active work" passive-scalar block of CLAUDE.md.

---

## A. Results — the host/device staleness audit (2026-08-07)

**Method.** Two passes, in this order.

1. **Static enumeration.** A script over `src/` that tracks subroutine scope,
   `!$omp target` regions and `declare target` routines, and reports every
   reference to a mapped array from HOST context (device kernels and
   declare-target helpers excluded). Run for `blk%q/qs/oldrhs`, `ibm%coef/mu`,
   `turb%nut/nut_sgs/fd`, all of `sst`'s transport arrays, `sc%wfYplus` and
   `bf%f` — i.e. every mapped array a device kernel WRITES. (Everything else
   that is mapped — the grid lines, the block metrics, the bc point lists,
   the turbulence metric tables, `sst%dwall/yeff/wallcell/wnorm/domwall`, the
   scalar config arrays — is host-filled before its map and never written on
   the device, so its host copy cannot go stale.) The enumeration came back
   SHORT: outside the init routines that fill these arrays before their map,
   the only host-context readers are the sites listed below.
2. **Probes.** An instrumented GPU build in a throwaway git worktree (never
   committed), selecting behaviour by a `MOBY_AUDIT` environment variable, so
   one build serves every probe. Each probe either DROPS an escape hatch or
   PERTURBS the host copy, and the run is compared with the same binary's
   baseline. A protected write site must MOVE when its `target update to` is
   dropped; a host-reading site must MOVE when the host copy is perturbed.

**Verdict table** (line numbers as of `fee53aa`, `src/moby_solve.f90` unless
stated):

| line | call | verdict | measurement |
|---|---|---|---|
| 137 | `read_field` | **PROTECTED** — `target update to(blk%q)` at `io.f90:890` | dropping it moves `channel_ibm` (`un` 2.08e+01) — the positive example the plan asked for |
| 139 | `zero_closed_halos` | **PROTECTED** — its own `target update to(blk%q)` | writing 3.0 instead of 0.0 into the closed halos moves `channel_ibm_refine` (`pn` 1.8e+27), so the host write does reach the device |
| 148–157 | `init_ibm` / `read_ibm_coeff_file` / `set_ibm_coeff` | **SAFE by construction** | file path fills the HOST before `enter_ibm_data`; analytic path fills the DEVICE after it |
| 162 | `init_ibm_band` | **SAFE by construction** — it reads `ibm%coef` INSIDE a target region, i.e. the device copy | analytic `wavy.ini` + `band_filter`: `nBand` = **6144 on GPU == 6144 on CPU**, although the host copy of `coef` is still all-zero at that point. A host read would have produced 0. |
| 183 | `target update from(blk%q)` | **PROTECTED / LOAD-BEARING** | dropping it moves cold-started `turb180` (`k` 2.90e-01, `omega` 5.26, `un` 2.28e-02) — the 2026-08-05 defect reproduced from the other side. Perturbing the host copy just after it ALSO moves the run (`k` 3.09e-02), i.e. `init_rans_transport` genuinely reads the host copy. |
| 194 | `target update from(ibm%coef)` | **PROTECTED, scope sufficient** | dropping it collapses the `wavy` ransgeom `wallcell` sum **2944 → 0** while `dwall`, `yeff`, `blocks`, `xc/yc/zc` stay bit-identical — exactly the arrays `classify_wall_cells` / `compute_wall_normals` need, and nothing more |
| 196 | `init_rans_geometry` | covered by 194 | (the `ibm%coef` reads live in `classify_wall_cells`) |
| 198 | `init_rans_transport` | covered by 183 | (see the 183 row) |
| 200 | `write_rans_geometry` | **SAFE — host-authoritative** | a host-only +0.125 on `sst%dwall` after `enter_rans_data` moves the dump by exactly 0.125 (so it does read the host copy), and no device kernel writes `dwall`/`yeff`, so the host copy IS the authority |
| 210 | `init_iddes_geometry` | **SAFE — host-authoritative** | a host-only +0.125 on `sst%yeff` moves `iddes180`'s `fd` by 9.99e-01 (the whole shielding function) — it reads the host copy, which nothing device-side writes. NOTE it uses `yeff`, not `dwall`: perturbing `dwall` alone is a no-op here. |
| 225 | `flow%setup_after_grid` | **SAFE — does not read `blk%q`** | a host-only +0.25 on `q` just before it: `max_abs` **0** on every dataset of `turb180` (7) and `iddes180` (8). With channel statistics ON (`min_channel`, 20 steps) the fields are still `max_abs` 0 and the stats file differs by 2.7e-12 — the same as a base-vs-base rerun (3.6e-12), i.e. the GPU atomics' accumulation order, not the perturbation. |
| 231 | `scalar_stats_setup` | **SAFE — does not read `blk%q`** | same probe; statically it reads only `[scalar]` config, the node lines and its own restart file |
| `moby_prepare.f90:154` | `target update from(ibm%coef)` before `write_case_file` | **PROTECTED** | the GPU-prepared case file equals the CPU-prepared one on EVERY dataset (`max_abs` 0), `coef_blocks` absmax 1e+28 — a missing update would have written zeros |
| `bodyforce.f90` `fill_trip` | the trip random-walk coefficients advance on the HOST each `trip_ts` | **SAFE — measured** | the suspicion was that `map(to: bf)` deep-copies `bf%trip_ak/bk` so the per-call `map(to: ak(1:nm))` would find them present and skip the refresh. It does not: with `trip_ts` shrunk to 0.005 so the walk advances 4× inside a 20-step run, CPU vs GPU is **3.9e-15** (`un`); control — the same case at `trip_ts` = 4.0 differs from it by `vn` **1.30e-01**, so the trip force is genuinely active and time-dependent |
| `bodyforce.f90` `update_bodyforce` | the `[force] type = custom` USER HOOK | **DOCUMENTATION GAP — FIXED** | the hook's comment invited reading `blk%q` "for controllers" with no warning that the host copy is stale throughout the time loop; a controller written to that comment would integrate the INITIAL field on GPU and the current one on CPU. Comment now says so and names the fix. |

Sites deliberately NOT probed, with the reason: `comm.f90`'s `sendbuf`/
`recvbuf` never have a host copy in play (GPU-aware MPI, `target data
use_device_addr`); `pressure_solver.f90`'s `phi`/`delta` are device-only
workspaces with no host reference outside their `allocate`; the `src/test_*`
drivers map nothing. `init_ibm_band` has **no committed gate case** — every
`band_filter` ini in the tree is an untracked `.`-prefixed experiment — so
its row above comes from an ini written for this audit, not from the suite.

**Conclusion: clean.** Every escape hatch in the solver is load-bearing and
every host-side consumer is either behind one or reading an array the host
owns. The one thing worth carrying forward is the METHOD, not the result:
the sharp test is to make the host-side change and check that the GPU output
moves, and it is cheap — one instrumented build with an environment switch
covers a dozen sites.

## B. Results — the verification debt (2026-08-07)

Every outstanding group was run. The per-gate numbers are in
`validation/scalar/README.md`'s "Re-gate 2026-08-07" table; the short form is
that **nothing moved because of the two 2026-08-05 fixes**, which is what §1
predicted, and that the re-run turned up ONE defect that had nothing to do
with them.

### B1. THE FINDING: the last step of a `t_final` run has a round-off `dt`, and its snapshot's `pn` is `O(1/dt)`

Seen first as a gate failure: `run_gates_s2.sh les` NaN'd 54 steps into the
statistics leg (`dt` → 0, `cfl` NaN). The cause is a solver-level trap, not a
driver quirk, and it is measured:

- `turbles.ini`'s relax leg runs from `t = 24.8` (step 49600) to
  `t_final = 30.0` at a fixed `dt = 5e-4` — exactly **10400** steps. It took
  **10401**.
- After step 60000, `t_current = 29.999999999975433`: the accumulated
  round-off of 10400 additions of 5e-4 is **2.46e-11**, while
  `run_should_continue`'s stopping tolerance is
  `max(1e-12, 100 eps |t_final|)` = **1e-12** (the floor wins: 100 eps 30 =
  6.66e-13). So the loop continues,
  `trim_dt_for_final_time` sets `dt = 2.46e-11`, and one more step runs.
- That step's projection solves for a pressure correction against
  `dt_gamma ~ 1e-11`, so the stored `pn` is amplified by `1/dt`:
  **|pn| = 1.5e6 at step 60001 against 9.1 at step 60000**, with the
  VELOCITY identical to 8.5e-6. The post-loop `write_field` stores exactly
  that field, so the FINAL snapshot of any `t_final`-terminated run is
  unusable as a restart — restarting from it blows the run up.

This is the true mechanism behind the landmine recorded in §1 (and in the S2
notes) as "the campaign's FINAL `*_50001.h5` carries |pn| ~ 1e6 of the
niter = 6 pn-drift mode". It is NOT the drift mode: the drift mode is the A2
velocity-neutral oscillation in `pn` on IBM runs at `niter = 6`, and it is
not what happens here (`turbles` has no body). Both make the final snapshot
a bad restart; only one of them is a two-line arithmetic problem.

**Fixed at the gate level** (this session): `run_gates_s2.sh` gained
`last_periodic`, and `les_legs` retargets from the last PERIODIC relax
snapshot rather than from `newest`, which picked the post-loop write.

**FIXED IN THE SOLVER 2026-08-07** (`step.f90`, `trim_dt_for_final_time`),
after the write-up above was corrected — the first version of the proposed
fix was an ABSOLUTE test and would not have caught this:

```fortran
remaining = dns%t_final - dns%t_current
if (remaining < FINAL_STEP_FRACTION*dns%dt) remaining = 0.0d0   ! 1.0d-6
dns%dt = min(dns%dt, max(0.0d0, remaining))
```

The loop's existing `if (dns%dt <= 0.0d0) exit` then fires. The test is
RELATIVE to the step because the round-off in `t_current` grows with the step
count (~`N eps t_final`, 6.9e-11 here) and outruns any fixed floor, while
`remaining/dt = 4.9e-08` identifies it scale-free. It is also start-value
dependent — each addition rounds at ~`eps t`, so summing 10400 x 5e-4 from
`t = 24.8` reproduces the shortfall while summing from 0 does not.

**Gated three ways.**

- **(A) Inert where it must be.** Both bit-exactness suites, CPU and GPU, vs
  `~/s5c_ref_binaries`: `run_bitexact.sh` 7/7 and `run_bitexact_s3.sh` 9/9,
  **max_abs 0** on all 32 case-runs. Every case in both suites is
  `nsteps`-terminated (`short_ini` sets `t_final = 0.0`), so the new branch
  is never taken there — which is exactly the point.
- **(B) Effective where it must be.** The `turbles` relax leg, the exact run
  that produced the finding, re-run with the fix: **10400 steps, not 10401**;
  the final snapshot is step **60000**, not 60001; `max|pn| = 9.093` against
  **1.5064e+06** before. And the step-60000 snapshot is **max_abs 0** against
  the pre-fix run's own step-60000 snapshot on `un vn wn pn nut theta` — the
  fix changes nothing up to the last real step, it just declines to take the
  fake one. The final snapshot of a `t_final` run is now a normal restart.
- **(C) Not over-triggering.** A GENUINE final partial step still runs:
  `wave.ini` (fixed `dt = 1e-3`) with `t_final = 0.0105` takes 11 steps, the
  last of length 5.0e-4 — half a step, 5e5 times the threshold — and ends at
  `t_current = 0.0105` exactly.

### B3. The archived GPU reference binary was itself stale (found 2026-08-07, while gating the phase timer)

A tail of the same defect class, one level up. `~/s5b_ref_binaries/`'s
PROVENANCE says the set carries both 2026-08-05 fixes. That is true of its
CPU binary and **false of its GPU binary**: it was archived at 16:42 on
2026-08-05, and the tree's own `build_gpu` is 16:46 — the four minutes in
which the `!$omp target update from(blk%q)` was added. So the archived GPU
reference is exactly the "correct on CPU, pure no-op on GPU" intermediate
state that §A's site 183 is about.

Proven three ways, all on `turb180`, 20 steps, GPU:

- the `s5a` and `s5b` GPU binaries produce IDENTICAL fields (`max_abs` 0 on
  all seven datasets) although the two FILES differ in size and md5 — so the
  s5b build predates the fix even though it is a distinct build;
- a GPU nofma binary built fresh from commit `8f60944` differs from the s5b
  one by `un` 2.2814559502704945e-02, `k` 2.8996297594464926e-01, `omega`
  5.2618740388259901e+00;
- those are the same numbers, to 13 digits, that §A's `no_q_from` probe
  produced by deleting that one directive.

It cost nothing this time only because the failure was loud: the phase
timer's GPU suite failed on exactly the three cold-started RANS cases and
nowhere else, which is the fingerprint. Against a quieter change it would
have been read as a regression in the change.

**Fixed:** `~/s5c_ref_binaries/` is the reference set now — CPU+GPU nofma
solve+prepare, built from a clean `git worktree` at commit `8f60944`, hash
recorded in its PROVENANCE, and cross-checked by the 7-case suite passing
`max_abs` 0 against it on CPU **and** GPU. `~/s5b_ref_binaries/PROVENANCE.txt`
carries a correction notice. **LESSON: cut a reference set from a COMMIT,
never from a working tree mid-edit, and record the hash.**

### B2. What the groups measured

| runner / group | outcome |
|---|---|
| `run_gates.sh uniform conserve conduction wave pr restart det` | ALL PASS, every number the recorded one (`cond_16/32/64` L2 = 2.814767e-03 / 7.054808e-04 / 1.764823e-04 and `wave` L2 = 1.554700e-02 / 3.909091e-03 / 9.786703e-04, both to every digit; `det` 1 == 4 ranks AND CPU == GPU `max_abs` **0** on `un vn wn pn s1`) |
| `run_gates_s2.sh les` | PASS — `theta_tau` 0.050841 (recorded 0.050860), `theta+/U+` 0.7204 / 0.8563 (0.7233 / 0.8558), Kader mean 0.0683 (0.0682); see B1 for what had to be fixed first |
| `run_gates_s2.sh band` | PASS — excess **−0.0265** (recorded +0.0012). Run down (`validation/scalar/README.md`, "Where the S2 `band` excess comes from"): the metric is a RAW rms ratio, and the core offset **1.0365 IS the two runs' wall-flux ratio 1.0360**; normalised by each run's own `theta_tau` the core matches to **1.0005** (the recorded 1.0002) and what remains is a localized **−2.6 % DEFICIT** at the interface — the const-1/2 restriction's known dissipation, and the opposite sign from a spurious band. The flux gap is the CONTROL's residual drift (`turbles` `theta_tau` 0.051343 → 0.050341 across the window; `turbslab` steady to 2e-5) |
| `run_gates_s2.sh det` | 1 == 4 ranks and CPU == GPU **max_abs 0** on `un vn wn pn nut theta` |
| `run_gates_s3.sh solid conserve balance prep refine missing det cyl` | ALL PASS; `det` **max_abs 0** on all six datasets for 1 vs 4 ranks, CPU vs GPU, the LES variant (7 datasets) and the file-IBM variant; `cyl` Nu 3.3655 with the CV cross-check at 0.00 % |
| `run_gates_s4.sh stats accum plane levels restart det noeffect cyl tools` | ALL PASS; `levels` (2.563e-14 / 2.122e-13) and `cyl` (Nu 3.3653, Q/Lz 3.722690e-01) reproduce to every digit. NOTE `s4 cyl` must run AFTER `s3 cyl` — it needs that group's snapshots and silently SKIPS otherwise (it did, on the first pass) |
| `run_gates_s5.sh det` | 1 == 4 ranks **max_abs 0** on nine datasets; CPU vs GPU 0.0 on eight and **5.55e-17** on `theta_kc` — the recorded 5.6e-17 |

## 0. Why this session exists

Two solver fixes landed on 2026-08-05:

1. `src/moby_solve.f90` — `apply_bc` + `exchange_halos` (+ a
   `target update from(blk%q)`) before the `[rans]` init, because
   `init_rans_transport`'s `k = 1.5 (tu/100 |u|)^2` initial condition read an
   unfilled velocity halo and came out a factor 4 low on the last plane of
   every block.
2. `src/modules/scalar_stats.f90` — solid cells excluded from the
   penalization accumulator, because the body-heat diagnostic double counted
   at `ibm_value = 0`.

Fix 1 CHANGES the k/omega initial condition of every cold-started RANS run,
so a lot of recorded gate numbers had to be re-measured. Most were, and were
unmoved. **What was not re-measured is listed in §1**, together with why each
is expected to be safe — the point of the session is to convert those
arguments into measurements, because the whole reason this branch is in good
shape is that it has never accepted an argument where a number was available.

And the second half of the session is §2: fix 1 was *correct on CPU and a
pure no-op on GPU* for a week's worth of reasoning, because host-side init
code reads a stale host copy of a device-mapped array. That is a CLASS of
defect, and exactly one instance of it has been found. Nobody has looked for
the others.

---

## 1. The verification debt

Every runner lives in `validation/scalar/` and takes a group argument. What
was re-run on 2026-08-05/06 after the fixes, and what was not:

| runner | groups | re-run? |
|---|---|---|
| `run_gates.sh` (S1) | `uniform conserve conduction wave pr det` | **none** |
| `run_gates_s2.sh` | `kays wferr sst les band det` | kays, wferr, sst ✔ — **les, band, det outstanding** |
| `run_gates_s3.sh` | `solid conserve balance prep refine det missing cyl` | **none** |
| `run_gates_s4.sh` | `stats accum plane levels restart det noeffect heat adia cyl tools` | heat, adia ✔ — **the rest outstanding** |
| `run_gates_s5.sh` | `unit ref sweep det` | unit, ref, sweep ✔ — **det outstanding** (an equivalent was run by hand: 1 vs 4 ranks incl. the x split, and CPU vs GPU, both from `wfs180_y30.ini`) |
| `run_bitexact.sh` | 7-case suite | ✔ CPU **and** GPU |
| `run_bitexact_s3.sh` | 9 scalar cases | ✔ CPU |

**Why each outstanding group is expected to pass, i.e. what a surprise would
mean.** Read these before running, so a failure is interpretable:

- **S1 (`run_gates.sh`)** — every case is turbulence-free and body-free, so
  neither fix can reach it. `run_bitexact_s3.sh` already proves the stronger
  statement (max_abs 0 against the pre-fix binaries) on overlapping inis. A
  failure here would mean one of the two fixes touches a path it has no
  business touching.
- **S2 `les` / `band`** — LES campaigns, restart-based, no `[rans]` section,
  so the IC fix cannot reach them; the heat fix is not called at all
  (`heat_interval` off). These are the two most expensive groups in the tree;
  budget them first or last deliberately, not by accident. LANDMINE (cost the
  S2 session an hour): restart from `*_49600.h5`, never from the campaign's
  FINAL `*_50001.h5` — those carry `|pn| ~ 1e6` of the niter = 6 pn-drift mode
  and blow up in one step.
  *(2026-08-07: the attribution in that last sentence is WRONG, and §B1
  above has the measured mechanism — the final snapshot is the post-loop
  write of a step whose `dt` was trimmed to the accumulated round-off in
  `t_current`, so its `pn` is amplified by `1/dt`. Same symptom, different
  cause, and it bites every `t_final`-terminated run, body or not.)*
- **S3 (all groups)** — analytic/file IBM without RANS, so the IC fix cannot
  reach them. The heat fix is a measured no-op at `ibm_value = 1`, which is
  what every S3 case uses; `balance` and `surface` are Python-side and
  unaffected. NOTE `run_gates_s3.sh prep` regenerates case files with
  `moby_prepare` — use the CPU build (canonical; the GPU build computes coef
  on device and differs by libm ulps).
- **S4 (the rest)** — statistics only; `stats`/`accum`/`plane`/`levels`
  compare the solver's rows against the snapshot's, which is a self-consistency
  identity that neither fix perturbs.
- **S5a `det`** — the hand-run equivalent passed at max_abs 0 including the x
  split, and CPU vs GPU came out at 5.6e-17. Running the group closes the
  bookkeeping; it needs the nofma binaries.

**Conventions.** `~/s5c_ref_binaries/` is the reference set (CPU+GPU nofma,
solve+prepare, commit-pinned to `8f60944`, PROVENANCE inside).
`~/s5a_ref_binaries/` is STALE for cold-started RANS cases, and so is the
GPU binary of `~/s5b_ref_binaries/` — see B3. Build both paths with the module loaded; always
`mpirun`. The re-run landmines (RANKS pinning, cetus' missing h5py, `pgrep`
self-matching, the `det` tolerance-0 nofma requirement) are listed in §12 of
`docs/next_session_scalar.md` — read them, they cost the last session real
time.

**Recording.** Numbers that reproduce need one line each in
`validation/scalar/README.md`'s re-gate table; a number that MOVES is the
finding of the session and needs the full treatment (what moved, by how much,
which fix, whether the converged answer or only a transient).

---

## 2. The host/device staleness audit — the part that matters

**The class.** `enter_*_data` maps a derived type's arrays to the device.
After that, the host copy and the device copy are independent. Any host-side
code that READS a mapped array sees whatever the host last wrote — device
kernels do not update it — and any host-side code that WRITES one is invisible
to the device until an explicit update. There are exactly two escapes, both
already used in `moby_solve.f90`:

```fortran
!$omp target update from(x)   ! device -> host, before host code reads x
!$omp target update to(x)     ! host -> device, after host code writes x
```

**The one known instance** (fixed 2026-08-05, `src/moby_solve.f90`):
`init_rans_transport` is host code reading `blk%q`, which
`enter_block_data` mapped ~60 lines earlier. The `apply_bc` +
`exchange_halos` inserted to fix the k IC are DEVICE kernels, so on GPU they
wrote the device copy and the host copy stayed zero-halo'd. The fix looked
complete on CPU and was a pure no-op on GPU: the "fixed" GPU binary
reproduced the pre-fix result BIT-FOR-BIT, and `k` kept its x-dependence
(x-spread 1.49e-02, against 0.0 on the CPU). Only running the GPU
bit-exactness suite exposed it.

**The audit.** Walk `src/moby_solve.f90`'s init block in order and, for every
host-side call after a map, establish which of the three states it is in:
protected (an explicit update precedes/follows it), safe by construction (the
array is host-filled and nothing device-side has touched it yet), or WRONG.
The call sites, with their line numbers as of `a69a615` and what makes each
one interesting:

| line | call | why it is on the list |
|---|---|---|
| 137 | `read_field` | writes `blk%q` on the HOST **after** `enter_block_data` (133) — the converse direction. GPU restarts are gated and work, so it must already update to the device; CONFIRM that and cite it as the positive example |
| 139 | `zero_closed_halos` | writes halos; establish whether it is host or device code |
| 162 | `init_ibm_band` | reads `ibm%coef` after `enter_ibm_data` (151/153); its comment says "from the device coefficients" |
| 194 | `target update from(ibm%coef)` | the EXISTING precedent, protecting `init_rans_geometry` — check its scope is actually sufficient |
| 200 | `write_rans_geometry` | runs AFTER `enter_rans_data` (199) and dumps host `dwall`/`yeff`/`wallcell` |
| 210 | `init_iddes_geometry` | reads `sst%dwall`/`sst%yeff`, mapped at 199 |
| 225, 231 | `setup_after_grid`, `scalar_stats_setup` | run after the init exchanges; if either reads `blk%q` on the host it is reading a copy that is stale again (the `target update from(blk%q)` at 183 is BEFORE the later exchange) |

Then widen it: `src/moby_prepare.f90` (its GPU build computes coefficients on
the device, then host code writes the case file), and any `enter_*_data` in
`gpu_runtime.f90` / `blocks.f90` whose arrays are touched host-side
afterwards. `grep -n "target update"` over `src/` gives the current set of
escapes; anything reading a mapped array without one nearby is a candidate.

**How to test a candidate rather than argue about it.** The argument is what
failed last time. Build CPU and GPU, run the case both ways, and compare —
a host/device staleness bug shows up as *the GPU result not moving when the
host code changes*, which is why the sharp test is: make the host-side change,
and check that the GPU output CHANGES. If a candidate has no case that
exercises it, say so explicitly rather than marking it safe.

**Deliverable.** A table of the call sites with a verdict each, any fixes with
their gates, and — if the audit finds nothing — that stated plainly, because
"we looked and it is clean" is a result worth recording. Add the class itself
to CLAUDE.md's conventions if it is not already sharp enough there: *host-side
code reading a device-mapped array after `enter_*_data` sees a stale copy.*

---

## 3. Next-session prompt

> Two jobs on branch `scalar`, in this order — read
> `docs/next_session_verification.md` in full first, then the STATUS header
> of `docs/next_session_scalar.md` and the "Active work" passive-scalar block
> of CLAUDE.md.
>
> **(a) The host/device staleness audit (§2) — do this FIRST, it is the half
> with a real chance of finding something.** On 2026-08-05 a fix to the RANS
> cold-start initial condition was correct on CPU and a PURE NO-OP on GPU,
> because `init_rans_transport` is host code reading `blk%q` while the
> `apply_bc`/`exchange_halos` inserted to fix it are device kernels — the
> "fixed" GPU binary reproduced the pre-fix result bit-for-bit. That is a
> class of defect and exactly one instance has been found. §2 lists the call
> sites in `moby_solve.f90`'s init block with line numbers and why each is on
> the list; classify every one as protected / safe-by-construction / wrong,
> test candidates by making a host-side change and checking the GPU output
> MOVES (arguing is what failed last time), then widen to `moby_prepare.f90`
> and the other `enter_*_data` consumers. Report a verdict table; "clean" is
> a result, and so is "this call site has no case that exercises it".
>
> **(b) The verification debt (§1).** Fix 1 changed the k/omega IC of every
> cold-started RANS run, so the gate groups in the table there were never
> re-measured after it. Run them: `run_gates.sh` (all), `run_gates_s2.sh
> les band det`, `run_gates_s3.sh` (all), `run_gates_s4.sh` minus
> heat/adia, `run_gates_s5.sh det`. §1 records why each is EXPECTED to pass —
> read that first so a failure is interpretable, and treat any moved number
> as the finding of the session rather than as noise to be tolerated.
>
> Conventions that are not negotiable: the reference set is
> `~/s5b_ref_binaries/`, NOT `~/s5a_ref_binaries/` (stale for cold-started
> RANS);   *(2026-08-07: use `~/s5c_ref_binaries/` — s5b's GPU binary turned
> out to be stale too, see B3)* build both paths with `module load toolkits/nvhpc/25.9`; always
> `mpirun`; the `det` groups compare at TOLERANCE 0 and need the nofma
> binaries. The re-run landmines — `turbles`/`turbslab`/`turbsst` pin
> `dims = 1 1 1` so they need `RANKS=1`; istmcetus has no h5py (run solver
> legs there via `validation/scalar/env_cetus.sh`, compare locally);
> `ps -C moby_solve` does not match the archived reference binaries and
> `pgrep -f`/`pkill -f` match the shell running them, so kill by PID — are in
> §12 of `docs/next_session_scalar.md`. Never declare a group done with a
> failing build or an unrecorded number.
>
> Record reproduced numbers in `validation/scalar/README.md`'s re-gate table
> and update this document's STATUS header. Do NOT start S5b, S5c or
> conjugate C1 in the same session.
