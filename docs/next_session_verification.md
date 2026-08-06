# Verification debt + the host/device staleness audit (branch `scalar`)

STATUS: **NOT STARTED.** Written 2026-08-06 at the end of the session that
created the debt (commits `2c90aea`, `7d54556`, `a69a615`). One session, two
halves, in this order — the audit is the part with a real chance of finding
something, so do not let the gate re-runs eat the whole slot.

Prerequisite reading: the STATUS header of `docs/next_session_scalar.md`
(what the two 2026-08-05 fixes were), `validation/rans_sst/README.md` (the
cold-start IC write-up, which is where the landmine is documented in full),
and the "Active work" passive-scalar block of CLAUDE.md.

---

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

**Conventions.** `~/s5b_ref_binaries/` is the reference set (CPU+GPU nofma,
solve+prepare, PROVENANCE inside); `~/s5a_ref_binaries/` is STALE for
cold-started RANS cases. Build both paths with the module loaded; always
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
> RANS); build both paths with `module load toolkits/nvhpc/25.9`; always
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
