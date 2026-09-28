# Next session: numerics review, steps 0–5

STATUS UPDATE (2026-09-27, second session): **STEP 6 DONE** on top of 0–5 —
the boundary rows are ONE affine write from tables resolved at init
(`bcKind/bcW/bcC`, velocity + pressure + scalar columns), `apply_bc` is
generic over a variable list, `apply_scalar_bc_q` is deleted, and the two
remaining kernels carry no geometry and no type branch. **REDUCED SCOPE,
recorded in the review (section 10 step 6 MEASURED + the section 4 header):
the rows did NOT become exchange entries**, because the same-level copy
entries read physical ghosts and pack runs first, so a BC launch before the
exchange is needed in any design — the "zero kernels" half of the plan was
unreachable and the abstraction half is what shipped. 61 comparisons at
`max_abs 0` (7-case CPU 1+4 ranks + GPU, 9-case CPU + GPU, freestream
outlets/parabola + red-black 2:1 + RANS inlet CPU 1+4 ranks + GPU, all
nofma; plus the 7-case suite at production flags). The pipe re-run below
FINISHED 2026-09-28 17:28 and step 2d is CLOSED: the production tables were
re-measured (review section 10 step 2 MEASURED (d)); the fluid-side thermal
statistics and the interface heat balance are within run-to-run scatter, the
solid-side decay moved toward the reference (+9 → +3 % deep excess), and the
pipe tutorial's `asset/` now ships the re-run's accumulators and figures with
a provenance note (README) and an addendum (report). Steps 7–12 additionally inherit: `apply_bc(blk, bc,
vars, outflow_copy)` is the one entry point for q-variable ghosts, so a new
transported q variable gets its BCs by filling a column
(`set_scalar_bc_rows`) — no new kernel; and NEVER run two `run_bitexact*.sh`
suites with the same `MODE` label concurrently in one tree (their output
prefixes collide and one deletes the other's snapshot).

STATUS: **STEPS 0–5 DONE (2026-09-27, one session; commits e126d1a → the
step-5 commit on `main`).** One item is still in flight and one outcome
overturned a prediction:

- **Step 0** — `~/numrev_ref_binaries/` (8 binaries, PROVENANCE, hash
  `6db60ef`). Use it for every comparison on `main` from here on.
- **Step 1 (F3) — CONTRADICTED, retracted.** On the B0-era binary
  (`5d44eb4`) the recorded instability reproduces (e-fold 33–37 t.u.) and
  `cheb_lmax = 2.2` grows FASTER (27–32 t.u.); `main` is stable at niter 6
  with either bound; the min_channel control cannot tell the bounds apart.
  No code change. The one physics change in the `5d44eb4 → main` span is the
  momentum convection form (divergence then, skew now) — a correlation, not a
  proof. `docs/numerical-methods.md` was NOT touched.
- **Step 2 (F2) — confirmed, default changed.** Conjugate scalars run the
  advective form (unset key resolves; explicit `divergence` refused); drift
  4.09e-5 → 0.0 on the wavy case. Suites 9/9 CPU + GPU and 7/7 CPU max_abs 0;
  C1 moved exactly its conservation gate (1.22e-16 → 2.235e-11, re-baselined
  as a leak band, 1e-9); C2, C3 all PASS. **IN FLIGHT: the pipe Nusselt
  re-run**, leg D from the 2026-09-20 settled state with the new binary, on
  istmcorax in `~/pipe_rerun_f2/` (`p_f2_stat.*`, 168 000 steps at 0.56
  s/step, started 2026-09-27 17:10, expected to end ~2026-09-28 19:00). When
  it ends: `tutorials/cht/pipe/pipe_stats.py` on the `p_f2_stat_*.h5`
  snapshots, `asset/compare_neuhauser.py thermal … --scalar c0` etc. against
  the old `p_pr_statsD.npz` and the Neuhauser reduction, then refresh
  `asset/` (make_caches.py) and the report. The run keeps `niter = 6` so
  the binary is the only change; the niter-12 re-measurement is the
  one-by-one pass's job.
- **Step 3 (F4)** — sailplane on `x_max_patch = outlet`; Neumann normal
  velocity is a config error. Prepared-vs-legacy 1-step max_abs 0; 200 steps
  bounded on plain Jacobi (sor 0.8) AND Chebyshev niter 12 on the GPU —
  red-black is no longer required (kept in the ini as the faster choice);
  freestream unchanged.
- **Step 4 (F5)** — one `apply_bc` per projection in both solvers; max_abs
  0 at PRODUCTION flags on the 7-case suite (CPU 1 + 4 ranks, GPU),
  freestream and redblack_interface (1 + 4 ranks); `proj_timing: apply_bc`
  −91.6 % on min_channel and pois_io.
- **Step 5 (F8, F6)** — `cflmax` documented as per-direction with an
  init-time `cfl:` print of the worst-case directional sum (inert: 7-case
  suite max_abs 0 at production flags, CPU and GPU). KMM180: the archived run
  HAD stepped at the Péclet cap `dt = 7.696e-6` (650 k steps for 5 h/u_τ);
  the ini is now `natural_dyw_plus = 0.5` and steps at `dtmax = 3.125e-4`
  (cfl 0.13, Péclet 0.23; 100-step GPU runs of both spacings, README
  written). **5b's verdict: the momentum half of step 11 has no customer.**

What steps 6–12 inherit: (i) the pipe result above (finish step 2d before
touching the CHT tutorials); (ii) `main`'s Chebyshev projection has NO
measured instability to fix, so any `lmax` change is a numerics preference
without a gate — leave it; (iii) the generated inputs the suites need
(`IC_turbles.h5`, `IC_turbslab.h5`, `IC_refine.h5`, `ibm_coeff_blocks.h5`)
are symlinked from the `mobydiff.scalar` / `mobydiff.bl` checkouts and are
not committed — `run_bitexact*.sh` fails their cases with an HDF5 "no such
file" otherwise; (iv) step 6 (one BC mechanism) now starts from a projection
loop that calls `apply_bc` once, which is the shape it wanted; (v) the
`tools/h5maxdiff` landmine still holds (no dataset args = no passive scalars).

Every measured number is in `docs/numerics_review_2026-09-26.md` section 10
under its step (the MEASURED blocks) and section 9's table carries the
per-finding verdicts. The original prompt follows, for the record.

---

Read `CLAUDE.md`, then `docs/numerics_review_2026-09-26.md` sections 9 and 10
in full before touching anything. Execute steps 0 to 5 of section 10, in
order, one at a time, each closed by its gate before the next starts. Do not
start steps 6–12. Record every measured number in the review document under
the step it belongs to (a "MEASURED" line per step) and update the STATUS
header of this file when you stop.

Build both `./compile.sh cpu` and `./compile.sh gpu` first, plus the
`cpu_nofma`/`gpu_nofma` modes when a step says nofma. Run every executable
through `mpirun`. Remote GPU hosts and their setup are in CLAUDE.md "Build
and run"; pin `CUDA_VISIBLE_DEVICES` and check `nvidia-smi` on istmcetus.

**Step 0 — reference set.** From THIS commit (record `git rev-parse HEAD`),
build CPU and GPU `nofma` binaries of `moby_solve` and `moby_prepare` into
`~/numrev_ref_binaries/` with a PROVENANCE file (hash, date, flags). The
2026-08-07 lesson applies: never cut a reference from a working tree
mid-edit.

**Step 1 — F3, Chebyshev `lmax` (experiment, no code first).** Case: the
boundary-layer configuration that produced the B0 instability
(`tutorials/turbulentBoundaryLayer`, the Blasius laminar precursor
`blasius2d.ini` path described in its README) with `[pressure] accel =
chebyshev`, `niter = 6`, then `12`, `dtmax = 0.5`; the same with an
additional `cheb_lmax = 2.2`. Log `max|div|` (add a temporary print if
needed, or use the existing runtime file) and `max|pn|` per 100 steps over
~100 time units. Prediction: default e-folds at ~36 t.u.; 2.2 does not.
Control: `tutorials/min_channel/input_gpu.ini` with both bounds, 200 steps,
divergence residual only. If confirmed: `pressure_solver.f90:185` auto bound
`2.0 -> 2.2` (1.1 x Gershgorin), one sentence in
`docs/numerical-methods.md`, and re-gate every Chebyshev case in
`validation/README.md` at its recorded physics band (the fields move by
construction, bit-exactness is not the gate). If NOT confirmed, record the
histories and stop the step; the review's causal claim is then retracted in
section 9.

**Step 2 — F2, conjugate convective mask.** (a) Test: the cheapest curved
conjugate case in `validation/conjugate/` (see its README; an oblique or
cylinder-type ini, not the pipe), initial scalar `θ ≡ 1` in fluid and solid,
flow on, `ibm_wall = conjugate`; run `[scalar] convection = divergence` and
`= advective`, 1000 steps, and report `max|θ − 1|` over FLUID cut cells
(cells with `phi > 0` and a neighbour with `phi < 0`; `tools/h5maxdiff` or a
short h5py script). Prediction: O(1e-2) drift in divergence, round-off in
advective. (b) If confirmed: in `validate_conjugate_config`
(`scalar.f90:~1351`) resolve an UNSET `[scalar] convection` to `advective`
when any scalar is conjugate, and make an explicit `divergence` with a
conjugate scalar an `error stop` naming review finding F2. (c) Gates: the
9-case scalar suite bit-exact vs step 0 (nofma, nothing non-conjugate may
move); `validation/conjugate/run_gates_c{1,2,3}.sh` — expect exactly ONE
gate to move, the `Σ C θ dV` drift with the flow on, which becomes an O(h)
cut-cell leak; re-baseline it in the README with the reason. (d) Re-run the
pipe Nusselt comparison (`tutorials/cht/pipe`, remote GPU) and update its
asset report. Steps (a)–(c) are one session; (d) is a background run.

**Step 3 — F4, Neumann normal velocity.** (a)
`tutorials/sailplane/input.ini`: replace the `x_max_u/v/w_type = neumann`
rows by `x_max_patch = outlet` (the pressure there becomes Dirichlet;
update the README). (b) `boundary.f90 resolve_face_bcs`/validation: a
`BC_NEUMANN` on the NORMAL component of a non-periodic face is a config
error whose message points at the outlet patch. (c) Gates: the sailplane
leg of `validation/prepare/run_gates_big.sh` is re-baselined (the BC changed
on purpose); the sailplane runs 200 steps on the DEFAULT Jacobi/Chebyshev
solver without the `solver = redblack` row (record whether it now can);
`validation/freestream/run_gates.sh` unchanged. Half a session.

**Step 4 — F5, `apply_bc` out of the projection loop.** Depends on step 3.
`pressure_projection` (`pressure_solver.f90:297-343`): remove the
per-iteration `apply_bc`, call it once after the last
`jacobi_apply`/`interface_correct` and BEFORE the final full exchange; same
in `redblack_projection` (once after the last colour). Gate: bit-exact
`max_abs 0` at PRODUCTION flags vs step 0's production build (the removed
writes were idempotent, so anything non-zero means a face kind was missed —
stop and find it) on the 7-case suite, `validation/freestream` and
`validation/redblack_interface` (1 and 4 ranks). Measure `proj_timing: bc`
before/after with `[output] profile = true`. Half a session.

**Step 5 — F8 and F6.** (a) `docs/configuration.md`: `cflmax` is compared
against the largest single-component Courant number; add an init-time print
of the worst-case SUM over directions from the same reduction in
`get_timestep_rates`. (b) In the sibling checkout `mobydiff.scalar` (or
`mobydiff.bl`) find the KMM180 run logs and read the `cfl`/`dt` the run
actually used; if the Péclet limiter cut `dt` to ~8e-6, set
`tutorials/channel_kmm180/input.ini` `natural_dyw_plus` to a wall-resolved
value (0.5) and the caps accordingly, verify with a 100-step print, and
record the decision in the tutorial. An hour each.

Stop conditions: a gate that fails and is not understood; a numerics
prediction (steps 1, 2) that is contradicted — record, do not force it.
Commit after each step with the step number in the subject. At the end,
write the MEASURED lines into the review, update this file's STATUS, and
list what steps 6–12 inherit.
