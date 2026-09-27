# Next session: numerics review, steps 0–5

STATUS: NOT STARTED (written 2026-09-27 at `ac7b925`). This is the prompt for
the session that executes the correctness tier of
`docs/numerics_review_2026-09-26.md` section 10. Copy the block below into
the session verbatim.

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
