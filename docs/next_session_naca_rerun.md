# Next session: re-run tutorials/naca/rans from scratch, stopping on `steady_tol`

## STATUS: DONE (2026-08-08)

The converged state is regenerated (`tutorials/naca/rans/c11_nose_660000.h5`,
t = 34.0, plus the level-11 baseline `c11_aoa5_460000.h5` at t = 29) and the
README table is refreshed with the current tooling. Cost: ~33 h stage 1 on
the RTX 3060 + ~46 h stages 2-3 on one istmcetus A6000 (0.46 s/step vs the
3060's 1.16 — worth moving hosts for).

Results: runtime CV budget C_L 0.5199 +- 0.0001, C_D 0.01294 +- 0.00017
(OF 0.5142 / 0.0134); Cp_min -1.7797, EXACTLY OpenFOAM's, where the
destroyed state read -1.7686. The regeneration reproduces that state to a
Cf median ratio of 0.999 (5-95 % 0.995-1.003) and a Cp rms of 0.0024, so
the whole Cf/OpenFOAM discussion in the README carried over unchanged.

Four things this doc got wrong or could not know, all now in the README's
"Convergence and stopping" section and in `run_case.sh`:

1. **The tolerance could not fire at the inis' 20-step sampling.** The CV
   momentum carries a coherent f ~ 35 U/c oscillation (lag-1 autocorrelation
   +0.94) that cancels in C_L/C_D but not in a one-interval difference, so
   the measure floors at ~2A/dt_sample = 0.19 in the transient and ~8e-3
   when converged. The fix is `force_sample_interval`, not the threshold:
   sampling ~0.25 t.u. apart drops the floor two decades. The doc's
   "steady_tol = 1e-3 reproduces the shipped convergence level" read the
   README's "+- 0.005 is the CV d/dt term" as something the 20-step measure
   could see; that +- IS this wobble at 0.5 t.u. differencing.
2. **Stages 1-2 are NOT the place to save time.** "5e-3 is plenty; they are
   re-converged after each interpolation anyway" is wrong for stage 2: C_L
   relaxes with tau ~ 9 t.u. after the L10 -> L11 interpolation, so its
   18 t.u. is ~2 time constants. Cutting it to 6 was tried and, on the same
   nose layout run in parallel, under-reads C_L by 2 % (0.5105 and still
   drifting, vs 0.5199).
3. **A steady certification is not a converged suction peak.** Stage 3
   certified at 8.5e-4 after 2 t.u.; the measure is a momentum budget the
   converged circulation already satisfies, while Cp needs ~3 t.u. Run the
   published 5 t.u. out and keep the tolerance as a floor.
4. **The final snapshot was dt-clipped anyway** (step 660002 after 2
   micro-steps onto t_final), so the README's regular-cadence rule still
   applies — 660000 is the state to post-process. It would NOT apply had
   the stage stopped on the tolerance, as the doc says.

Solver change made for this (committed with the campaign): the forces trace
gained a `dmomdt` column carrying |2 dmom/dt|/qref, the quantity
`steady_tol` thresholds — without it the tolerance cannot be sized from a
real trace, since the reported C_L/C_D already have -dmdt folded in.

## The prompt

> Reproduce the converged state of `tutorials/naca/rans` from scratch
> (`./run_case.sh scratch`), and drive each stage to convergence with the new
> `[case.airfoil] steady_tol` instead of its hand-tuned `t_final`. The shipped
> converged state was destroyed, so `restart` is not available and the staged
> L10 -> L11 -> nose protocol has to run in full; it is a multi-day GPU job, so
> plan it as a detached, monitored pipeline rather than a foreground run.
> Read `docs/next_session_naca_rerun.md` first — it lists what survives on
> disk, where to put the tolerance, how to size it, and what to verify at the
> end. Deliverables: the regenerated converged state, and the tutorial's
> published C_L/C_D/Cp/Cf table refreshed with the CURRENT tooling (the
> runtime control-volume budget and the fixed `cv_forces.py`), with the
> README's provenance note updated.

## Why this run is worth the GPU time

Two independent reasons, both from the 2026-08-03 session (commits `c0cd1bb`,
`57f5665`):

1. **The published numbers are stale.** The table in the tutorial README
   (C_L 0.520 / C_D 0.0128 production; 0.514 / 0.0130 baseline) came from
   `postProcess/cv_forces.py` BEFORE its half-cell tangential-velocity
   collocation fix, so it shifts by ~0.2 %. The converged state that produced
   them no longer exists, so there is no way to recompute without re-running.
2. **`steady_tol` needs a real production case.** It has only been exercised
   on the Re 40 cylinder (stops itself at t = 105.4) and negatively on the
   shedding Re 100 cylinder. This case is the intended use: a steady RANS run
   whose `t_final` values (12 / 30 / 35) are hand-tuned guesses from a
   previous campaign. Replacing them with a measured criterion is the point.

## What is on disk, and what is not

- `assets/geometry/ibm_coeff_c11_nose.h5` (647 MB) SURVIVES — stage 3's
  `moby_prepare` (~40 min at 20 ranks) is already paid for.
- `assets/geometry/ibm_coeff_c11.h5` (L11) and `.ibm_coeff_l10.h5` are ABSENT;
  `run_case.sh scratch` regenerates both, each costing a prepare.
- NO state files (`c11_nose_640000.h5`, `c11_aoa5_*.h5`, the L11 baseline) —
  all destroyed. `./run_case.sh restart` will refuse.
- `assets/geometry/n0012_b11.stl` and all OpenFOAM reference data survive, so
  the comparison target is intact.

## Cost, and why the tolerance changes it

RTX 3060 numbers from the README: the nose case runs 1.09 s/step at 8.69M
cells (~12.7 h per time unit); stage 3 alone (t = 30 -> 35) is ~42 h. Stage 2
is the expensive one as configured — L11 at dt 5e-5 from t = 12 to 30 is
~360k steps, order 100 h — and stage 1 (L10, dt 1e-4, t -> 12) adds ~15-20 h.
So the configured protocol is roughly a week of GPU time.

That is exactly what `steady_tol` is for here: stages 1 and 2 only have to be
converged enough to hand off to the next refinement level, and their `t_final`
values were chosen with a safety margin. Let the tolerance cut them short.
Keep each stage's `t_final` as a SAFETY NET so a stage that never converges
still terminates.

## Where to put `steady_tol`, and how to size it

The stages are built by `run_case.sh` with `sed` from `c11_aoa5.ini`
(stages 1 and 2) and run `c11_aoa5_nose.ini` directly (stage 3). Add the key
per stage in `run_case.sh` rather than to the committed inis, so each stage
gets its own tolerance — and so the committed tutorial inis keep reproducing
the documented `t_final` behaviour for anyone who wants it.

**Sizing the tolerance.** The criterion is the budget's unsteady term in
coefficient units, `|2 dmom/dt| / (U_inf^2 c L_span)`, larger component. The
README already tells you the achieved steadiness of the shipped state: its
quoted C_L +- 0.005 / C_D +- 0.0008 "is the CV d/dt term, not sampling
error". So at t = 35 this case was still moving at order 1e-3 in coefficient
units. Start from that:

- stage 3 (production): `steady_tol = 1e-3` reproduces the shipped
  convergence level. Do NOT reach for 1e-4 blind — it may cost days more, or
  never trigger. Watch the trace first (see below) and tighten only if it is
  clearly still falling.
- stages 1-2 (handoffs): `5e-3` is plenty; they are re-converged after each
  interpolation anyway.

**Sizing the WINDOW — this is the part that needs care.** The ini samples
forces every 20 steps, i.e. every 5e-4 time units at dt 2.5e-5. The default
`steady_samples = 3` therefore spans 1.5e-3 t.u. — about a minute of wall
clock against a ~20 t.u. convergence. That is far too short to mean
"converged"; a transient plateau would trip it. Set `steady_samples` so the
window covers a physically meaningful stretch:

    steady_samples = 2000     # 2000 x 20 steps x 2.5e-5 = 1.0 time unit

Keep `force_sample_interval = 20` (the dense force trace is useful, and a
denser difference only helps the noise floor, which sits several decades below
1e-3). Scale `steady_samples` per stage for its own dt: stage 1 at dt 1e-4
needs 500 samples for the same 1 t.u., stage 2 at dt 5e-5 needs 1000.

**Before committing to a tolerance, measure.** Run stage 1 with `steady_tol`
unset for a few thousand steps, plot the measure from the forces trace
(`|2 dmom/dt|/qref` is reconstructible as the sample-to-sample change of the
reported coefficients' unsteady part, or just instrument it temporarily), and
pick the tolerance from where the curve flattens. Guessing costs days here.

## Two properties of the criterion that matter for this case

- **It is pressure-free.** The measure is built from the box momentum, i.e.
  velocities only, so it is immune to the stored-`pn` drift that afflicts the
  reported C_L/C_D on long runs. The stopping decision is trustworthy even if
  the coefficients printed alongside it are not.
- **It avoids the dt-clipped final snapshot.** The README's standing warning —
  "ALWAYS post-process a REGULAR-CADENCE snapshot", because a run ending on
  `t_final` writes its last field on a dt-clipped micro-step that inflates the
  stored incremental `pn` — does not apply when the run stops on `steady_tol`.
  It stops on an ordinary step, so the final write is directly usable. Verify
  this rather than assuming it (check the last two `dt` values in the log).

## What to verify at the end

1. **The physics still matches OpenFOAM**: Cp/Cf overlay via
   `./postprocess.sh`, and the C_L/C_D table. Expect the ~0.2 % shift from the
   `cv_forces.py` fix; anything larger is a real change and needs explaining.
2. **Runtime budget vs offline tool**: the solver's final sampled C_L/C_D
   against `cv_forces.py --boxes 1.5 --aoa 5` on the final snapshot. They
   should agree to ~0.01 % (the tool's residual is its block-edge one-sided
   gradients). A larger gap points at `pn` drift in the runtime series — this
   case runs `niter = 18` specifically to keep the stored pressure clean
   enough for the border integrals, and this is the check that it worked.
   If it did not, use the clean-p protocol from `validation/cylinder/README.md`.
3. **Box independence**: `cv_forces.py --boxes 1.5 2.5` should bracket the
   runtime value; the spread IS the error bar the README quotes.
4. Update the README table, its provenance note (currently says the numbers
   predate the fix), and `CLAUDE.md`'s naca entry.

## Landmines (all previously paid for)

- `module load toolkits/nvhpc/25.9` does NOT work on this machine; set PATH
  explicitly to `/opt/Nvidia/nvhpc/Linux_x86_64/25.9/{compilers,comm_libs/12.9/hpcx/latest/ompi}/bin`.
  A detached job that misses this dies at MPI_INIT.
- `moby_prepare` is MPI-parallel ONLY (no OpenMP in the CPU build): give it
  `-n 20 --bind-to none`, not threads. ~40 min for the nose case.
- Run detached (`setsid nohup`), and kill by PID (`pgrep -x`) — `pkill -f
  <pattern>` matches the invoking shell's own command line.
- Do NOT run `tools/compare_fields.py` on these snapshots: it reassembles onto
  the finest lattice (69 GB at this refinement) and gets OOM-killed, taking
  the session with it. Use `validation/prepare/compare_snapshots.py`.
- Guard the deletion of any `c11_*` file. The previous converged state was
  lost to a careless `rm -f c11_nose_6400[0-9]*.h5`.
