# Next session: read the HoreKa pair, then a second cleaning and profiling pass

STATUS: handout written 2026-10-06 (head `b4861e1` + this file), after the
session that fixed the body-on-outlet blow-up (`docs/next_session_rough_outlet.md`
section 4) and queued the final rough case on HoreKa. Nothing below is started.
The user's order (2026-10-06): **(1) read the HoreKa results; (2) a SECOND
cleaning and profiling pass; (3) THEN the section-6 simplifications, the
performance items and the `config.f90` `case default`; (4) the review's
feature tail stays for the future.**

## 1. Reading the HoreKa results

Two chains on `dev_accelerated` (one job per user at a time, so they
alternate; each chunk is 12 500 steps = 250 t.u., ~40 min on 4 A100):

| run dir (`…/xt8786-roughBoundaryLayer/optimiseBlockRefinement/`) | sides | code | start | end step |
|---|---|---|---|---|
| `rough_tbl_run` | `smooth` (production_stats.ini) and `rough` (roughness OUT at x = 580, the fallback) | `moby-tbl-rough` (b27936e-era, BEFORE the fix) | smooth from `smooth_687500.h5`, rough from `start_field.h5` (step 675000) | 775000 |
| `rough_io_run` | `rough` only: roughness THROUGH the outlet, `cflmax 1.0` | `moby-tbl-rough-io` (branch `rough-outlet`, c7b2773) | `smooth_737500.h5` of the pair | 837500 |

Check first (`status.txt`, one line per chunk start/end; `provenance.txt`
lists every job of the chain): no `exit` code other than 0, the chain not
"BROKEN" (then resubmit by hand with the command it prints), and on the
`rough_io_run` side that `partition:` and `boundary: … zeroed` lines in the
chunk logs read the same as the first chunk's. The dimple faces are the
thing the fix was for: a chunk's `rough_runtime.txt` must show L2 div of
order 1e-5 and Linf ~1.1 throughout (the 2026-10-05 failure read Linf 2 at
+1000 steps and 2e29 at +2000).

Copy to the workstation (no h5py on HoreKa; 200 MB each):
`stats_smooth_775000.h5`, `stats_smooth_712500.h5` from `rough_tbl_run`;
`stats_rough_837500.h5`, `stats_rough_775000.h5` from `rough_io_run`
(the fallback's `stats_rough_775000.h5` / `_712500.h5` from `rough_tbl_run`
only if the through-outlet run failed). Then

    cd tutorials/turbulentBoundaryLayer
    python3 assets/postpro/rough_compare.py \
        --smooth stats_smooth_775000.h5 stats_smooth_712500.h5 \
        --rough  stats_rough_837500.h5  stats_rough_775000.h5

Windows are differences of two cumulative files: smooth = the last 1250
t.u. of the pair; rough = steps 775000..837500 (its first flow-through after
the restart, ~650 t.u. = 32 500 steps, is discarded that way). Deliverables
(handout section 2 of `next_session_rough_outlet.md`, unchanged):

- the von Karman balance per x band on BOTH sides (`momentum_integral.py`
  conventions; c_f on the roughness = the balance's d theta/dx, validated to
  ~1 % on smooth data): smooth within ~5 % to L_x − 100, within 20 % to two
  units from the outlet; rough holding the balance THROUGH the rough region
  and — new with the through-outlet run — to the outlet: the x = 580..650
  bands are now rough, there is no rough-to-smooth relaxation to report;
- k+ along the plate (expected 9..11 on the rough u_tau) and ΔU+ above the
  crests at x = 300..550 against 3.72 (MacDonald et al. 2016; report both the
  mean-plane and the crest origin, the two differ by ~0.5);
- the last column: stored p rms of the outlet column, the dimple faces'
  statistics (anything above a few 1e-2 mean inflow there is worth a look —
  the cut-down read ±0.1 instantaneous once turbulent);
- the smooth statistics replace `assets/mobydiff/xyz_3200_384_136/data.nc`
  and `code_comparison.png` via `reproduce.py` (README paragraph "The table
  is the 10 000-time-unit statistics of 2026-10-02" changes; say the window
  is 1250 t.u.);
- `tests_record.md` section 14 gets the numbers; the CLAUDE.md bullet of
  2026-10-05 ("the HoreKa pair is RUNNING") gets replaced; the handout's
  STATUS updated; then push (commits `a33cdac c7b2773 da46a99 3927bc8
  b4861e1` + this file are local).

If the through-outlet run blew up: the chunk logs name the step; dump that
chunk's last snapshot locally (`pocket_history.py` style, block of the
trough column) BEFORE touching the scheme — the fix was gated on a laminar
and a 4 M-cell turbulent case, the production case has Re_theta ~ 700 and
k+ ~ 10 at the plane, so a new mode there is a new finding, not a retry.

## 2. The second cleaning and profiling pass (BEFORE section 3)

Same shape as 2026-09-26: a thorough read of the code as it is NOW (the
outlet work, step 7, F1/F7, the backflow datum, the trip port, the rough
wall keys all landed since the review), then a fresh profile of the
production step. Inputs:

- `docs/numerics_review_2026-09-26.md` is the previous pass; its sections
  1–5 and 9–10 are closed, section 6 is the open list (see 3 below). Write
  the new pass as its own dated document, do not append.
- Profile on HoreKa with the campaign matrix
  (`overheadTest/horeka/exchange/run_matrix.sh`, results
  `results_horeka_2026-09-25.md`/`_2026-10-01.md` as the last baselines),
  at 4/8/16 ranks, PLUS the rough case (`rough_jacobi.ini` runs again since
  2026-10-05) and `min_channel`; ncu per kernel where a bucket moved
  (`run_ncu.sh`; read `launch__registers_per_thread` and
  `sm__warps_active`, never a step time alone). Every published ratio has a
  moved denominator since 2026-09-14 — re-measure, do not quote.
- Candidates for the cleaning half: dead keys/branches left by the outlet
  work (`BC_OUTFLOW` as a row TYPE with no row; the `outflow_copy`
  references in comments), the `interface_correct` outlet arguments that are
  gone, duplicated local copies of constants, drivers in `validation/` that
  are not gates any more (`prepare/run_gates_step7.sh`), the C3 transient
  CONTROL rows (an August binary on today's inputs — rebuild the control or
  retire the rows), the `docs/next_session_*.md` set (34 files: fold the
  DONE ones into a short index or move them under `docs/history/`).
- One rule from this session for the pass: a gate's bound is taken from the
  state the gate READS (the pocket divergence bound was set from t = 150
  and failed at t = 50).

## 3. After the pass: simplifications, performance, `case default`

In this order unless the profile says otherwise.

1. **`config.f90` `case default`**: an unknown key in a known section is
   silently ignored — it has produced a wrong-but-successful run at least
   three times (2026-09-25 campaign, `stats_start_step`, the August control
   rows). One hour: error-stop with the section and key, then run every ini
   in `tutorials/` and `validation/` through `moby_prepare` to catch the
   keys that only existed by accident.
2. **Section 6 of the review, still open:** (5) the three
   `interface_correct` kernels → one launch over a face list (~1 ms/step at
   16 ranks); (6) split `scalar_transport_kernel` by mode (conjugate vs
   plain; the dormant branches cost occupancy, the register lesson of
   `jacobi_apply`), bit-exact if expressions move verbatim; (7) the
   per-substage `qs → q` copy, 6 doubles/cell — measure before deciding
   (double-buffering fights the flat-array mapping); (8) fuse
   `cheb_combine` into `jacobi_compute_phi`.
3. **Performance items left by the campaign:** the fused predictor split
   into per-component kernels (the register lever is exhausted on it; ncu
   before and after); the fp64 divide in `compute_rdenom`, now only inside
   body blocks; a reduced divergence round for red-black (needs its own
   read-set proof: a cell at i = 0 reads the low halo of all three
   components); the unexplained `phi_exchange` locality shift of the
   divergence-list reordering; plan entry 7-5 (per-leaf weight column +
   weighted Morton split) only if a case is load-imbalanced.
4. Every item: nofma suites + outlet suite + `body_outlet` (incl. `pocket`)
   at `max_abs 0` against `~/backflow_ref_binaries` (c7b2773), production
   flags only when the arithmetic source is byte-identical.

## 4. For the future (not this session, not the next): the review's feature tail

- Step 10, implicit solid conduction for conjugate scalars (Chebyshev-Jacobi
  Helmholtz on solid + cut cells, Crank–Nicolson per substage): only if the
  CHT campaigns are rerun at niter 12.
- Step 11, line-implicit wall-normal momentum diffusion: only if a
  viscous-limited production case appears.
- Step 12, second halo layer → TVD/van Leer for the RANS scalars (the
  SD7003 front) → optional bounded passive-scalar convection and a
  4th-order central experiment.
- Accepted residuals, reopen only with a case that needs them: the
  first-order restriction ghost at a level jump on a physical wall (6.8 %
  coarse-side wall-cell nut), the coarse-owns-face v' excess at 2:1
  interfaces (the const-1/2 price), C4 corners (needs an `s_t` estimate C2
  showed does not exist), question 9 of the body-at-outlet work (a level
  jump on an outlet next to a body, never measured).
- The one-by-one re-measurement pass of `validation/README.md` (the long
  tier still carries niter-6 / pre-skew numbers) resumes whenever a session
  has a GPU and no plan step.
