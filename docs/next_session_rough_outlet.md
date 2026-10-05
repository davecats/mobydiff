# Next session: the body on the outlet plane (open defect) and the rough-wall pair

STATUS: handout written 2026-10-05 evening, after the session that merged the
`turbulentBoundaryLayer` branch (`bf95993`) and launched the smooth/rough
pair on HoreKa. Nothing here is started.

## Where things stand

- `main` is pushed through `ce9eed4`. Reference set for body-free and
  outlet-free gating is still `~/body_ref_binaries` (`d5cbda8`); the merge,
  the restart fix and the wall keys are `max_abs 0` against it on the 7-case,
  9-case and outlet suites (CPU and GPU, nofma), and the wavy P0 geometry is
  bit-identical. **Cut a new reference set from `ce9eed4` before touching
  the outlet again** (`validation/scalar/compile_nofma.sh` or the four
  `compile.sh` modes, ONE PER CALL, banners checked).
- The HoreKa pair: `optimiseBlockRefinement/rough_tbl_run`, self-chaining
  12 500-step chunks (job 5183680 and its successors; the driver re-submits
  itself with a retry loop). Sides `smooth` (production_stats.ini, the
  predicted outlet, from `smooth_687500.h5`) and `rough`
  (`production_stats_rough.ini`, roughness x = 120..160 in, 580..620 out,
  from the developed smooth field at step 675000). END_STEP 775000 = 2000
  t.u. per side. Chunk-end copies `stats_<side>_<step>.h5` make any window a
  difference of two files (`rough_compare.py` takes two per side).
- The rough case was meant to run its roughness THROUGH the outlet. It
  cannot: see the defect below. It ends at x = 580 instead; the usable
  region (x <= 550) is unaffected, the rough-to-smooth step at ~600 is not.

## 1. The defect: an immersed body on the outlet plane is unstable in a turbulent trough-cut geometry

What was measured (2026-10-05, `tests_record.md` section 14):

- Full case (3200 x 384 x 136, roughness through the outlet, restarted from
  the developed field): L2 div 1.3e-4 at step +1000 with Linf 2.02 (smooth
  side 1.23), 2e29 at +2000.
- Cut-down `~/outlet_runs/rough_cans/cut_skew.ini` (492 x 384 x 22 on
  100 x 100 x 4.333, one spanwise wave, cold Blasius start, trip at 10,
  roughness in over 30..60 and through the outlet at 100): blows up between
  t = 48 and 52, Linf 1.06 -> 3.7 -> 3e16 in 200 steps. `cut_divergence.ini`
  the same at t = 52: NOT the convection form.
- `cut_end.ini` (identical, `wall_x_end = 60`, the surface below the wall by
  x ~ 80): 300 t.u. stable, L2 div 6e-6, Linf 1.07-1.10. The one-variable
  A/B: it is the body on the outlet plane.
- Dumps `dump_2400/2500/2600.h5` (same directory, nb = 82 48 22, 48 leaves)
  locate it: the LAST column i = 490 (of 492), rows j = 0..2 (y < 0.03,
  the trough bottoms, which are tangent to the domain wall because
  `wall_offset = amp_x`), k = 11..12 in block 34 (origin 410 0 0). There the
  OUTLET face coefficient is graded 139 / 0.074 / 0.013 (rows 1..3) while
  the face one cell behind it (i = nb) is SOLID 2.2e27 (rows 1..2) or graded
  13.7 / 0.094 (rows 3..4). The solid-cell pressure of that column goes
  0.019 -> 0.056 -> 13.3 over 200 steps while the stored outlet plane
  `un_xmax` goes 1.01 -> 7.4; the interior is quiet until then.

So the configuration is B1's "n solid, f fluid" limit (`r = mu_f/mu_n ~
1e27`), which `validation/body_outlet/` gated only with laminar Poiseuille
between slabs, plus graded (partially penalized) outlet faces, which no gate
had at all. The convection form, the time step (dt 0.02, CFL 0.35 on the
sum) and the inlet (zeroed data: 0) are excluded.

Where to look first (investigate before changing anything, the outlet
sessions' rule):

1. **Stage by stage for the last fluid cell whose upstream face is solid.**
   Write the table of `docs/next_session_outlet.md` ("One substage, stage
   by stage") for the cell (nb, j, k) with `q(nb)` solid and `q(nb+1)`
   fluid or graded: the predictor's `I_n` of a SOLID face (its stencil reads
   the outlet face value in the convective product and the diffusion),
   `predict_outlet_faces` with `r ~ 1e27` (what `r*qs(n)` carries once
   `qs(n)` holds only round-off of the solid DOF), the projection's
   diagonal of a cell with ONE open face (the outlet, `dnHigh*mu(f)*d1P`),
   the stored pressure `p += phi/dt_gamma` of that dead-end cell, and how
   the next predictor reads it back through `dp(f)` (the Dirichlet ghost
   `2 p_out - p(nb)`). The linearised loop face <-> dead-end pressure has
   a marginal eigenvalue on paper; the red-black projection's six
   iterations on a cell whose diagonal is one face, and the convective
   `q(f)^2` term, decide the sign. Measure the growth rate from the dumps
   (per-substage, not per-step) before theorising.
2. **The graded outlet face.** A face with `0 < coef < SOLID` on the plane
   is a partially penalized outlet DOF: it sits in the Jacobi denominator
   with its own `mu_f`, is predicted with `mu_f` and corrected with
   `cfHigh*mu_f`. Check whether the SPD pair still holds for it (the gate
   for that is the Poiseuille slab case with the slab edge placed INSIDE a
   cell on the outlet plane, not 0.3 dy off a face).
3. **A reproducer small enough for the CPU**: shorten `cut_skew.ini` to the
   last 20 units (inlet = the developed profile is not available; a
   periodic-in-x rough channel with an outlet is not the same case) — or
   keep the 10-minute A6000 run, it is cheap enough. `[output]
   field_interval = 50` from step 2400 catches the onset.
4. Candidate fixes to WEIGH, not to apply first: treating an outlet face
   whose upstream face is solid as a solid face (it is a dead-end pocket
   open to the outlet, the physical answer is "no flux"); a solid-cell
   pressure anchor at the plane; evaluating the face's own `N` from the
   neighbour's FLUID faces only. Each needs the Poiseuille slab gates, the
   outlet suites (`run_bitexact_outlet.sh`, `validation/freestream/
   run_gates.sh` low/refined/restart) and `validation/body_outlet/` at
   `max_abs 0` for every configuration that is NOT this one, and the
   cut-down A/B as the gate for this one.

## 2. The pair, when it is in

Wait for `rough_775000.h5` and `smooth_775000.h5` (eight chunks each; the
rough side lags by one chunk and the chain shares the user's one-job-at-
a-time dev queue, so count on a day or two). Then, on the workstation
(no h5py on HoreKa; copy the stats files, 200 MB each):

    python3 tutorials/turbulentBoundaryLayer/assets/postpro/rough_compare.py \
        --smooth stats_smooth_775000.h5 stats_smooth_712500.h5 \
        --rough  stats_rough_775000.h5  stats_rough_712500.h5

(window = the last 1250 t.u.; the rough side's first flow-through after
its restart ends at step ~707500). Deliverables:

- the von Karman balance per band for both sides: the smooth side must
  reproduce the predicted-face outlet zone of the earlier grid (balance
  within ~5 % to L_x - 100); the rough side must hold the balance through
  the rough region, and the x = 580..650 bands show the rough-to-smooth
  relaxation (expected: c_f drops below the smooth value for a few δ,
  Antonia & Luxton 1972 in reverse);
- k+ along the plate (expected 9..11 on the rough u_tau) and ΔU+ above the
  crests at 300..550 against 3.72 (MacDonald et al. 2016, full-span channel
  at Re_tau 180; the comparison is loose: k+ is not exactly 10, the layer
  is not a channel, and the virtual origin is worth 0.5 — report ΔU+ for the
  mean-plane origin and for a crest origin);
- the smooth statistics replace `assets/mobydiff/xyz_3200_384_136/data.nc`
  and `code_comparison.png` via `reproduce.py` (README paragraph "The table
  is the 10 000-time-unit statistics of 2026-10-02" then changes; note the
  window is 1250 t.u., not 10 000, unless the chain is extended);
- `tests_record.md` section 14 gets the numbers; the CLAUDE.md bullet of
  2026-10-05 gets its RUNNING replaced.

If a chunk fails: `status.txt` has one line per chunk end; the driver stops
the chain after two consecutive failures of a side and prints the hand
command; `attempt1_no_body/` and `attempt2_body_on_outlet/` hold the two
failed attempts of this session.

## 3. Decisions (taken by the user, 2026-10-05)

1. **The defect is FIXED, not documented away.** "There is no reason why this
   should not work." Section 1 is the next session's whole first part.
2. **Order: 1 then 2.** The rough case exists to have a rough boundary layer
   AND to validate the body on the outlet. So: fix the defect, gate it on
   the cut-down A/B and the suites, THEN rerun the rough case with the
   roughness THROUGH the outlet (`wall_x_end` removed from
   `production_stats_rough.ini`). The pair running now (roughness ending at
   580) is kept as the smooth rerun plus a rough fallback; its rough side is
   not the final one.
3. Sampling length 2000 t.u. per side: ok as is.
4. **Both sides stay CONS** for the CaNS comparison — recorded as the
   EXCEPTION to the production form (skew); no skew side.
5. **`cflmax = 1.0` whenever a case is rerun** (it was 0.8 = 46 % of the
   RK3 limit since cflmax became the sum over directions). Apply it to the
   rough rerun of item 2 and to any other ini touched from now on; do not
   sweep the inis for its own sake.
6. HoreKa housekeeping: done 2026-10-05 (the two failed-attempt directories,
   `tbl_stats_run`, `outlet_tbl_run`, the worktrees `moby-outlet-new`,
   `moby-outlet-ref`, `moby-tbl-stats`).
