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

## 4. Investigation write-up (2026-10-06, main `0c9c8ba`; reference set `~/rough_ref_binaries`)

Reference set cut first (four builds, banners `0c9c8ba`, provenance inside).
Everything below is measured on the CPU with `build_cpu` variants; the
turbulent cut-down A/B is on the local RTX 3060.

### The history of the mode, from the dumps

`~/outlet_runs/rough_cans/analysis/{onset,history}.py` read the dumps of the
cut-down (`dump_100 .. dump_2700`, block 34 = the last column, rows 1..6,
k = 11/12 = the two cells straddling the spanwise trough). Coefficients at
that column (the saddle of the surface sits ON the plane: in x the surface
descends from 0.065 two cells upstream to 0 at the plane, in z it rises on
both sides, so the last column is a dimple open to the outlet):

| face | row 1 | row 2 | row 3 | row 4 | row 5 |
|---|---|---|---|---|---|
| u, outlet plane (`nb+1`) | 139 | 0.074 | 0.013 | 0 | 0 |
| u, one cell behind (`nb`) | solid | solid | 13.7 | 0.094 | 0.029 |
| u, two cells behind (`nb-1`) | solid | solid | solid | solid | solid |

The row-3 outlet face is the one that goes: its neighbour `q(nb)` is held
near zero by its coefficient (13.7: `u(nb)` reads 4e-8 .. 8e-4 through the
run) and the face behind that is solid.

| t | 2 | 6 | 12 | 20 | 26 | 33 | 38 | 44 | 48 | 50 | 52 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `u_out` row 3 | +0.07 | **+0.29** | +0.25 | +0.15 | +0.07 | 0 | −0.04 | −0.11 | −0.22 | −0.42 | −7.4 |
| `p(nb)` row 3 | −2.0e-3 | −2.7e-3 | −2.5e-3 | −1.7e-3 | −1.0e-3 | 0 | +2.7e-4 | +1.2e-3 | +2.7e-3 | +5.8e-3 | 0.43 |

The physical velocity at that face (y ≈ 0.03 in a laminar layer at
x = 100) is ~0.02: the face is 10x off from the first dump, an OUTFLOW that
follows the pocket pressure (`u_out ∝ −p(nb)`: outflow while the pocket is
below the outlet pressure, inflow once above it), crosses zero at t ≈ 33,
then grows as INFLOW at an accelerating rate (per 2 t.u.: x1.5 at t = 40,
x1.9 at t = 50) and blows up in finite time. The rows above (4..6) read the
physical 0.01–0.03 throughout. The fluid entering at the plane leaves
through the ceiling of the cell: `v(nb, 4) dx = −u_out dy` to three digits
— the ceiling is 19x the plane face, so the escape velocity, and the
pressure that drives it, stay small.

### A reproducer small enough for the CPU

`~/outlet_runs/pocket/pocket.ini`: the same wall line (`cans`, dy_1 0.0105),
the same dx (0.205), one spanwise wave, FIVE streamwise waves on
lx = 21.74, roughness in over x = 5..10 and through the outlet, NO trip,
Blasius cold start, red-black 6, `cflmax 1.0`. 450k cells, 4 CPU ranks,
0.13 s/step, 2500 steps = 50 t.u. in 5 minutes. It reproduces the whole
history with the same coefficient table (row 3: 0.012 on the plane, 12.8
behind, solid behind that): `u_out` row 3 +0.35 at t = 6, through zero at
t ≈ 25, exponential inflow, blow-up at **t = 34** (`pocket_history.py`).
Laminar: the instability needs neither turbulence nor the trip.

### Stage by stage: what the row-3 face is given

The face obeys `q_f(new) = sf q_f + mu_f (dt_a rhs_n + dt_b oldrhs_n)
− mu_f dt_g dp_f` with `rhs_n` the NEIGHBOUR'S non-pressure right-hand side
(the zero-normal-gradient-of-R closure of 2026-10-01). With `q_n ≈ 0`,
`q_nn = 0` (solid) and the tangential velocities at the column ≈ 0, the
neighbour's stencil reduces to

    rhs_n = −(q_n + q_f)²/(4 dx) + nu (q_f − 2 q_n + q_nn)/dx² + ...
          ≈ −q_f²/(4 dx)        + nu q_f/dx²

Two things handed to the face with the wrong sign, because they are the
neighbour's interaction WITH the face: (a) `+nu q_f/dx²` is anti-diffusion
on the face (rate 0.054/t.u. at this dx and Re); (b) `−q_f²/(4 dx)` is a
sign-definite quadratic — a drag for outflow, a SOURCE for inflow: the
momentum entering through the plane is dumped at the wall half a cell
behind the face's own centre (`u_c(nb) = q_f/2`) and nothing in the face's
equation balances it, because the pressure that would (a stagnation
pressure in the pocket) is set by the projection from continuity alone, and
continuity is satisfied by the tiny escape velocity through the 19x larger
ceiling. The energy of the face gains `|q_f|³/(4 dx)` on inflow: finite-time
blow-up, as measured. This is the open-boundary backflow problem of the
literature (a fixed-pressure plane with inflow has no energy bound;
Heywood, Rannacher & Turek 1996; Bruneau & Fabrie 1996; Bazilevs et al.
2009), reached here through a geometry — a dimple on the plane — that
makes inflow at an outlet face the steady state rather than a transient.

### Closures tried on the reproducer (all in `predict_outlet_faces`, 4 CPU ranks)

| closure | result |
|---|---|
| as on main: `R_f = R_n` | blow-up at t = 34 |
| C1: the face's own normal fluxes with zero gradient BEYOND the plane (ghost `q_g = q_f`): convection `−(4 q_f² − (q_n+q_f)²)/(4 dx)`, diffusion `−nu (q_f − q_n)/dx²`; tangential terms and forces still the neighbour's; own RK memory | blow-up at **t = 31**: fixes the viscous sign, triples the quadratic (`−3 q_f²/(4 dx)`); the ram term is the driver, not the anti-diffusion |
| Z: symmetric ghost `q_g = q_n` (zero gradient AT the face): normal convective flux difference ≡ 0, diffusion `2 nu (q_n − q_f)/dx²` | blow-up at **t = 10**, the other way: with no ram drag the OUTFLOW runs away under a pocket pressure below the outlet value (1.06 on the plane at t = 8) — the old pressure mode's shape. The borrowed quadratic is load-bearing for outflow |
| BV: the neighbour's convection kept; the face's own viscous flux (zero viscous flux through the plane, `−nu (q_f − q_n)/dx²`); and on INFLOW (`s q_f < 0`) the outlet datum is the TOTAL pressure `p_face = p_out − q_f²/2` | **stable** to t = 150 (`bvlong`): the pocket inflow saturates at −0.083 (t = 38), −0.066 at t = 50, −0.032 at t = 150; the early +0.3 outflow is unchanged |
| **B: the datum alone (the shipped scheme)** | **stable**: −0.094 at t = 38, −0.077 at t = 50, the same history otherwise; the turbulent cut-down 300 t.u. (below) |

The total-pressure datum is the Bernoulli statement for fluid entering from
a reservoir at `p_out`: the kinetic energy it brings is paid by the
pressure at the plane, not fed into the face; the boundary energy flux
becomes `p_out |u_n|` and the cubic source is gone. In the face equation it
is `dp_f −= s q_f² d1f` on inflow only (the half-cell gradient of a datum
lowered by `q_f²/2`): with the borrowed `−q_f²/(4 dx)` the inflow sees a
net RESTORING `+3 q_f²/(4 dx)`, outflow keeps its drag. Low and high faces
are one statement through the outward sign `s`.

### The scheme, as shipped

ONE amendment to the 2026-10-01 face equation, in `predict_outlet_faces`
(step.f90), symmetric in the side through the outward sign `s`:

    s  = -1 on a low face, +1 on a high one
    if (s q_f < 0) dp_f -= s q_f^2 d1f

i.e. the held datum is the STATIC pressure `p_out` where fluid leaves and
the TOTAL pressure `p_out − q_f²/2` where it enters (the ghost moves by
`−q_f²`, the half-cell gradient by `−s q_f² d1f`). The closure of R
(zero normal gradient of the neighbour's non-pressure right-hand side)
stays as ratified. Three lines; dormant on every outlet without backflow,
so the body-free suites and every outlet case without inflow are
bit-exact against `~/rough_ref_binaries` (see the gates below).

Tried and dropped: a second piece that replaced the borrowed normal
viscous term by the face's own `−nu (q_f − q_n)/h²` (zero viscous flux
through the plane) with its own RK memory. It removes the anti-diffusion
(a), but (a) is bounded by the datum anyway (pocket inflow −0.094 without
it, −0.083 with it), in smooth outflow the borrowed second difference is
the better approximation of `nu u_xx` (the do-nothing flux is a spurious
`nu u_x/h`; Lamb-Oseen moved 9.80e-5 → 9.84e-5), and it moves every outlet
case whose flow varies along the normal at all (the Poiseuille restart
field, 1e-6 after 50 steps), i.e. it costs bit-exactness on the whole
outlet suite for no measured gain. The clean statement is: the closure is
one thing, the boundary condition is another, and only the latter was
wrong.

### Gates (2026-10-06, the datum-only tree; reference `~/rough_ref_binaries`, `0c9c8ba`)

Bit-exactness at nofma, `max_abs 0` on every dataset incl. the RANS
scalars and the `un_xmax` plane:

| suite | CPU | GPU |
|---|---|---|
| 7-case (`validation/scalar/run_bitexact.sh`) | 7/7, 4 ranks | 7/7 |
| 9-case (`run_bitexact_s3.sh`) | 9/9 | 9/9 |
| outlet (`validation/freestream/run_bitexact_outlet.sh`: oblique, pois_io, lamboseen, outlet_high, outlet_low, blasius2d, cyl_re100) | 7/7 at 1 rank AND at 4 ranks | 7/7 |

(The datum is dormant on every one of them: no outlet face sees backflow,
the cylinder wake at 10 D included.) Physics gates, all PASS:

- `validation/freestream/run_gates.sh`, all groups (oblique exact 0.0;
  Poiseuille 1.637e-3 / 1.521e-3, last-cell p 1.247e-3, drift 5.6e-16;
  Lamb-Oseen reflected fraction 9.84e-5; ranks exact; config; the mirrored
  low-face pair 2.0e-15; the −x Poiseuille and vortex to every digit;
  refined-from-rest 1.5–2.7e-13, 1 == 4 ranks exact; restart == continuous
  on five cases) — the numbers of 2026-10-01 reproduced.
- `validation/body_outlet/run_gates.sh` profile, mirror, permute, limits,
  ranks, restart, twin (CPU, 4 ranks, SHORT) and gpu (nofma pair): ALL
  PASS; the new **pocket** group PASS (dimple face −0.077 at t = 50, rms
  divergence 1.0e-4; the reference binary blows up at t = 34).
- The pocket reproducer CPU == GPU `max_abs 0` (nofma pair, same case
  file, 1500 steps = 30 t.u., the datum live from t ≈ 25): the one gate on
  which the new line is NOT dormant.
- `validation/blasius/run_blasius.sh` (CPU): theta 1.65 %, H 0.82 %,
  du/Ue 1.48e-3, dv 1.29e-2 — the record to every digit.
- The 4 M-cell turbulent cut-down (`~/outlet_runs/rough_cans/bv/
  cut_fin300.ini`, RTX 3060, 0.089 s/step) with the roughness THROUGH the
  outlet: 300 t.u., L2 divergence 1.9e-4 → 5.7e-6 (the `wall_x_end = 60`
  twin: 6.0e-6), Linf 1.10, the dimple face between −0.08 and +0.13 once
  the layer is turbulent (max|p| 0.16–0.18 at the plane from t ≈ 200
  against 0.04 before; the reference blew up at t = 52).

Reference set for anything on top of this: cut one from the commit
(`~/rough_ref_binaries` is the pre-fix set; the pocket case and any outlet
with backflow differ from it by design).
