# The calibration box (sailplane junction, step 2)

Section 4 of `docs/next_session_sailplane_junction.md`: the truncated model
at the production topology, run four times with SST + wall functions to
confirm the set-up, measure the sensitivities that decide whether the 0.3 deg
uncertainty of alpha_box matters, and provide the initial field of the
production run. Everything in CHORD units (root chord 0.831 m = 1), Re_c
3e5, U = 1.

## Files

- `calib_a1.5_y3.5.ini` -- the master ini (alpha 1.5 deg, the cut at 3.5 c).
  The other three are derived from it:

      for a in 1.0 2.0; do sed -e "s/a1\.5/a$a/g" -e "s/pitched 1.5 deg/pitched $a deg/" \
          calib_a1.5_y3.5.ini > calib_a${a}_y3.5.ini; done
      sed -e "s/y3\.5/y4.5/g" -e "s/^ny = 56/ny = 72/" -e "s/^ly = 3.5/ly = 4.5/" \
          calib_a1.5_y3.5.ini > calib_a1.5_y4.5.ini      # (plus two comment lines)

- `frue_a{1.0,1.5,2.0}.stl` (git-ignored, 14 MB each): the CAD STL pitched
  nose-up about the root quarter chord, (x, z) = (2.84, 0.47) m:

      cd tutorials/sailplane/junction
      for a in 1.0 1.5 2.0; do
        python3 stl_sections.py rotate "../FRUE V0 ohneRundung.stl" $a calib/frue_a$a.stl
      done

- `*.case.h5` (git-ignored): prepared on the workstation CPU build on the
  rank count of the solve, `mpirun -n 4 ../../../../build_cpu/moby_prepare
  calib_aA_yY.ini` (logs `*.prep.log`, kept beside them, git-ignored).
- `../horeka/submit_calib.sh`: one 4-A100 `dev_accelerated` job per case
  (account `hk-project-p0024333`, D5), resumable chunks; the read-out is
  `../strip_loads.py` on the snapshots.

## The strip-load call for these cases

    python3 ../strip_loads.py calib_a1.5_y3.5_NNNNN.h5 --case calib_a1.5_y3.5.case.h5 \
        --stl frue_a1.5.stl --scale 1.20337e-3 --translate 4 0 9.434 \
        --wing-y0 0.57 --wing-y1 3.5 --wing-xmax 8.5 --strips 0.1 \
        --cp-station 0.722 1.203 --out a1.5_y3.5_NNNNN

The fuselage side at the wing is y = 0.467 m = 0.562 c (`stl_sections.py`
slices at x = 2.8-3.4 m), so the strips start at 0.57 c; the stations are
y = 0.6 and 1.0 m. `--alpha` stays 0: the STL is pitched, the inflow is
axial.

## Record

(filled as the steps complete -- RE-VALIDATED convention: binary hash, ini,
every number)

### Prepare scaling (2026-10-26, build_cpu at a2e0d9b)

The leaf counts, and why they are 7x the handout's estimate:

| case | leaves | refined | case file |
|---|---|---|---|
| calib_a1.0_y3.5 | 372 070 | 366 382 | 12.1 GB |
| calib_a1.5_y3.5 | 371 586 | 365 896 | 12.1 GB |
| calib_a2.0_y3.5 | 372 453 | 366 766 | 12.1 GB |
| calib_a1.5_y4.5 | 420 115 | 412 725 | 13.7 GB |

Data point 2 (istmhydra, EPYC 7662 at 2.1 GHz, 48 ranks per case, two
cases at a time, `prep_all.sh`): wall 3125-3243 s for the 3.5 c cases,
3630 s for the 4.5 c case; max RSS 2.1-2.5 GB per rank (100-120 GB per
case over 48 ranks). The rank count is free with an explicit `[blocks]
nb`: the sailplane tutorial case prepared on 2 and on 4 ranks is
bit-identical in every dataset (checked the same day); moby_prepare is
MPI-parallel only (no host OpenMP), so use many ranks. For step 5 (3-5 M
leaves): 10x the leaves at the same per-leaf cost is 8-10 rank-hours of
Ryzen per 100 000 leaves, i.e. a HoreKa CPU node for a few hours.

Data point 1 (workstation, Ryzen 5900X at 4.95 GHz, 4 ranks):

`calib_a1.5_y3.5.ini` as first written, `refine_body_levels = 5` (c/512 on
the surface), 4 ranks: **371 586 leaves (365 896 refined), 190 M cells**
-- the handout expected 55 000 leaves / 28 M cells. The body band is the
reason: the prepare refines every block whose span dilated by 2 cells
touches the surface (`classify_block_geometry`), so at c/512 the band is
2-3 blocks thick on each side of about 16 c^2 of wetted surface, not the
one-block band the estimate assumed. Cost: about 60 min from the STL load
to "writing case file" (classification + IBM coefficients + wall
distance), 50 GB resident over the 4 ranks (62 GB machine); the run was
stopped during the case-file write by the session's memory guard, so no
file and no wall time. Consequences: the RTX 3060 (12 GB) cannot hold the
solve, and 4 A100 would step at about 0.7 s (NACA scaling): one
`accelerated` job of 5 h per case (`../horeka/submit_calib.sh`) instead
of a `dev_accelerated` hour.

### Smoke test (2026-10-26, istmhydra CPUs, build_cpu at a2e0d9b)

`smoke_hydra.sh`: calib_a1.5_y3.5 for 20 steps on 16 and on 32 ranks
(`.smoke16.ini` / `.smoke32.ini` = the case ini with nsteps 20,
field_interval 0). Both ran to "main loop ended" with the boom through
the outlet (x_max declared outlet, 656 192 Dirichlet velocity rows, 0
zeroed inside the body; 3 299 546 IBM wall cells, 18 527 775 solid cells;
`report_partition` 2 007 744 cross-rank face cells on 16 ranks). Step
time 44.3 s on 16 ranks, 24.8 s on 32 (2.1 GHz EPYC, 190 M cells), RSS
5.7 / 2.9 GB per rank. Rank identity 16 vs 32 on `un vn wn pn` + the RANS
scalars: `max_abs 0` on all seven datasets (un vn wn pn nut k omega, 190 252 032 cells each; `smoke_diff.txt`).

Found by this test: without `x_min_patch = inlet` the RANS inferred a
WALL on the inflow face ("RANS domain face x_min: wall (inferred)": the
Dirichlet tangential velocities of the inflow look like a no-slip wall to
the inference), which would have put the inlet plane into the SST wall
distance and the wall-cell classification. All four inis now declare it;
the case file carries only the BODY distance, so no re-prepare was needed.
The smoke runs themselves were made before the fix (identical inis on
both rank counts, so the identity result stands).

### The runs (HoreKa, 2026-10-27, build_gpu of `$W/moby-junction` at a2e0d9b, nvhpc 25.3)

The `accelerated` queue estimated a 21 h wait (priority from age only;
the account's fair share is 0.40), so the main case runs as a chain of
1-hour `dev_accelerated` chunks of 3000 steps (job 5213577 + successors,
`submit_calib.sh` without NO_CHAIN) and the three sensitivity cases wait
on `accelerated` (5 h each, jobs 5213563/5213564/5213565, `afterany` the
first chunk, which built the worktree). First chunk, node hkn0403, 4 x
A100-SXM4-40GB: 22.0-22.1 GB per GPU at 98-99 % utilisation (190 M
cells: 116 B per cell per GPU), `report_partition` 920 320 cross-rank
face cells (2 x 2 x 1), "RANS domain face x_min: inlet (declared)", no
warning in the set-up. Step time 0.75 s (2000 steps in
25 min, snapshot `calib_a1.5_y3.5_2000.h5` at 01:03).

**The cold start DIVERGED** (found on that snapshot, 01:10): at step
2000, t = 0.164, dt had collapsed to 2.9e-5 with the Courant number at
its limit 1.0 (the snapshot's `cfl` attribute is (Courant, Peclet) =
(1.0, 0.027)); max |u| 25.5, |v| 28.3, |w| 26.3, |p| 2731; 25.7 M cells
with |u| > 3, 22.3 M of them on level 5, 23 M in y < 0.5 (the fuselage
side), concentrated at x 4.5-5.5 (forebody), 6.5-8 (the junction) and
9-10 (the boom), only 20 000 near the wing cut at y > 3. The largest
values sit on the wing near the cut (7.94, 3.35, 10.07) and behind the
root trailing edge (8.3, 0.34, 9.89). Chain cancelled (5213577,
5213655), the three sensitivity jobs HELD. The eddy viscosity in the
snapshot (max 349 nu on level 5) bounds dt only at 5.4e-4, so the
collapse of dt is the Courant number of the diverged velocity, not a
diffusion limit. Diagnostics: 50 steps with a snapshot every 10 on
istmhydra (`diag_hydra.sh`, `.diag50.ini`), 400 steps with a snapshot
every 50 on 4 A100 (`../horeka/submit_diag.sh`, `.diag_base.ini`, job
5213711).

**Root cause (01:40): the run started FROM REST.** `[flow] initial_u`
defaults to 0 and the ini did not set it, so the inflow at u = 1 was
switched on against a fluid at rest: an impulsive, GLOBAL transient in
which every column of the 190 M-cell domain must accelerate at once
through the pressure. The approximate projection (Chebyshev-Jacobi, 12
sweeps per step) cannot transmit that: after one step (`proj_cheb12_1.h5`,
from rest) u = 1 only within a few cells of the inlet, p = 424 at level 0
next to the inlet, |div| h rms 3.8e-2 on level 0 (the converged NACA
field has 1e-7 to 1e-8, `snapshot_divergence.py`). The front then reaches
the fuselage with a divergent field: step 10 max |u| 3.4 just upstream of
the nose, step 20 max |u| 6.4 and p 2117 AT the nose on level 5, 2 M
cells above 3 U spread over the whole forebody (x 4-5.25, every y row
equally -- NOT the symmetry plane, NOT the wing cut), with a 2-cell
alternation in y. Fix: `initial_u = 1.0` in all four inis (the uniform
stream; the body's penalization then makes only the local impulsive
start every immersed-boundary run makes). The from-rest numbers above
are the record of the failure; the snapshots were deleted.
Diagnostics relaunched with the uniform field: istmhydra 50 steps every
10 (`diag_hydra.sh`) and six one-step projection variants (cheb 12/24/
48/96, rb 6/20, `diag_proj_hydra.sh`), HoreKa 400 steps every 50 (job
5213722) and the same with dt capped at 1e-4 (5213723).

**Uniform start (02:30): stable, and a second finding.** HoreKa 400 steps
(job 5213722, 0.85 s/step): no warning, no cell above 3 U after step 150;
the trailing-edge peak of the starting vortex is 1.9 U at step 10, 3.4 at
step 100 (25 cells), 2.4 at steps 200 and 400; dt settles at 6.9e-4
(Courant 1.0); the body-band mass error (level 5 rms) decays 6.6e-3
(step 1) -> 1.2e-3 (10) -> 6.8e-4 (400), levels 0-3 at 1e-4 or below.
With dt capped at 1e-4 (job 5213723) the same transient: peak 4.9 U at
t 0.01, 2.35 at t 0.04. One step from the uniform field, level-5 rms
mass error vs sweeps: Chebyshev 12 / 24 / 48 = 6.6e-3 / 4.4e-3 / 2.7e-3
(`proj_*_1.h5`), i.e. only as 1/sqrt(niter) at a cost linear in niter.

But the FAR-FIELD PRESSURE carries a global inlet-to-outlet mode
(`farfield_pressure.py`, level-0 slab means; the velocity stays at
1.000 +- 0.005): baseline <p> in the inlet chord +2.9 (step 100), +5.1
(200), -4.7 (400) against 0 at the outlet; with dt 1e-4: +14.8 (step
100), +7.2 (200), -10.5 (300), -13.3 (400). The period is about 400
STEPS in both (0.25 c/U at dt 6.5e-4, 0.045 at dt 1e-4) and the
amplitude scales as 1/dt: the unconverged smooth component of phi,
turned into pressure by the 1/dt of the projection -- a numerical mode
of the 12-sweep Chebyshev-Jacobi projection on a 176 x 56 x 320 base
grid (the lowest mode is damped by about pi^2/320^2 per Jacobi sweep),
not a flow. The converged NACA field (thin domain, 660 000 steps at dt
2.5e-5) has p in [-0.9, 0.5]. A linear-in-x mode cancels in the lift of
a closed body but puts +-0.9 across a chord into Cp. Variants running:
red-black sor 1.5 niter 6 (job 5213794), niter 20 (5213795, the
sailplane tutorial's setting), Chebyshev 48 (5213796), 400 steps each
with a snapshot every 50; read with `farfield_pressure.py` ON HoreKa
(h5py 3.14 is available on the login node).

Variant results (03:48; inlet-chord <p> at steps 50, 100, ..., 400; the
outlet chord stays within +-0.3, <u> within +-0.5 %):

| projection | s/step (4 A100) | inlet-chord <p> over steps 50..400 | verdict |
|---|---|---|---|
| Chebyshev 12 (base) | 0.85 | +0.3, +2.9, -4.7, +5.1, (200-400 not copied), -4.7 | sustained, about 5 |
| red-black sor 1.5, 6 | 0.97 | +3.5, -1.2, +0.2, -1.8, +3.7, -0.9, +2.1, -3.6 | sustained, 2-4 |
| red-black sor 1.5, 20 | 2.35 | -5.5, -3.6, -4.3, +2.4, +0.5, +4.7, +0.6, +2.6 | sustained, 3-5 |
| Chebyshev 48 | 2.24 | +12.4, -3.4, -2.3, +2.7, -0.8, -0.9, +0.9, -0.5 | DECAYS to below 1 by step 250 |

The mode flips sign every 50-100 steps; neither red-black setting damps
it, 48 Chebyshev sweeps do (four times the sweeps of the base at 2.6x the
step time: the projection is 0.038 s per sweep, the rest of the step
0.4 s). The main-case chain (job 5213797, Chebyshev 12, started 03:47)
is left running as the long-time test: its step-2000 snapshot (t about
1.3) says whether the mode decays with the starting transient at 12
sweeps; if not, the chain switches to niter 48 at the next chunk (the
ini is the authority on restart) and a case costs about 13 h instead of
5 h on 4 A100.

Chain, chunk 1 (Chebyshev 12), step 2000, t = 1.40 (04:18): inlet-chord
<p> +1.3 (from about 5 during the transient), |p| > 3 in 28 M cells (160 M
at step 100), max |u| 1.59 on the wing, no cell above 3 U, level-5 mass
error rms 7.7e-4, dt 7.6e-4; nut peaks at 178 nu at the wing cut (y =
3.5). The mode DECAYS with the starting transient at 12 sweeps, so the
chain stays at the handout's niter 12 and the three sensitivity jobs
were released (5213563/4/5 on `accelerated`). Checked again at every
chunk end (`farfield_pressure.py` on the newest snapshot, on HoreKa).

Chunk 1 end (04:27): 3000 steps in 2357 s (0.79 s/step), t = 2.12, dt
6.4e-4 now at the PECLET limit 0.5 (Courant 0.82: the eddy viscosity has
grown), max |u| 1.57, no cell above 3 U, inlet-chord <p> -1.8 (+1.3 at t
1.40): the mode settles at an amplitude of 1.5-1.8 at 12 sweeps instead
of vanishing. A smooth global mode cancels in the lift of a closed body
but offsets Cp by about 1 at the wing against an inlet reference and by
+-0.3 across a chord. PLAN: run the transient at niter 12 (cheap) and
switch the ini to niter 48 for the chunks that write the read-out
snapshots (the last 3 c/U), where the mode decays within 250 steps; the
same "polish" chunk follows each sensitivity run. The strip-load Cp
reference becomes the far field at the wing's x (`--pinf-box`), not the
inlet column. Solver-level consequence for the production run: the
approximate projection needs a coarse-grid correction for domains of this
size -- for the user to decide.

Chunk 2 end (05:10): step 6000, t = 3.57, 0.79 s/step, dt 4.2e-4 (Peclet
0.5, Courant 0.54), max |u| 1.54, inlet-chord <p> +0.31, NO cell with |p|
above 3, p max 1.44 at the junction leading edge (7.17, 0.55, 9.98): the
mode keeps decaying at 12 sweeps (5 -> 1.5 -> 0.3 over t 0.2 / 2.1 /
3.6), so the niter-48 polish may not be needed; decided at the end on
the inlet-chord <p> of the read-out snapshots (target below 0.05).

Chunk 3 end (06:30): step 9000, t = 4.78, dt 4.0e-4 (Peclet-bound), max
|u| 1.52, inlet-chord <p> -0.16, p max 1.40 at the junction leading edge;
the mode decays further (0.31 -> 0.16 over t 3.6 -> 4.8).

### Dry run of the read-out on the step-9000 snapshot (t = 4.78, NOT converged)

`strip_loads.py` on `calib_a1.5_y3.5_9000.h5` with the README call plus
`--pinf-box 7.0 8.5 1.0 3.0 14.0 16.0` (the far field above the wing):
608 s, 6.0 GB RSS on the workstation; 1.93 M pieces, 1.24 M inside the
domain, all on level 5; Cp from 99 % of them, Cf from 84 %. Geometry
check at the y = 1.0 m station from the pieces: chord-line angle +3.83
deg to the flow (= 1.5 pitch + 2.3 incidence), positive camber, chord
0.924 c (0.768 m): the pitched STL is right.

Preliminary flow at t 4.8 (the lift of an impulsive start is at about 85 %
by the Wagner function here): strip c_l 0.28 at the root rising to 0.41
at the cut; upper-surface Cf negative from x/c 0.60 to the trailing edge
(a trailing-edge separation over the rear 40 % of the flapped section,
f_sep 0.35 on every strip), lower surface attached; station y = 1.0 m
upper Cp plateau -1.22 at x/c 0.2-0.3, both surfaces at -0.48 at the TE
and the stagnation at +0.43: a uniform Cp offset of about -0.55 from the
far-field reference (the global pressure mode is still 0.16 at the inlet
chord; the lift is insensitive to it, Cp curves are not -- reference
them to the TE or to the local far field of the final snapshots).

XFOIL on the same section (`xfoil_polars.py` sequence, 240 panels, Re
2.77e5, alpha 3.8 to the chord line): free transition Ncrit 9 c_l 1.036
with the laminar bubble 0.53-0.88 (the step-1 basis, lifting line C_L
1.0); FORCED transition at 0.02/0.02 (the basis for a fully turbulent SST
box) c_l 0.793, c_d 0.0244, turbulent separation from x/c 0.75 to the
TE, Cp at the TE +0.02 (`liftingline/polars/sec_y1.0_turb.csv`). The box
at t 4.8 separates 0.15 c earlier and carries half that lift; judged on
the converged snapshots, with the fully turbulent polar as the reference
for the SST + wall-function model and the free-transition one for the
transfer of alpha (section 4.4).

Chunk 4 end (07:13): step 12000, t = 5.97, dt 4.0e-4, max |u| 1.66,
inlet-chord <p> +0.08, p max 0.91 at the junction leading edge (the
physical stagnation level); the mode is below 0.1 from t 6 on.

Chunk 5 end (08:14): step 15000, t = 7.16, dt 4.0e-4, max |u| 1.86 and p
max 1.25 both at the junction leading edge (7.23, 0.6, 10.04: the
horseshoe region, growing with time), inlet-chord <p> -0.03.

Chunk 6 end (09:16): step 18000, t = 8.34, max |u| 1.55, p max 1.10 at
the junction leading edge, inlet-chord <p> +0.015. Chunks start within
minutes of each other now (one per hour).

Chunk 7 end (09:57): step 21000, t = 9.52, max |u| 1.55, p max 1.16 at
the junction leading edge, inlet-chord <p> -0.003: the mode is gone.
Chunk 8 waits behind the user's own dev job (one dev job per user).

Chunk 8 end (10:42): step 24000, t = 10.70, max |u| 1.65, p max 1.02 at
the junction leading edge, inlet-chord <p> +0.004. Chunk 9 queues behind
the user's own CHT dev chain (alternating, one dev job per user).

### Read-out at t = 10.7 (step 24000; the lift has converged)

`strip_loads.py` on `calib_a1.5_y3.5_24000.h5` (same call as the dry
run): area-weighted mean strip c_l 0.381 (0.371 at t 4.8: +3 %), c_dp
0.021, c_df 0.005, c_m -0.025; c_l(y) 0.31 at the root strip rising to
0.41 at the cut; upper-surface separation (Cf < 0) from x/c 0.68 at the
y = 0.6 m station and about 0.65 at y = 1.0 m, f_sep 0.30-0.38 on every
strip, lower surface attached. The Cp is clean now (stagnation 1.03, TE
0.00 on both sides, the far-field box at (7-8.5, 1-3, 14-16) is a valid
reference once the global mode is gone). Station y = 1.0 m, upper Cp at
x/c 0.05 / 0.15 / 0.3 / 0.5 / 0.7 / 0.9: -0.53 / -0.67 / -0.75 / -0.51 /
-0.18 / -0.15; lower +0.02 / +0.01 / -0.02 / -0.08 / -0.14 / -0.19. XFOIL
fully turbulent at the same 3.8 deg: upper -1.24 / -1.03 / -1.06 / -0.80
/ -0.15 / -0.02, lower +0.51 / +0.31 / +0.21 / +0.12 / +0.04 / -0.06,
separation from 0.75, c_l 0.79. Both surfaces of the box carry less load
than XFOIL, as if the section saw about 2 deg less incidence; the
geometry is right (chord line +3.83 deg), and the free-slip cut mirrors
the wing (no tip downwash in the box, by design). The prime suspect is
the near-wall resolution: at c/512 the first cell centre sits at y+
about 15 (Cf 0.003-0.005), right above the laminar/log switch of the
two-layer wall function (`WF_YPLUS_LAM` 11.5, Weber Eqs. 4.39-4.42), so
the wall shear and the separation of the flap's pressure recovery are
computed in the buffer layer. A finer body band (c/2048, y+ 3-4 with
`wall_treatment = resolved`, about 4x the cells) would decide it -- a
cost question for the user. The sensitivities of section 4.3 remain
measurable within the model.

Shape of c_l(y) against the lifting line (section 4.3, first row), both
normalised by the value at the cut: box 0.74 / 0.80 / 0.84 / 0.90 / 0.93
/ 0.95 / 0.96 / 0.97 / 0.99 at y = 0.62 / 0.92 / 1.22 / 1.52 / 1.82 /
2.12 / 2.42 / 2.72 / 3.32 c; lifting line (`strips_cl1.0_ncrit12.csv`,
no fuselage) 1.08 at 0.61, 1.02 at 0.9-1.2, 1.00 at the cut. The box
loses a quarter of the lift at the root strip and recovers over 2 c of
span: the fuselage junction of the wall-function RANS (a separated
corner), not the lifting line's 5 % root bump. The level (0.38 vs 1.04,
expected 0.03-0.06 below) and the shape both fail the handout's
expectation; the finding of step 2 is that the SST + wall-function box
at c/512 does not reproduce the attached, lightly loaded flap flow the
lifting line assumes, so alpha cannot be transferred on its basis
without the near-wall resolution question settled first.

Chunk 9 end (12:23): step 27000, t = 11.88, max |u| 1.57 (now at the
wing near the cut, 7.1, 3.2), p max 1.17 at the junction leading edge,
inlet-chord <p> +0.001. Snapshots from 27000 on are in the read-out
window (the last 5 c/U).

Chain end (2026-10-27 14:17): twelve chunks, the last trimmed to 1925
steps by t_final; final snapshot `calib_a1.5_y3.5_34925.h5` at t = 15.00,
kept on HoreKa and in `~/horeka_stage/junction/`. Read-out window
snapshots: 27000 (t 11.88), 30000 (13.1), 33000 (14.2), 34925 (15.0).

The three sensitivity cases ran their `accelerated` job on 2026-10-27
from 15:50 / 16:04 / 16:06 (nodes hkn08xx, 0.79 / 0.80 / 0.88 s per step
for a1.0 / a2.0 / y4.5): the 20 000-step chunk ended at step 20000, t
about 9, after 4.4-4.9 h, inside the 5 h limit but with no time left for
a second chunk; the jobs were killed at the limit before writing their
exit line (the step-20000 snapshots are complete). Resubmitted 2026-10-29
09:20 as jobs 5216676 / 5216677 / 5216678 (restart from 20000, to
t_final).

Final-state checks (`farfield_pressure.py`, `snapshot_extrema.py` on
HoreKa, 2026-10-29): main case at t 15.00: inlet-chord <p> +0.002, max
|u| 1.55 and p max 1.20 at the junction, no cell above 3 U, <u> outlet
1.0014. Sensitivity cases at step 20000: a1.0 at t 9.78 (dt 4.3e-4, p max
1.00), a2.0 at t 8.75 (dt 3.7e-4, p max 1.17), 4.5 c cut at t 9.13 (dt
3.9e-4, p max 1.05); all mode-free (|<p>| at the inlet below 0.006) and
healthy. The outlet through the boom: <u> on the outlet chord 1.0011-
1.0015 in every case, no backflow event in any log.

## The read-out (section 4.3 of the handout), 2026-10-29

Binary: `$W/moby-junction/build_gpu/moby_solve` at a2e0d9b, nvhpc 25.3,
4 A100-40GB; inis as committed (initial_u 1.0, x_min_patch inlet,
Chebyshev niter 12, cflmax 1.0, pecletmax 0.5); `strip_loads.py` with the
README call plus `--pinf-box 7.0 8.5 1.0 3.0 14.0 16.0`. Main case: the
four snapshots 27000 / 30000 / 33000 / 34925 (t 11.88 / 13.10 / 14.20 /
15.00). Sensitivity cases: their step-20000 snapshots (t 9.78 / 8.75 /
9.13 for alpha 1.0 / alpha 2.0 / the 4.5 c cut; continuations to t 15
are queued, the main case moved by 0.3 % between t 10.7 and 15).

| measurement | result | against the handout's expectation |
|---|---|---|
| c_l(y), alpha 1.5, 3.5 c | area-weighted mean 0.381 (0.3806 / 0.3805 / 0.3803 / 0.3816 over the four snapshots: scatter 0.001, the flow is steady); c_l(y) 0.31 at y 0.62 c, 0.34 at 1.0, 0.37 at 1.4, 0.39 at 1.8, 0.40 at 2.4, 0.42 at the cut; c_dp 0.0213, c_df 0.0049, c_m -0.025 | level 0.38 vs 1.04 (expected 0.98-1.01); shape 0.74 of the cut value at the root strip vs 1.08 (lifting line): FAILS both |
| dc_l/dalpha, root region y 0.6-1.0 m | c_l 0.280 / 0.310 / 0.344 at y 0.62 c and 0.320 / 0.340 / 0.384 at 1.22 c for alpha 1.0 / 1.5 / 2.0; mean over the wing 0.348 / 0.381 / 0.414: 0.066 per deg, linear (0.0328 and 0.0325 per half degree) | about 0.1 per deg expected; two thirds of it, the separated flap |
| junction Cp, y = 0.6 m, alpha 1.0 / 1.5 / 2.0 | upper plateau Cp(0.3 c) -0.79 / -0.81 / -0.81, Cp(0.5) -0.48 / -0.47 / -0.45; lower Cp(0.3) -0.15 / -0.12 / -0.08; suction peak -1.28 / -1.44 / -1.15 at 0.16 / 0.10 / 0.15 c (the peak values are the noisy leading-edge band) | the upper surface barely answers to 0.5 deg; the load change comes from the lower surface (+0.035 per 0.5 deg) |
| junction Cp, y = 1.0 m | upper Cp(0.3) -0.73 / -0.75 / -0.76, lower Cp(0.3) -0.05 / -0.03 / +0.01; TE Cp -0.10 / -0.14 / -0.11 (no recovery: separated) | same |
| separation line (upper Cf < 0 from, x/c, x > 0.3), alpha 1.0 / 1.5 / 2.0 | y = 0.6 m: 0.73 / 0.68 / 0.62; y = 1.0 m: 0.64 / 0.62 / 0.62; f_sep 0.30-0.38 on every strip | "a turbulent separation near the TE at most" expected: it is a 35 % chord separation on the whole span, moving 0.05 c forward per 0.5 deg at the junction |
| 3 c vs 4 c cut, alpha 1.5 | root-region c_l(y) 0.312 / 0.334 / 0.350 / 0.378 (4.5 c cut, t 9.1) vs 0.311 / 0.333 / 0.350 / 0.377 (3.5 c cut at t 10.7) at y 0.62 / 0.92 / 1.22 / 1.52 c: within 0.3 % | the 3 c cut STANDS for production |
| the outlet through the boom | every log clean, no backflow event, outlet <u> 1.0011-1.0015, final |div| h rms 7e-4 on the body band (1e-4 on levels 0-3) | the 1e-5..1e-4 of the NACA is not reached on the band with 12 sweeps (see the projection note above) |
| prepare scaling | see "Prepare scaling" above: 372-420 k leaves, 52-60 min at 48 EPYC ranks, 2.1-2.5 GB per rank | 7x the handout's leaf estimate |

Snapshots kept: `calib_a1.5_y3.5_34925.h5` (the initial field of step 7)
on HoreKa and in `~/horeka_stage/junction/`, plus the read-out window
and the step-20000 restart files on HoreKa; strip tables and pieces in
`~/horeka_stage/junction/*_strips.csv`, `*_pieces.npz`, `*_cp_y*.csv`.

