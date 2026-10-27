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
CHAIN_PROGRESS_PLACEHOLDER
