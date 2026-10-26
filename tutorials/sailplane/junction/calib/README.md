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
