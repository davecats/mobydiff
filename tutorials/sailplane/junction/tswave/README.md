# The TS-wave test (sailplane junction, step 3a)

Section 5.1 of `docs/next_session_sailplane_junction.md`: does a
Tollmien-Schlichting wave grow through the immersed plate at the production
resolution as Orr-Sommerfeld says it should? The go / no-go of the DNS core:
the penalized wall (second order at fine resolution) must neither damp nor
seed the wave beyond a few percent of the linear growth. Everything in inlet
momentum thicknesses theta_in = 1, U = 1, `[flow] re` = Re_theta,in.

## Files

- `orr_sommerfeld.py` -- the reference: the spatial Orr-Sommerfeld
  eigenvalue of the Blasius layer (Chebyshev collocation, quartic pencil in
  alpha, parallel flow). Check: Jordinson (1970) Re_delta* 998, omega 0.1122
  -> alpha = 0.30859 - 0.00571 i here against 0.30866 - 0.00573 i (0.02 % /
  0.4 %), converged over N 100-140 and y_max 30-60.
  `--F 0.9e-4 --re-range 500 1800 100`: branch I at Re_delta* 770, branch II
  1340, peak -alpha_i 0.0078 at 1100, N = 2.5 across (amplification e^2.5).
- `tsbase.ini` -- the base flow: the immersed flat plate (the analytic wavy
  wall at zero amplitude, `wall_offset = 10.8`) on a uniform h = 0.45 grid,
  5560 x 424 x 8 cells (18 M, 35 445 blocks of 8^3 after the 1390 buried
  ones; the prepare takes 68 s and 5 GB on one CPU rank), Re_theta,in 250,
  the Blasius u + v inlet raised to the plate (`blasius_y0`), the pressure
  outlet at x_max and at the top as in `validation/blasius`, impulsive
  start (`initial_u = 1`), two flow-throughs (12 500 steps of dt 0.4);
  `tsbase2.ini` continues it by 15 000 steps to t 11 000 (see Results:
  at t 5000 the downstream half was still in the start-up transient).
  Resolution: 6.8 cells per delta* at the forcing station (Re_delta* 753),
  8.9 at 1000, 12 at 1340 -- the production core's delta*/7.
- `tsw.ini` -- the forced run from the base field: the vibrating ribbon
  (`[force] type = trip` with the two keys added for this test,
  `trip_omega` = harmonic spanwise-uniform forcing and `trip_y0` = the
  plate height; `blasius_y0` is the third) at x = 200, F = 0.9e-4
  (omega 0.0225, period 279 = 698 steps), amplitude 1e-5, 16 000 steps
  (6400 theta/U: the wave front at c_g ~ 0.35 U passes the end of the
  read-out range, x 2127, plus three periods); `tsw_phases.ini` then
  samples one period in 8 snapshots (every 87 steps).
- `ts_analysis.py` -- the read-out: delta*, theta, H of the base flow
  against Blasius; the disturbance u' = u - u_base projected on omega over
  the 8 phases, its maximum over y as A(x), d ln A/dx against -alpha_i at
  the local Re_delta*, and the N factor from the forcing station against
  the integrated Orr-Sommerfeld curve (`viz/ts_growth.png`).

## Run

```bash
module load toolkits/nvhpc/25.9; export OMPI_MCA_hwloc_base_binding_policy=none
mpirun -n 1 ../../../../build_cpu/moby_prepare tsbase.ini          # tswave.case.h5
mpirun -n 1 ../../../../build_gpu/moby_solve tsbase.ini > tsbase.log
./chain.sh            # tsbase2, tsw, tsw_phases, then the two analysis calls below
python3 ts_analysis.py --base tsbase_27500.h5 --ini tsw.ini                    # base check
python3 ts_analysis.py --base tsbase_27500.h5 --forced tswp_*.h5 --ini tsw.ini # growth
```

Step times on the RTX 3060 (18 M cells): 0.69 s/step alone (tsw), 1.04
(tsbase2, shared part of the time), 1.47 with two other runs on the card
(tsbase, 2026-11-15); the analysis takes 4 min for the 9 snapshots.

## Results

- BASE FLOW AT t 5000 (two flow-throughs, 2026-11-15): the layer follows
  Blasius within 0.2 % in delta* and 0.7 % in theta up to x 800 (Re_delta*
  1000, H 2.577..2.586 against 2.591), then runs THIN downstream: delta*
  -2 % at x 1200, -6 % at 1600, -12 % at x 2127 with H falling to 2.45.
  The deviation is a monotone function of the station's age ratio
  U t / (x + x_v) (3.7 at x 800, 1.9 at x 2127), i.e. the impulsively
  started plate still relaxing from its thinner, fuller Rayleigh state --
  Stewartson's steady state needs a ratio of about 4 -- not a wall or
  pressure-gradient effect (the inlet station, ratio 8.8, is within 0.1 %).
  So the base runs on to t 11 000 (`tsbase2.ini`; ratio 4.1 at x 2127)
  before the forcing; the forced stage was lengthened to 16 000 steps at
  the same time (the earlier 8400 only carried the wave front to x ~ 1400).
  A drifting base would otherwise leak into the omega projection: a ramp
  Delta u over one period projects to 0.33 Delta u with 8 phases, and the
  t 5000 drift at x 2127 gives ~1e-3, ten times the wave.
- BASE FLOW AT t 11 000 (`tsbase2.ini`, 2026-11-16): Blasius within 0.36 %
  in delta* and 1.0 % in theta over x 100..2127, H 2.574..2.586 (2.591).
  The immersed plate carries the laminar layer; the base is converged.
- TS GROWTH (forced 16 000 steps to t 17 400, 8 phases, `analysis.log`,
  `viz/ts_growth.png`, `viz/amax.csv`): the wave is linear (u'/U 1.2e-5 at
  the reference station x 262 behind the ribbon's near field, 7.8e-5 at
  its peak). Against parallel Orr-Sommerfeld at the local Re_delta*:

  | Re_delta* | x | -alpha_i OS | d lnA/dx measured | N_OS | N_meas |
  |---|---|---|---|---|---|
  | 900 | 527 | 0.00143 | 0.00117 | 0.24 | 0.31 |
  | 1000 | 784 | 0.00184 | 0.00166 | 0.67 | 0.68 |
  | 1100 | 1068 | 0.00177 | 0.00164 | 1.20 | 1.16 |
  | 1200 | 1378 | 0.00128 | 0.00119 | 1.68 | 1.60 |
  | 1300 | 1716 | 0.00041 | 0.00030 | 1.98 | 1.88 |

  Branch II: the amplitude peaks at Re_delta* 1334 (OS neutral at 1337).
  The peak N is 1.90 against 2.00 (the OS maximum 2.00). The eigenfunction
  |u'|(y) at Re_delta* 1000 lies on the OS one within 0.06 of the peak,
  its inner maximum at y/delta* 0.51 (OS 0.55), the node at 2.0. The
  growth rate over Re_delta* 850..1300 is 0.89 of parallel OS -- an 11 %
  deficit, uniform over the band; non-parallel corrections to the u'_max
  growth go the other way at this F, so the deficit is numerical: the
  grid (delta*/6.8..12) or the wall. A ripple of +- 5 % in A at the TS
  wavelength appears beyond Re_delta* 1250, in the last 1000 theta before
  the outlet (a potential-field reflection at omega); it is outside the
  averaging band. VERDICT so far: the wall neither seeds nor damps the
  wave qualitatively (shape, branch II, N within 5 %), the 11 % growth
  deficit is referred to the h/2 check below.
- THE h/2 CHECK (`h2/`, istmcetus A6000, 2026-11-16 04:43 -> 09:42; 0.41
  s/step): the same layer at h = 0.2 theta_in (delta*/14 at the ribbon) on
  the sub-range Re_delta* 700..1150, inlet Re_theta,in 270 so that the
  ribbon at x 100 sits at Re_delta* 755 as in the parent, lx 1100.8,
  5504 x 568 x 8 = 25 M cells (45 408 leaves, prepare 67 s / 6.5 GB), dt
  0.18, base 32 000 steps (t 5760, age ratio 3.7 at the end of the
  read-out), forced 17 800 steps (omega 0.0243), 8 phases every 180 steps
  (`chain_h2.sh`; the analysis runs on the workstation, istmcetus has no
  h5py). The base at t 5760 is Blasius within 0.11 % (delta*) / 0.28 %
  (theta), H 2.585..2.590 -- three times closer than the parent's 0.36 /
  1.0 %. The wave (reference x 157, Re_delta* 784, u'/U 1.3e-5):

  | Re_delta* | -alpha_i OS | d lnA/dx at h/2 | ratio h/2 | ratio h (parent) |
  |---|---|---|---|---|
  | 900 | 0.00155 | 0.00135 | 0.87 | 0.82 |
  | 950 | 0.00183 | 0.00177 | 0.97 | 0.86 |
  | 1000 | 0.00199 | 0.00190 | 0.95 | 0.90 |
  | 1050 | 0.00201 | 0.00192 | 0.96 | 0.92 |
  | 1100 | 0.00191 | 0.00214 | 1.12 | 0.93 |

  (the OS rates are per theta_in of each case; the parent's column is read
  at the same Re_delta* from `viz/amax.csv`). The growth rate over
  Re_delta* 850..1092 averages 0.947 of OS against 0.889 on the parent
  grid (a +-3 % ripple at the TS wavelength rides on A over the whole of
  this shorter domain, the outlet being 2.5 wavelengths past the read-out
  end); the eigenfunction at Re_delta* 1000 is within 0.015 of its peak
  (parent 0.058), the inner maximum at y/delta* 0.57 (OS 0.55; parent
  0.51). The last rows (Re_delta* > 1050) sit within 200 theta of the
  outlet and carry its ripple. A ninth phase file (the run's final step,
  off the 180-step lattice because the restart step 49 800 is not a
  multiple of it) had to be dropped: `ts_analysis.py` now does that itself
  and prints the phase spacing.
- VERDICT 2026-11-16: GO. Halving the spacing halves the growth deficit
  (11 % -> 5.3 %) and closes the eigenfunction, with the same immersed plate
  under both grids: the deficit is the grid, not the wall. The penalized
  wall neither seeds nor damps the Tollmien-Schlichting wave within the
  resolution of the test. The design number for the DNS band: at
  delta*/7 (the production core as planned) TS growth rates come out
  ~10 % low, at delta*/14 ~5 % low (Richardson, second order: ~2.6 % at
  the limit); a 10 % low growth rate moves an N = 9 transition estimate
  downstream by ~10 % of the amplification length. The production core's
  TS band should therefore be refined to delta*/10..14, or the 10 % carried
  as a known bias.
