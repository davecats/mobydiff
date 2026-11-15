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

Step times on the RTX 3060: 0.9 s/step alone, 1.47 s/step with two other
runs on the card (tsbase, 2026-11-15).

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
- (the TS growth read-out: pending, chain running since 2026-11-15 20:52)
