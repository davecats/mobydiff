# Step 3b: the projection-residual floor on the unforced SD7003 bubble

Handout `docs/next_session_sailplane_junction.md` section 5.1. The question:
does the residual the iterative projection leaves behind act as a disturbance
floor that moves the transition front of a laminar separation bubble? The
test is the SD7003 at Re_c = 6e4, alpha = 4 deg, with no forcing and no model,
run from one developed state under three projections -- Chebyshev `niter = 12`
(the production setting), Chebyshev 24, red-black SOR 1.5 with 12 sweeps --
and read for the separation, transition and reattachment positions. A moved
front is a finding (more iterations in the DNS band, never a compensating
trip).

No scale-resolving SD7003 case existed (`validation/sd7003`, removed
2026-09-26, was the gamma-Re_theta RANS benchmark; `tutorials/naca/LES/` is
an unbuilt plan for the NACA 0012). This directory is the case; the
coordinates `sd7003.dat` come from the removed directory (git 6f2c6c1^).

## The case

| item | value |
|---|---|
| geometry | SD7003, chord 1 along x, leading edge at (4.5, 6.0), not rotated; incidence 4 deg in the free stream (`[case.airfoil]`), span z periodic, 0.25 c (Galbraith & Visbal 2010: 0.2 c) |
| domain | 12 c x 12 c x 0.25 c, base c/32 isotropic (384 x 384 x 8), `nb = 8`, octree (`refine_dims = xyz`, as the production grid) |
| finest level | 5 = c/1024 on the whole surface and in the box x/c 0.09..1.06, y 5.969..6.125: the bubble, the detached shear layer and the turbulent layer lie inside one level, no 2:1 interface between the wall and 0.07 c above it; Delta+ ~ 4 behind the reattachment (the production core's spacing in wall units), delta*/3 in the laminar layer at x/c 0.2 |
| coarser boxes | c/512 to 0.5 c behind the trailing edge, c/256 to 1.5 c, c/128 to 3.5 c, c/64 to the outlet, following the wake at +4 deg |
| cells | 101 116 leaves = 51.8 M (34.3 M at c/1024, 9.0 M at c/512); case file 2.5 GB, prepared in 12 min on 8 CPU ranks |
| model, projection | none; Chebyshev 12 unless the branch says otherwise |
| time step | 4e-4 (`cflmax = 1`, which binds at times: Courant 0.75..1.0), 2500 steps per c/U |
| 3D seed | the geometry, grid and start are exactly span-uniform, and `[flow] initial_noise` is a dead key (read, never used), so stage 0 runs a weak random wall-normal trip at x/c 0.3 (`trip_amp = 0.02`, velocity gain ~ 1e-3 U) for the first 3 c/U; it is OFF in every later stage |

Stages: `sd_seed.ini` (7500 steps, t 0 -> 3, trip on), `sd_dev.ini` (22 500,
t 3 -> 12, the common developed state), then the branches `sd_c12.ini`,
`sd_c24.ini`, `sd_rb.ini` (25 000 each, t 12 -> 22, a snapshot every 0.1 c/U).
`sd4_*.ini` is the same chain one level coarser (c/512, 22.0 M cells, dt 8e-4)
-- small enough for the workstation's RTX 3060, run first as the shakedown
and as the Delta+ ~ 8 data point (`chain4.sh`).

Read-out: `sd_reduce.py` turns every snapshot into span-averaged first and
second moments on the finest lattice (2.5 GB -> a few MB; the chain deletes
the snapshot), `sd_bubble.py` averages them in time and prints separation,
transition (max_y -u'v' = 0.001 U^2, the criterion of Galbraith & Visbal) and
reattachment per run, with the amplification curves max_y u'u', w'w' (x/c),
whose level ahead of the bubble IS the disturbance floor.

```bash
~/ibmc/bin/python ../../../../tools/make_airfoil_stl.py selig --file sd7003.dat \
    --xc 4.5 --yc 6.0 --chord 1.0 --lz 0.25 --resample 2880 --out sd7003.stl
mpirun -n 8 ../../../../build_cpu/moby_prepare sd_seed.ini sd.case.h5      # sd4_seed.ini sd4.case.h5
setsid nohup ./chain4.sh > chain4.out 2>&1 < /dev/null &
python3 sd_bubble.py --run c12 'sd4_c12_*.red.npz' --run c24 'sd4_c24_*.red.npz' --run rb 'sd4_rb_*.red.npz'
```

## Timing (2026-11-16, binary of be57d25, Chebyshev 12)

| grid | device | s/step | device memory | note |
|---|---|---|---|---|
| c/512, 22.0 M | RTX 3060, alone | 0.86 | 7.0 GB | 39 ns per cell and step; projection 71 % of the step in the profiled 300-step run |
| c/1024, 51.8 M | A6000 (istmcetus), SHARED with another user's job at 99 % load | 1.73 | 16.2 GB | an upper bound; projection 84 % |
| c/1024, 51.8 M | A6000 alone | ~0.9 (estimate) | | from the step-3a pair on the same two cards (A6000 = 0.43 x the RTX 3060 per cell) |
| c/1024, 51.8 M | 1 H100 (HoreKa 2, jobs 133191 / 133253) | 0.396 | 18.6 GB host | Chebyshev 24: 0.722; red-black 12: 0.868 |
| c/1024, 51.8 M | 2 H100, one node | 0.208 | | 130 k face cells across ranks |
| c/1024, 51.8 M | 4 H100, one node | 0.111 | | Chebyshev 24: 0.202 (x 1.81); red-black 12: 0.244 (x 2.19); 232 k face cells across ranks; 3.55 x one GPU |

The whole c/1024 chain is 105 000 steps (30 000 to the common state, 3 x
25 000): 4.8 h on 4 H100, about 20 GPU-hours with the snapshot output. The
regular `gpu-h100` queue estimated a week on 2026-11-16 and the dev node an
hour, so the chain runs in one-hour jobs of `dev-gpu-h100`
(`../horeka/chain_sd7003_hk2.sh`, first job 133264, run dir
`$(ws_find junction)/sd7003_run`). The first timing job lost its 4-rank
runs to `--map-by numa --bind-to core` (two ranks failed to start);
`--bind-to none` since.

## Results

(none yet: the c/512 chain started 2026-11-16 12:33 on the RTX 3060)
