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

### c/512 twin (2026-11-17, RTX 3060, binary of be57d25)

Three branches from the same developed field (step 15 000, t 11.8), 12 500
steps each, 78 snapshots in the window t 14.0..21.5 (7.5 c/U);
`sd_bubble.py --run c12 'sd4_c12_*.red.npz' --run c24 ... --run rb ... --out viz4`
(`viz4/bubble.png`, `viz4/bubble_*.csv`). Positions in x/c on the suction
side; the two half windows (t 14.0..17.8 and 17.8..21.5) in brackets.

| projection | separation | transition (-u'v' 1e-3) | reattachment | c_l | c_d |
|---|---|---|---|---|---|
| Chebyshev 12 | 0.256 (0.258, 0.253) | 0.534 (0.532, 0.544) | 0.726 (0.712, 0.737) | 0.588 | 0.0300 |
| Chebyshev 24 | 0.257 (0.261, 0.253) | 0.538 (0.540, 0.552) | 0.728 (0.723, 0.739) | 0.590 | 0.0299 |
| red-black 12 | 0.264 (0.265, 0.262) | 0.532 (0.528, 0.552) | 0.705 (0.698, 0.708) | 0.582 | 0.0288 |
| literature | 0.22..0.30 | 0.53..0.58 | 0.65..0.70 | | |

- **Chebyshev 12 against 24: nothing moves.** Transition 0.004 c and
  reattachment 0.002 c apart on the full window, against a half-window
  scatter of 0.01..0.025 c in every run (the bubble is still lengthening
  slowly: the second half is later in all three). Doubling the iterations
  buys nothing at this resolution.
- **Red-black 12 reattaches 0.02..0.03 c earlier** than either Chebyshev run,
  in both half windows, with the same transition position; separation 0.007 c
  later, c_l 1 % and c_d 4 % lower. That is about the size of the scatter and
  wants the c/1024 run before it is called a difference between the solvers.
- **The upstream level of w'w' does not set the front.** At x/c 0.1..0.2
  red-black carries w'w' 1.7e-9 (flat from the leading edge), Chebyshev 12
  4e-11 and Chebyshev 24 1e-10; from x/c 0.25 on the three curves are one, and
  they saturate at the same place. u'u' at x/c 0.1..0.2 is 1e-6..1e-5 in all
  three and changes by a factor 3 between the half windows: it is the
  two-dimensional unsteadiness the shedding bubble radiates upstream, not a
  projection residual.
- The bubble is 0.03..0.06 c longer than the literature band at reattachment
  (c/512, Delta+ about 8 behind reattachment); separation and transition are
  inside it.

### c/1024 (2026-11-18, HoreKa 2, 4 H100, binary of 0c0c01f)

52 M cells, three branches from step 30 000 (t 12.0), 25 000 steps each, 78
snapshots in t 14.0..21.5; `viz5/bubble.png`, `viz5/bubble_*.csv`. Half
windows in brackets. Step times 0.111 (Chebyshev 12), 0.201 (Chebyshev 24),
0.243 s/step (red-black 12); about 19 GPU-hours of stepping for the whole chain.

| projection | separation | transition (-u'v' 1e-3) | reattachment | c_l | c_d | u'u' at x/c 0.1..0.2 |
|---|---|---|---|---|---|---|
| Chebyshev 12 | 0.257 (0.260, 0.254) | 0.500 (0.552, 0.495) | 0.731 (0.740, 0.721) | 0.595 | 0.0256 | 9.0e-5 |
| Chebyshev 24 | 0.262 (0.263, 0.261) | 0.507 (0.506, 0.508) | 0.721 (0.720, 0.721) | 0.592 | 0.0250 | 1.8e-5 |
| red-black 12 | 0.269 (0.270, 0.269) | 0.512 (0.516, 0.509) | 0.700 (0.700, 0.700) | 0.587 | 0.0243 | 3.2e-6 |
| literature | 0.22..0.30 | 0.53..0.58 | 0.65..0.70 | | | |

- **The three-dimensional breakdown does not move with the projection.**
  w'w' of the three runs is one curve from x/c 0.3 to saturation at x/c
  0.65, although upstream red-black carries 4e-8 (flat from the leading
  edge) against 1e-11 for the Chebyshev runs. As at c/512, the upstream
  spanwise level is not what seeds the bubble.
- **Chebyshev 12 is not quiet at this resolution.** It carries a
  SPANWISE-COHERENT streamwise fluctuation over the laminar and separated
  layer (max u'u' 9e-5 at x/c 0.1..0.2, 1e-2 at x/c 0.4) that Chebyshev 24
  cuts by 3..5 and red-black 12 by 10..30, and its bubble wanders: transition
  0.552 in the first half window and 0.495 in the second, reattachment
  0.740 and 0.721. Chebyshev 24 and red-black 12 are stationary to 0.007 c
  between the halves. At c/512 the same comparison showed no ordering
  (Chebyshev 24 above 12): twelve fixed iterations lose ground as the grid is
  refined.
- **Red-black 12 reattaches 0.02 c ahead of Chebyshev 24** (0.700 against
  0.721, both stationary) and separates 0.007 c later, the same sign and
  size as at c/512: a difference between the two solvers that the grid does
  not remove. Red-black sits on the edge of the literature band; which of the
  two is the converged answer is NOT decided by this test (it wants a
  Chebyshev 48 or red-black 24 branch, 25 000 steps, about 3 GPU-hours).
- c/512 to c/1024: separation and reattachment within 0.01 c for each
  solver; transition by the -u'v' threshold 0.02..0.03 c earlier on the fine
  grid; c_d 15 % lower (0.025 against 0.030), c_l 1 % higher.

**Verdict for the production core.** The projection does not set where the
flow turns three-dimensional, so no trip and no compensation is called for.
Chebyshev 12 is too loose for a c/1024-class DNS band (coherent unsteadiness,
a wandering bubble); Chebyshev 24 (x 1.81 per step) or red-black 12 (x 2.19)
is the price there, and the 0.02 c between those two is the uncertainty this
test leaves on a reattachment line.

Chain notes: 9 one-hour jobs on `dev-gpu-h100` between 2026-11-16 17:29 and
2026-11-18 05:20. The chain stopped twice without a successor (QoS `dev`
allows 3 submitted jobs per user and other jobs of the same user held the
slots; the retry loop dies with the job) -- a watcher on the workstation
that resubmits when no chain job is queued is the remedy. The reduced
snapshots and force traces are in the run directory on HoreKa 2 and (not in
git) here; the end fields `sd_{c12,c24,rb}_55000.h5` on HoreKa 2 only.
