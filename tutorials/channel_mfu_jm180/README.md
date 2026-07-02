# Minimal flow unit channel at Re_tau = 180

The smallest periodic box (2 x 2 x 2, Jimenez & Moin 1991) that still sustains
near-wall turbulence, on a 48 x 128 x 80 grid. Because it is small it runs
quickly, so it is the right case for **experimenting** with the solver.

## Files

- `input.ini` — the run configuration (6-iteration Chebyshev-Jacobi pressure
  projection).
- `restart.h5` — a turbulent field to restart from (~16 MB, shipped with the
  repo so the case runs out of the box).

## Run it

```bash
# build once (from the repo root)
./compile.sh cpu

# run (single MPI rank is fine)
mpirun -n 1 ./build_cpu/main tutorials/channel_mfu_jm180/input.ini
```

It writes wall-normal statistics to `channel_stats.h5`, which you can plot with:

```bash
python3 tools/plot_stats.py channel_stats.h5
```

## Experiment: your own wall boundary conditions

This case is ideal for trying time- and space-varying wall boundary
conditions. There is exactly ONE place to edit: the subroutine
`wall_velocity` in `src/modules/wall_bc.f90`. It returns the velocity of the
bottom and top walls as a function of the position on the wall `(x, z)`, the
time `t`, and `v_sensed` (the wall-normal velocity a short distance into the
flow, for feedback control). After editing, rebuild with `./compile.sh cpu`
and rerun. Examples you can implement in that one function:

- **Spanwise oscillating wall:** `w = A*sin(omega*t)`
- **Streamwise travelling wave of spanwise wall velocity:**
  `w = A*sin(kx*x - omega*t)`
- **Blowing / suction:** `v = v0` (keep the net mass flux near zero)
- **Opposition control:** `v = -v_sensed`. Here the driver already senses the
  wall-normal velocity on a detection plane `SENSOR_OFFSET` cells into the flow
  (y+ ~ 11 by default, tunable at the top of `wall_bc.f90`) and passes it in;
  you just oppose it. Add a gain if you like: `v = -gain*v_sensed`.

## Prescribing the start time

Set `t_start` under `[time]` in `input.ini` to choose the initial time of the
run; it overrides the time stored in `restart.h5`.
