# Turbulent channel flow at Re_tau = 180

A full plane channel (Kim, Moin & Moser 1987) on a 256 x 128 x 256 grid:
periodic in x (streamwise) and z (spanwise), no-slip walls in y, driven by a
constant streamwise body force. This is the classic DNS validation case.

## Files

- `input.ini` — the run configuration (already set to the production numerics:
  6-iteration Chebyshev-Jacobi pressure projection).
- `restart.h5` — a developed turbulent field to restart from (256^3, ~270 MB).
  It is **not** shipped inside the git clone because of its size; ask the
  instructor for it, or generate one yourself.
- `channel_stats.h5` — a **sample** statistics file so you can try the plotting
  tool straight away, without running the (large) simulation first.

## Run it

```bash
# build once (from the repo root)
./compile.sh cpu

# run (single MPI rank is fine)
mpirun -n 1 ./build_cpu/main tutorials/channel_kmm180/input.ini
```

The run samples wall-normal statistics and writes them to `channel_stats.h5`.

## Look at the statistics

```bash
python3 tools/plot_stats.py channel_stats.h5
```

This shows the mean velocity profile (also in wall units), the Reynolds
stresses, and the total-stress budget split into its viscous and turbulent
parts, and it prints u_tau and Re_tau.

## Experiment: your own wall boundary conditions

There are two single-edit-point control hooks:

- **Wall boundary conditions:** `src/modules/wall_bc.f90`, subroutine
  `wall_velocity` — a moving wall, an oscillating wall, blowing/suction, a
  travelling wave, or opposition control (it receives the sensed wall-normal
  velocity).
- **Steady volume force:** `src/modules/volume_force.f90`, subroutine
  `body_force` — e.g. the Schlatter & Canton streamwise vortices. Off unless
  `input.ini` sets `[force] enabled = true` and `type = steady`.

Edit the hook, rebuild, and rerun. The `channel_mfu_jm180` tutorial is a much
smaller, faster case for such experiments (see its README for details).

## Prescribing the start time

`restart.h5` stores the time at which it was saved. Set `t_start` under
`[time]` in `input.ini` to override it — for example, to restart a developed
field and count time from zero.
