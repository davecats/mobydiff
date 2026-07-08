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

## Visualise a field in ParaView

By default the case writes only statistics. To also write 3D field snapshots,
set a positive `field_interval` under `[output]` in `input.ini` (e.g.
`field_interval = 500`); the solver then writes a `channel_field_<step>.h5`
every that many steps.

The field HDF5 stores the solver's block-refinement layout, which ParaView
cannot open portably on its own. Convert a snapshot to a native VTK file first:

```bash
python3 tools/field_to_vtr.py channel_field_500.h5      # -> channel_field_500.vtr
```

Open the resulting `channel_field_500.vtr` in ParaView (File > Open) and colour
by `un`, `vn`, `wn`, or `pn`. The `.vtr` is a self-contained native VTK
rectilinear grid: it needs no reader choice and opens the same on Linux and
Windows, including files on a WSL drive. Convert a whole run at once and
ParaView loads the numbered series as a time animation:

```bash
python3 tools/field_to_vtr.py channel_field_*.h5
```

`field_to_vtr.py` reassembles the blocks onto the finest grid and never
modifies the input `.h5`; it needs the `vtk` Python module (or run it with
ParaView's `pvpython`). The solver also drops a `.xdmf` sidecar next to each
snapshot — that works with ParaView's legacy "XDMF Reader" on Linux, but the
`.vtr` is the reliable cross-platform route.

For a quick static look without ParaView, `tools/plot_field_slices.py` writes a
PNG of the three centre-plane cross-sections (x-y, x-z, y-z) for every variable:

```bash
python3 tools/plot_field_slices.py restart.h5      # -> restart.png
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

## Experiment: your own volume force

A different control strategy applies a steady (time-independent) volumetric
force inside the flow — e.g. the streamwise vortices of Schlatter & Canton
(doi:10.1007/s10494-016-9723-8). This has its own single edit point: the
subroutine `body_force` in `src/modules/volume_force.f90`, which returns the
force `(fx, fy, fz)` at a point `(x, y, z)`. The Schlatter & Canton form is
written out there ready to paste in.

Unlike the wall condition, the volume force is OFF by default. Turn it on by
adding to `input.ini`:

```ini
[force]
enabled = true
type    = steady
```

Then edit `body_force`, rebuild, and rerun. (Remember `beta` must be a multiple
of `2*pi/Lz` so the force is periodic in the spanwise direction — here
`Lz = 2`, so the fundamental is `beta = pi`.)

## Prescribing the start time

Set `t_start` under `[time]` in `input.ini` to choose the initial time of the
run; it overrides the time stored in `restart.h5`.
