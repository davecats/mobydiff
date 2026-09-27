# Turbulent channel, Re_tau 180 (Kim, Moin & Moser 1987)

`[case] name = channel`, 288 x 136 x 288 on 4pi x 2 x 2pi, natural-stretched
y, constant mean pressure gradient (`forcing_x = 1`, friction units), cold
start from the channel initialiser's mean profile + disturbance, red-black SOR
projection. Statistics go to `channel_kmm180_stats.h5`
(`tools/check_parabolic_channel.py`, `tools/channel_stats_profile.py`,
`tools/channel_loglaw.py`), the per-step runtime line to `runtimedata.txt`.

```bash
mpirun -n 4 ../../build_gpu/moby_solve input.ini
```

## Wall spacing and the time step (decided 2026-09-27, numerics review F6)

The ini carried `natural_dyw_plus = 0.05` -- a first cell 0.05 wall units
off the wall, KMM's own first collocation point. For a second-order
finite-difference scheme with explicit diffusion that buys nothing and costs
the time step: the Peclet limiter `pecletmax / ((1/Re) / dy_wall^2)` caps
`dt` at 7.7e-6 against the ini's nominal 3.125e-4. **The archived run
confirmed it** -- the restart written on 2026-08-03
(`mobydiff.scalar/tutorials/channel_kmm180/channel_kmm180_restart.h5`,
metadata) carries `dt = 7.696034892725926e-06`, `cfl = [0.0037, 0.5]` (the
Peclet limit binding, the Courant number 0.4 %) and `step = 650000` at
`t = 5.0`: 650 k steps for 5 h/u_tau, i.e. 26 M steps for the ini's
`t_final = 200`. An interpolated twin at dy_wall+ = 0.5 sits beside it
(`channel_kmm180_restart_dywplus0p5.h5`, the same field on the coarser
line).

The ini now sets `natural_dyw_plus = 0.5`, a wall-resolved first cell
(dy_wall = 2.77e-3, ny still 136; the stretched line's centre spacing is
unchanged at 0.021). Limits at this spacing: Peclet rate 726 -> `dt <=
6.9e-4`; Courant rate ~ 20/0.0436 -> `dt <= 1.7e-3`; so `dtmax = 3.125e-4`
binds and the run steps 40.6x faster. Verified 2026-09-27 with 100-step GPU
runs of this ini (`dt` and `cfl` attributes of the snapshot):

```
natural_dyw_plus = 0.5   dt = 3.125000e-04  cfl = 0.1306  Peclet = 0.2268  dy_wall+ = 0.498
natural_dyw_plus = 0.05  dt = 7.696035e-06  cfl = 0.0032  Peclet = 0.5000  dy_wall+ = 0.053
```

-- the 0.05 line reproduces the archived run's `dt` to every digit, and the
0.5 line is dtmax-bound with both limiters well below their caps.

Any KMM statistic recorded before this date was measured on the 0.05 line
at `dt = 7.7e-6`; the case is on the one-by-one re-validation list
(`validation/README.md`).
