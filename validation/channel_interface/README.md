# Turbulent-channel validation of the 2:1 interface treatment

Re_tau = 180 channel (u_tau = 1: one time unit = one eddy turnover).
Three cases:

> **Solver (defaults set in the inis).** The projection is
> Chebyshev-accelerated damped Jacobi (`accel = chebyshev`) at `niter = 12`
> (6 until 2026-09-26) with `sor = 0.8` -- `sor < 1` is REQUIRED, simple
> Jacobi DIVERGES for `sor > 0.8` (the old red-black SOR used 1.5; the
> red-black twin of `developed/` runs at 1.5 with `solver = redblack`). On
> the uniform reference (1000 steps) Chebyshev held ~3.4x lower divergence
> than plain Jacobi at equal niter for ~1.6 % more cost -- see
> `divergence_comparison.png`. `cflmax = 1.0` since the 2026-10-21 rerun
> (dtmax binds).
>
> **Note (2026-09-26): the reflux setup originally described here no longer
> exists.** The refined cases were first run with `[blocks] momentum_reflux =
> true` (Berger-Colella reflux of the interface advective momentum), and
> `run_validation.sh` had a `NOREFLUX=1` mode for the on/off comparison. That
> comparison (`developed/README.md`) showed the reflux to be the u'/v'
> coarse-cell band artifact, so the key was REMOVED (2026-07-01) together with
> `interface_constant_half` / `interface_skew`: the solver now always runs the
> validated const-1/2, reflux-off interface, and the `NOREFLUX` mode is gone.
> The uniform-grid conservation gates it cited (`validation/momentum_interface`,
> `tools/divsum.py` / `momsum.py`) were removed with the `MOBY_*` diagnostic
> hooks they needed; recover them from git history.

| case          | grid                          | interfaces            |
|---------------|-------------------------------|-----------------------|
| reference     | uniform 256 x 128 x 256       | none                  |
| refined_y110  | 128 x 64 x 128, wall bands of 24 base cells refined to level 1 | y+ = 112 |
| refined_y55   | 128 x 64 x 128, wall bands of 16 base cells refined to level 1 | y+ = 55  |

The y line is `natural` stretching (blend 16, dyw+ = 0.5 shape) at the
base resolution; the reference uses `[grid.*] subdivided = true`, which
builds each of its lines as the midpoint subdivision of the base line —
bitwise the refined cases' level-1 lines. The refined near-wall region
therefore has *exactly* the reference resolution (fine wall spacing
y+ = 0.77); deviations are attributable to the interface treatment and
to the coarser (2x) far field, not to near-wall resolution.

## Re-validated 2026-10-22: the first recorded numbers

RE-VALIDATED 2026-10-22, commit 532ef54 (banner), production flags, istmcetus
A6000, 1 rank, `run_validation.sh gpu 1` 6.8 h for the three cases (transient
5 + statistics 20 turnovers each; 0.09-0.12 s/step), `tools/channel_interface_validation.py`
against the uniform 256 x 128 x 256 reference (u_tau 0.99696). No result had
been recorded for these cases before (the README's solver box was the only
record). Settings changed for the rerun: `cflmax = 1.0` (dtmax 3.125e-4
binds), `sor = 0.8` explicit (the ICs come from the KMM180 red-black run),
`niter = 12` since 2026-09-26.

| case | u_tau (dev) | max mean-U deviation (of centreline) | divergence rms: interface / fine interior / coarse interior |
|---|---|---|---|
| refined_y110 (interface y+ 112) | 1.00439 (+0.75 %) | 1.22 % at y+ 98 | 2.0e-4 / 7.9e-4 / 3.4e-5 |
| refined_y55 (interface y+ 55) | 1.00309 (+0.61 %) | 1.58 % at y+ 50 | 6.6e-4 / 1.0e-3 / 4.8e-5 |

Reported as they come out, nothing tuned: the mean-velocity deviation
peaks just below the interface (the coarse side's first rows) at 1.2-1.6 %
of the centreline velocity, the wall shear sits within 0.8 % of the
reference, and the interface rows' divergence residual is BELOW the fine
interior's (the projection's residual on the finest cells), i.e. the 2:1
transfer adds no divergence of its own. The y+ 55 interface sits in the
buffer layer and costs more than the y+ 112 one, as expected. The profiles
and spectra are in `plots_y110/` and `plots_y55/` (untracked).

## Running (sized for a real GPU node, not a laptop)

```bash
./compile.sh gpu
cd validation/channel_interface
./run_validation.sh gpu 1          # or: cpu <nranks>
```

The script first runs the quick interface-decay gate, then generates
initial conditions by interpolating
`tutorials/channel_kmm180/channel_kmm180_restart.h5` (288x136x288) onto
the case grids (`tools/make_channel_restart.py`), and for each case runs

- a transient leg, t = 0..5 (discarded; the interpolated IC is not
  divergence-free and the flow needs to re-equilibrate), then
- a statistics leg, t = 5..25 (20 eddy turnovers), with channel stats
  sampled every 50 steps (per-level files: `channel_stats.h5` +
  `channel_stats_l1.h5`) and field snapshots every 4000 steps for
  spectra (~20 snapshots; ~0.5 GB each for the reference — adjust
  `field_interval` in the .ini files if disk is tight).

At dtmax = 3.125e-4 the statistics leg is ~64k steps per case
(reference: roughly 0.1-0.2 s/step on a data-centre GPU).

## Post-processing

```bash
python3 ../../tools/channel_interface_validation.py \
    --reference runs/reference/stats --refined runs/refined_y110/stats \
    --label "interface y+112" --out plots_y110
python3 ../../tools/channel_interface_validation.py \
    --reference runs/reference/stats --refined runs/refined_y55/stats \
    --label "interface y+55" --out plots_y55
```

Produces `profiles.png` (U+, rms, Reynolds stress, mean-U deviation,
interface height marked), `spectra.png` (streamwise/spanwise E_uu at
y+ ~ 15, just below the interface on the fine side, and just above it
on the coarse side), and prints u_tau deviation, the maximum mean-U
deviation, and divergence residuals (rms/max) in the interface band vs
the fine and coarse interiors.

Report the deviations as they come out — nothing here is tuned.
