# Turbulent-channel validation of the 2:1 interface treatment

Re_tau = 180 channel (u_tau = 1: one time unit = one eddy turnover).
Three cases:

> **Solver + reflux (defaults set in the inis).** The projection is
> Chebyshev-accelerated damped Jacobi (`accel = chebyshev`) at `sor` (omega)
> `= 0.8`, `niter = 6`. `sor < 1` is REQUIRED -- simple Jacobi DIVERGES for
> `sor > 0.8` (the old red-black SOR used 1.5). At this niter the projection is
> under-converged for both accelerators (it plateaus on the large scales; a
> multi-level Schwarz is the planned fix), but on the uniform reference (1000
> steps) Chebyshev holds ~3.4x lower divergence than plain Jacobi at the same
> niter=6 for ~1.6% more cost -- see `divergence_comparison.png`. The refined
> cases run with `[blocks] momentum_reflux = true` -- the Berger-Colella reflux of
> the interface advective momentum (conserves the 2:1 interface momentum flux, the
> localized `-<u'v'>` / mean-shear defect; see `docs/interface_review.md` ii-iii
> and `validation/momentum_interface`, where it is gated to round-off
> conservation). The refine bands are full-extent planes (no corners), so the
> reflux conserves the interface momentum exactly there.
>
> Run `./run_validation.sh gpu <n>` (reflux ON by default). For a reflux-on vs
> reflux-off comparison add `NOREFLUX=1` -- those runs land in
> `runs/<name>_noreflux/`. (`tools/divsum.py` / `momsum.py` assume uniform cell
> volumes, so they are not meaningful on this stretched grid; conservation is
> gated on the uniform `momentum_interface` cases. Set `[pressure] accel = jacobi`
> to compare against plain damped Jacobi.)

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
