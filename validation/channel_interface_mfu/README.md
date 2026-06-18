# Minimal-flow-unit 2:1-interface channel validation

A fast variant of `../channel_interface`. Same Re_tau = 180 turbulent channel
and the same 2:1-interface treatment, but in a **minimal flow unit**: a small
2 x 2 x 2 box instead of the full 4π x 2 x 2π domain. ~17x fewer cells, so it
runs much faster and is convenient for quick stability / qualitative checks and
iteration.

| case          | grid (base)   | fine level   | interface |
|---------------|---------------|--------------|-----------|
| reference     | 48 x 128 x 80 (uniform) | —  | none      |
| refined_y110  | 24 x 64 x 40, both wall bands -> level 1 | 48 x 128 x 80 | y+ = 112 |
| refined_y55   | 24 x 64 x 40, both wall bands -> level 1 | 48 x 128 x 80 | y+ = 55  |

The wall-normal (`y`) line is **identical** to the full `channel_interface`
case (natural stretching, blend 16, dyw+ = 0.5), so the interfaces sit at the
same y+ and the near-wall resolution matches. The streamwise/spanwise counts
(24, 40) are the multiples of the block size nb = 8 closest to the full
channel's base spacing (dx = 0.083, dz = 0.050 vs the full 0.098 / 0.049). The
reference uses `[grid.*] subdivided = true`, so its grid is bitwise the refined
cases' level-1 lines.

Initial conditions are interpolated from the full channel restart
(`tutorials/channel_kmm180/channel_kmm180_restart.h5`) by sampling its 2 x 2
corner window in x and z (`tools/make_channel_restart.py --lx 2 --lz 2
--nx 24 --nz 40`). That snippet is not periodic over the small box, so the
transient leg lets the MFU re-equilibrate before statistics are sampled.

## Running

```bash
./compile.sh gpu
cd validation/channel_interface_mfu
./run_validation.sh gpu 1          # or: cpu <nranks>
```

Each case runs a transient leg (t = 0..5, discarded) and a statistics leg
(t = 5..25). At dtmax this is far cheaper than the full case.

**Caveat on statistics:** the MFU has ~40x less homogeneous plane area than the
full domain, so converged turbulence statistics need correspondingly longer
time-averaging. Treat this case as a fast stability and sanity check (does the
interface stay stable, do the mean profile and near-wall structure look right),
not as a source of fully converged Reynolds-stress profiles. Use the full
`channel_interface` case for converged statistics.

## Post-processing

```bash
# Profiles / spectra / divergence (the two channel halves are kept separate):
python3 ../../tools/channel_interface_validation.py \
    --reference runs/reference/stats --refined runs/refined_y110/stats \
    --label "MFU y+112" --out plots_y110

# x-y cross-section of u, v, w, p from any field snapshot:
python3 ../../tools/plot_field_section.py \
    runs/refined_y110/stats/channel_field_XXXX.h5 --z 1.0
```
