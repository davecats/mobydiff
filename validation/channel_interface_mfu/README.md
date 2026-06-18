# Minimal-flow-unit 2:1-interface channel validation

A faster variant of `../channel_interface`. Same Re_tau = 180 turbulent channel
and the same 2:1-interface treatment, but in a **minimal flow unit**: a small
2 x 2 x 2 box instead of the full 4π x 2 x 2π domain. Roughly half the cells of
the full case (and a much smaller homogeneous footprint), so it runs faster and
is convenient for quick stability / qualitative checks and iteration.

| case          | grid           | resolution vs reference | interface |
|---------------|----------------|-------------------------|-----------|
| reference     | 48 x 128 x 80 (uniform) | — (baseline) | none |
| refined_y110  | 48 x 128 x 80 centre, both wall bands -> level 1 (96 x 256 x 160) | = in centre, 2x finer at walls | y+ = 112 |
| refined_y55   | 48 x 128 x 80 centre, both wall bands -> level 1 (96 x 256 x 160) | = in centre, 2x finer at walls | y+ = 55  |

The block-refined case has the **same resolution as the reference in the
centre** and is **2x finer at the walls** -- so it is **never coarser than the
reference**. The coarse side of each 2:1 interface therefore carries the same
spectral content as the reference, so an interface artifact is not confused
with under-resolution aliasing (which would appear if the refined coarse zone
truncated more of the spectrum than the reference).

Grids are anchored on the full `channel_interface` case's native 24x64x40
lines: the reference and the refined **centre** are both the midpoint
subdivision of those lines (48x128x80, `[grid.*] subdivided = true`), bitwise
identical; the refined **walls** subdivide once more (96x256x160). The y-line
matches the full case's fine level, so the interfaces sit at the same y+.

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
