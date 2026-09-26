# Tools reference

The `tools/` directory holds the utilities that support `mobydiff`: geometry helpers,
verification checks against exact solutions, field comparison, restart generation, and
post-processing/plotting. The Python tools need Python 3 with `numpy`, `h5py`, and
`matplotlib`; each script's docstring is its full usage reference.

The field readers understand the solver's block-table snapshot layout (and the legacy
global-3D layout of old files).

> The preprocessor for immersed-body cases is the Fortran executable `moby_prepare`, not a
> tool in this directory — see [Running](running.md#preprocessing-immersed-bodies-moby_prepare).

---

## Geometry

| Tool | Purpose |
|------|---------|
| `make_airfoil_stl.py` | Extruded cylinder / NACA 4-digit / Selig-file airfoil STLs for quasi-2D cases (`--span y` for the `refine_dims = xz` orientation). |
| `make_geometry_stl.py` | ASCII STLs for the conjugate-heat-transfer gates: tilted plane, cylinder, annular pipe shell. |
| `check_annulus.py` | Checks an annular pipe body in a prepared case file against the analytic geometry. |
| `mobygeom.py` | **Retired** Python STL-to-IBM preprocessor, kept only as the independent cross-implementation reference for the `validation/prepare/` gates (plus its STL generators and watertightness tests). See [`tools/README_mobygeom.md`](../tools/README_mobygeom.md); it reads the grid from a `moby_prepare` case file (`--grid-file case.h5`). |

---

## Verification / comparison

| Tool | Checks | Invocation |
|------|--------|------------|
| `check_beltrami.py` | Error vs the exact 3D Beltrami / ABC flow (`[flow] initial = beltrami`) | `python3 tools/check_beltrami.py FIELD.h5` |
| `check_parabolic_channel.py` | Poiseuille profile of a forced laminar channel | `python3 tools/check_parabolic_channel.py FIELD.h5 [--tolerance T]` |
| `check_interface_decay.py` | The `validation/interface_decay/` gate: noise on a refined patch must decay | see the case README |
| `compare_fields.py` | Two solver outputs (bit-exact refactor checks) | `python3 tools/compare_fields.py REF.h5 CAND.h5 [DATASETS...] [--tolerance T] [--export-global FILE]` |
| `h5maxdiff.c` | Max \|a−b\| per dataset, in C (for machines without `h5py`/`h5diff`); build line in its header | `h5maxdiff A.h5 B.h5 [DATASET...]` |

`compare_fields.py` reassembles block-table fields onto the finest lattice (which can
exhaust memory at deep refinement); `--export-global` writes the reassembled global field for
visualization. With no dataset arguments it discovers the datasets present. `h5maxdiff`'s
default list is `un vn wn pn`, `nut` and the RANS scalars — name passive scalars (e.g. `s1`)
explicitly.

Interface diagnostics used by the 2:1-interface validation: `interface_diagnostics.py`,
`interface_coarse_gate.py` (Beltrami slab fields), `channel_band_profile.py`,
`channel_interface_validation.py`, `patch_interface_diff.py`, `patch_interface_stats.py`
(channel wall bands and embedded patches).

---

## Restart utilities

### `make_channel_restart.py`

Generates initial-condition / restart files for the 2:1-interface channel validation by
interpolating an existing channel restart (`--source`) onto a target grid (uniform
reference, wall-band-refined, base, or embedded patch; `--refine-dims xz` supported):

```bash
python3 tools/make_channel_restart.py --mode {reference,refined,base,patch} \
    --source SRC.h5 --out OUT.h5 [--band-cells 24] \
    [--refine-box x0 x1 y0 y1 z0 z1] [--dyw-plus 0.5] [--nx 128 --ny 64 --nz 128]
```

---

## Post-processing / plotting

The channel plotters write a PNG and accept one or more runs as `FIELD.h5:LABEL` (the label
is used in the legend).

| Tool | Produces |
|------|----------|
| `channel_loglaw.py` | Mean streamwise velocity in wall units ($U^+$ vs $y^+$, semilog), both walls folded. |
| `channel_stats_profile.py` | Single-snapshot mean and rms fluctuation profiles ($x,z$-averaged per wall-normal row). |
| `plot_channel_stats.py` | Time-averaged channel statistics from a run's `stats_file` (combines all refinement levels; `_l1`, `_l2`, … are found automatically). |
| `slice_channel.py` | An $x$–$y$ mid-span cross-section of a wall-band-refined channel field ($u,v,w,p$), reassembled onto the finest grid. |
| `plot_patch_slice.py` | Cross-section of a refined-patch field vs a base-grid control, patch outline overlaid. |
| `plot_beltrami_fields.py`, `plot_beltrami_slices.py` | Beltrami solver/exact/error slices, and single-block vs refined cross-sections over time. |
| `scalar_stats.py` | Reads the passive-scalar statistics files (`[scalar] stats_file`). |

Example:

```bash
python3 tools/channel_loglaw.py loglaw.png run.h5:mycase reference.h5:ref
python3 tools/plot_channel_stats.py stats.png channel_stats.h5:mycase
```

---

## Performance

| Tool | Purpose |
|------|---------|
| `moby_tune.sh` | Finds a machine's best rank-to-GPU mapping (`MOBY_GPU_ORDER`) by measurement: `tools/moby_tune.sh <moby_solve> <case.ini> [outdir]`. |
| `partition_analysis.py` | Offline analysis of how a block-to-rank partition cuts the 2:1 interface. |
