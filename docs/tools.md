# Tools reference

The `tools/` directory holds the Python utilities that support `mobydiff`: geometry
preprocessing for the immersed boundary method, verification checks against exact solutions,
restart generation, and post-processing/plotting. They require Python 3 with `numpy`,
`h5py`, and `matplotlib`.

All of these read the solver's HDF5 output directly and understand both the legacy global-3D
layout and the block-table (refined) layout.

---

## Geometry / preprocessing

### `mobygeom.py`

Converts STL triangle meshes into IBM coefficient fields. Grid geometry is always imported
from a Fortran-generated `mobygrid` HDF5 file, so the preprocessor reuses the solver's exact
node lines, stretching, periodicity, and staggered coordinates.

```bash
python3 tools/mobygeom.py <subcommand> [options]
```

Main subcommands:

| Subcommand | Purpose |
|------------|---------|
| `stl-ibm-coeff` | Compute IBM coefficients from one or more STL meshes (the main path). |
| `block-active` | Write per-block solid-removal flags into a coefficient file. |
| `block-table` | Write block-table IBM coefficients for refined (AMR) runs. |
| `make-sphere-stl` / `make-bent-pipe-stl` / `make-lidded-bent-pipe-stl` | Generate STL test bodies. |
| `check-stl-geometry` | Run generic STL mesh / probe classification checks. |
| `test-stl-sphere`, `test-stl-bent-pipe`, … | Validate coefficients against analytic geometries. |
| `stress-stl-watertightness` | Damaged-STL robustness stress tests. |

Typical `stl-ibm-coeff` invocation (see the [`sailplane` tutorial](tutorials.md#sailplane-external-aerodynamics-with-ibm)):

```bash
python3 tools/mobygeom.py stl-ibm-coeff \
    --geometry body.stl --output ibm_coeff.h5 --grid-file grid.h5 --re 1.0e5 \
    [--scale S] [--translate X Y Z] [--tile-size NX NY NZ] [--jobs N] \
    [--check-fluid-points file.txt]
```

The full workflow — exporting the grid with `./build_cpu/mobygrid input.ini grid.h5` and the
complete option list — is documented in [`tools/README_mobygeom.md`](../tools/README_mobygeom.md).

---

## Verification / validation checks

Each of these compares a solver field against a known solution and reports L2/L∞ error. See
[Validation](validation.md) for how they fit into the verification workflow.

| Tool | Checks against | Invocation |
|------|----------------|------------|
| `check_beltrami.py` | Exact 3D Beltrami / ABC flow (2π-periodic cube) | `python3 tools/check_beltrami.py FIELD.h5` |
| `check_tgv.py` | Exact 2D Taylor–Green vortex | `python3 tools/check_tgv.py FIELD.h5 [--error-map out.png]` |
| `check_parabolic_channel.py` | Poiseuille profile of a forced channel | `python3 tools/check_parabolic_channel.py FIELD.h5 [--tolerance T]` |
| `compare_fields.py` | Another solver output (bit-exact refactor checks) | `python3 tools/compare_fields.py REF.h5 CAND.h5 [DATASETS...] [--tolerance T] [--export-global FILE]` |

`compare_fields.py` handles both output layouts, reassembling refined fields onto the finest
lattice; `--export-global` writes a single reassembled global field for visualization.
Datasets default to `un vn wn pn` (the three velocity components and pressure).

---

## Restart utilities

### `make_channel_restart.py`

Generates initial-condition / restart files for the 2:1-interface channel validation by
interpolating an existing channel restart onto a target grid (uniform reference or
wall-band-refined block-table):

```bash
python3 tools/make_channel_restart.py --mode {reference,refined,base,patch} \
    --out OUT.h5 [--source SRC.h5] [--band-cells 24] \
    [--refine-box x0 x1 y0 y1 z0 z1] [--dyw-plus 0.5] [--nx 128 --ny 64 --nz 128]
```

---

## Post-processing / plotting

These write a PNG and accept one or more runs as `FIELD.h5:LABEL` (the label is used in the
legend).

| Tool | Produces |
|------|----------|
| `channel_loglaw.py` | Mean streamwise velocity in wall units ($U^+$ vs $y^+$, semilog), both walls folded. |
| `channel_stats_profile.py` | Single-snapshot mean and rms fluctuation profiles ($x,z$-averaged per wall-normal row). |
| `plot_channel_stats.py` | Time-averaged channel statistics from a developed run's `stats_file` (combines all refinement levels; `_l1`, `_l2`, … are found automatically). |
| `slice_channel.py` | An $x$–$y$ mid-span cross-section of a wall-band-refined channel field ($u,v,w,p$), reassembled onto the finest grid. |

Example:

```bash
python3 tools/channel_loglaw.py loglaw.png run.h5:mycase reference.h5:ref
python3 tools/plot_channel_stats.py stats.png channel_stats.h5:mycase
```
