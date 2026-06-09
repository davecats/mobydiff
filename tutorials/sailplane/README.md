# Sailplane Generic-Flow Tutorial

The source STL is in millimetres and is symmetric about its original `y = 0`
plane. The tutorial uses metres, keeps the symmetry plane at computational
`y = 0`, and simulates only the `y >= 0` half-domain.

The original STL bounding box is:

```text
min  = (3.226455e-10, -9.744634, -0.250183) m
max  = (8.644325,      9.744634,  1.700000) m
size = (8.644325,     19.489268,  1.950183) m
```

The computational domain in `input.ini` is:

```text
Lx = 5.0 * length = 43.22162499838677
Ly = 2.5 * span   = 48.72317
Lz = 5.0 * height = 9.750915
```

The coefficient generation command scales the STL by `0.001` and translates it
by `(17.288649999032064, 0.0, 4.150549)`. This places the full mirrored STL at
the centre of the corresponding full domain, while the solver only uses the
positive-`y` half.

Generate the IBM coefficients with:

```bash
cd /home/davide/Codes/FDM/tutorials/sailplane
../../build_cpu/mobygrid input.ini sailplane_grid.h5
/home/davide/ibmc/bin/python ../../tools/mobygeom.py stl-ibm-coeff \
  --geometry "FRUE V0 ohneRundung.stl" \
  --output sailplane_ibm_coeff.h5 \
  --grid-file sailplane_grid.h5 \
  --re 1.0e5 \
  --scale 0.001 \
  --translate 17.288649999032064 0.0 4.150549 \
  --check-fluid-points fluid_points.txt \
  --jobs 2 \
  --tile-size 32 36 16
```

If the grid changes in `input.ini`, rerun `mobygrid` and regenerate
`sailplane_ibm_coeff.h5`. The coefficient generator reads `nx`, `ny`, `nz`,
`lx`, `ly`, `lz`, periodicity, and the exact solver node coordinates from
`sailplane_grid.h5`.

`mobygeom.py` uses padded bounding-box culling for external STL bodies by
default, so only the grid points near the transformed sailplane are classified.
The fluid probe file is used as a fail-fast validation check. Use
`--fluid-points fluid_points.txt` only for an additional expensive ambiguous
winding repair pass, and add `--fluid-ray-scope all` only when deliberately
debugging the older full ray-vote path. Use `--no-bbox-cull` only for debugging
against the older full-grid path.

The tutorial is currently sized as a large one-step smoke case for the 6 GB
Quadro RTX 3000 in this workstation. It is intentionally kept a little below
the estimated device-memory ceiling so both the solver and the Python STL
preprocessing remain stable under WSL.
