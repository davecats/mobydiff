# A3 INCREMENT 0 — multi-level (3-level) refine_body prerequisite gates

Levels > 1 were lightly exercised before the airfoil physics (les_ibm used
levels 2 = l1 only): these gates run `mobygeom.py block-table --levels 3` +
`refine_body` with `refine_levels = 2` on the committed cylinder geometry
(`validation/cylinder/cylinder.stl`, D = 1 at (6.0, 8.02)) on a small
128x128x8 / lx = ly = 16 / lz = 0.25 grid, nb = 8 (D/h = 32 at level 2 —
the committed single-level cyl_re40 resolution, at 2.9x fewer cells).
Leaf table: 912 leaves (224 l0 + 192 l1 + 496 l2; 16 l2 leaves buried).
The solver's Fortran leaf builder cross-checks the file's blocks table
row-by-row at read, so every run below also gates Python==Fortran at 3
levels.

- `uniform.ini`   gate (a): uniform oblique flow (u, v, w all nonzero)
                  across every l0-l1 and l1-l2 interface/edge/corner,
                  preserved EXACTLY. Runs on the ZERO-FORCE TWIN
                  (`make_uniform_twin.py`): same touch-driven refinement,
                  buried masks zeroed + blocks table rebuilt without
                  burial removal (no closed faces in the flow), coef = 0
                  (no IBM force) — the body shapes the block layout but
                  cannot touch the constant field.
- `dwall.ini`     gate (b): per-level dwall/yeff on the REAL 3-level file
                  ([rans] dump_geometry). All x/y faces are declared
                  patches (inlet/outlet), so no domain wall is min'ed in;
                  `check_dwall_cylinder.py` checks the dump against an
                  independent exact distance to the as-built float32 STL
                  prism (2D ring distance + cap planes for inside-xy
                  points, cross-validated against a brute-force
                  point-to-all-triangles evaluation).

## Results (2026-07-13, all PASS)

- uniform: max|u-u0| = max|v-v0| = max|w-w0| = 0.0 and pn spread = 0.0
  after 50 steps, on CPU 1 rank, CPU 4 ranks AND GPU; 1 == 4 ranks
  byte-identical (compare_fields --tolerance 0).
- dwall: max|dwall-ref| = 3.6e-15 / 1.6e-15 / 1.6e-15 at levels 0/1/2
  (yeff identical); ransgeom dump 1 == 4 ranks identical.
- FOUND + FIXED while gating: `mobygeom.py` computed the dwall_blocks
  tiles with trimesh `proximity.on_surface`, which returns distances
  O(1e-6) too large near the quantized surface (measured 1.117e-6 at a
  level-2 cell 2.5e-3 from the cylinder; the flat T1 slabs were immune —
  axis-aligned planes). Now `igl.point_mesh_squared_distance` (exact
  AABB, matches a brute-force all-triangles reference to 1e-15).
  `validation/rans_geometry` flat gates re-run with the new generator:
  still exactly 0.0, and the regenerated les_ibm file still differs from
  the committed one only by the added dwall_blocks.

## Workflow

```bash
./setup.sh          # grid.h5 + 3-level block-table file + zero-force twin
./run_gates.sh      # uniform (1/4 ranks + GPU if built) and dwall gates
```
