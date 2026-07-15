# 2D block refinement (refine_dims = xz) — R2D gates

Gates for the anisotropic 2:1 refinement phases R2D-0..3
(docs/next_session_refine2d.md): blocks refine in x and z only, the y
direction keeps ONE global (possibly stretched) node line at all levels.
`[blocks] refine_dims = xyz` (default octree) | `xz` (quadtree).

## R2D-0 — builder + lines + io (PASS 2026-07-15)

Infrastructure only: quadtree leaf table (2x1x2 children), mixed Morton
ids (CANONICAL xz form: y tile index in key bits 42+ above the 2D (x,z)
Morton curve of the finest-lattice coords, x even / z odd bits), shared
y node line in every level column, per-direction level scalings
(`2**(l*refMask(d))` everywhere), the `refine_dims` file attribute
(written only in xz mode, xyz files byte-identical; restart cross-checks
it), mobygeom `--refine-dims xz` mirror, compare_fields.py per-direction
reassembly. 2:1 interfaces ERROR OUT (comm.f90 resolve_neighbors) until
R2D-1.

- **Gate (a)** — 7-case nofma suite bit-exact (xyz untouched), CPU AND
  GPU: PASS (max_abs 0 on all datasets incl. nut/k/omega/gamma/rethetat).
- **Gate (b)** — `bash gate_leaftable.sh`: mobygeom block-table
  (`--refine-dims xz --refine-box`, out-of-domain sphere STL so body
  classification runs with all-fluid masks) vs the solver builder
  (`leaftable_test`, src/test_leaftable.f90) row-by-row on a non-cubic
  6x4x8 lattice, 3 levels: PASS — xz 576 leaves (levels 112/272/192)
  and the xyz control 1424 leaves both identical row-by-row.
- **Gate (c)** — all-refined xz box (allref_xz.ini: 32^3 nb=8, every
  block split once in x,z -> 256 level-1 leaves, no interfaces,
  Beltrami, 5 steps) vs the doubled-x,z single-level twin
  (twin_xz.ini: 64x32x64 with `[grid.x/z] subdivided = true`,
  bitwise the level-1 lines): un/vn/wn/pn max_abs = 0.0 on CPU AND GPU
  (cetus A6000); 1 == 4 CPU ranks max_abs = 0.0 in xz mode.
- **Barrier** — box_xz.ini (mixed-level xz table) error-stops at
  exchange setup with "[blocks] refine_dims = xz: 2:1 interfaces are
  not implemented yet (R2D-1)". Never runs.

Run gate (b): `PY=$HOME/ibmc/bin/python bash gate_leaftable.sh`
Run gate (c): solver on allref_xz.ini and twin_xz.ini, then
`compare_fields.py twin_xz_5.h5 ar_xz_5.h5 un vn wn pn --tolerance 0`.
