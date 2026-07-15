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

## R2D-1 — exchange entries (PASS 2026-07-15)

The xz-quadtree interface transfers, entirely in the entry GENERATION
(comm.f90) — the exchange kernels are untouched (the per-dim affine
gather-map form covers the new types with no branches):

- resolve_neighbors: parent/child lattice coords via
  parent_coord/child_origin (per-direction), PROLONG parity tq masked to
  0 in unrefined dims, the child adjacency parity filter applied only to
  refined dims — an x/z-normal 2:1 face is fed by 2 fine sub-entries, a
  y-normal face by 4 (2x2 in-plane; the genuinely new type: the y line
  CONFORMS, so the normal direction is a plain copy row and only the
  tangential x,z resolution differs 2:1).
- interface_boxes / entry_gather_map: every UNREFINED direction of an
  interface entry uses the same-level copy form (full nb tangential
  range; the adjacent source row normally). RESTRICT dst quarter offsets
  ride the masked tq.
- entry_blend: the covering-cell index is per-direction; at a conforming
  y face aHalf == bHalf makes the pressure-ghost blend degenerate to 1
  (injection at the geometrically correct location) with no special case.
- iface_restrict_normal returns 0 for a conforming normal (the gather
  already reads the single face-adjacent row).

NOTE: the projection y-type face_grad (conforming normal metric, SPD
pair) and the conservative y-face reconciliation arrive in R2D-2 —
R2D-1 is gated on uniform flow ONLY (phi = 0 makes the interface metric
unused; non-uniform xz runs before R2D-2 are NOT validated).

- **Gate (a)** — 7-case nofma suite bit-exact vs R2D-0 becbda3, CPU AND
  GPU: PASS.
- **Uniform oblique flow** (uniform_xz.ini: 64x32x64, nb=8, central box
  refine_levels=2 -> 472 leaves [192, 248, 32], all face orientations /
  edges / corners at l0-l1 AND l1-l2; 32 fine-side cross-level y-face
  adjacencies + 80 x / 80 z): (u,v,w) = (0.9397, 0.3420, 0.2) preserved
  EXACTLY after 50 steps (max dev 0.0, pn spread 0.0) — mass residual
  exactly zero by the same token. 1 == 4 CPU ranks max_abs 0.0;
  CPU == GPU max_abs 0.0.
- **Zero-force body twin** (unibody_xz.ini, the multilevel_body pattern
  in xz mode): mobygeom block-table --refine-dims xz --levels 3 on the
  cylinder STL (508 leaves [220, 96, 192]) + make_uniform_twin.py (now
  refine_dims-aware); inlet/outlet BCs; uniform flow EXACT (0.0, pn
  spread 0.0) on 1 CPU rank, == 4 ranks and == GPU to max_abs 0.0. Also
  exercises the file-based refine_body xz read path (anisotropic mask
  shapes + solver builder cross-check against the xz blocks table).
