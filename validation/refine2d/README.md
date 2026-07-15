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

## R2D-2 — projection y-type + physics gates (2026-07-15)

Solver change: face_grad/face_grad_denom/face_grad_corr gain a
per-direction `refined` flag ([blocks] refine_dims) — a 2:1 face normal
to an UNREFINED direction shares its node line on both sides (the
xz-quadtree y face), so its pressure gradient is the UNIFORM stencil
d1f in BOTH the Jacobi denominator and the velocity-face correction
(the SPD pair); refined directions keep the locked composite (2/3, 4/3)
branches, xyz mode is arithmetically identical. The y-face
reconciliation (sync + owned-face corrections + conservative restrict/
inject ghosts) was already direction-agnostic after R2D-1. The restart
refine_dims cross-check now applies only to BLOCK-layout files (legacy
global-3D restarts carry no layout to clash).

- Uniform gates (uniform_xz + unibody_xz) re-PASS exact.
- **Laminar order (Beltrami patch)**: the classic laminar channel is
  BLIND to xz-only refinement (the parabola has no x/z variation), so
  the decaying-Beltrami exact error takes its place (bp_xz_32/64.ini:
  central xz patch, dt halved with h). L2(vel) error 6.63e-3 (base 32)
  -> 9.30e-4 (base 64): order 2.83. Controls on the SAME setup: uniform
  grids 1.39e-5 / 3.47e-6 (the patch error is interface-dominated) and
  the VALIDATED xyz octree patch 9.20e-3 -> 1.23e-3 (order 2.90): the
  xz interface error is ~30% SMALLER than the octree's at equal base
  resolution (one direction conforms), converging at the same order.
  1 == 4 CPU ranks bit-exact (max_abs 0.0) with nonzero interface
  transfers.
- **RANS xz smoke** (rans_xz.ini: turb180 column, xz-refined wall
  bands, k/omega scalar halos across tangentially-2:1 y faces): 50
  steps finite and sane; 1 == 4 ranks EXACT (max_abs 0 incl. k/omega/
  nut); CPU vs GPU round-off class on default-FMA builds (u 2.5e-14,
  omega 9.1e-12 on its ~1e4 wall scale).
- **Turbulent xz wall-band channel** (PASS, corax RTX 5090, ~4 h):
  run_developed.py --refine-dims xz (IC_xz.h5 via make_channel_restart.py
  --refine-dims xz; 6656 leaves = 512 core l0 + 6144 band l1; wall bands
  x,z at 256-equivalent, y the SHARED stretched 64-line — the classic
  anisotropic wall-AMR layout), transient t 0..5 + developed stats
  t 5..25 (80k steps, niter = 6 Chebyshev throughout = the convergence-
  unchanged evidence), analyzed by analyze_xz_channel.py vs the
  committed campaign stats (figure xz_channel_stats.png):
  - NO interface band: u'/v'/w'/-u'v' jump ratios at both interface
    rows 0.995-1.039 (the y-refined campaign's reflux band was 1.56
    excess / 0.31 kink; its VALIDATED reflux-off case sits in the same
    <=1.04 class).
  - Core (0.7<y<1.3, resolution identical to the uniform128 control):
    fluctuation ratios u'/v'/w'/-u'v' = 0.976/0.972/0.956/0.978 —
    matching the VALIDATED xyz reflux-off signature (0.978/0.949/
    0.955/0.977, the known ~5% const-1/2 small-scale under-transmission;
    v' is actually CLOSER to 1 in xz) with mean-U shifted toward the
    better-resolved reference exactly like the validated case (bulk U
    15.90 vs uniform128 15.26, reference 15.68).

