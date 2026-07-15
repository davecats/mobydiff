# Next session(s) — 2D block refinement (refine x,z only; y fixed)

User decision 2026-07-15: before the NACA AoA sweep in the new strategy,
build ANISOTROPIC 2:1 refinement — blocks refine in x and z only, the y
direction keeps ONE global (possibly stretched) line at all levels. The
unrefined direction is FIXED as y by convention:

- quasi-2D airfoil: span along y (the case re-orients, see below) — the
  span stops paying 2^L cells per level: the L5 airfoil drops from 65k
  to ~4-6k leaves, i.e. L5 "ground-truth" LE quality (the R1 fan result)
  at less than today's L4 cost;
- turbulent channel: y = wall-normal keeps the stretched global line
  while the WALL-PARALLEL resolution follows the distance from the wall
  — classic anisotropic AMR for wall turbulence.

Config: `[blocks] refine_dims = xyz` (default, today's octree — must
stay bit-exact through every phase) `| xz` (the new quadtree mode).

## What exists (inventory — read these first)

- Phases 0-3d (docs/block_refinement_strategy.md, CLAUDE.md): cubic
  nb^3 blocks, octree leaf table along the finest-lattice Morton curve,
  per-level midpoint-subdivided node lines, 26-neighbour 2:1 smoothing,
  exchange entries with per-dim affine gather maps (entry_gather_map),
  RESTRICT/PROLONG sampling on the source side, 4-sub-face enumeration
  per coarse face, const-1/2 velocity transfer (hardwired ON), the
  low-block-owns-face momentum stance (3c), face_grad composite
  projection stencils (SPD denominator/correction pairs), conservative
  copy reconciliation. ALL of it validated bit-exact/EXACT — the xyz
  path must not move.
- mobygeom block-table: per-level classification (bbox-windowed since
  79473ce — deep levels no longer allocate full-level rasters),
  chunked block-mask counts (750b205), coef/dwall tiles at each leaf's
  level, --keep-buried (cf68225), --refine-box COMBINED with body
  classification (9af9bcc). All per-direction machinery already exists
  (grid_axis_nodes/subdivided_args are per-direction).
- validation/multilevel_body (increment 0 gates), validation/naca0012
  (R1 fan tables = the reference metrics), the 7-case bit-exact suite
  (memory: bit-exact-gates-short).

## Design

### Leaf table and lattices (blocks.f90)

- `split_leaf` in xz mode: 2x2x1 children (x,z halve; the child keeps
  the parent's nb y-cells at the SAME global y-spacing). The per-level
  block lattice is nTiles_x*2^l x nTiles_y x nTiles_z*2^l — the y
  count is level-independent.
- Leaf ids along a mixed Morton curve: interleave (x,z) bits at the
  finest lattice, append the y tile index (or interleave with y at
  weight 1 — pick ONE canonical form, document it; both ends of every
  exchange derive entry order from it). lidOf/levelOffset sizing
  follows (4^l per level instead of 8^l).
- Node lines: lineX/lineZ midpoint-subdivided per level as today;
  lineY(:, l) = the single global line for every l. The grid_type
  subdivided-line identity (uniform fine line == subdivided coarse
  line, bitwise) continues to hold per refined direction.
- 2:1 smoothing stays the full 26-neighbour rule on the mixed lattice
  (y-offset neighbours included): a coarse block may not meet a
  grandchild ANYWHERE, or a y-face would carry a 4:1 tangential jump.
- Buried removal / refine_body: classify_block_geometry and the
  file-mask path work per level with the anisotropic spacing — the
  dilated-window logic is already per-direction.

### Interfaces (comm.f90, the core of the work)

Two topologies in xz mode:

1. **x/z-normal 2:1 faces** — like today but ONE refined tangential
   dim: a coarse face is fed by 2 fine sub-faces (parity in the single
   refined tangential direction; y is a plain copy dim). The
   entry_gather_map per-dim affine form expresses this directly (the y
   map becomes identity); sub-entry enumeration drops from 4/2/1
   (face/edge/corner) to 2/1 in-plane.
2. **y-normal faces between DIFFERENT levels — the genuinely new
   object.** The y-spacing is EQUAL on both sides (conforming normal
   direction: no deep-halo asymmetry, no fine-owns-face problem), but
   the tangential (x,z) resolution differs 2:1: a coarse block's y-face
   footprint is tiled by the y-faces of 4 fine blocks (2x2 in x,z).
   Transfers are tangential-only:
   - cell-centred p: RESTRICT = 4-sample tangential average; PROLONG =
     inject the covering coarse value (the const-1/2 stance carries
     over verbatim — injection, no tangential interpolation);
   - v (face-staggered in y, the CONFORMING dim): fine and coarse DOFs
     lie on the same y-plane; RESTRICT = 4-sample tangential average of
     matching-plane values, PROLONG = injection;
   - u/w (face-staggered in a refined dim): matching-face samples in
     their own staggered dim as today (1 matching sample), tangential
     average in the other.
   Momentum stance: both sides predict their y-face DOFs (spacing
   conforms — no storage obstruction), the final exchange reconciles
   conservatively to one owner (reuse the low-owns rule for
   determinism). The 3c/interface-normal asymmetry analysis does NOT
   apply here (it was a consequence of non-conforming normal spacing).
- Corner/edge extension: "combined neighbour occupied at ANY level"
  generalizes; the tangential halo extension rules re-derive per
  topology. The gather-map refactor (all ops as per-dim affine maps +
  weights) is the right substrate — no new kernel branches.

### Projection (pressure_solver.f90)

- x/z 2:1 faces: today's composite face_grad, unchanged in form.
- y 2:1 faces: the normal metric conforms, so the face gradient is the
  UNIFORM stencil; only the flux-area bookkeeping is 2:1 (4 fine fluxes
  sum against one coarse flux). The denominator must count the same
  sensitivity the correction applies (the SPD pair discipline — this is
  where instability hides; add the y-type branch strictly additively,
  do not touch the existing branches).
- Chebyshev bounds: Gershgorin lmax = 2 still holds (row sums of the
  weighted operator unchanged in form); lmin auto-derivation unchanged.

### io / tools

- blocks table: origin x,z in level-l cells, origin y in GLOBAL cells;
  level as today. The (nBlocksGlobal, nb^3) dataset layout, hyperslab
  writes, restart, named scalars: unchanged. xyz-mode files stay
  readable bit-identically; xz-mode is a new file variant (store
  refine_dims as an attribute; the solver cross-checks it).
- compare_fields.py: reassembly replicates coarse cells in x,z only
  when the attribute says xz (factor 2^(lmax-l) per refined dim, 1 in
  y); surviving-block masks likewise.
- mobygeom block-table: `--refine-dims xz` — subdivided_args refines
  x/z only; level_masks windows/counts, coef and dwall tiles work
  unchanged on the anisotropic level grids (verify the staggered-point
  helpers take per-direction lines — they do).
- slice_field.py / the fan metrics: already per-direction via the
  attrs; add the dims attribute handling.

### The y-span airfoil case (phase 3 of this work)

`[case.airfoil] span = z (default) | y`. With span = y: chord along x,
LIFT direction z; freestream (U cos a, 0, U sin a); inlets x_min,
z_min, z_max; outlet x_max; y periodic; drag/lift unit vectors rotate
to the x-z plane; make_airfoil_stl gains --span y (extrude along y).
Penalization-force reduction and the A2 exactness bookkeeping are
direction-agnostic. STL/mobygeom: extrusion past both y faces,
length-periodic in y.

## Phased plan (each gated before the next; the xyz default must be
## bit-exact vs HEAD at EVERY phase — nofma 7-case suite, CPU AND GPU)

- **R2D-0: builder + lines + io.** refine_dims config; quadtree leaf
  table + mixed Morton ids + per-level lattices; y-line sharing; blocks
  table attribute; xz-mode INTERFACES ERROR OUT (the 3a pattern).
  Gates: (a) suite bit-exact (xyz untouched); (b) xz box refinement
  leaf tables Python (mobygeom --refine-dims xz --refine-box) ==
  Fortran row-by-row; (c) an all-refined xz case (every block level 1)
  runs and is bit-exact vs the double-resolution-in-x,z single-level
  twin (the dyadic-line identity per refined dim, the 3a gate form).
- **R2D-1: exchange entries.** x/z faces (2 sub-entries) + the y
  tangential-2:1 type + edges/corners; copy-first entry ordering and
  the canonical wire order re-derived for the mixed Morton form.
  Gates: uniform oblique (u,v,w nonzero) flow through an xz-refined
  patch EXACT (0.0 incl. pn spread) across every face orientation,
  edge and corner, 1==4 ranks byte-identical, CPU AND GPU; global mass
  residual round-off; the multilevel_body zero-force-twin pattern at 3
  levels in xz mode.
- **R2D-2: projection + momentum stance.** face_grad y-type branch
  (conforming normal metric, 2:1 flux bookkeeping), SPD
  denominator/correction pair, conservative y-face reconciliation;
  const-1/2 tangential transfers.
  Gates: uniform still EXACT; laminar channel with an xz-refined band
  vs the uniform-fine reference (order ~2 away from interfaces, the 3c
  gate form); Re_tau 180 channel with wall-distance-adaptive xz
  resolution (coarse core, fine near-wall bands) vs the uniform
  reference: developed stats, NO interface bands (the
  channel_interface gate machinery + reflux-off lessons apply);
  Chebyshev convergence unchanged.
- **R2D-3: refine_body xz + the y-span airfoil.** mobygeom
  --refine-dims xz body classification; dwall tiles; keep-buried;
  [case.airfoil] span = y + STL --span y; fan bench at L5-xz (expect
  the R1a collapse: 0.0015-class near strips) and cost measurement vs
  L4-3D.
  THEN: the NACA AoA sweep (0/4/8 at L5-xz, keep-buried, forces) and
  the SD7003 rerun on the 2D-refined grid — the original goals, now at
  ground-truth LE quality.

## GPUs for test runs (status 2026-07-15; memory: remote-hosts,
## machine-nvhpc-module)

- **Local workstation**: RTX 3060 12 GB (cc86). `module load
  toolkits/nvhpc/25.9`; `./compile.sh cpu && ./compile.sh gpu`;
  nofma reference builds via the scratchpad build_nofma.sh pattern
  (CMAKE_Fortran_FLAGS=-Mnofma, OPENMP_OFFLOAD_FLAGS="-mp=gpu
  -gpu=nofma"; save main_ref copies BEFORE touching solver code).
  Suite ~35 min. Good for: the bit-exact suite, smokes, small gates.
- **istmcetus**: 2x RTX A6000 49 GB (cc86) — the LOCAL build_gpu
  binary runs directly (shared filesystem). SHARED MACHINE: check
  `nvidia-smi` first, pin `CUDA_VISIBLE_DEVICES` to the free GPU
  (GPU0 is frequently occupied by another user — NEVER touch it).
  module load works normally over `ssh istmcetus 'bash -lc "..."'`.
  Good for: big-memory cases (the 121 M-cell L6 fit only here).
- **istmcorax**: RTX 5090 32 GB (cc120) — needs its OWN build dir
  `build_gpu_corax` (no 25.9 modulefile: `export
  NV=/opt/Nvidia/nvhpc/Linux_x86_64/25.9; PATH=$NV/compilers/bin:
  $NV/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH`, cmake with
  `-DOPENMP_OFFLOAD_FLAGS="-mp=gpu -gpu=cc120"`). ~2.4x the 3060.
  Good for: the long physics runs (fan benches ~25 min, aoa cases ~2 h
  at L4).
- **Parallel etiquette and mechanics**: one job per GPU; fan
  independent runs out across all free GPUs (the naca sweep + SD7003
  ran on three hosts at once); launch detached (`setsid ... < /dev/null
  & disown`, logs on the shared FS); babysitter scripts (until-loop on
  a completion line, then ssh-launch) hand jobs between hosts without
  supervision. **REBUILD build_gpu_corax after EVERY solver change**
  (`ssh istmcorax '... cmake --build build_gpu_corax -j'`) — a stale
  remote binary silently ignores new ini keys and reproduces the old
  run bitwise (cost one bench cycle on 2026-07-15). Health checks with
  `ps aux | grep ... | grep -v grep` (never pgrep -f: it self-matches
  the invoking shell); watch forces files for FINITE garbage too, not
  just NaN (the L6 blow-up printed 1e115 finite values).

## Watch for

- The xyz default is sacred: every phase ends with the 7-case nofma
  suite bit-exact, CPU AND GPU. New branches must be unreachable when
  refine_dims = xyz.
- The mixed Morton form changes leaf ids for xz mode only — never for
  xyz files/restarts.
- The y-face SPD pair (denominator counts exactly what the correction
  applies) — the interface work's core lesson; add branches strictly
  additively.
- Smoothing must bound the level jump across y-faces too (26-neighbour
  rule on the mixed lattice), or a y-face sees 4:1 tangentially.
- exchange_scalar_halos consumers (RANS k/omega/gamma/ret, the phi
  loop, the band-filter dilation) inherit the new y-type entries
  automatically via the gather maps — gate one RANS xz case early
  (rans_inlet-style smoke) to catch scalar-path surprises.
- The band filter and keep-buried are per-level generic — re-run their
  gates once in xz mode (les_ibm xz-band + a keep-buried forces case).
- mobygeom: the leaf builder mirror (build_leaf_table_py) must stay
  rule-for-rule identical to blocks.f90 — the solver's row-by-row
  cross-check is the gate.

## NEXT-SESSION PROMPT

> Read docs/next_session_refine2d.md (the whole file) and CLAUDE.md.
> Branch claude/jacobi-interface; HEAD has the A3 phases + LE-fan
> resolution (refinement = reference answer, [ibm] band_filter =
> production option) committed. Build fresh nofma reference binaries
> (main_ref, CPU+GPU) at HEAD before touching code; rebuild
> build_gpu_corax after every solver change (stale-remote-binary
> landmine). Implement 2D block refinement ([blocks] refine_dims = xz,
> y never refined) in gated phases R2D-0..3 per the doc, one solver
> job at a time, parallelizing INDEPENDENT gate runs across the free
> GPUs (local 3060 / istmcetus A6000s CUDA_VISIBLE_DEVICES-pinned,
> GPU0 usually foreign / istmcorax 5090 via build_gpu_corax + PATH
> exports; nvidia-smi first, one job per GPU, setsid+shared-FS logs,
> ps-not-pgrep, watch forces for finite garbage not just NaN).
> R2D-0: quadtree builder + mixed Morton + shared y-line + io attr,
> interfaces error out; gates: suite bit-exact (xyz untouched),
> Python==Fortran xz leaf tables, all-refined-xz == doubled-xz-
> resolution twin bit-exact. R2D-1: exchange entries (x/z faces 2
> sub-entries; NEW y tangential-2:1 type, conforming normal spacing);
> gates: uniform oblique u,v,w flow EXACT through an xz patch at every
> orientation/edge/corner, 1==4 ranks, CPU==GPU, mass round-off.
> R2D-2: projection y-type face_grad (SPD pair!) + const-1/2
> tangential transfers + conservative y-face reconciliation; gates:
> uniform EXACT, laminar xz-band channel order check, Re_tau 180
> wall-distance-adaptive xz channel vs uniform reference (developed
> stats, no interface bands), a RANS xz smoke for the scalar halos.
> R2D-3: mobygeom --refine-dims xz + refine_body + dwall +
> keep-buried; [case.airfoil] span = y (+ make_airfoil_stl --span y);
> fan bench at L5-xz vs the R1 tables (validation/naca0012 README) +
> cost vs L4-3D. THEN the NACA AoA sweep (0/4/8, keep-buried forces)
> and the SD7003 rerun on 2D-refined grids. Deferred, do NOT start:
> TVD/van-Leer scalar upwind (measurement-justified, separate), ghost
> convection, smoothed-mask/Brinkman, convective outlet, GPU
> profiling. git: stage explicit paths only; commit per gated phase
> with results in the message; record findings in the validation
> READMEs as done throughout A3.
