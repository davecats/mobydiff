# Block-Based Local Refinement and Solid-Block Removal Strategy

Design notes for extending mobydiff with (1) local mesh refinement following the
Building-Cube Method (Nakahashi & Kim, AIAA 2004-434) and its HPC formulation in
CUBE (Jansson et al., IJHPCA 2019), and (2) removal of decomposition blocks that
lie entirely inside an immersed boundary. No code is changed by this document;
it records the agreed strategy and a phased, verifiable implementation plan.

Decisions taken up front:

- **Static refinement**: the block layout is decided at mesh-generation time
  (geometry-adaptive, in `mobygrid`) and is immutable during a run. Dynamic
  adaptation and load re-balancing are designed for, but deferred (Section 12).
- **The mesh stays fully Cartesian**: blocks are index-space boxes over global
  per-direction node lines. All numerics inside a block are exactly the current
  second-order staggered stencils; nothing changes in the per-point formulas.

## 1. Core idea

Replace "one structured box per MPI rank" with "many equal-sized structured
blocks per rank". Every block has the same number of cells `nb^3` (input
parameter, e.g. 16 or 32, must be even for red-black). Local resolution is set
by the *block size in physical space*, i.e. by its refinement level `l`:
a level-`l` block covers `nb` cells of the level-`l` global grid, where the
level-`l` grid is obtained from the base grid by `l` rounds of cell bisection.
Neighbouring blocks differ by at most one level (2:1 rule, as in both papers).

This is the BCM trade: a small, fixed amount of geometric metadata per block,
in exchange for completely regular, branch-free inner loops — which is exactly
what the GPU wants. Equal `nb` makes every block an identical unit of work,
storage, communication and (later) load balancing.

## 2. The block set: one flat object, GPU-resident

Following Jansson et al. (Section 2.3: no tree at solve time, only linear
arrays with neighbour adjacency), all per-block state lives in flat arrays with
the block as the **last** index. One derived type owns everything, mapped once
with `target enter data` like the existing `field_type`/`grid_type`:

```fortran
integer(C_INT), parameter :: FACE_SAME=1, FACE_COARSE=2, FACE_FINE=3, &
                             FACE_PHYS=4, FACE_CLOSED=5

type :: block_set_type
    integer(C_INT) :: nb = 16          ! cells per direction, all blocks
    integer(C_INT) :: nBlocks = 0      ! blocks owned by this rank
    integer(C_INT) :: nLevels = 1

    ! --- metadata (small, host + device mirror) ---
    integer(C_INT), allocatable :: globalId(:)      ! (nBlocks), Z-order id
    integer(C_INT), allocatable :: level(:)         ! (nBlocks)
    integer(C_INT), allocatable :: origin(:,:)      ! (3,nBlocks) cell origin in
                                                    ! level-l global index space
    integer(C_INT), allocatable :: nbrKind(:,:)     ! (26,nBlocks) FACE_* above
    integer(C_INT), allocatable :: nbrRank(:,:)     ! (26,nBlocks)
    integer(C_INT), allocatable :: nbrSlot(:,:)     ! (26,nBlocks) local slot or
                                                    ! remote global id
    ! --- per-block 1D metrics, sliced from per-level node lines ---
    real(C_DOUBLE), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)      ! (-1:nb+2,NVAR,nBlocks)
    real(C_DOUBLE), allocatable :: d1x(:,:,:), d1y(:,:,:), d1z(:,:,:)! (0:nb+1,NVAR,nBlocks)
    real(C_DOUBLE), allocatable :: lapXm(:,:,:), ...                 ! (0:nb+1,NVAR,nBlocks)

    ! --- fields ---
    real(C_DOUBLE), allocatable :: q(:,:,:,:,:)      ! (0:nb+1,0:nb+1,0:nb+1,NVAR,nBlocks)
    real(C_DOUBLE), allocatable :: qs(:,:,:,:,:)     ! (...,NVEL,nBlocks)
    real(C_DOUBLE), allocatable :: oldrhs(:,:,:,:,:) ! (1:nb,...,NVEL,nBlocks)
end type block_set_type
```

Notes:

- `ibm%coef`/`ibm%mu` and `les%nut` get the same trailing block dimension.
- A face whose neighbour is refined (`FACE_FINE`) has 4 fine neighbours; store
  them in a small side table indexed from `nbrSlot` rather than widening the
  26-entry table.
- No allocatable-inside-array-of-derived-types: OpenMP offload maps each flat
  array once; kernels receive plain contiguous arrays. This is the single most
  important GPU-clarity decision.

Kernels gain one outer loop and keep their bodies verbatim:

```fortran
!$omp target teams distribute parallel do collapse(4)
do b = 1, blk%nBlocks
  do k = 1, nb
    do j = 1, nb
      do i = 1, nb
        ! current momentum / sweep body, with g%d1x(i,var) -> blk%d1x(i,var,b)
```

Because all blocks are identical in shape, `collapse(4)` produces one large,
perfectly regular iteration space: occupancy and coalescing are independent of
how many blocks a rank owns, and removed blocks (Section 6) simply do not
appear — zero masking in inner loops, zero wasted memory.

## 3. Per-level grid lines: stretched grids survive refinement

`init_grid_direction` is reused unchanged. The only new ingredient is the
family of global node lines:

- Level 0: the current `xNode/yNode/zNode` built from the configured
  distribution (uniform/cosine/tanh/natural).
- Level l+1: midpoint bisection of every level-l cell.

Midpoint subdivision guarantees that a level-l cell is exactly the union of its
2^3 children — required for conservative interface transfer — while preserving
arbitrary stretching of the base line. Each block slices its `(-1:nb+2)`
coordinate window from its level's line (using the existing `face_at` /
`cell_center_at` halo extension at physical boundaries) and computes `d1`/`lap`
metrics exactly as today. Inside a block nothing is new.

## 4. Block-layout generation (extends `mobygrid`, serial, offline)

1. **Root tiling**: cover the domain with level-0 blocks, `globalSize/nb` per
   direction (the configured grid is the level-0 line; `nb` must divide
   `globalSize`).
2. **Geometry refinement**: recursively split blocks whose dilated bounding box
   intersects the immersed surface (classification via `mobygeom` / `isInBody`)
   until the finest level is reached. This is the octree pass of Nakahashi &
   Kim — the tree exists only here, never in the solver.
3. **2:1 smoothing**: split any block with a neighbour more than one level
   finer, iterate to fixed point.
4. **Finest-level buffer at the wall**: enforce that every block whose dilated
   box intersects the body surface is at the finest level, plus ≥1 block of
   buffer. Consequence: *level interfaces never touch the IBM region* (the BCM
   paper makes the same recommendation since coarse-fine transfer is the least
   accurate operation). The penalization, `set_ibm_coeff`, and the wall physics
   only ever see uniform finest-level resolution.
5. **Solid-block removal**: drop every block whose *dilated* region (all cells
   plus one halo layer, sampled at cell centres and at the three staggered face
   locations) lies inside the body. Dilation guarantees no surviving fluid
   stencil ever reads a removed block. Mark the faces of surviving neighbours
   as `FACE_CLOSED`.
6. **Ordering and output**: number surviving blocks along a Z-order
   space-filling curve (Jansson Section 6.1) and write the block table
   (id, level, origin, neighbour table, per-face kinds) plus per-level node
   lines to the grid HDF5 file. Rank assignment at solver start is the linear
   distribution `n = floor((N + P - p - 1)/P)` over Z-order — contiguous in
   space, computable without communication, and restart works on any rank
   count.

As implemented (Phase 3d), the offline tool is `mobygeom.py block-table`:
it classifies per-level `block_touch_l{l}`/`block_buried_l{l}` rasters
(dilated block windows sampled at cell centres and the three staggered
locations, on midpoint-subdivided node lines), mirrors the solver's leaf
builder in Python to enumerate the surviving leaves, and writes the
`blocks` table plus per-leaf IBM coefficient tiles (`coef_blocks`,
ghost-inclusive, evaluated at each leaf's own level) into the
coefficient file. The solver verifies its independently built leaf table
row-by-row against the file's and reads its slot range directly; with
`refine_body` it reads the per-level masks instead of classifying
analytically. Legacy single-level global-grid coefficient files remain
readable for unrefined runs.

## 5. Halo exchange: from rank-boxes to block-pair lists

`comm.f90` already has the right skeleton: precomputed send/recv boxes, flat
packed buffers, one flat GPU loop over all points for pack/unpack. The
generalization is to make the *exchange entry*, not the rank, the unit:

- At init, build a list of exchange entries, one per (block face/edge/corner,
  neighbour) pair, each carrying: source slot, destination slot + box,
  destination rank, and a per-dimension *gather map* precomputed from the
  node lines (`entry_gather_map`). Every transfer is the same weighted
  gather: for destination index i along a dimension, the source rows are
  `base..base+cnt-1` with `base = ishft(ga*i + gb, -gs)`, and variables
  staggered along that dimension sample the single matching face row.
  Same-level copies (ga=1, cnt=1), restrictions (the two covering fine
  rows: ga=2 tangentially, the row pair across the face) and
  prolongations (the covering coarse row: a halving shift tangentially,
  the tq-dependent row across the face) are all instances — there is no
  per-operation branch in the kernels, and the pressure ghost blend of
  Section 6 is just a destination-completion weight pair (wp, wpDst)
  applied on the receiving side.
- Entries are ordered **same-level copies first** (per peer), with prefix
  point counts, so the copy-only exchange the projection uses per colour
  (`interp=.false.`) is literally a prefix of the full one: shorter
  messages and shorter flat loops, no runtime filtering.
- **Local entries** (neighbour on same rank, the common case after Z-order
  distribution) are executed as a single device kernel — a flat loop over the
  concatenated entry points, exactly the current `pack_q_boxes` indexing
  pattern, but writing directly into the destination halo. No MPI, no host.
- **Remote entries** are grouped per neighbour rank into one message; pack and
  unpack kernels are the current ones with the entry-offset table replacing the
  26-neighbour offset table. The gather is applied *inside* pack so the wire
  format stays a flat real array (destination-point values).
- The nonblocking `start/finish_halo_exchange` split survives unchanged, and
  the internal/external block zoning of Jansson (overlap local exchange and
  interior kernels with MPI) composes naturally with
  `docs/nonblocking_overlap_strategy.md` later.

Halo width stays 1 (second-order stencils, including the diagonal terms in the
momentum cross-fluxes — hence the full 26-direction adjacency, as today).

## 6. 2:1 level interfaces

Operations (Jansson Fig. 3, adapted to the staggered grid):

- **Fine → coarse (`RESTRICT`)**: coarse halo value = average of the 2^d fine
  values it covers. For the face-normal staggered velocity on the shared face
  this is the *conservative* choice: with midpoint subdivision the four fine
  sub-faces tile the coarse face exactly, so the coarse face flux equals the
  sum of fine fluxes.
- **Coarse → fine (`PROLONG`)**: v1 = injection (copy of the covering coarse
  value), as in BCM; v2 = trilinear interpolation for second-order interface
  accuracy. Make the operator a switch so accuracy can be assessed by diff.

Interface ownership rule, to keep the projection conservative and simple.
The original draft of this section said "the fine side owns the shared
face" unconditionally; the implementation (Phase 3c) revises this to
**the low-side block owns the shared face**, for the following storage
reason. A block owns the staggered faces it computes, `u(1..nb)`; its
`u(nb+1)` is halo. A fine block EAST of an interface holds the fine face
DOFs as its `u(1)` — interior storage, momentum-predictable (the stencil
reads the prolonged halo) and sweep-correctable: there fine-owns-face
works exactly as drafted. But when the fine blocks sit WEST, the
fine-level face DOFs are their `u(nb+1)` halos: fine momentum cannot
predict them (the stencil would need `q(nb+2)`, beyond the one-cell
halo), and restricting them back would write the coarse neighbour's
INTERIOR `u(1)` plane with variable-discriminating entries while the
prolong of the same DOFs reads it — an intra-kernel cycle.

Ownership defines the *reconciled* value: after the projection, the full
exchange overwrites the non-owner's halo copy with the owner's face
(RESTRICT 4-sub-face average when the high side is coarse, PROLONG
injection when the high side is fine). Conservation is exact in both
orientations: the restricted coarse flux is the equal-area average of
the fine sub-fluxes (midpoint subdivision halves cells exactly), and the
injected fine sub-fluxes sum to the coarse flux by construction.

**Relaxation at interfaces (revised after Phase 3d, the hard-won part).**
The original draft relaxed interface faces one-sidedly: the owner
corrected the face, the other side masked it (out of `denom` and the
corrections, like a wall) and only consumed the exchanged value in its
divergence. Both this and the variant that keeps the full denominator
are *unconditionally unstable*: the non-owner cannot push back on a
live flux that keeps moving, the resulting block iteration is a
non-contractive splitting, and a pressure-jump mode localized at the
interface grows exponentially (per-step gain independent of dt,
viscosity and `sor`; seeded by round-off, so short or symmetric gates
never see it). The cure is the BCM regime (Nakahashi), made conservative:

- **Symmetric corrections, own copies.** Every non-pinned face of a cell
  is in the denominator and is corrected by that cell's relaxation -
  including the block's own copy of a 2:1 interface face (`u(1)` on one
  side, the `u(nb+1)` halo on the other). Only walls and closed faces
  (`FACE_PHYS`/`FACE_CLOSED`) are pinned.
- **Stage-frozen ghosts.** The per-colour exchanges inside the
  projection run only the same-level copy prefix of the entry lists
  (`interp=.false.`, Section 5), so each side relaxes against
  stage-frozen ghosts and its own evolving copy of the interface faces.
  This is a regular damped block iteration; the one-sided
  overwrite-every-colour variant is what diverges.
- **Conservative reconciliation.** The final full exchange of the
  projection snaps the non-owner copies back to the owner's values, so
  the state entering the next stage is owner-consistent and mass
  telescopes exactly.
- **Consistent pressure ghosts.** Plain injection puts the coarse cell
  value at the fine halo centre, so the momentum pressure gradient at
  the interface is evaluated over a gap 1.5x larger than the fine metric
  assumes - a systematic overdriving of the interface velocity. PROLONG
  *face* entries for `p` therefore blend: ghost = wp*coarse +
  wpDst*first-interior with the weights from the node lines (uniform
  2:1: ghost = (2 p_C + p_f)/3). In the exchange this is the
  destination-completion weight pair of the gather (Section 5), applied
  on the receiving side. Edges and corners stay plain injection; only
  face halos enter the pressure gradient.
- **Covering-cell prolong sources.** A PROLONG entry's source row in an
  off-dimension is the *covering* coarse cell of the halo layer, which
  depends on the fine block's parity tq (an edge/corner's coarse
  neighbour spans past the shared boundary). The off-based row (1 or
  `nb`) is correct only for faces, where 2:1 smoothing forces tq to the
  touching side.

Accuracy: a coarse-owned face carries a uniform (coarse-resolution)
flux across its four sub-faces until the trilinear prolong upgrade;
fine-owned faces carry full fine resolution. With the finest-level
wall buffer (Section 4) interfaces sit in smooth flow, where this is a
second-order-consistent approximation.

A manufactured-field halo audit (`MOBY_HALO_AUDIT=1`, see `main.f90`)
checks every exchange-written halo cell against the design semantics
above on the actual block layout; transfers should be verified with it
before debugging the physics.

Red-black parity: with `nb` even, define each block's `colorOffset` from the
parity of its global level-l origin (`modulo(sum(origin),2)`), the direct
generalization of the current `modulo(sum(localSize(:,0)-1),2)`. Across a 2:1
jump checkerboard parity is not meaningful; the interface relaxes in the
block-Gauss-Seidel sense via the per-colour exchanges, the same regime the
solver is already in at rank boundaries.

## 7. Closed faces: removal inside the immersed boundary

A `FACE_CLOSED` face (neighbour removed) is an exact zero-flux face. As
implemented (Phase 2), the mechanism is the per-block face-mask machinery
already used for physical walls, NOT the `mu = 0` device this section
originally proposed (`mu` keeps its penalized value; making it exactly zero
would have required special-casing `update_ibm_mu`):

- the halo layer and the pinned interface velocity are zeroed once at init
  and never written again: momentum skips the pinned face (the same
  per-block start masks as walls), the SOR sweep window starts inside, and
  both the `denom` merge() and — new in Phase 2 — the face-correction
  merge() carry the same no-flux condition (on physical walls the masked
  correction was a dead value that apply_bc overwrote, so masking it
  changes nothing there);
- no exchange entries point at removed blocks, and the tangential halo
  extension generalizes from "physical wall" to "combined edge/corner
  neighbour absent", which is identical when nothing is removed;
- apply_bc serves `FACE_PHYS` faces only.

The closed face sits at least one cell inside the solid (the dilation
margin; a full block once the Phase-3 finest-level wall buffer exists), so
it carries no physics. Measured on the wavy-wall and sphere cases: fluid
cells agree EXACTLY with the no-removal run; velocities everywhere differ
at most by the `SOLID*mu` penalization residual (~1e-26); the only O(1)
differences are the decoupled penalized pressure inside surviving solid
cells, which nothing reads back.

Safety criterion recap: removable ⇔ block dilated by one cell is solid at all
four variable locations (cell centres and the three staggered grids). The
analytic IBM classifies at solver init; the file-based IBM reads the
`block_active` table written by `mobygeom.py block-active` into the
coefficient file. Everything else (partially solid blocks) stays and is
handled by the IBM exactly as now.

Payoff: for a typical immersed body the interior of the geometry plus the
penalized far-solid cells disappear from memory, FLOPs, *and* the SOR iteration
count pressure (no huge-`coef` cells iterating uselessly).

## 8. What changes in each module (survey)

- `init.f90`: add `block_set_type`; `init_grid` becomes per-level line
  generation + per-block metric slicing (reusing `init_grid_direction`);
  `init_field` allocates the 5D arrays.
- `comm.f90`: exchange-entry lists (Section 5); allreduces unchanged.
- `step.f90`: outer block loop; `uStartX/vStartY/wStartZ` and the boundary
  start indices become per-block per-face data derived from `nbrKind`;
  CFL/Peclet limits loop over blocks then allreduce as today.
- `pressure_solver.f90`: outer block loop; `sweepLo` / `pressureNeumann*`
  generalize to per-block per-face masks (one `integer(C_INT)` face-kind read
  per cell-face replaces the current scalar flags).
- `ibm.f90`: trailing block dimension; `set_ibm_coeff` loops blocks (wall
  blocks are all finest-level by construction); `mobygeom` coefficient file
  gains the block table layout.
- `boundary.f90`: `apply_bc` applies to blocks with `FACE_PHYS` faces only.
- `io.f90` / `field_hdf5.c`: datasets become `(nb,nb,nb,NVAR,nBlocks_global)`
  plus the block table; any-rank-count restart via Z-order linear distribution.
  (The run-length compression of Nakahashi & Kim is unnecessary with HDF5
  chunking, but the block table is the same compact metadata.)
- `mobygrid.f90`: Section 4. The solver never refines, never searches, never
  walks a tree.

## 9. Time stepping

One global `dt`, limited by the finest level (the per-block `d1` metrics feed
the existing CFL/Peclet reduction). No sub-cycling: BCM also advances all
cubes synchronously, and the RK3 + projection structure is untouched. Local
time stepping would break the single-pressure projection and is explicitly out
of scope.

## 10. Choosing `nb`

Halo memory/compute overhead per block is `(1+2/nb)^3 - 1`: ~40% at `nb=16`,
~20% at `nb=32`, ~10% at `nb=64`. Small `nb` ⇒ tighter geometry adaptation and
better load granularity; large `nb` ⇒ less halo traffic and fewer exchange
entries. Start with `nb=32` in 3D (CUBE uses 16^3, BCM 64^2 in 2D); make it an
input and measure. `nb` must be even (red-black) and ≥ 4 (so restriction
stencils never span more than one neighbour).

## 11. Phased implementation and verification

Each phase leaves the code releasable and is verified before the next.

1. **Phase 0 — block container, 1 block/rank.** Introduce `block_set_type`
   with `nBlocks=1`, `nb = localSize`. Pure refactor; results must match the
   current code to machine precision (`tools/compare_fields.py`).
2. **Phase 1 — many same-level blocks per rank.** Z-order distribution,
   block-pair halo entries (COPY only), per-block face masks replacing the
   per-rank flags. Verify: channel DNS statistics vs Phase 0; divergence norm
   identical behaviour; single-rank multi-block == multi-rank single-block.
3. **Phase 2 — solid-block removal.** `mobygrid` classification + `FACE_CLOSED`
   handling. Verify: wavy-wall IBM case gives identical fields in fluid cells
   with and without removal (removal must be a pure no-op on the physics);
   measure memory/runtime gains.
4. **Phase 3 — 2:1 refinement.** Per-level lines, RESTRICT/PROLONG in
   pack/unpack, fine-owns-face rule, finest-level wall buffer in `mobygrid`.
   Verify: (a) uniform-flow preservation across interfaces to round-off;
   (b) global mass conservation (sum of divergence) to round-off;
   (c) Taylor-Green / channel with an artificial refinement patch vs uniform
   fine reference; (d) IBM case (sailplane tutorial) vs uniform-fine result.
5. **Phase 4 — performance.** Internal/external block zoning + nonblocking
   overlap (merges with `nonblocking_overlap_strategy.md`); profile pack/unpack
   vs sweep kernels.

## 12. Deferred: dynamic adaptation and load balancing

The block object is deliberately sufficient for both, so nothing in Phases 0-4
needs rework:

- **Sensor**: per-block Laplacian-filter sum (Nakahashi & Kim eq. 15) over the
  velocity field; split/merge decisions on whole blocks only.
- **Re-mesh**: regenerate the block table (the octree pass is cheap and
  serial-per-rank on metadata), conservatively transfer fields
  (restriction/injection already exist), re-run the Z-order distribution.
- **Load balancing**: weighted dual-graph partitioning with per-block weights
  (cells + IBM cost), as Jansson Section 6.5; only relevant once block counts
  per rank are large and IBM cost is unbalanced.

## References

- K. Nakahashi, L. Kim, *Building-Cube Method for Large-Scale, High Resolution
  Flow Computations*, AIAA 2004-434. (Cube generation, 2:1 smoothing, removal
  of cubes inside the body, ghost-cell transfer, adaptive cube refinement.)
- N. Jansson, R. Bale, K. Onishi, M. Tsubokura, *CUBE: A scalable framework
  for large-scale industrial simulations*, IJHPCA 33(4), 2019. (Flat non-tree
  block arrays, halo interpolation at 2:1 faces, Z-order linear distribution,
  overlapped exchange, dynamic load balancing.)
