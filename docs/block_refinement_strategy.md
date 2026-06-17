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

**As shipped (default scheme).** For `nLevels > 1` the production path is the
*composite projection* — see the companion `docs/composite_projection_strategy.md`
for the full derivation. The block ABOVE each 2:1 shared face owns it at its
interior low face `v(1)` (the block below holds it as a `v(nb+1)` halo), and the
red-black sweep reconstructs that owned face from the in-projection pressure
change `Δp = p − p_start`, `v = v* − dt·mu·ifGrad·(Δp_above − Δp_below)`, so the
coupled coarse-fine SPD system is relaxed in situ. The cross-level coupling is
carried by the per-colour `[u,v,w,p]` exchange (RESTRICT/PROLONG); single-level
grids never enter the interface path and stay bit-exact. The transfer operators
and the ownership rule described below are unchanged — but the *relaxation* it
describes ("symmetric BCM relaxation with stage-frozen ghosts") was the
predecessor scheme, and §6a a revision that was explored but not adopted; the
composite projection superseded both. The rest of this section is kept for the
design rationale and the conservation/transfer details that still apply.

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

### 6a. Explored (not adopted): fine-authoritative normal velocity + momentum reflux

**Status: superseded.** This subsection records the turbulence-validation
diagnosis that motivated revisiting the interface treatment and a candidate
fix (fine-authoritative normal velocity + Berger-Colella reflux). That
candidate was *not* implemented; the shipped response was the composite
projection (above-block-owns, faces reconstructed from the in-projection
pressure change — see the lead of §6 and `composite_projection_strategy.md`).
The diagnosis below remains the reference for what an interface scheme must
get right in energetic turbulence.

Turbulence validation (`validation/channel_interface`, interfaces at
y+=55 and y+=112, in the buffer/log layer) shows the "low-side owns"
rule is not accurate enough when an interface sits in energetic
turbulence: a spurious Reynolds stress appears localized at the
interface (the −⟨u'v'⟩ excess is almost entirely a *correlation*
overshoot, ρ_uv up to +5-6 %, not an rms change), the streamwise
spectrum piles up at the fine grid-scale just below the interface and is
deficient just above, the mean shear-stress balance carries a constant
excess in the fine band with a jump across the interface (an
un-refluxed momentum flux), and a spurious spanwise mean flow grows in
time. The error scales with the local turbulence intensity (y+=55 worse
than y+=112) and is wall-asymmetric, tracking the two opposite-handed
interface orientations.

Root cause: the interface treatment is conservative for **mass** but not
for **momentum**, and the interface-normal velocity is handled
asymmetrically. In this storage convention a block computes the normal
faces `v(1..nb)` (its low boundary + interior), and its high boundary
face `v(nb+1)` is a halo. So at a *fine-low* interface (fine band below,
coarse above — the bottom wall) the shared normal face is the coarse
block's interior `v(1)` and the fine block's halo `v(nb+1)`: the **coarse
side computes the interface-normal velocity at coarse resolution and the
four fine sub-faces are injected copies** (the fine side never carries a
fine-scale wall-normal velocity at the interface, exactly where the
turbulent flux lives). At a *fine-high* interface (top wall) it is the
reverse and the fine side is already authoritative. `interface_boxes`
encodes the same asymmetry: a high (`off=+1`) restrict writes the
boundary normal face at index `nb+1`, but a low (`off=-1`) restrict
writes only the halo at index 0, never the boundary face at index 1.

The fix has two parts.

**Part 1 — fine-authoritative normal velocity via uniform redundant
top-face computation (uniform Route B).** Rather than special-casing the
interface, every block momentum-computes *its own* top normal face —
`u(nb+1)`, `v(nb+1)`, `w(nb+1)` — redundantly with the block above/
right/in front, exactly as the SOR sweep already recomputes its open
halo layer redundantly to make results independent of `nb` and rank
count. With every block holding its own high-side normal face the
compute/own asymmetry of the "low-side owns" rule simply dissolves: at a
fine-low interface (bottom wall) the four fine sub-faces are computed by
the fine block at fine resolution — there is no injected coarse predictor
— and the coarse interface face is the area-restriction of those four,
identical to the fine-high (top wall) case. The two walls become
symmetric and the normal-velocity transfer direction is *not*
phase-dependent: the normal velocity always reconciles fine→coarse
(restrict), and coarse→fine prolongation only fills the fine block's
stencil halo, never overwrites a computed face.

What this requires:

- **Two-cell high-side halo.** The top-face stencil for `v(nb+1)` reaches
  `v(nb+2)`, so `q`/`qs` carry a two-cell halo on the *high* side of each
  direction (`0:nb+2`); the low side stays one cell (`0`), because the
  top-face stencil never reads below `v(0)`. Use `0:nb+2`, not
  `-1:nb+2` — the low-side ghost would be cosmetic halo for nothing.
- **`oldrhs` high-side layer, recomputed.** `oldrhs` (RK history) grows
  to `1:nb+1` on the high side and is recomputed redundantly each
  substage. It is a pure function of exchanged state, so no extra message
  and it cannot desync.
- **Two-deep normal halo fill.** The exchange fills the normal-velocity
  halo two layers deep: a same-level COPY writes the neighbour's
  `v(1),v(2)` into `v(nb+1),v(nb+2)`; a coarse-above interface PROLONG
  fills both fine layers. This is what makes the redundant top-face
  computation bit-identical to the neighbour's for same-level blocks —
  the Inc-2 canary (channel nb=4 still bit-exact vs Phase 2).
- **Metric at the top face.** `lapYp(nb+1)` (the coefficient of `v(nb+2)`)
  must carry the correct one-sided spacing.
- **Symmetric, ownership-based sweep masks.** A block corrects a normal
  interface face iff it owns it (the neighbour across is *coarser*,
  `FACE_COARSE`) or the face is open; it does not correct `FACE_FINE`
  (the finer neighbour owns) or pinned faces. Interface faces stay in
  `denom` on both sides (the Phase-3d stability requirement). The coarse
  interface face becomes the restriction of the four fine sub-faces — a
  slave used only in the coarse cell's divergence and tangential-momentum
  flux, never an independent DOF nor re-injected. Structurally the low
  (`off=-1`) restrict of the normal component must write the boundary
  face at index 1, which `interface_boxes` currently omits.

Computing the top face on *every* block (including same-level) is more
than the interface strictly needs — the asymmetry is only at 2:1 faces —
but it buys branch-free, uniform inner-loop logic at the price of a wider
halo and some redundant FLOPs. That overhead (the `0:nb+2` halo and the
universal redundant top-face computation) is a deliberate, known cost and
stays on the Phase-4 optimisation list; the simpler version is to compute
the top normal face only where a block is fine at a `FACE_COARSE`
interface.

Gated increments (each independently falsifiable). As implemented, the
momentum change was targeted to `FACE_COARSE` high faces (the fine-low
case) rather than computed on every block: a same-level top face is
overwritten by the halo refresh anyway, so computing it everywhere is
dead work, and targeting keeps the same-level path untouched.

- **Inc 1** (done, validated 2026-06-15) — `q`/`qs` carry a `0:nb+2`
  high-side halo, `oldrhs` `1:nb+1`; `lap*(nb+1)` filled; the same-level
  COPY fills two high-side layers. No momentum change. channel nb=4
  bit-exact (CPU+GPU), uniform-flow exact, `MOBY_HALO_AUDIT` clean.
- **Inc 2a** (done) — the momentum predictor and the copy-back compute
  the top normal face `v(nb+1)` at a `FACE_COARSE` high face (the fine
  block of a fine-low interface). Exchange unchanged ⇒ the face is still
  overwritten by the existing prolong ⇒ bit-exact vs Inc 1 everywhere.
- **Inc 2b** (done, validated 2026-06-15) — interface ownership flip,
  entirely in the exchange (`comm.f90` `faceNrm`): the fine block's
  PROLONG no longer injects its owned `v(nb+1)` (writes only the `nb+2`
  stencil and the tangential/pressure halo); the coarse block's RESTRICT
  writes `v(1) = ¼Σ v^f(nb+1)` to its interior face (index 1) **and** its
  outward halo `v(0)` (index 0). Two hard-won points:
  - **Sweep corrections stay SYMMETRIC** (`face_pinned` only). Making the
    coarse skip its interface face (`no_correct` on `FACE_FINE`) is the
    one-sided splitting that blows up — re-confirmed here, reverted.
  - **The RESTRICT must keep writing the coarse's normal *outward halo*
    `v(0)`**, not only the interface face `v(1)`. The coarse still
    *predicts* `v(1)` in its own momentum, whose stencil reads `v(0)`;
    a dropped `v(0)` leaves it stale and breaks the constant-field
    cancellation, injecting a divergence-free spurious tangential mode at
    fine-low interfaces (uniform flow no longer exact). With `v(0)`
    restored: uniform-flow spread **exactly 0** at any SOR count,
    divergence 0, `interface_decay` contracts (200 steps), channel nb=4
    bit-exact, halo audit clean — CPU and GPU.
- **Inc 4** (convective y-reflux done, validated 2026-06-15;
  `src/modules/reflux.f90`) — tangential-momentum reflux (Part 2). The
  coarse cell's interface flux is replaced by the area-restriction of the
  four fine fluxes (Berger-Colella); the correction vanishes identically
  for uniform flow. Done: the convective flux at y-interfaces (the channel
  bands), both tangential components, same-rank fine neighbour, CPU+GPU.
  Gates: uniform-flow spread 0 (CPU+GPU), reflux non-trivially active on
  the noisy patch with `interface_decay` still contracting, channel nb=4
  bit-exact (no-op without refinement). Pending (all vanish for uniform,
  so the exact gates are unaffected): the viscous flux, x/z interfaces
  (body refinement only), and the cross-rank flux exchange (a y-interface
  whose fine side is off-rank currently aborts). The turbulent
  −⟨u'v'⟩ closure is validated on the big machine.

**Part 2 — tangential-momentum reflux.** After the predictor, the
   convective+viscous flux of the tangential momentum (u, w) through the
   interface differs between the coarse computation (one face) and the
   fine computation (sum of four). Accumulate both, restrict the fine
   flux, and correct the interface-adjacent cells by the mismatch
   (Berger-Colella reflux). The correction is constructed to vanish
   identically for uniform/constant flow so the exact-uniform-flow gate
   and single-level bit-exactness are preserved.

Gates for the revision: the existing exact gates must still hold
(uniform-flow spread 0, channel nb=4 bit-exact vs Phase 2, mass to
round-off, `MOBY_HALO_AUDIT` clean, `interface_decay` contraction over
200 steps — run at every step because this is the instability-prone
region); then the turbulence validation should show the −⟨u'v'⟩
correlation overshoot, the spectral pile-up and the wall asymmetry
collapse.

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
   pack/unpack, the interface ownership/relaxation rule, finest-level wall
   buffer in `mobygrid`. Verify: (a) uniform-flow preservation across
   interfaces to round-off; (b) global mass conservation (sum of divergence)
   to round-off; (c) Taylor-Green / channel with an artificial refinement
   patch vs uniform fine reference; (d) IBM case (sailplane tutorial) vs
   uniform-fine result.
5. **Phase 4 — performance.** Internal/external block zoning + nonblocking
   overlap (merges with `nonblocking_overlap_strategy.md`); profile pack/unpack
   vs sweep kernels.

**As-built status** (the per-phase log with commit hashes and exact gate
results lives in `CLAUDE.md`; this is the summary):

- Phases 0–2 completed and bit-exact as planned.
- Phase 3 shipped the **composite projection with above-block-owns** (§6 lead),
  not the original fine-owns-face rule; geometry-driven refinement
  (`refine_body`) works for both analytic and file-based (`mobygeom`) IBM.
- Phase 4 (performance) so far, all bit-exact: the halo exchange is one
  weighted gather with a precomputed per-point→entry lookup (no per-point
  binary search); the redundant projection-entry velocity exchange was
  removed; the intermediate composite exchanges ship pressure only on the
  cross-level (interface) entries; and the IBM penalization update is skipped
  when no immersed boundary is present. Nonblocking compute/comm overlap is
  still open.
- I/O: fields are written per block as `(nBlocksGlobal, nb_z, nb_y, nb_x)`
  datasets plus the `blocks` table; each write also emits a ParaView **XDMF**
  sidecar (`<prefix>_<step>.xmf`) — one structured grid per block referencing
  the HDF5 hyperslabs — so block-decomposed and 2:1-refined fields load
  directly without reassembly.

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
