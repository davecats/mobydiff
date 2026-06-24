# 2:1 edge/corner reconstruction — strategy for a consistent (≥O(h)) corner

Branch `claude/jacobi-interface`. Status after commit `e3a4a7a`.

## The problem

The momentum predictor needs accurate deep-halo ghosts across a 2:1 interface
(see `validation/momentum_interface`). The ghost is a one-sided extrapolation
from the fine interior, `q(0) = c1·q(1) + c2·q(2) + c3·q(3)`:

- **Cubic** `(3,-3,1)` — 2nd-order accurate (the planar slab fine band is 2.0),
  but its L1 norm is 7, so it **amplifies a high-k mode ×7**.
- On a *single* (planar) interface the ×7 is harmless (planar `interface_decay`
  is stable). At a **2:1 edge/corner** a block reconstructs in 2–3 directions and
  those amplified halos feed the *same* corner cell's cross-advection
  (`uv`/`uw`/`vw`); the combined gain **blows up** (3D-patch `interface_decay`).
  Decisive: the blow-up is **dt-scaled** (smaller dt → less growth) ⇒ the
  PREDICTOR reading amplified halos, not a projection null mode; and
  `MOBY_NORECON` is fully stable.

## What is already done (commit e3a4a7a)

Lower the extrapolation **order** at edge/corner blocks (`nIf ≥ 2`, where `nIf` =
number of 2:1 interface faces of the block) to the **linear** ghost `2q1-q2`
`(2,-1,0)`, L1 norm 3 — bounded enough to stay stable under the corner coupling.
Planar blocks (`nIf < 2`, including all channel wall bands) keep the cubic.

Result: stable; planar slab still 2nd order; the corner operator is now
**consistent/converging** (patch fine-band rms 1.4e-3→5.3e-4, ~14× smaller than
the previous O(1)/diverging skip).

## Why it is not yet clean O(h)

The linear ghost makes the **normal diffusion vanish**: with `q(0)=2q1-q2`,
`d²v/dy² ≈ (q(0)-2q(1)+q(2))/h² = 0`, an **O(1)** error wherever the true
curvature ≠ 0. So the corner fine-band order is still sub-O(h) (~0.5–0.86,
decreasing). This is fundamental: a **consistent diffusion needs the ghost to
carry curvature ⇒ at least a quadratic extrapolation**, which is exactly the
(3,-3,1)-class that amplifies. One cannot get both from the ghost alone.

## Strategy: decouple advection and diffusion at the interface

The fix is to stop using the (amplifying) deep-halo ghost for the **diffusion**,
and treat advection and diffusion separately at edge/corner cells:

1. **Advection** — keep the bounded **linear** deep-halo ghost (`2q1-q2`). It is
   O(h)-consistent for the advective flux and does not amplify. (Already in.)

2. **Diffusion (the new piece)** — replace the deep-halo-based second derivative
   at the interface-adjacent cell with a **composite one-sided stencil** that
   does NOT read an extrapolated deep halo, so nothing amplifies:
   - Use the **owned/restricted interface FACE value** (already exchanged,
     conservative) plus the fine interior `q(1),q(2),q(3)` to build `d²/dn²` at
     the interface cell with one-sided coefficients on the *fine* node line.
   - This is the generalisation of the reverted **increment 1**
     (`blocks.f90 correct_interface_diffusion`/`fix_lap_interface`, which rebuilt
     the interface-row Laplacian coefficients with the correct 1.5·h metric) —
     revive it, but (a) for the NORMAL component too and (b) gated to edge/corner
     blocks (`nIf ≥ 2`) so planar blocks keep the validated cubic path.
   - Because the stencil reads only the owned face + interior (no ×7 ghost), it
     is bounded ⇒ stable, and one-sided 2nd-order coefficients give a consistent
     (O(h) or better) curvature.

With (1)+(2) nothing in the corner momentum stencil amplifies, and both the
advective and the diffusive operators are consistent ⇒ the corner fine band
should rise from ~0.5 to ≥1 (O(h)), ideally toward 2.

### Alternatives considered (lower priority)

- **Limited (TVD/minmod) ghost everywhere** — bounded, but minmod degenerates to
  the linear slope for smooth fields, so it has the same `dif=0` problem; it does
  not give consistent diffusion. (Verified: linear-everywhere drops slab v to ~0.9.)
- **Filtered reconstruction source** — pre-smooth the fine interior to remove the
  amplifiable tangential high-k before the cubic. Heuristic (filter strength), and
  it dissipates physical content; only if the composite stencil proves hard.

## Gates (must all hold)

- `tutorials/interface_decay` 3D-patch **stable** (the must-not-regress gate).
- Planar slab fine band still **2.0** (`run_gate.sh`, slab Axis-1).
- **Patch fine band order rises from ~0.5 to ≥1** (the new target;
  `run_gate.sh`, patch Axis-1).
- Axis-2 mass round-off; reflux still conserves; **bit-exact with no interface**;
  **CPU==GPU** bit-identical.

## Pointers

- Reconstruction: `reconstruct_interface_halos` in `src/modules/step.f90`
  (the `c1/c2/c3` coefficient switch on `nIf`).
- Diffusion coefficients: `blk%lapY*` etc. built in `src/modules/blocks.f90`;
  the momentum diffusion is in `momentum` (`src/modules/step.f90`).
- Reverted inc-1 reference: git history around commit `4960507`.
- Memory: `corner-reconstruction-todo`.
