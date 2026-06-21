# Cubic tangential prolong — design, partial result, open residual

Status: **attempted, partial, reverted** (2026-06-21, unsupervised). The codebase
is at the clean committed state; this note captures the work so we can finish it
together.

## Goal

After the normal-velocity diffusion-Laplacian fix (`9b5234f`), the 2:1 interface
scheme is order-2 in the *normal* velocity (u at an x-interface) but the
*tangential* velocity floors at order ~1.5 (v at an x-interface). Root cause
(established earlier): the momentum reads the tangential velocity halo `v(0)` into
`d2v/dx2`, a SECOND derivative across the interface, which amplifies the halo's
O(hc^2) prolongation error by 1/hf^2 to O(1). The centred stencil is stable but
needs an O(hf^4) halo; a one-sided fine-only stencil is consistent but
anti-diffusive (it blew up the 3D patch — see interface_improvement_log.md). The
remaining stable route is to make the **prolong cubic** (O(hc^4) halo).

## Implementation (what was tried)

In `comm.f90`, extend the per-dim gather:

- `gather_taps`: `idx(0:3)`, `w(0:3)`, n up to 4. For `lin==1` (linProlong, the
  final exchange) use CUBIC on the fixed 4-point window `[base-1, base+2]` for
  every interpolated dim. The fixed window avoids the halo-depth asymmetry
  (1 low / 2 high halos): a covering coarse cell `base` is interior `[1,nb]`, so
  `base-1>=0` and `base+2<=nb+2` always hold; fall back to the 2-tap linear only
  if `base` is outside `[1,nb]` (block corner) so the read stays in bounds.
  - cell-staggered (`var/=dim`): `(15,135,-27,5)/128` at par0 (fine at base-1/4),
    `(-7,105,35,-5)/128` at par1 (base+1/4). Verified exact on x^3.
  - face-staggered own dim (`var==dim`): coincident face (par0) stays 1 tap;
    midpoint (par1) is `(-1,9,9,-1)/16`.
- `gather_point`: bump `i1/i2/i3` and `w1/w2/w3` to `(0:3)`. The accumulation
  triple-loop already uses n1/n2/n3, so it needs no other change. Cubic only
  touches the interface PROLONG entries; the same-level COPY fast path and the
  injection relaxation are untouched, so single-level stays bit-exact and the SOR
  contraction (stability) is unchanged.

## Result (slab_x TGV, GPU)

| | committed (normal-fix) | + cubic |
|---|---|---|
| L2(u) order 32->64 / 64->128 | 2.29 / 2.11 | 2.29 / 2.11 |
| L2(v) order 32->64 / 64->128 | 1.95 / **1.51** | 2.07 / **1.66** |

Cubic IMPROVES v (1.51 -> 1.66) and is **stable** (interface-decay passes) and
**single-level bit-exact** (uniform 2.182e-5). But it does NOT reach full order 2.

## Open residual (the part to finish, supervised)

The cubic helped but capped at ~1.66. The obvious suspect was wrong:

- **Deep coarse halo (nb+2)**: the cubic's 4th tap (`base+2`) can land on the
  coarse block's deep high halo. `interface_boxes` fills only nb+1 on a RESTRICT
  (`ext=merge(2,1,op==OP_PROLONG)`), leaving nb+2 stale. Setting `ext=2` to fill
  nb+2 from the RESTRICT changed the result by **0** (bit-identical v=1.2750e-4) —
  so the cubic's base+2 tap does NOT read the cell that fix touched. My geometry
  mapping of which coarse cell `base+2` hits must be off. **First thing to check.**

Hypotheses to investigate next:
1. Which physical coarse cell does each cubic tap actually read? Instrument the
   gather for one interface point and print `base, idx(0:3)` and the coarse
   values vs the analytic field. Confirms whether nb+2 (or which cell) is read
   and whether it is current.
2. Is the *normal* dim of the PROLONG actually interpolated, or injected? The
   `lLin` flags from `entry_gather_map` mark "tangential dims of a PROLONG" as
   linear; if the interface-NORMAL dim (the one d2v/dx2 differentiates) is
   injected (lLin=0), the cubic never engages on the dim that matters and we are
   only improving the in-face dims. The early v stencil dump showed v(0)
   interpolated in x, which argues lLin(x)=1 — but verify per entry.
3. Multi-dim requirement: v(0)'s y/z interpolation error appears only in v(0)
   (not v(1)/v(2)), so it does NOT cancel in d2v/dx2 -> every interpolated dim
   must be O(hf^4). Confirm the face dim (var==dim) cubic is engaging at the
   midpoints and not silently staying linear.
4. Pre-asymptotic: the cubic leading constant is ~16x (hc^4 = 16 hf^4); 64->128
   may not be asymptotic. A 256 point would tell, but that is a slow run.

## Verification still owed before committing any cubic

- rank-independence (CPU 1 vs 2 ranks) of the cubic gather — it changes the
  SHARED gather_point, which is exactly what guarantees bit-identical
  local-copy vs off-rank-pack arithmetic. A pre-revert run was launched but had
  not returned.
- refined 3D patch stability + the interface_decay gate (the one-sided attempt
  passed slabs but blew up the 3D patch; re-check the cubic on the patch).
