# Composite Projection Strategy for 2:1 Interfaces

Design notes for replacing the block-by-block-then-reconcile interface treatment
with a **composite** (single-system) pressure projection — the formally-clean
route to a stable, conservative 2:1 staggered interface under the incompressible
projection. Companion to `block_refinement_strategy.md` (the BCM block plan,
whose §6a "uniform-B + reflux" revision this supersedes) and
`interface_review.md` (the diagnosis that motivates it; see §vi–§vii). No code is
changed by this document; it records the agreed strategy.

## 1. Why this document exists

`interface_review.md` §vii establishes an impossibility result: in the
block-SOR-then-reconcile framework, a symmetric fine-authoritative 2:1 normal
velocity and a contractive projection are mutually exclusive. The cause is a
**redundant interface DOF** — the coarse block keeps its own copy of the shared
normal face, relaxes it, and the two copies are reconciled *after* the solve.
That is a splitting of an operator that is not a projection, so it carries no
contractiveness guarantee, and in practice it injects an O(1), dt-independent
divergence at the coarse interface cell every substage (measured: `MOBY_DIV_AUDIT`,
O(1) at the interface, exactly 0 in the interior).

The composite projection removes the redundancy: the interface has ONE
authoritative DOF set, the discrete divergence and gradient are exact adjoints
across the interface, and the pressure-Poisson is solved as a single coupled
system. Contractiveness and exact (local) mass conservation then follow **by
construction, not by tuning** — which is the whole point.

## 2. The formal requirement

The projection step is `u^{n+1} = u* − G p` with `p` solving `L p = D u*`, where
`D` is the discrete divergence, `G` the discrete gradient, `L = D G`. The
operator `P = I − G L⁻¹ D` is an **orthogonal** projection (`‖P‖₂ = 1`, hence
non-amplifying and contractive) **iff** `D` and `G` are adjoint in the discrete
inner product:

    G = −Dᵀ    ⟺    ⟨D u, p⟩_cells = −⟨u, G p⟩_faces    (summation by parts)

with the cell inner product weighted by cell volumes and the face inner product
by face control-volumes (area × normal spacing). When `G = −Dᵀ`,

    L = D G = −D Dᵀ

is symmetric negative semidefinite (SPD up to sign and the constant-pressure
null space) for **any** `D`, because it has the form `−M Mᵀ`. So once `D` is
fixed conservatively, taking `G = −Dᵀ` makes the composite Poisson SPD
automatically, and any SPD-convergent solver (SOR, CG, multigrid) yields a
contractive projection. There are no free weights to tune.

## 3. Interface DOF layout: no redundancy

At a 2:1 y-interface (fine cells size `h` below, coarse cells size `2h` above;
the interface is an x–z plane) the four fine faces are the authoritative
normal-velocity DOFs. The coarse interface cell has **no** independent
normal-velocity DOF on its fine-facing side — it sees the area-sum of the four
fine faces. There is therefore nothing to reconcile; the jump that broke B2
cannot exist.

## 4. Discrete divergence `D` (conservative by construction)

`D` is the control-volume divergence (sum of outward face fluxes / cell volume).
Only the interface rows differ from the uniform stencil:

- **Fine interface cell** (volume `h³`): standard; its top face is one fine
  interface DOF `v_f`, flux `v_f·h²`.
- **Coarse interface cell** (volume `8h³`): its fine-facing flux is
  `Σ_{m=1..4} v_f^(m)·h²` (the four covering fine faces); the other five faces
  are standard coarse faces.

Conservation is exact: the coarse cell's interface flux is identically the sum of
the fine fluxes crossing the same surface (the fine faces tile the coarse face
under midpoint subdivision). Global mass telescopes, **and** the local
interface-cell divergence is well defined with no double counting.

## 5. Gradient `G = −Dᵀ`: the interface stencil is pinned, not tuned

Take the transpose of `D` in the (cell-volume, face-control-volume) inner
product. Worked in the normal (y) direction, reduced to 1D per unit tangential
area:

- Fine cell `F`: centre `y = −h/2`, volume (per area) `h`, pressure `p_F`.
- Coarse cell `C`: centre `y = +h`, volume (per area) `2h`, pressure `p_C`.
- Shared interface face velocity `v` at `y = 0`.

The `v`-coefficients in `D` are `+1/h` in `(Du)_F` and `−1/(2h)` in `(Du)_C`.
With cell weights `vol_F = h`, `vol_C = 2h`, the `v`-terms of `⟨D u, p⟩` are

    (v/h)·p_F·h + (−v/(2h))·p_C·(2h) = v (p_F − p_C).

Matching `⟨D u, p⟩ = −⟨u, G p⟩` with the face control-volume weight `S_v` gives
`(G p)_v = (p_C − p_F)/S_v`. The staggered `v`-cell spans the two centres, so
`S_v = h/2 + h = 3h/2` (per area), hence the **adjoint interface gradient**

    (G p)_v = (p_C − p_F)/(3h/2) = (2/3)(p_C − p_F)/h.

Two points worth stating plainly:

1. **No freedom.** Once `D` is the conservative control-volume divergence, the
   interface gradient is the pressure difference over the cell-centre gap `3h/2`,
   weighted by the staggered metric. No blend parameter, no choice.
2. **This is the gradient the existing blend was already reaching for.** The
   `(2 p_C + p_F)/3` pressure-ghost blend, evaluated at fine spacing `h`, gives
   `((2 p_C + p_F)/3 − p_F)/h = (2/3)(p_C − p_F)/h` — identical. That is why the
   Step-4 weight sweep found the blend "near-optimal," and why Step-4's "weighted
   adjoint = injection with the 1.5× dual-volume weight" is the **same operator**
   in a different bookkeeping convention. **The gradient weight was never the
   problem.** What was missing is solving the coupled system this gradient
   defines, instead of using it as a frozen ghost in a block-by-block relaxation
   followed by a reconcile. The composite projection finishes what Step 4 started,
   in the right solver framework — and the old "overdriving/consistency" worry
   evaporates, because an exact projection is exactly divergence-free, so there is
   no residual gradient to compensate.

Exact code coefficients depend on the inner-product convention; the **source of
truth is the harness check `D + Gᵀ = 0` to round-off** (Section 9), reusing the
Step-3/Step-4 small-matrix machinery. Instrument it; do not assume it.

## 6. Composite Poisson and the solve

`L = −D Dᵀ` is SPD up to the constant-pressure null space (fixed by pinning one
reference pressure, or enforcing the compatibility condition `Σ D u* = 0` to
round-off). Solve `L p = D u*` over the **composite** grid as one coupled system:

- The interface rows of `L` couple the coarse interface-cell pressure to the four
  fine interface-cell pressures (through the summed fine faces and their adjoint
  gradients) and vice versa. The relaxation **must** use these composite rows —
  the coarse interface pressure update reads the fine pressures and the fine
  updates read the coarse pressure, in the same sweep. This is the one structural
  change from the current solver: **no frozen cross-level ghost, no
  `interp=.false.` + final reconcile.**

Solver options:

- **Composite red-black / multicolor SOR** (minimal change). SPD-ness guarantees
  convergence. The 2:1 jump breaks global checkerboard 2-colorability at the
  interface; use a consistent multicolor ordering there, or a damped-Jacobi
  smoother on the interface band (both converge for SPD). Contractiveness is the
  SPD guarantee, not a tuning target.
- **Geometric multigrid** (optimal, larger change). The refinement levels are a
  natural coarse-grid hierarchy (the IAMR choice); gives mesh-independent
  convergence. Recommended long-term; not required for correctness.

Single global `dt` is a simplification, not a complication: because mobydiff does
NOT subcycle in time, there is ONE composite projection per RK substage — no
level/sync-projection split (the temporal machinery ABC/IAMR need only because
they subcycle). The composite operator here is purely spatial.

## 7. Predictor / momentum side

The fine block computes its interface faces (fine-authoritative, as B2 intended).
The coarse interface cell's interior advance uses the summed fine interface flux
for its interface term, consistent with the same `D` used in the projection. The
tangential-momentum reflux is the momentum-side analog of the conservative
composite divergence (it corrects the coarse interface-adjacent cell's tangential
convective+viscous flux to the area-sum of the fine fluxes). With the composite
projection there is no post-reconcile overwrite, so the predictor and the
projection see the same interface fluxes.

## 8. Exact vs approximate; the high-k caveat

mobydiff is staggered, where the **exact** MAC projection is the natural, clean
object — target the exact composite MAC projection above, not the cell-centered
**approximate** projection. The approximate projection is a workaround for the
non-solenoidal null space of the exact *cell-centered* projection, which
staggered discretizations do not have; it is formally weaker (only approximately
divergence-free) and so the wrong choice here.

Caveat (Olshanskii et al.): the exact composite staggered operator on graded
grids can admit weakly-controlled **high-k** harmonic interface modes
(checkerboard-type) in the joint null space of `D` and `Dᵀ` that the projection
does not remove. The current channel defect is **low-k**, which the exact
projection removes cleanly, so no filter is needed for it. If a high-k interface
mode appears, the principled cure is a consistent high-pass filter /
hyperviscosity that provably vanishes on smooth fields (NOT an ad-hoc blend or
sponge), added only if observed and verified to leave the uniform-flow and
single-grid gates exact.

## 9. Phased implementation and verification

Each phase leaves the code releasable and is verified before the next.

- **P-a — operators + adjointness harness (no solver change).** Implement the
  composite `D`, `G = −Dᵀ`, `L = −D Dᵀ` at interface cells (face-kind masks
  select the interface rows). Verify on the minimal 2-block matrix harness (reuse
  Step-3/Step-4): `D + Gᵀ = 0` to round-off; `L` symmetric and negative
  semidefinite. No dynamics yet. **This phase is cheap and de-risks everything —
  if `D + Gᵀ` is not zero on the real layout, stop and fix the operators before
  any solver work.**
- **P-b — composite solve on a static 2:1 patch.** Replace the block-reconcile
  pressure solve with the composite SOR (or MG) on a single static refined patch.
  Gates: uniform-flow spread 0; `L`-residual to round-off; the projected velocity
  divergence-free **locally** at the interface cell (the O(1) residual
  `MOBY_DIV_AUDIT` measured is now ~0 by construction); channel-nb4 (no interface)
  bit-exact vs S0.
- **P-c — turbulent channel.** Run `channel_interface` from the KMM restart.
  Because the projection is now orthogonal/contractive, `interface_decay`
  contraction and the low-k loop-gain ≤ 1 should pass **by construction**
  (verify, do not assume); global **and** local interface mass to round-off; then
  the −⟨u'v'⟩ correlation overshoot, spectral pile-up and wall asymmetry should
  collapse — the correctness payoff B2 was meant to deliver, now on a stable base.

Standing exact gates throughout: uniform-flow spread 0; channel-nb4 bit-exact;
`MOBY_HALO_AUDIT` clean; build BOTH `-Mnofma` / `-Mnofma -gpu=nofma`.

## 10. Effort and risk

Core effort: restructure the pressure solve so the interface rows of `L` couple
the levels in one system (composite SOR/MG), replacing `interp=.false.` +
reconcile. Fiddly parts: the multicolor ordering across the 2:1 jump, and the
loss of "every block identical" uniformity at interface cells (masked by
face-kind — a small fraction of cells). Convergence-rate risk: composite SOR may
be slow near interfaces (the 2:1 jump stiffens `L`); multigrid removes this but is
a larger change. Formal upside: contractiveness is a theorem (SPD), verification
is "is `L` SPD" rather than "did it blow up," and the result is exactly (locally)
divergence-free — no tuned weights, no dissipation, no reconcile.

## References

As in `interface_review.md` §v: Berger & Colella (1989); Almgren, Bell, Colella,
Howell & Welcome (1998); Martin & Colella (2000) / Martin, Colella & Graves
(2008); IAMR/AMReX; Olshanskii et al. (octree-MAC, staggered C/F spurious modes);
and the mimetic-finite-difference orthogonal-decomposition result `D = −Gᵀ ⟹`
discrete orthogonal projection.
