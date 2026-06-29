# Next session — edge/corner + LES validation, and a less-dissipative 2:1 interface

Branch `claude/jacobi-interface`. This is the CURRENT authoritative handoff for the
2:1 refinement work; it supersedes the deleted `next_session_*` / `interface_band`
handouts. Lineage (still valid background): `docs/jacobi_interface_handout.md`
(why the branch exists, the damped-Jacobi restart), `docs/interface_review.md`
(the decisive "collocated papers, staggered code" observation),
`docs/interface_projection_derivation.md` (SPD projection),
`docs/momentum_interface_handout.md` (the reflux), `docs/corner_reconstruction_strategy.md`
(edge/corner reconstruction), `docs/block_refinement_strategy.md` (master plan).
Memory: `interface-validation-suite`.

## State (resolved 2026-06-29) — the FLAT 2:1 interface works

- **Constant-1/2 interface is the DEFAULT** (`[blocks] interface_constant_half`):
  inject the velocity prolong + skip the cubic deep-halo reconstruction. Stable
  (the cubic/metric weights destabilise; const-1/2 = V&V's order-for-energy trade).
  The pressure projection is symmetric/SPD at the interface (`face_grad` composite
  stencil + conservative copy reconciliation) — what keeps Chebyshev-Jacobi stable.
- **`momentum_reflux` should be OFF.** The u'/v' bands were the reflux injecting the
  fine-side resolved Reynolds-stress flux (avg(u^2)-(avg u)^2) into the coarse cell
  — conservation-correct but the AMR coarsening band. Reflux OFF removes the band at
  ZERO stability/divergence cost.
- **Turbulent validation done** (developed channel Re_tau 180, t=5..25, 2:1 wall
  bands): reflux-off matches the uniform-256 reference — `-<u'v'>` to 0.2%, no
  spurious band (x-y/z-y cross-sections clean). Isolated against a matched
  **uniform-128** control, the only residual is a **~5% under-transmission of
  small-scale v'/w' energy** into the coarse core (the const-1/2 restriction is
  mildly dissipative; a loss, not a band; mean transport exact). Half of the raw
  core deficit vs uniform-256 is just coarse resolution; half is this interface loss.

Validation assets (`validation/channel_interface/`): `interface_benchmark/`
(~250-step stability + u'/v' band metrics, fast), `developed/` (run_reflux_study.sh,
run_reference.sh, run_uniform128.sh — two-leg stats; analysis scripts produce the
4-way, isolation, and x-y/z-y cross-section figures), `reference.ini` (uniform-256),
`uniform128.ini` (base-128 control). Beltrami slab energy gate: `MOBY_KEBAL` (now
per-component u/v/w) on `validation/beltrami/slab_y_diag.ini`.

## Open item 1 — VALIDATE edges and corners (single refined block)

The turbulent validation so far is FLAT-FACE only: the wall bands span the full
x-z plane, so the only 2:1 interfaces are y-faces — no edges (two interface faces
meeting) or corners (three). A **single refined block** (or a small embedded patch)
in the coarse interior exercises all 6 faces + 12 edges + 8 corners.

- Uniform-flow / Beltrami exactness on edges+corners already passes (Phase 3b/3c
  gates; `validation/beltrami/refined_fast.ini` is the 3D patch). The corner
  reconstruction is the slope-limited cubic (`docs/corner_reconstruction_strategy.md`,
  `lim_extrap`, corner-gated) — but note that is part of the cubic deep-halo
  reconstruction which is SKIPPED under const-1/2. Check what const-1/2 actually
  does at edges/corners (the velocity prolong there keeps plain injection — see
  `entry_blend` "Edges and corners keep plain injection").
- BUILD a turbulent test with an embedded refined patch (edges+corners in the flow,
  not just at the walls), e.g. a channel-core patch or a refine_body case, and run
  the same gates: interface_benchmark-style stability + band metrics, then developed
  stats vs a uniform-fine reference. Gate: stable, no band/spurious energy at the
  edges/corners, mean transport correct. The reflux is OFF, so the edge/corner
  question is purely the transfer + projection consistency there.
- Risk to watch: edges/corners were historically the hardest (corner diffusion is
  anti-diffusive — see corner-reconstruction-todo memory). With reflux off and
  const-1/2, re-confirm they are stable in TURBULENCE, not just Beltrami.

## Open item 2 — VALIDATE LES across the 2:1 interface

LES (`les.f90`, `les%nut`) has not been validated WITH refinement in turbulence.
- The eddy viscosity `nut` is a scalar carried per block with a trailing block index
  and exchanged (`exchange_scalar_halos(c, les%nut, blk)`); confirm its cross-level
  RESTRICT/PROLONG at the interface is consistent (it rides the same gather as phi).
- The Smagorinsky/WALE strain rate at the interface cells reads the (const-1/2)
  velocity halos — check the model contribution is sane across the resolution jump
  (the coarse side has a larger filter width; the model length scale must track the
  local cell size per level, which it should via the per-level metrics).
- Run a refined LES channel (higher Re, coarser base so LES is actually active) and
  compare to a uniform-fine LES reference: mean profile, stresses, and the `nut`
  profile across the interface. Gate: no `nut` discontinuity/band, stable, mean OK.
- Inert-path check first: a single-level LES run must be bit-exact with the
  interface code present (LES + const-1/2 are both inert without a 2:1 interface).

## Open item 3 — IMPROVE the interface formulation (the ~5% v'/w' loss)

The residual is the const-1/2 restriction being mildly dissipative (low-pass) for
the small-scale tangential components. Reducing it is an accuracy-vs-stability trade
(a less-dissipative transfer is exactly what reintroduced the band/instability that
const-1/2 cured). Candidate levers, in increasing cost/risk:

1. **Two-sided symmetric conservative momentum-flux reconciliation** (extend the
   3d-file pressure-side BCM symmetric relaxation to the momentum interface):
   distribute the interface flux mismatch across BOTH adjacent cells (each its own
   copy, reconciled to a conservative both-value) instead of the one-sided coarse-
   only reflux. Could de-asymmetrise and soften the loss WITHOUT the band — but must
   preserve the skew-symmetric energy property (symmetrising a skew operator is the
   delicate part) and the storage constraints below. This is the most promising,
   reuses existing pressure machinery.
2. **A higher-order / energy-consistent (adjoint, R = P^T) tangential transfer** for
   v'/w' only. Earlier shown to have ~no headroom on the ENERGY (KEBAL) metric, but
   its real value is operator SYMMETRY for Chebyshev (see the const-1/2 + Chebyshev
   analysis); it could reduce the dissipative loss if it keeps the projection SPD and
   the band off. Gate on the interface_benchmark (band NOT regrown) + KEBAL slab.
3. **Fine-owns interface velocity** — the principled staggered-AMR route (Almgren et
   al. 1998 composite projection). BLOCKED by storage, not the solver: the exchange
   writes HALOS only (never an interior plane), and the staggered normal face is an
   interior index for one block / a halo index for the other, so one orientation
   needs to slave a coarse INTERIOR plane the exchange cannot write. Needs either
   interior-write transfer entries or a 2-cell halo — a global change. The
   segregated Jacobi/Chebyshev solver makes the slaving tractable (covered-face flux
   as a BC) where the old red-black couldn't, but the storage is the real obstacle.
4. **Collocated interface** (or whole-solver) — would remove the staggered face-
   ownership problem entirely (interface_review.md: the papers are collocated, this
   code is staggered) but sacrifices the energy-conserving staggered framework and
   needs Rhie-Chow. Big; last resort.

Decide whether the ~5% even warrants this: for a refined run the core is
intentionally coarse, so a ~5% small-scale loss there may be acceptable. Lever (1)
is the natural next experiment if it does.

## Guardrails (every change)

`validation/channel_interface/VALIDATION_CASES.md`: interface_benchmark stability +
band metrics; KEBAL slab (per-component); bit-exact no-interface (single-level
inert); 1-rank==2-rank (dims N 1 1, x/z split, NOT y). Build `-Mnofma` (CPU) /
`-Mnofma -gpu=nofma` (GPU) for bit-exact compares. Run ONE GPU case at a time;
`pkill build_*/main` between runs. Module: `module load
/opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3`.
