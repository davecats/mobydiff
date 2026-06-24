# Momentum-predictor 2:1-interface gate

The 2:1-interface **projection** (pressure correction) is divergence-consistent
and mass-conserving to round-off (branch `claude/jacobi-interface`, commit
`d2e743a`). This gate instruments the remaining piece: the **momentum predictor**
(advection + diffusion, `src/modules/step.f90`) at the interface, which was both
unfixed and untested.

It is a manufactured **operator-truncation** test — a single RHS evaluation, not
a time-marched flow. One step is sufficient; this measures operator consistency,
not evolution. Two axes, both mandatory (a change can pass one and break the
other):

- **Axis 1 — accuracy / order.** `MOBY_PREDONLY` (no projection) + `MOBY_RHSDUMP`
  dumps the discrete momentum RHS `L_h(u) = -div(uu) + (1/Re)lap(u)` (the
  predictor's `blk%oldrhs`, captured on the pristine field) as `un/vn/wn`.
  `tools/rhsband.py` compares it to the analytic `L(u)` of the `tgv3d` field at
  the staggered points and reports the convergence order over nx=32/64/128, split
  interior / coarse-band / fine-band and per component. Interior must be ~2
  (sanity); the **fine-band order is the verdict**.
- **Axis 2 — continuity null-space.** The projection preserves the *mean* of
  `div(qs)`, so the global `Sum(vol*div(qs))` after predictor+sync (the net 2:1
  interface flux mismatch) is the one continuity error it can never remove. A
  real one-step run with `MOBY_DIVDUMP` + `tools/divsum.py` must keep it at
  round-off.
- **Axis 3 — momentum conservation** (`MOBY_TERMDUMP` → `tools/momsum.py`). The
  divergence-form advection and the viscous flux conserve total momentum: over a
  periodic domain every interior face flux telescopes, so `Sum(vol*term_i)` is
  round-off on a single-level grid. The only faces that do not cancel are the
  2:1 interface faces (coarse flux ≠ summed fine sub-face fluxes), so this global
  sum IS the momentum-conservation error — the momentum analogue of the mass gate.
  Run `MOBY_TERMDUMP=1/2/3` (u/v/w), feed the `_adv` (and `_dif`) dump to
  `momsum.py`. (Conservation is an integral property; per-cell `vol*term` is the
  physical flux divergence, nonzero everywhere — so the metric is the global sum,
  with the `uniform` case as the zero reference and a single-interface case to
  isolate one interface's local imbalance.)

### Momentum-conservation baseline (inc 5, no reflux)

`Sum(vol*adv)` per component, tgv3d, nx=64 (rel = vs `Sum(vol*|adv|)`):

| case | u | v | w |
| --- | --- | --- | --- |
| uniform | -1.5e-16 | -7.1e-17 | -2.9e-15 | (round-off — tool valid, operator conserves) |
| slab (y-iface) | +8.9e-7 (1.2e-8) | **+4.9e-4 (6.8e-6)** | +8.9e-7 (1.2e-8) |
| patch (x,y,z) | +1.2e-4 (1.7e-6) | +1.2e-4 | +1.2e-4 |

`Sum(vol*dif)` is round-off (~1e-17) for every component/case — **the viscous
momentum flux is already conserved across the interface**; only advection leaks.

Decomposition of the advective leak (key for any reflux):
- **Tangential momentum** (e.g. `uv`/`wv` through a y-interface, the component
  *parallel* to the interface): a CLEAN flux register — the flux lives ON the
  interface y-face — but small for this manufactured field (slab u,w ~9e-7). In
  real turbulence this is the Reynolds-stress flux and is the physically decisive
  one (interface_review §iii: the un-refluxed tangential flux is the −⟨u'v'⟩
  defect).
- **Normal momentum** (`vv` through a y-interface): DOMINATES here (slab v
  4.9e-4) but is the *hard* case — the v-momentum control volume is centred ON
  the interface and its `vv` flux is evaluated at cell centres that do **not**
  align across the 2:1 interface (coarse `y_int−h_c/2` vs fine `y_int−h_f/2`), so
  there is no single interface flux to register. This is the staggered
  normal-velocity-on-the-interface difficulty (interface_review §A, §vii) and
  needs a momentum-consistent interface treatment (the composite/approximate
  projection route), not a plain flux register.

The field is `initial = tgv3d` (added to the generic case):
`u=sin(kx)cos(ky)cos(kz)`, `v=cos(kx)sin(ky)cos(kz)`, `w=cos(kx)cos(ky)sin(kz)`.
Every component varies in every direction, so the wall-normal velocity varies in
the normal direction at all three interface orientations — the property Beltrami
lacks (`dv/dy=0`), and the term the interface gets wrong. It is NOT an NS
solution; it exists only for this operator gate. `lap(u_i) = -3k^2 u_i`, so the
analytic diffusion is trivial and the advection is closed-form from the 9
analytic first derivatives (see `tools/rhsband.py`).

## Run

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
cmake --build build_cpu -j            # -Mnofma reference
./validation/momentum_interface/run_gate.sh build_cpu/main /tmp/mi_gate
```

`rhsband.py --re 1e30` on a huge-Re dump isolates the **advection** operator
(diffusion ~0); the normal-Re dump is the combined operator.

## Baseline (the unfixed predictor, 2026-06-24, commit after `d2e743a`)

- **Interior: order ~2.0** everywhere (uniform, slab, patch) — instrument valid.
- **Fine band: order ~0 (O(1) inconsistent)**, rms *grows* with refinement
  (slab v: -0.27; patch all: -0.20→-0.45). **Dominantly the advection stencil**
  (huge-Re run ≈ full): a fine cell touching the interface reads a halo holding a
  plain-injected / prolonged *coarse-resolution* value, then applies a *fine*-metric
  advection/diffusion stencil across it → O(h_coarse)/h_fine = O(1) truncation.
  This is the momentum analog of the projection's `D·u_exact = O(1)` defect; it
  needs a composite advect/diffuse stencil (momentum-interface-todo piece 2),
  and the fine block predicting its own interface face (piece 1).
- **Coarse band, wall-normal velocity: order ~1** (slab v: 1.36→0.98); tangential
  u,w stay ~2.0. The coarse cell's normal-momentum flux reading the restricted
  fine data is the un-refluxed coarse-fine flux mismatch (Berger–Colella reflux,
  piece 3).
- **Axis 2: round-off** (slab/patch `Sum(vol*div(qs))` ≈ 1e-14) — the existing
  post-predictor sync already conserves global mass even for this
  normal-velocity-varying field. Any momentum-interface fix must keep this.

So: the predictor is **0th-order at the fine interface band, 1st-order on the
coarse-band normal velocity**, while global continuity is already clean. The fix
(step.f90) targets the fine-band advection/diffusion stencil first.

## Increment — tangential ghost blend (exchange, `comm.f90`)

A fine cell touching the interface reads its tangential-velocity halo from a
plain coarse **injection** placed at the coarse cell centre — 3/2 fine cells from
where the fine advection/diffusion stencil expects it — an O(1) interface
truncation (fine-band advection order ~0; diffusion rms diverged 1.33→2.6→5.1).
The pressure already solved the identical problem with a ghost **blend**
(`entry_blend`, ghost = (2·coarse + fine)/3 placed at the fine halo centre). This
increment extends that blend, in the exchange gather, to the velocity components
**tangential** to the face (`lDir(var)==0`): one elegant condition change, no
kernel edits, fixing advection AND diffusion together. The tangential halos feed
only the momentum stencil (never the divergence), so conservation is untouched;
the wall-**normal** component is left to the face-staggered increment.

(An earlier attempt, increment 1, corrected only the cell-centred *diffusion* via
interface-row Laplacian coefficients in `blocks.f90`; it was reverted in favour
of this blend, which subsumes it and also fixes advection. A wrong turn caught
along the way: correcting the FACE_FINE/coarse side regressed the coarse band —
the restricted halo already lands at the coarse centre.)

Result: tangential **u,w fine-band order 0 → ~0.9→0.6** (advection + diffusion;
held below a clean 1 by the still-broken wall-normal `v` contaminating the `uv`
cross-terms). Regressions all pass: non-interface bit-exact, coarse band
unchanged, Axis-2 round-off (2.8e-14), real-Beltrami conservation round-off,
projection div-free PROJONLY gate still 0.0, CPU==GPU bit-identical.

## Increment 3 — wall-normal deep-halo reconstruction (`reconstruct_normal_halo`)

The blend (inc 2) excluded the wall-**normal** velocity. Term-by-term
(`MOBY_TERMDUMP`) pinned its defect to exactly the wall-normal terms across the
interface — `adv_y = ∂(vv)/∂y` (0th order) and `dif_y = ∂²v/∂y²` (diverging
O(1/h)) — at an orientation-B fine block's **owned** interface face (j=1), which
reads the deep halo `v(0)` below it. Every tangential term (incl. the tangential
Laplacian) is already 2nd order. An exact-`v(0)` test confirmed an accurate ghost
restores both to 2nd order.

The velocity prolong fills `v(0)` with one coarse face value — O(h) inaccurate
*tangentially*, which the wall-normal stencils amplify (`/h`, `/h²`). Fix:
`step.f90` reconstructs `v(0)` by **quadratic extrapolation from the fine side**,
`v(0) = 3v(1) − 3v(2) + v(3)` (per direction, at `physLow == FACE_COARSE` faces),
called right before the predictor. It is tangentially accurate, **purely local**
(no cross-block coarse reads / exchange-ordering race), and the deep halo never
enters the divergence, so conservation is untouched.

Result: `adv_y` 0→**2.0**, `dif_y` diverging→**2.0** (term-by-term); the full
wall-normal `v` fine band goes from ~0/diverging to **~0.6**, at **parity with
u,w** — all three components consistent. Slab and 3D patch (corners/edges, all
orientations) symmetric. Regressions all pass: non-interface bit-exact (inert
without an interface), Axis-2 round-off (2.8e-14 slab, 1.1e-14 patch),
real-Beltrami conservation round-off, projection PROJONLY gate unaffected (the
reconstruction runs only in the predictor path), CPU==GPU bit-identical.

## Increment 4 — tangential deep-halo reconstruction (`reconstruct_interface_halos`)

Inc 3 left the tangential cross-advection terms (`adv_x = ∂(vu)/∂x`,
`adv_z = ∂(vw)/∂z`, ~1st order) carrying the *tangential*-velocity halos'
tangential-injection residual: the blend (inc 2) places them correctly NORMAL to
the face but piecewise-constant TANGENTIALLY (one coarse value across the covered
fine cells = O(h)), so a tangential derivative of that ghost converges only ~0.6.

Fix: generalise inc 3's local cubic fine-side extrapolation
(`q(0)=3q(1)-3q(2)+q(3)`) from the wall-normal component to ALL THREE velocity
components in every 2:1 interface deep-halo row, at BOTH orientations, over the
full halo PLANE including the in-plane halo ring (`0..n+1`). The ring matters: the
interface-face cross-advection reaches the neighbouring tangential halo column
(`v`'s `∂(vu)/∂x` at `i` reads `u(i+1,0,k)`, so the edge cell `i=nx` needs
`u(nx+1,0,k)` reconstructed — without the ring, slab `v` stuck at ~0.5 from that
single edge cell per row). The three orientations run in order x,y,z so an
edge/corner reads the already-reconstructed column of the earlier plane. The
reconstruction reads the fine column at the same in-plane index, so it is
tangentially accurate; all reconstructed cells are deep halos that never enter
the divergence, so conservation is untouched (the OWNED interface normal face is
left alone). It supersedes the tangential ghost blend's role at these halos
(the blend still fills cells the reconstruction's per-orientation loops skip).

Result (per-term, v-momentum, slab fine band): `adv_x`/`adv_z`
`6.9e-3 → 1.7e-4`, order **0.6 → 2.00**; every term (`adv_x/y/z`, `dif_x/y/z`)
now 2nd order. **Slab fine-band order 0.6 → ~2.0 for all three components**
(u,w 1.97→1.86; v 1.99→2.00). Patch (corners/edges) fine band `0.6 → 1.76→1.36`
(the residual is the corner double-extrapolation; the coarse band there is also
only ~1.5). Regressions all pass: non-interface bit-exact (0.0), Axis-2 round-off
(slab +2.7e-14, patch +1.1e-14), PROJONLY div-free Beltrami 0.0 everywhere,
real-Beltrami conservation round-off (divpre +5.3e-17), coarse band unchanged,
CPU==GPU bit-identical (slab + patch).

## Increment 5 — coarse-side normal deep-halo reconstruction (`reconstruct_interface_halos`)

Layer 2: lift the coarse band from ~1st to 2nd order. The defect (pinned
term-by-term to the coarse cell adjacent to a `physLow==FACE_FINE` interface,
i.e. coarse-above-fine) is `adv_y` (order 0.99) and `dif_y` (order ~0, constant
error): the coarse cell advects / diffuses its interface face by reading the deep
halo `q(i,0,k)` one coarse cell INTO the fine region, which the exchange fills
with the RESTRICTION — a 4-point average of the covering fine faces. A
face-average differs from the point value the coarse stencil wants by O(h²) (the
tangential curvature), so the pointwise `∂(vv)/∂y` is O(h) (~1st) and `∂²v/∂y²`
is O(1) (~0th). This is the un-refluxed coarse-fine flux mismatch (Berger–Colella
"the fine flux is the accurate one").

Fix: the SAME local cubic extrapolation from the coarse interior
(`q(0)=3q(1)-3q(2)+q(3)`), applied to the normal component's deep halo at
`physLow(d)==FACE_FINE` faces. It gives a point-accurate, tangentially-accurate
ghost for the momentum stencil; the face-average stays in the OWNED interface
face `q(1)` (in the divergence) for mass conservation, and only the deep halo
`q(0)` — never in the divergence — is reconstructed. It is the Berger–Colella
fine-authoritative idea realized as a LOCAL reconstruction (codebase principle:
prefer local reconstructions over the racing gather), conservation-neutral and
vanishing for uniform/linear flow. The other orientation's coarse cell reads the
interface face directly and is already 2nd order, so only the FACE_FINE low face
is treated.

Result (per-term, v-momentum, coarse cell adjacent to interface): `adv_y`
0.99 → **3.00**, `dif_y` -0.02 → **2.00**; all coarse-band terms 2nd–3rd order.
**Slab coarse-band v order 0.98 → 2.98**; patch coarse band `1.47 → 2.04` (all
components). Fine band unchanged (slab ~2.0, patch 1.76→1.36). Regressions all
pass: non-interface bit-exact (0.0); Axis-2 round-off (slab +2.5e-14, patch
+1.2e-14); PROJONLY div-free Beltrami 0.0 everywhere; real-Beltrami conservation
round-off (divpre +5.5e-17); fine band + interior unchanged; CPU==GPU
bit-identical (slab + patch).

**Status:** the 2:1 momentum predictor is now **2nd order at the interface** for
all components, both bands, both orientations (clean on the slab; ~1.4–2.0 on the
patch, the residual being corner double-extrapolation). Mass conservation stays
round-off. The remaining momentum-CONSERVATION (flux-register equality, vs the
accuracy this increment delivers) is not gated here and is a separate property.
