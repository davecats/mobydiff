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

**Still open:** the wall-**normal** velocity (its interface face is fine-owned and
its deep halo is face-staggered — needs a face-staggered-aware blend / the fine
predicting its own face, and tangential interpolation to reach a clean 2nd order
rather than the current ~1st-order injection residual).
