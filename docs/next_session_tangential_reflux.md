# Next session — finish the 2:1 momentum interface to clean 2nd order

Branch **`claude/jacobi-interface`**, tag **`jacobi-momentum-interface`** (HEAD
`a5fc325`, pushed). The momentum predictor interface is now **consistent
(converging) everywhere** — no 0th-order or diverging terms. Two accuracy layers
remain, each lifting a residual from ~1st to 2nd order. They are independent;
do them in either order, gating each.

## Read first
- `docs/momentum_interface_handout.md` — full state + the GATE SUITE (the
  per-term test `MOBY_TERMDUMP`/`rhsterms.py`, the order gate
  `MOBY_RHSDUMP`/`rhsband.py`, the continuity gate `MOBY_DIVDUMP`/`divsum.py`,
  the regressions). BUILD/RUN commands are there.
- memory `momentum-interface-todo` (the three-piece map + this session's results)
  and `beltrami-stability-scope` (single-step only).
- `docs/interface_review.md` §iii and §v (uniform-B + Berger–Colella reflux).

## Settled (do not re-litigate)
- The fine-band defect is now ONLY the tangential cross-advection terms
  (`adv_x`/`adv_z` ~0.6); `adv_y`, all `dif_*` are 2nd order. Verified
  term-by-term (`tools/rhsterms.py --var 2 --rows`).
- The coarse-band defect is ONLY the wall-normal velocity (~1st); coarse-band
  tangential is 2nd order.
- The wall-normal v fix is a LOCAL fine-side reconstruction
  (`reconstruct_normal_halo`, step.f90, `q(0)=3q(1)-3q(2)+q(3)`), chosen over a
  coarse-bilinear-in-the-gather because the latter races on coarse halos in the
  shared exchange. Prefer local reconstructions for the same reason.
- Conservation: the deep halo `q(0)` is not in the divergence (reconstructing it
  is conservation-neutral). The owned interface face is — leave it alone.

## Task — pick a layer, gate it, then the other

### Layer 1: tangential interpolation of the tangential velocity halos (fine band → 2nd)
The tangential cross-advection `adv_x = ∂(vu)/∂x`, `adv_z = ∂(vw)/∂z` read the
tangential velocity deep halos (e.g. `u(i,0,k)` at a y-interface). Increment 2's
blend places them correctly NORMAL-to-the-face but leaves them piecewise-constant
TANGENTIALLY (one coarse value across the fine cells = O(h)); a tangential
derivative of that converges ~0.6. Make those halos tangentially smooth.
- Preferred: a LOCAL reconstruction (avoid the cross-block gather race — see the
  handout). The tangential halo at a fine cell could be reconstructed from the
  injected coarse value + the fine interior tangential structure, or via a
  carefully-ordered second gather pass that reads only filled coarse halos.
- Gate: `rhsterms.py --var 1/2/3` — `adv_x`/`adv_z` must climb toward 2nd; the
  `rhsband.py` fine band must rise from ~0.6 toward ~2. Keep Axis-2 round-off,
  non-interface bit-exact, projection PROJONLY 0.0, CPU==GPU.

### Layer 2: Berger–Colella tangential-momentum reflux (coarse band → 2nd)
The coarse cell adjacent to the interface advects its NORMAL velocity with the
restricted fine flux; the coarse-minus-summed-fine flux mismatch is un-refluxed
(coarse-band normal v ~1st). Add the flux-register reflux correction to the
adjoining coarse cell. This is the only step that restores MOMENTUM conservation
across the interface (mass already conserved). Construct it to VANISH for uniform
flow (preserves the exact uniform-flow/bit-exact gates).
- Gate: `rhsband.py` coarse-band normal velocity must rise from ~1st to ~2nd;
  uniform-flow exactness and all regressions preserved.

## Workflow
1. Build CPU+GPU. Run `validation/momentum_interface/run_gate.sh build_cpu/main
   /tmp/mi_gate` for the baseline order (confirm the numbers in the handout).
2. Implement one layer. Re-run the per-term gate to confirm the targeted terms
   improve; re-run ALL regressions (handout list).
3. Build GPU, confirm CPU==GPU bit-identical. Commit with the gate numbers.
4. Repeat for the other layer. Update the README + memory + this doc.

## Gotchas
- CPU `-Mnofma` reference; GPU bit-identical. One run at a time; never
  `pkill -f build_cpu/main` live. `/tmp` is swept — regenerate inputs.
- Single-step operator tests only; no long-time Beltrami.
- The shared gather feeds the projection's round-off conservation — touching it
  is high-risk; prefer local reconstructions and re-verify Axis-2 every time.
