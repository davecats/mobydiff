# Momentum predictor at the 2:1 interface — handout

Branch **`claude/jacobi-interface`**, tag **`jacobi-momentum-interface`** (HEAD
`e261df2`; the code milestone is `a5fc325`, inc 3). Read this first, then the
memory `momentum-interface-todo`, then
`docs/jacobi_interface_handout.md` (the projection work this builds on) and
`docs/interface_review.md` §iii (the uniform-B / reflux analysis).

## Where we are

The **projection** (pressure correction) was already divergence-consistent and
mass-conserving to round-off (damped Jacobi, both orientations, sync,
interface-row phi restrict — see the projection handout). This work fixed the
**momentum predictor** (advection + diffusion, `src/modules/step.f90`) at the
interface, which was both unfixed and untested.

**The momentum interface operator is now CONSISTENT (converging) everywhere** —
no 0th-order or diverging terms remain. Two accuracy layers are left (both push
a residual from ~1st to 2nd order); they are NOT consistency defects.

### Commits (oldest → newest)
- `05182d2` **operator-truncation gate**: `initial = tgv3d` manufactured field +
  `MOBY_RHSDUMP` + `tools/rhsband.py` + `tools/divsum.py` +
  `validation/momentum_interface/{run_gate.sh,README.md}`.
- `4960507` inc 1 cell-centred interface diffusion via Laplacian coeffs —
  **superseded** by inc 2 (kept in history; blocks.f90 is back to pre-inc-1).
- `88b0fae` inc 2 **tangential ghost blend** (comm.f90): extends the pressure
  `entry_blend` to the velocity components tangential to the face. Fixes u,w.
- `c83c3c3` **per-term operator dump** `MOBY_TERMDUMP` + `tools/rhsterms.py`.
- `a5fc325` inc 3 **wall-normal deep-halo reconstruction** (step.f90). Fixes v.

### Current order (slab, 32/64/128, with all increments incl. inc 4)
```
interior     2.0   (all components)            <- 2nd order
coarse_band  u,w 2.0   v 1.0                    <- v NOT 2nd (reflux = Layer 2)
fine_band    u,v,w ~2.0                         <- 2nd order (inc 4 done)
```
Per-term (v-momentum, fine band): ALL terms (`adv_x/y/z`, `dif_x/y/z`) = 2.0.
(Inc 4 lifted `adv_x`/`adv_z` from ~0.6 to 2.00.) Patch (corners/edges) fine
band is ~1.4–1.8 — the residual is the corner double-extrapolation, expected.

## The gate suite (how to validate)

All gates are single-step operator tests (NOT time evolution — Beltrami/tgv3d
are run for at most ONE step; see memory `beltrami-stability-scope`). The field
`initial = tgv3d` is a manufactured field where every velocity component varies
in every direction, so the wall-normal velocity varies in the normal direction
at all three interface orientations (Beltrami's `∂v/∂y=0` lacks this).

Build (CPU is the `-Mnofma` reference; GPU must stay bit-identical):
```
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
cmake --build build_cpu -j && cmake --build build_gpu -j
```
One-shot driver (generates inputs, runs both axes over 32/64/128):
```
./validation/momentum_interface/run_gate.sh build_cpu/main /tmp/mi_gate
```
NOTE the `/tmp` scratchpad is swept periodically — regenerate inputs with
`run_gate.sh` rather than relying on leftover `.ini` files.

### Axis 1 — accuracy / convergence order  (`MOBY_RHSDUMP` → `rhsband.py`)
Dumps the predictor's discrete RHS `L_h(u) = -div(uu) + (1/Re)lap(u)` (the
`blk%oldrhs` captured at the first substage on the pristine field) as `un/vn/wn`
of `<prefix>_rhs`. `rhsband.py` compares to the analytic `L(u)` at the staggered
points and prints the order over 32/64/128, split interior / coarse-band /
fine-band, per component.
```
MOBY_PREDONLY=1 MOBY_RHSDUMP=1 mpirun -x MOBY_PREDONLY -x MOBY_RHSDUMP -n 1 \
    build_cpu/main slab_<n>.ini          # PREDONLY: predictor only, no projection
python3 tools/rhsband.py r32_rhs.h5 r64_rhs.h5 r128_rhs.h5
```
- interior must be ~2 (instrument sanity); the **fine-band order is the verdict**.
- `--re 1e30` on a huge-Re dump isolates the **advection** operator (diffusion ~0).

### THE PER-TERM TEST — attribute a defect to one term  (`MOBY_TERMDUMP` → `rhsterms.py`)
This is the decisive diagnostic: it recomputes, from the post-IC-exchanged field
exactly as the predictor would, **each momentum term of a chosen component
SEPARATELY** — the three advection flux-divergence terms (`adv_x/adv_y/adv_z`)
as `un/vn/wn` of `<prefix>_adv`, and the three Laplacian terms (`dif_x/dif_y/
dif_z`, times 1/Re) as `un/vn/wn` of `<prefix>_dif`. The files are written at
step `_0` (the hook runs before the loop, then `stop`s).
```
MOBY_TERMDUMP=2 mpirun -x MOBY_TERMDUMP -n 1 build_cpu/main slab_<n>.ini   # 2 = v
python3 tools/rhsterms.py P32_adv_0.h5 P32_dif_0.h5 P64_adv_0.h5 P64_dif_0.h5 \
        P128_adv_0.h5 P128_dif_0.h5 --var 2          # order per term
python3 tools/rhsterms.py P64_adv_0.h5 P64_dif_0.h5 --var 2 --rows          # per y-row
```
`--rows` prints each term's error per y-row of the orientation-B fine blocks,
which localizes a defect to the interface FACE (j=1) vs interior rows. This is
how the v defect was pinned to exactly `adv_y` (0th) and `dif_y` (diverging)
while every tangential term was clean — and how to verify a fix term-by-term.
The hook reconstructs the wall-normal halo first (matches the predictor).

### Axis 2 — continuity null-space  (`MOBY_DIVDUMP` → `divsum.py`)
The projection preserves the MEAN of `div(qs)`, so the global `Sum(vol*div(qs))`
after predictor+sync (the net 2:1 interface flux mismatch) is the one continuity
error it can never remove — it MUST stay at round-off. Run a REAL step (no
PREDONLY/PROJONLY) with `MOBY_DIVDUMP` and sum the `_divpre` dump:
```
MOBY_DIVDUMP=1 mpirun -x MOBY_DIVDUMP -n 1 build_cpu/main slab_64.ini
python3 tools/divsum.py <prefix>_divpre_1.h5     # must be ~1e-14
```

### Mandatory regressions (run after every change)
- **non-interface bit-exact**: a single-level (no `[blocks]`/`refine`) run must be
  byte-identical (interface fixes must be inert without `FACE_COARSE`).
- **Axis-2 round-off** on slab and 3D patch (`refine = 1.6 4.7 ...` all dirs).
- **projection div-free PROJONLY gate**: `MOBY_PROJONLY=1 MOBY_DIVDUMP=1` on a
  Beltrami slab — `divpost` fine-band must stay 0.0 (the projection consistency
  must not regress). Reconstruction runs only under `.not. projOnly`, so it is
  inert here by construction, but check it.
- **real-Beltrami 1-step**: conservation round-off + no blow-up + finite fields.
- **CPU == GPU bit-identical** on the slab `_rhs` dump.

## What is left (the next session's work)

Layer 1 (tangential interp) is **DONE** (inc 4 below); the fine band is 2nd
order. Only Layer 2 (reflux) remains. The operator is consistent everywhere.
Detailed scoping is in `docs/next_session_tangential_reflux.md`.

1. **Tangential reconstruction of the tangential velocity halos** — **DONE**
   (increment 4, `reconstruct_interface_halos` in `step.f90`). Generalised inc
   3's local cubic fine-side extrapolation `q(0)=3q(1)-3q(2)+q(3)` from the
   wall-normal component to ALL THREE velocity components in every interface
   deep-halo row, BOTH orientations, over the full halo PLANE incl. the in-plane
   ring `0..n+1` (the ring is load-bearing: `v`'s `∂(vu)/∂x` at the interface
   face reaches `u(nx+1,0,k)`; without it slab `v` stuck at ~0.5 from one edge
   cell per row). Three orientations run x,y,z so an edge/corner reads the
   earlier plane's reconstructed column. Purely local (no gather race), deep
   halos only (conservation-neutral). RESULT: slab fine-band order 0.6 → ~2.0
   all components (`adv_x`/`adv_z` 0.6 → 2.00 per-term); patch 0.6 → 1.4–1.8
   (corner residual). All regressions pass (see README inc 4).

2. **Berger–Colella tangential-momentum reflux** → lifts the coarse band to 2nd
   order. The coarse cell adjacent to the interface advects its NORMAL velocity
   using the restricted fine flux; the coarse-minus-summed-fine flux mismatch is
   un-refluxed (coarse-band normal v ~1st order; tangential u,w already 2nd).
   This is the only piece that restores MOMENTUM conservation across the
   interface (mass is already conserved). "Constructed to vanish for uniform
   flow" is the constraint to keep the exact gates. See `docs/interface_review.md`
   §v (Berger & Colella 1989; Almgren et al. 1998).

## Gotchas (carried over)
- CPU `build_cpu` (`-Mnofma`) is the reference; GPU `build_gpu`
  (`-Mnofma -gpu=nofma`) must stay bit-identical.
- Run cases ONE AT A TIME; never `pkill -f build_cpu/main` while a run is live
  (kills your own mpirun — use the `[b]uild` bracket trick if you must).
- `/tmp` scratchpad is swept — regenerate `.ini` with `run_gate.sh`.
- Single-block / single-level must stay bit-exact.
- Beltrami/tgv3d are SINGLE-STEP operator tests; long-time stability is NOT a
  target (`beltrami-stability-scope`).
- The deep halo `q(0)` (below an owned interface face) never enters the
  divergence — reconstructing it is conservation-neutral. The OWNED interface
  face `q(1)`/`q(nb+1)` does — leave its fine-owns handling alone.
