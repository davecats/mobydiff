# Next session — the fine-owns/coarse-owns interface ORIENTATION asymmetry

Branch `claude/jacobi-interface`. Supersedes `docs/next_session_v_normal.md` (its
adjoint-normal-transfer premise is DISPROVEN below). Read
`docs/interface_band_handout.md` and memory `interface-validation-suite` for the
full lineage.

## What this session settled (2026-06-26c)

Davide ran the converged developed stats (t=5..25, 2 GPU); results in
`validation/channel_interface/developed/runs/{default,skew}`. Combined with a new
per-component KEBAL split, the picture is now decisive and REDIRECTS the work.

### 1. The residual interface ENERGY defect is TANGENTIAL u, not normal v.

Added a per-component band-SKEW split to `print_step_ke_balance` (main.f90,
`MOBY_KEBAL`, diagnostic-only, production bit-exact). It prints
`SKEW band u= .. v= .. w= ..`.

- **Clean Beltrami slab** (`slab_y_diag.ini`, y-band, V = interface-normal):
  | config | u (tang) | v (NORMAL) | w (tang) | total |
  |---|---|---|---|---|
  | old baseline (metric/cubic) | +0.149 | +0.073 | -0.095 | 0.127 |
  | **const-1/2 default** | +0.150 | **+0.0095** | -0.097 | 0.062 |
  The const-1/2 default already cut the interface-NORMAL v band-SKEW by 87%
  (0.073->0.0095). The residual 0.062 is the TANGENTIAL u (+0.15) + w (-0.097).
- **Turbulent channel** (1 step, GPU): u = **-150**, v = +1.55, w = +0.34. The
  interface-normal v is ~1% of the streamwise u. (Channel SKEW is stretching-
  contaminated: interior +145 -- the per-component RATIO is the signal, not the
  absolute.)

CONSEQUENCE: the planned "volume-weighted ADJOINT transfer for the interface-NORMAL
component" (next_session_v_normal.md) has ~no energy headroom. Mathematically the
current normal pair (injection-prolong + average-restrict) IS already the
volume-weighted adjoint pair for face areas, and the const-1/2 default made it
nearly energy-consistent. Do NOT build it. (The dominant u-band is the documented
low-k, predictor-advection-driven, transfer-IMMUNE intrinsic resolution-jump band;
4+ prior sessions ruled out transfer/filter fixes for it.)

### 2. The v' kink is REAL (survives time-averaging) and there is a clean
### ORIENTATION asymmetry -- the new, sharper lead.

Converged time-averaged rms at the two wall-band interfaces (default run,
`runs/default/stats/channel_stats.h5` merged with `_l1`):

| interface | who OWNS the shared face | u' (across) | w' (across) |
|---|---|---|---|
| **bottom** y~0.643 | FINE (low side) | 1.61 -> **2.37** -> 1.55 | 0.859 -> 0.828 -> 0.784  (**smooth, monotone**) |
| **top** y~1.357 | COARSE (low side) | 1.75 -> **2.68** -> 1.75 | 0.874 -> **0.942** -> 0.942  (**kink UP**) |

w' has NO mean shear, so a top-only w' kink is a pure ORIENTATION-dependent
artifact (the two walls are physically symmetric; only the interface treatment
breaks top/bottom symmetry). The single orientation-dependent treatment is the
Phase-3c **low-block-owns-face**: at the bottom (fine below coarse) the FINE block
owns the shared 2:1 face (good); at the top (fine above coarse) the COARSE block
owns it (the defect). v' carries the same asymmetry (its kink 0.24-0.31; the
benchmark's tracked metric).

This is exactly the documented Phase-3c DEVIATION (CLAUDE.md, §6): the original
doc-6 design wanted UNCONDITIONAL fine-owns-face, but it was "unrealizable in this
storage convention (fine face DOFs in unpredictable halos; restriction would write
the coarse INTERIOR u(1) plane against the prolong reading it)", so they settled for
low-owns. We now SEE its cost: the coarse-owns orientation pollutes w' and v'.

## Mechanism (working model -- VERIFY before fixing)

Under the const-1/2 default, `reconstruct_interface_halos` is OFF (main.f90:340).
So the asymmetry is in the exchange transfer + reflux + the post-predictor
`exchange_halos(..., syncface=.true.)` (main.f90:394) that writes the low-owned
shared face. `interface_normal_dim` (comm.f90:741) skips the cross-level
PROLONG/RESTRICT only for off(d)==+1 (a block keeping its own HIGH face q(nb+1)),
making the LOW block the owner. At a coarse-owns face the shared normal v is
predicted at COARSE resolution and the fine side is slaved (PROLONG injection) ->
the fine cells adjacent cannot be div-free at fine resolution -> the projection
compensates with spurious tangential (w') motion. At a fine-owns face both sides
see a fine-resolution face and stay consistent.

## Do this, in order

1. CONFIRM the mechanism cheaply before any big change: localize whether the top
   (coarse-owns) w'/v' kink comes from (a) the normal-v PROLONG injection to the
   fine halo, (b) the syncFace, or (c) the reflux. The slab has BOTH orientations
   (middle refined -> bottom = coarse-owns, top = fine-owns -- OPPOSITE of the
   channel); per-component KEBAL + a per-orientation roughness probe on the slab is
   the fast (~5s) testbed. The benchmark (~2.5min GPU) gives stability + early u'/v'
   banding but NOT the converged w' kink (only 14h stats show that) -- weigh this.
2. FIX APPROACH (Davide's call -- two very different efforts):
   - **Targeted**: make the normal-v PROLONG to the fine halo at a COARSE-owns face
     a fine-resolution reconstruction (not coarse injection), so the fine side near
     a coarse-owned face is not slaved to one coarse value. Smaller, testable on the
     benchmark for stability; converged w' only via stats.
   - **Architectural**: restore UNCONDITIONAL fine-owns-face (reverse the 3c
     deviation). Principled and symmetric, but the documented hard problem ("fine
     face DOFs in unpredictable halos"). High risk (many prior interface changes
     blew up); needs the full gate suite + stats.
3. GATE every change: interface_benchmark stability PASS + banding not worse; KEBAL
   slab (now per-component); bit-exact no-interface; 1-rank==2-rank (dims N 1 1).
   The decisive metric (converged w'/v' asymmetry) needs the 14h developed stats
   (run_developed.py --ranks 2) + the existing uniform-fine reference.

## State / artifacts

- Per-component KEBAL committed (main.f90). `MOBY_KEBAL=1` on slab_y_diag.ini.
- Converged stats in `validation/channel_interface/developed/runs/{default,skew}`
  (default_stats.png, skew_stats.png). default==skew statistically (KESKEW inert
  on the asymmetry, as expected -- it's not an energy defect).
- Benchmark current default: STABLE, u' excess x1.56, v' kink 0.306.
- Transfer code map: comm.f90 `copy_local_entries` (1281, the weighted gather;
  velInject toggle 1293; owned-face skip 1393-1395), `entry_gather_map` (684),
  `interface_normal_dim` (741), `iface_restrict_normal` (756), `entry_blend` (776).
  step.f90 `reconstruct_interface_halos` (301, OFF under const-1/2),
  `reflux_compute_flux`/`reflux_accumulate`/`reflux_apply`,
  `skew_interface_correction`. main.f90 predictor 332-403, syncFace 394.
