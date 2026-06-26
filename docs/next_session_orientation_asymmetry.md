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

## MECHANISM CONFIRMED (2026-06-26c): the MOMENTUM REFLUX creates the bands.

Decisive toggle on the benchmark (250 steps, const-1/2, GPU): momentum_reflux
OFF vs ON, per-interface bands from the final field (rms_rows):

| interface | u' excess ON->OFF | v' kink ON->OFF |
|---|---|---|
| bottom (fine-owns) | x1.56 -> **x1.00** | 0.238 -> **0.025** |
| top (coarse-owns)  | x1.45 -> **x1.02** | 0.306 -> **0.040** |

Reflux OFF essentially ERASES both the u' streak spike and the v' normal step,
at NO cost to stability or divergence (div_l2 ~2.3e-3, div_max ~0.089, mass
round-off -- IDENTICAL to reflux ON). The bands are a momentum-reflux artifact,
NOT the interface transfer/ownership and NOT an energy (KEBAL) defect (KEBAL is
inert on the reflux because the band is momentum structure, not skew-energy).

WHY (step.f90 reflux_accumulate:634, reflux_compute_flux:478): the reflux replaces
the coarse interface-cell flux F_coarse by the restricted fine flux avg(F_fine).
For the normal component F=(q_halo+q_int)^2, so avg(F_fine)=avg(of squares) >>
(avg)^2 ~ F_coarse (Jensen gap) whenever the fine side carries resolved
fluctuations. The reflux thus INJECTS the fine-side resolved Reynolds-stress flux
into the under-resolved coarse interface cell -> it piles up as the u'/v' band
(exactly Cevheri & Stoesser's "energy accumulation at grid coarsening"). The
reflux conserves the MEAN interface flux (its design purpose, the -<u'v'> defect,
validation/momentum_interface) but pumps the FLUCTUATING stress onto the coarse
cell. The reflux was built/validated with the OLD recon-based interface (memory
momentum-interface-todo); under const-1/2 recon is OFF and the reflux's
band-creating side-effect dominates with no offsetting benefit measured.

The orientation asymmetry (w' clean at fine-owns / kink at coarse-owns) is the
reflux acting only on the coarse side (physLow/High==FACE_FINE) -- the coarse cell
that receives avg(F_fine) is the low block at the top interface, the high at the
bottom; the asymmetry is which side carries the resolved-stress injection.

## Do this, in order

1. THE TRADE (decisive, needs stats): run the developed stats with reflux OFF
   (run_developed.py with momentum_reflux=false) vs ON vs the uniform-fine
   reference. Reflux exists for the MEAN -<u'v'> interface flux; confirm whether
   reflux-OFF degrades the mean profile / Reynolds shear at the interface (the
   thing reflux conserves) or whether it is simply better here. ~14h GPU.
2. FIX (if reflux-off harms the mean): make the reflux transfer the MEAN flux only,
   not the fluctuating part -- e.g. reflux the x,z-averaged (homogeneous-direction)
   flux mismatch, or a single-valued constant-1/2 flux, so it conserves -<u'v'>
   without injecting the fine-side resolved stress into the coarse cell. Gate on
   the benchmark (u'/v' band NOT regrown, stable) + the mean profile via stats.
   Investigate the 0.25 scaling / Jensen gap in reflux_accumulate for correctness.
3. GATE every change: interface_benchmark stability PASS + band metrics (u'/v')
   tracked; KEBAL slab (per-component); bit-exact no-interface; 1-rank==2-rank.
   The decisive metric (converged mean + Reynolds stresses) needs the 14h stats.

The earlier adjoint-normal-transfer and fine-owns-ownership leads are now
SECONDARY: with reflux off the residual bands are ~1.0-1.02 excess / 0.025-0.040
kink (near zero). Confirm the reflux is the whole story (converged stats) before
chasing them.

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
