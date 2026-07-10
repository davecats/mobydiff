# IDDES T5 gates — full IDDES hybrid (SST + WALE)

Physics gates for IDDES phase T5 (docs/next_session_iddes.md).
Increment 1 was the DDES shielding; increment 2 (current) is the full
IDDES elevating/WMLES branch (Gritskevich et al. 2012 SST-IDDES). The
blend lives in turbulence.f90, in the RANS-retention convention:
fd = max(fd_dt, f_B) with fd_dt = tanh((C_dt1 r_dt)³) on
r_dt = ν_t/(κ² y_eff² |∇u|) and the geometric f_B; the k-destruction
length l_hyb = fd·(1+f_e)·l_RANS + (1−fd)·C_DES·Δ enters D_k
point-implicitly (Δ = the IDDES wall-aware mesh length by default);
nut = fd·nut_rans + (1−fd)·nut_sgs (our validated WALE blend — textbook
SST-IDDES has no separate SGS model; the SST limiter plays that role).
Evaluation toggles: `[turbulence] iddes_cdt1` (20 vs the DDES 8),
`iddes_clip` (Spalart's max(0, l_RANS − l_LES) clipping), `iddes_delta`
(iddes | cbrt).

```bash
cd validation/iddes
./run_gates.sh                  # all gates, one job at a time
BIN=../../build_cpu/main RANKS=4 ./run_gates.sh iddes180
python3 check_gates.py          # needs python3 + h5py + numpy
```

Gate groups (see run_gates.sh):

- `iddes180`  gates (a)/(b): WMLES-style developed channel (Re_tau 180 on
              the LES-validation 64x48x64 grid, dx+~35/dz+~18/y+_1~0.5).
              Transient to t=5, statistics leg to t=25. Checks: f_d -> 1
              through the RANS wall layer / -> 0 in the core (final
              snapshot's fd dataset); mean profile has no gross log-layer
              mismatch vs the pure-WALE stats reference
              (../channel_interface/les/runs/uniform/stats) and the T2
              RANS turb180 profile.
- `fd0`       gate (c) SGS limit: [turbulence] fd_force = 0 makes the
              blend nut = nut_sgs EXACTLY (0*x + 1*y), so 20 steps must be
              bit-exact vs the pure-WALE twin on the same IC.
- `fd1`       gate (c) RANS limit: fd_force = 1 on the turb180 grid,
              restarting from the converged T2 field — the k-sink becomes
              sqrt(k)/l_RANS, the same fixed point as beta* omega k through
              different arithmetic, so the profile must hold the turb180
              answer to round-off-class drift.
- `ibm`       gate (d): the les_ibm off-grid IBM channel (ibm180
              conditions) runs model = iddes stably through the IBM wall
              treatment. Needs ../rans_geometry/ibm_coeff_blocks_l1.h5.
- `ranks`     gate (e): iddes-on 20 steps, 1 rank == 4 ranks EXACT.
- `toggles`   evaluation only (no pass/fail): t=5..10 stats legs from the
              SAME default transient for default vs iddes_cdt1=8 vs
              iddes_clip=true vs iddes_delta=cbrt; check_gates.py prints
              log-layer deviations + fd band means side by side.

The model /= iddes full bit-exactness gate (e) runs OUTSIDE this
directory: rebuild both sides -Mnofma / -gpu=nofma and compare the
standard short list (min_channel, les_ibm ± refine_body, Beltrami
y-slab, turb180, wf180_y30, lam30t) with tools/compare_fields.py.

## RESULTS 2026-07-10, increment 2 (full IDDES elevating branch —
## all gates PASS; local GPU runs, nofma builds)

- (a) fd(y+) on the final t = 25 snapshot: 1.000 EXACTLY through y+ ≤ 15
  (the geometric f_B guarantees retention out to d_w ≈ 0.53 h_max = y+ 18.6
  on this grid), band means 0.872 at y+ 5–25 / 0.034 at 25–60 (DDES
  increment: 0.67 / 0.10 — the handover moved outward as required), core
  max 1.8e-4.
- (b) log-layer mean U (full t = 5..25 stats): max dev 0.5% vs the
  pure-WALE stats reference and 0.7% vs the T2 RANS turb180 profile
  (DDES increment: 3.0%/2.8% — the elevating branch bought a ~5x
  log-layer improvement on top of the wider RANS coverage). Centreline
  U = 17.78.
- toggle evaluation (t = 5..10 stats legs from the same transient;
  log-layer max dev vs WALE / vs RANS-T2; fd means y+ 5–25 / 25–60):
  default (C_dt1 = 20, no clip, IDDES Δ): 0.004/0.008, 0.871/0.035;
  iddes_cdt1 = 8: 0.005/0.010, 0.870/0.023 (slightly narrower shield,
  marginally worse); iddes_clip = true: 0.004/0.008, 0.871/0.035
  (indistinguishable from default — the clipping never binds on this
  case); iddes_delta = cbrt: 0.013/0.005, 0.906/0.066 (thicker fd tail,
  3x worse vs WALE). VERDICT: keep the Gritskevich defaults (C_dt1 = 20,
  plain convex blend, IDDES Δ).
- (c) fd_force = 0: bit-exact vs pure WALE (max_abs 0 on u/v/w/p/nut,
  20 steps; fe is zeroed with the force, so the identity still holds).
  fd_force = 1 from the converged turb180 field: drift after 2000 steps
  max|du|/u_max = 9.8e-13, dk/k_max = 2.2e-12 (fe = 0 under force makes
  l_hyb = l_RANS exactly).
- (d) iddes_ibm 2000 steps: finite, bounded (max|u| 25.4, nut ≤ 1.9e-2),
  fd in [0, 1].
- (e) 1 rank == 4 ranks EXACT (max_abs 0 incl. k/omega/nut/fd); CPU vs
  GPU 20 steps EXACT (max_abs 0 on every field — the increment-1 tanh
  ulp gap does not reappear; the exp() geometric fields are host-side).
  model /= iddes bit-exact vs the pre-increment baseline (nofma,
  max_abs 0 incl. k/omega/nut/gamma/rethetat) on min_channel 4-rank,
  les_ibm ± refine_body, Beltrami y-slab, turb180, wf180_y30, lam30t —
  CPU AND GPU.
- NOTE: the tog_clip leg took 2.8 h wall for a ~5 min job with healthy
  fields throughout — host suspend/GPU throttle during the run, not a
  solver issue (the sibling legs ran at 0.029 s/step).

## RESULTS 2026-07-10, increment 1 (DDES shielding — HISTORICAL; the
## elevating branch above supersedes these fd numbers)

CONVENTION (found the hard way): the stored `fd` is the RANS-RETENTION
weight tanh((8 r_d)^3) = 1 − f_d^Spalart. Implementing Spalart's
1 − tanh(...) verbatim with the doc's blend hands the WALL layer to the
SGS model: measured fd(wall) = 0 and a +16% log-layer error.

- (a) fd(y+): 1.000 below y+ 5, handover through the buffer (mean 0.67
  at y+ 5–25, 0.10 at 25–60), < 0.001 in the core. The early handover is
  the known DDES-with-resolved-content behaviour; extending RANS coverage
  outward is the elevating branch's purpose (increment 2).
- (b) log-layer mean U (full t = 5..25 stats): max dev 3.0% vs the
  pure-WALE stats reference, 2.8% vs the T2 RANS turb180 profile
  (bar: 15%; the intermediate t = 5..10 window read 1.8%/1.9%).
- (c) fd_force = 0: bit-exact vs pure WALE (max_abs 0 on u/v/w/p/nut,
  20 steps) — the blend identity 0·nut_rans + 1·nut_sgs is exact.
  fd_force = 1 from the converged turb180 field: drift after 2000 steps
  max|du|/u_max = 9.8e-13, dk/k_max = 2.2e-12 (round-off-class only).
- (d) iddes_ibm 2000 steps: finite, bounded (max|u| 24.9, nut ≤ 2.1e-2),
  fd in [0, 1].
- (e) 1 rank == 4 ranks EXACT (max_abs 0 incl. k/omega/nut/fd); CPU vs
  GPU 20 steps ≤ 1.8e-14 (u), 1.4e-14 (fd) — intrinsic-ulp class.
  model /= iddes bit-exact vs the post-STEP-0 baseline (nofma, max_abs 0)
  on min_channel, les_ibm ± refine_body, Beltrami y-slab, turb180,
  wf180_y30, lam30t — CPU AND GPU.
