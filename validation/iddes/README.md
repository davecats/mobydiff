# IDDES T5 gates — DDES-shielding hybrid (SST + WALE)

Physics gates for IDDES phase T5, first increment: the DDES blend
(docs/next_session_iddes.md). The blend lives in turbulence.f90:
f_d = 1 − tanh((8 r_d)³) shields the RANS near-wall layer; the
k-destruction length l_hyb = f_d·l_RANS + (1−f_d)·C_DES·Δ enters D_k
point-implicitly; nut = f_d·nut_rans + (1−f_d)·nut_sgs.

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

The model /= iddes full bit-exactness gate (e) runs OUTSIDE this
directory: rebuild both sides -Mnofma / -gpu=nofma and compare the
standard short list (min_channel, les_ibm ± refine_body, Beltrami
y-slab, turb180, wf180_y30, lam30t) with tools/compare_fields.py.

The full IDDES f_B/f_e elevating branch is a SEPARATE second increment,
gated on log-layer-mismatch reduction on the same iddes180 case.

## RESULTS 2026-07-10 (all gates PASS; local GPU runs, nofma builds)

CONVENTION (found the hard way): the stored `fd` is the RANS-RETENTION
weight tanh((8 r_d)^3) = 1 − f_d^Spalart. Implementing Spalart's
1 − tanh(...) verbatim with the doc's blend hands the WALL layer to the
SGS model: measured fd(wall) = 0 and a +16% log-layer error.

- (a) fd(y+): 1.000 below y+ 5, handover through the buffer (mean 0.67
  at y+ 5–25, 0.10 at 25–60), < 0.001 in the core. The early handover is
  the known DDES-with-resolved-content behaviour; extending RANS coverage
  outward is the elevating branch's purpose (increment 2).
- (b) log-layer mean U (t = 5..25 stats): max dev 1.8% vs the pure-WALE
  stats reference, 1.9% vs the T2 RANS turb180 profile (bar: 15%).
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
