# RANS T2 gates — k-omega SST transport, resolved wall mode

Physics gates for IDDES phase T2 (docs/next_session_iddes.md). These are
the LONG runs, packaged to run on a remote machine:

```bash
# on the big machine (after rsync + ./compile.sh cpu at the repo root)
cd validation/rans_sst
RANKS=8 ./run_gates.sh            # all cases, sequentially; or one name:
RANKS=16 ./run_gates.sh ibm180

# back home (or anywhere with python3 + h5py + numpy)
./check_gates.sh
```

Cases (each writes `<prefix>_<laststep>.h5` + `<prefix>.log`). RESULTS
2026-07-08 (runs on the remote machine, checks local) — all PASS:

- `laminar.ini`   gate (a): Re_tau 10 channel, SST on, tu = 1%. Parabola
                  to 2.4e-4 of u_max; k decayed 4e-3 -> 3e-47. (At
                  Re_tau 30 / tu 5% the no-transition SST finds its
                  textbook weakly-turbulent branch — k self-sustains at
                  nut/nu ~ 0.5 and the parabola flattens 12% — so the
                  gate probes the genuinely subcritical regime.)
- `turb180.ini`   gate (b): developed channel Re_tau 180 (natural y,
                  y+_1 ~ 0.5). Log line to 4.9%; U+ centreline 18.16 vs
                  DNS 18.20 (0.2%); u_tau = 1.0008.
- `turb395.ini`   gate (b): Re_tau 395. Log line to 6.5% (the pure
                  kappa/B line deviates from real profiles by several %
                  over the overlap — hence the DNS centreline anchor):
                  U+ centreline 20.16 vs DNS 20.13 (0.15%);
                  u_tau = 1.0016.
- `ibm180.ini`    gate (c), THE key IBM gate: the les_ibm off-grid
                  file-IBM plane channel (walls mid-cell, uniform y,
                  y+_1 ~ 2-3) through the IBM wall treatment
                  (dwall_blocks + wall-cell omega pinning). Log line to
                  4.3% with the walls at y = 0.259375/2.259375, u_tau
                  0.966. Needs ../rans_geometry/ibm_coeff_blocks_l1.h5
                  (../rans_geometry/setup.sh regenerates it).
- `refine180.ini` gate (d): uniform-y Re_tau 180 with both wall bands
                  refined to level 1 (2:1 interfaces at y = 0.25/1.75).
                  NO interface band: jump/local-variation ratios 0.58
                  (k), 1.11 (omega), 0.90 (nut), 0.73 (u); core profiles
                  match the RESOLVED turb180 reference to 0.5% (u/k),
                  2.8% (omega), 1.7% (nut).
- `base180u.ini`  informational coarse-wall control (uniform y+_1 ~ 2.8,
                  no refinement): its under-resolved sublayer feeds a
                  spurious core-k plateau (~1.9 vs the resolved ~0.9) —
                  the demonstration of WHY the wall bands get refined.
                  Not a pass/fail gate.

Gate (e) (CPU==GPU + LES/no-model bit-exactness, nofma) is NOT here — it
is short and runs locally (see the phase notes in CLAUDE.md; cap nsteps
~20 per the bit-exact-gates-short rule).

`rans_channel_check.py` holds the pass criteria (tolerances overridable
per invocation in check_gates.sh).
