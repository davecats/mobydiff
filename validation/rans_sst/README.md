# RANS T2+T3 gates — k-omega SST transport (resolved walls + wall functions)

Physics gates for IDDES phases T2 and T3 (docs/next_session_iddes.md).
These are the LONG runs, packaged to run on a remote machine (the T3
wall-function sweep is tiny and also runs fine on a laptop; ibm180wf ran
on the local GPU):

```bash
# on the big machine (after rsync + ./compile.sh cpu at the repo root)
cd validation/rans_sst
RANKS=8 ./run_gates.sh            # all cases, sequentially; or one name:
RANKS=16 ./run_gates.sh ibm180
./run_gates.sh t3                 # just the T3 wall-function set

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

T3 wall-function cases (`[rans] wall_treatment = wall_function`; all
checked with `--mode wallfn` against the RESOLVED turb180 field, which
carries the DNS centreline anchor to 0.2%, so the anchor is transitive).
RESULTS 2026-07-08 (runs local: sweep on 4-rank CPU, ibm180wf on the
local GPU) — all PASS:

- `wf180_y30/45.ini` T3 gate (a): turb180 coarsened to uniform ny = 6/4
                  (y+_1 = 30/45; ny = 6 is not nb-divisible, so it runs
                  with [blocks] nb unset — rank-box blocks). Implied U+
                  centreline vs DNS 18.20: 1.2% / 0.7%; u_tau from the
                  delivered wall stress (nu + nut_1) U_1/y_1 = 1.0000.
- `wf180_y05/y15/y22.ini` T3 gate (b): graceful degradation across the
                  buffer range (y+_1 = 5/15/22.5). Implied centrelines
                  -3.1% / +2.8% / +3.0% — a mild overshoot, NO
                  double-counting dip (that would be a deficit). First
                  cells below y+ 30 sit 10-19% above the resolved
                  profile: the log-line error at the anchor cell,
                  informational. base180u (resolved, y+_1 2.8) is the
                  marginal control.
- `ibm180wf.ini`  T3 gate (c): ibm180 through the wall-function blend
                  (200k steps, local GPU, 18.6 ms/step). The k-based y+
                  keeps every wall cell on the viscous branch, whose
                  arithmetic is exactly the resolved treatment: the
                  steady profile matches the T2 resolved ibm180 field to
                  ROUND-OFF (u 4.5e-16, nut 1.1e-16) — zero regression.

Gates (d)/(e) (resolved-mode bit-exactness vs T2 8991192 + wall-function
determinism, nofma) are NOT here — they are short and run locally (cap
nsteps ~20 per the bit-exact-gates-short rule). T3 results: resolved
bit-exact (max_abs 0, CPU AND GPU) on min_channel/les_ibm ± refine_body/
Beltrami y-slab/turb180; wall-function 1-rank == 4-rank exactly, CPU vs
GPU <= 2e-13 (the `log()` intrinsic in the wall-function branch differs
by an ulp between host and device libm — resolved mode remains exactly
CPU == GPU).

`rans_channel_check.py` holds the pass criteria (tolerances overridable
per invocation in check_gates.sh).
