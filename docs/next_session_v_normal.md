# Next session — quantify + fix the v' interface-NORMAL step

Read `docs/interface_band_handout.md` first (full state; "Session 2026-06-26" and
"-26b" are the latest). Branch `claude/jacobi-interface`. Memory:
`interface-validation-suite`.

## Where we are (done, committed)

- **Constant-1/2 energy-conserving interface is the DEFAULT** (`[blocks]
  interface_constant_half`, commit 1428641): inject the velocity prolong + skip the
  cubic deep-halo reconstruction. It makes the refined channel STABLE (div_max
  settles ~0.056 over 500 steps) where the old metric/cubic interface BLEW UP
  ~step 200. Single-level runs bit-exact; 1-rank == 2-rank bit-exact (x-split);
  rank→GPU auto-pinned, so 2-GPU runs work.
- **MOBY_KEBAL** convective KE-balance gate (47961fe) + optional **MOBY_KESKEW /
  interface_skew** skew-symmetric band correction (4c43741). The V&V levers cut the
  band energy 0.127 → 0.054 (no cubic) → 0.045 (const-1/2 prolong) → 0.036 (local
  skew); the residual 0.036 is the cross-interface area-mismatch term.
- Validation cases curated: `validation/channel_interface/VALIDATION_CASES.md`
  (interface benchmark (stability + banding), KEBAL gate, bit-exact, multi-rank, developed stats).
- Developed-flow stats setup: `validation/channel_interface/developed/`
  (`run_developed.py`, two-leg, 2-GPU). `tools/plot_channel_stats.py`.

## THE LEAD: the v' interface-NORMAL step (handout -26b)

Single-snapshot profiles (t=0.078, EARLY) show the interface signatures sort by
component role:
- u' (tangential, mean-shear): big rms SPIKE at both interfaces (the streak band).
- w' (tangential, no shear): almost NO spike.
- **v' (wall-NORMAL = interface-normal): a SOLVER-CREATED STEP** — IC smooth across
  the interface (0.620→0.616→0.609), step 250 fine side SUPPRESSED (→0.491) /
  coarse side RAISED (→0.741). v is the component with the special 2:1-face
  treatment (restriction / shared face / fine-owns-face), so this is the prime
  suspect for a REAL interface-normal defect, distinct from the u' band and
  consistent with the residual KEBAL 0.036.

## Do this, in order

1. **Get CONVERGED statistics** (the snapshot at t≈0.08 is not converged; the
   two-wall asymmetry there is just single-realization, already in the IC ~5%).
   Run the developed case on 2 GPU:
   `cd validation/channel_interface/developed && python3 run_developed.py --ranks 2`
   and the `--skew` variant. Also run a **uniform-fine reference** (256x128x256, no
   interface — `validation/channel_interface/reference.ini` style) for the same
   t-window so the interface profiles have a ground truth.
   Quantify: is the v' step still there in the time-average? How big vs the
   uniform-fine v' at the same y? Does interface_skew reduce it?
2. **If the v' step persists** (likely): build the **volume-weighted ADJOINT
   transfer** — restriction = transpose of prolongation (V&V §2.1.3, the one untried
   piece). The current pair is interp/inject prolong + simple-average restrict
   (non-adjoint). The interface-NORMAL velocity is the component the restrict/
   shared-face touches, so an adjoint (energy-consistent) transfer is the principled
   cure for the v' step and the residual KEBAL 0.036. Gate it; verify the KEBAL band
   SKEW → round-off on the slab, the v' step shrinks in the developed stats, bit-
   exact no-interface, 1==2 rank.
3. **Optional accuracy**: replace the disabled cubic reconstruction with a
   constant-1/2 LINEAR-AVERAGE ghost (keeps some interface accuracy while satisfying
   V&V constant-1/2) instead of skipping it entirely.

## Guardrails (every change)

Re-run the validation cases (`VALIDATION_CASES.md`): interface benchmark PASS (stability) + band metrics (u' excess, v' kink) tracked,
KEBAL slab, bit-exact no-interface, 1==2 rank. Don't chase interface truncation
order — chase the KE balance / the v' step. CPU `-Mnofma`, GPU `-Mnofma -gpu=nofma`
for bit-exact; one GPU case at a time; pkill stray `build_*/main` between runs.
