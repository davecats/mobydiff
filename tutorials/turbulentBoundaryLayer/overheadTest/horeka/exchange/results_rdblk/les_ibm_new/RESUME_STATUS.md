# LES <-> IBM coupling validation — STATUS

Branch `claude/jacobi-interface`. See `README.md` for the full case + results table.

## COMPLETE — converged developed-statistics campaign run (t=5..25, 51 snapshots)

The campaign ran on a faster GPU (a_wale / b_none / c_refine to ~25 600 steps).
ALL gates PASS converged (`ibm_les_profiles.png`):
- Gate 1 solid nut==0 exact (51 snapshots). Gate 2 band/core nut ratio 0.05.
- Gate 3 a_wale + b_none + grid-aligned `../les/` collapse on the log law; WALE
  bulk U=15.09 > no-LES 14.69; resolved stresses symmetric across both IBM walls.
- Gate 4 nut steps 2.03x up into the coarse core across the y=0.75/1.75 2:1
  interface — smooth step, NO band.
- Gate 5 stable ~25 600 steps. Gate 6 N/A. Gate 7 CPU==GPU 4.6e-14.

Nothing outstanding. (Optional follow-ups below.)

## Earlier dev-GPU confirmation (superseded by the converged campaign)

The LES<->IBM coupling (`ibm_aware` solid-cell nut masking, `les.f90`) is validated
on an off-grid IBM plane-wall channel. **No solver code changed** — the existing
mask + WALE `sd2` already give the physical `nut->0` into the wall, no band.

- **Gate 1 PASS** solid-cell `nut==0` exactly (every snapshot).
- **Gate 2 PASS** no spurious wall-nut spike: band/core nut ratio 0.05 (physical
  `nut->0`). The prime suspect does NOT materialise; no band-aware damping needed.
- **Gate 3 PASS (sanity)** law of the wall recovered (sublayer `U+~y+`, log layer
  `2.44 ln y+ +5`, bulk U=15.5). Converged stats = the campaign below.
- **Gate 4 PREVIEW PASS** case (c) `nut(y)` smooth across the y=0.75 2:1 interface,
  no band; nut steps UP into the coarse core (physical delta^2 direction).
- **Gate 5 PASS** stability: case (a) 4000 steps + case (c) 400 steps, div bounded,
  mass ~1e-15, no NaN. The 2:1-interface x IBM x LES triple does not blow up.
- **Gate 6 N/A** no solver code changed.
- **Gate 7 PASS** CPU vs GPU case (a) 10 steps agree to 4.6e-14 (FMA round-off) —
  the offloaded `ibm_aware` masking branch is correct on GPU.

Evidence runs (gitignored): `runs/a_transient/`, `runs/a_stats/` (12 nut snapshots),
`runs/c_smoke/`.

## Reproduce the campaign

```bash
cd validation/channel_interface/les_ibm
[ -f IC_refine.h5 ] || ./setup.sh            # regenerate the gitignored case-c files (needs geometry venv) — or rsync them
python3 run_ibm_les.py --arch gpu --case all --mpirun "<mpirun>"   # a_wale, b_none, c_refine; t=5..25
python3 measure_nut.py   --run runs/a_wale/stats              # gates 1-2
python3 ibm_les_stats.py                                       # gates 3-4 -> ibm_les_profiles.png
```

## Optional follow-ups (not blocking)
- Strict CPU==GPU bit-exact (nofma builds) -> 0.0, vs the 4.6e-14 FMA result here.
- box-filter the DNS reference for a stricter resolved-stress gate (as in `../les/`).
