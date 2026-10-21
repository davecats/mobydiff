# Minimal-span channel: the suite's blocks + 2:1 + Chebyshev case

A Re_tau = 180 channel (natural-stretched y, `forcing_x = 1`, the channel
initialiser's mean profile + transition disturbance) shrunk in span to a
minimal-flow-unit slab, with both near-wall bands refined to level 1 (flat
y-interfaces at y+ = 112, as in `validation/channel_interface`). It keeps the
channel velocity (~18) that makes the 2:1 interface projection a real test
-- a laminar Poiseuille base never exercised it -- at ~1/200th of the
validation's cost. It is the first case of the 9-case bit-exactness suite
(`validation/scalar/run_bitexact.sh`, 20 steps, 1 == 4 ranks), i.e. the
suite's blocks + 2:1 interface + Chebyshev-Jacobi (`niter = 12`) case.

| ini | grid | steps | role |
|---|---|---|---|
| `input.ini` | 128 x 64 x 8 | 400 | the SUITE case (the driver caps it at 20 steps) |
| `input_gpu.ini` | 128 x 64 x 16 | 4000 | the refined case at full length, a stability run |
| `input_gpu_noref.ini` | 128 x 64 x 16 | 4000 | its uniform control (`[blocks] refine` removed) |
| `input_gpu_short.ini` | 128 x 64 x 16 | short | a quick GPU smoke test |

```bash
module load toolkits/nvhpc/25.9
mpirun -n 1 build_cpu/moby_prepare tutorials/min_channel/input_gpu.ini   # the case file, once
mpirun -n 1 build_gpu/moby_solve tutorials/min_channel/input_gpu.ini
```

Inspect max|u|, max|v| over the `channel_field_*.h5` snapshots (a few lines
of h5py over `un`/`vn`): the refined run must stay bounded, turbulent O(20),
like the control.

History: this directory was the Step-1 confirmation gate of the original 2:1
interface projection blow-up (`docs/history/interface_review.md` section vi),
run with the old red-black `sor = 1.5 / niter = 3` before the interface
treatment was reworked; those numbers no longer reproduce and were removed
on 2026-10-21 (`git show 532ef54:tutorials/min_channel/README.md` has them).
