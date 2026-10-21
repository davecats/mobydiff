# Interface decay gate

Quick stability gate for the 2:1 refinement interface treatment: white
noise (`[flow] initial_noise`) on a 3D refined patch in a periodic box,
no forcing. All velocity and pressure extrema must contract under
viscosity; exponential growth localized at the interfaces is exactly how
the original one-sided interface relaxation failed, and smooth-flow
gates (uniform flow, band-refined channels) are blind to it.

Run and check (about a minute on a GPU):

```bash
mpirun -n 1 ../../build_cpu/moby_prepare input.ini   # the case file, once
mpirun -n 1 ../../build_gpu/main input.ini
python3 ../../tools/check_interface_decay.py .
```

The checker compares max|u|,|v|,|w|,|p| between step 20 (after the
initial projection has settled the noise) and step 200 and fails unless
everything decayed and stayed finite.

RE-VALIDATED 2026-10-21, commit 532ef54 (banner), production flags, RTX 3060,
1 rank, 14 s (`niter = 3`, as the ini declares): PASS -- step 20 -> 200
max|u| 2.07e-2 -> 2.68e-3, |v| 2.23e-2 -> 3.01e-3, |w| 2.00e-2 -> 2.88e-3,
|p| 7.14e-2 -> 6.55e-3, a contraction of 7.7x (velocity) / 10.9x
(pressure). No earlier number was recorded for this gate. Settings changed
for the rerun: none.
