# Interface decay gate

Quick stability gate for the 2:1 refinement interface treatment: white
noise (`[flow] initial_noise`) on a 3D refined patch in a periodic box,
no forcing. All velocity and pressure extrema must contract under
viscosity; exponential growth localized at the interfaces is exactly how
the original one-sided interface relaxation failed, and smooth-flow
gates (uniform flow, band-refined channels) are blind to it.

Run and check (about a minute on a GPU):

```bash
mpirun -n 1 ../../build_gpu/main input.ini
python3 ../../tools/check_interface_decay.py .
```

The checker compares max|u|,|v|,|w|,|p| between step 20 (after the
initial projection has settled the noise) and step 200 and fails unless
everything decayed and stayed finite.
