# A3 INCREMENT 3 — SD7003 transition benchmark (Re_c = 6e4, aoa = 4)

The gamma-Re_thetat transition model on the standard low-Re LSB benchmark,
through the full airfoil stack (freestream composition, file IBM from the
UIUC Selig coordinates spline-resampled to 720 points, 5-level
refine_body with --keep-buried, SST + transition, scalar inlets holding
k/omega/gamma/Re_thetat~ freestream values, penalization forces).

- tu = 0.1 %, nut_ratio = 1 => omega_inf = 0.09: ambient decay time
  1/(beta* omega_inf) = 123 c/U — the inlet tu reaches the airfoil
  unchanged (the naca0012 tu = 5 ambient-decay lesson does not bite).
- Grid = validation/naca0012 (Delta = 1.465e-3 c at the surface; laminar
  y+_1 ~ 1.8 at Re 6e4, the docs' 3-4-level arithmetic).

Gates (check_sd7003.py; published gamma-Re_thetat RANS-LM scatter —
Windte/Radespiel, LES reference Galbraith & Visbal):

- LSB present (reversed-flow patch on the suction side);
- transition location x_t/c = 0.5 +- 0.1 (near-wall gamma crossing 0.5);
- C_L within +-15 % of ~0.60, C_D within +-15 % of ~0.022;
- the gamma-front chordwise smearing (stations from gamma 0.1 to 0.9, in
  level-4 cells) is MEASURED AND REPORTED ONLY — it is the number that
  decides the separate TVD/van-Leer + second-scalar-halo increment.

## Workflow

```bash
./setup.sh          # grid.h5 + resampled STL + keep-buried block-table (~1 h)
mpirun -n 1 ../../build_gpu/main aoa4.ini     # ~7 h on the RTX 3060
python3 check_sd7003.py sd7003_aoa4_*.h5 forces_aoa4.txt
```
