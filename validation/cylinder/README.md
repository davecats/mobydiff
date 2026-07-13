# A1/A2 cylinder gates — airfoil flow case + penalization-integral forces

Validation for phases A1 (the `[case] name = airfoil` module) and A2 (the
runtime C_L/C_D statistic F = int coef*u dV) of `docs/next_session_airfoil.md`.

Quasi-2D cylinder, D = 1 at (6.0, 8.02) in a 16 x 16 x 0.25 box, uniform
h = 1/32 (D/h = 32), z periodic (nz = nb = 8). The 0.02 D vertical offset
seeds the Re = 100 shedding; it is inconsequential at Re = 40. The force is
sampled every `force_sample_interval` steps in `after_step` (per-block device
sums scattered into the global block table -> the allreduce is exact and the
final sum runs in global-id order, so the sampled force is independent of the
rank count BY CONSTRUCTION).

Setup (STL + per-Re coefficient files; coef = SOLID/re so one file per Re):

    ./setup.sh          # needs build_cpu/mobygrid + the ibmc venv
                        # (trimesh + shapely + mapbox_earcut + h5py)

Runs (GPU recommended; one job at a time):

    mpirun -n 1 ../../build_gpu/main cyl_re40.ini     # steady drag
    mpirun -n 1 ../../build_gpu/main cyl_re100.ini    # vortex shedding
    mpirun -n 1 ../../build_cpu/main empty.ini        # empty-domain zero force

Gates (python3 check_cylinder.py ...), RESULTS 2026-07-13:

- `steady forces_re40.txt` — PASS: C_D = 1.6924 +- 2.7e-5 over the last 20%
  (t = 80..100), |C_L| = 4.0e-4. The value sits ABOVE the unbounded-flow
  1.5-1.6 band as expected for this setup: ~6% Dirichlet-far-field blockage
  at 16D plus the first-order penalization's effective diameter (~D + h =
  1.03D); the hard gate band is 1.4-1.7.
- `strouhal forces_re100.txt` — PASS: St = 0.168 (spectral FUNDAMENTAL of
  C_L over t = 100..200), mean C_D = 1.448, mean C_L = 2.4e-4 (symmetric),
  C_L amplitude 0.51. Shedding self-starts at t ~ 40 from the 0.02 D offset.
  NOTE: the confined/penalization C_L carries a 3rd harmonic of comparable
  power (period exactly T0/3), so zero-crossing counting reports 3x St —
  the metric takes the lowest spectral peak within 35% of the largest.
- `empty forces_empty.txt` — PASS (exact): C_L = C_D = 0.0, and the
  aoa = 5 deg freestream preserved exactly (case-composed twin of
  validation/freestream/oblique.ini).
- force determinism — PASS: forces.txt for 1 vs 4 CPU ranks BYTE-IDENTICAL
  (the ordered global-id reduction); CPU vs GPU sampled C_L/C_D difference
  0.0 on this case.
- `cv <clean-p snapshot> --re 40 --cd-pen <C_D>` — the Gauss/CV cross-check:
  outer-box momentum + pressure + viscous fluxes on the converged steady
  snapshot must reproduce the penalization C_D to discretization error.
  Independent validation of BOTH the force statistic and the in/outflow
  faces (the inlet/outlet fluxes dominate the border balance).
  RESULT 2026-07-13: PASS — box [4,9]x[6,10]: C_D = 1.6906 (0.1% off the
  penalization 1.6924!), box [3,11]x[4,12]: 1.733 (2.4%), box
  [2,14]x[2,14]: 1.803 (6.5%; longer faces accumulate more collocation/
  gradient error).
  CAVEAT (found 2026-07-13): long IBM runs at production niter = 6
  accumulate a large VELOCITY-NEUTRAL oscillating mode in the stored
  pressure (std ~4e2 here; the same family as the known channel pn drift,
  memory: pressure-volume-average-drift). u, the forces (u-only) and the
  dynamics are untouched, but an instantaneous pn snapshot is useless for
  border fluxes. Getting the clean-p snapshot: copy the converged restart,
  ZERO its pn (h5py), restart with niter = 60 for ~300 steps (forces hold
  1.69 throughout; the converged projection rebuilds the physical p in a
  few substages — stagnation +0.56, outlet column pinned to 2e-5).
  Do NOT restart with the polluted p at niter = 60 (the accumulated
  spurious grad-p loses its self-consistent sloppy-projection compensation
  -> violent transient, C_L ~ O(100), dt collapse), and do NOT run the
  whole case at niter = 60 (~15x production cost here).

Note: the penalization force is the TOTAL force (pressure + friction
combined; modeled turbulent stress included automatically through the
mu-masked eddy-viscosity correction) — the split needs surface integration,
which is not planned.
