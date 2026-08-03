# A1/A2 cylinder gates — airfoil flow case + control-volume forces

Validation for phase A1 (the `[case] name = airfoil` module) of
`docs/next_session_airfoil.md` and for the runtime C_L/C_D statistic, which
is the CONTROL-VOLUME momentum budget over `[case.airfoil] cv_box`
(`docs/configuration.md`, `docs/next_session_cv_forces.md`).

Quasi-2D cylinder, D = 1 at (6.0, 8.02) in a 16 x 16 x 0.25 box, uniform
h = 1/32 (D/h = 32), z periodic (nz = nb = 8). The 0.02 D vertical offset
seeds the Re = 100 shedding; it is inconsequential at Re = 40. The force is
sampled every `force_sample_interval` steps in `after_step` (per-block device
sums scattered into the global block table -> the allreduce is exact and the
final sum runs in global-id order, so the sampled force is independent of the
rank count BY CONSTRUCTION).

Setup (STL + per-Re coefficient files; coef = SOLID/re so one file per Re):

    ./setup.sh          # needs moby_prepare + the ibmc venv
                        # (trimesh + shapely + mapbox_earcut + h5py)

Runs (GPU recommended; one job at a time):

    mpirun -n 1 ../../build_gpu/main cyl_re40.ini     # steady drag
    mpirun -n 1 ../../build_gpu/main cyl_re100.ini    # vortex shedding
    mpirun -n 1 ../../build_cpu/main empty.ini        # empty-domain zero force

## The force statistic changed (2026-08-03)

These gates were originally built on the PENALIZATION integral
`F = int coef*u dV`, which has been REMOVED from the solver: it is exact
bookkeeping only while the solid interior is present, and production runs
remove the buried core (`[blocks] remove_solid`). The runtime now reports the
control-volume budget, and every ini here carries a `cv_box`.

CONSEQUENCE, not yet done: `forces_re40.txt` and `forces_re100.txt` are still
the OLD penalization series, so the `steady` and `strouhal` gates below quote
penalization numbers. They must be REGENERATED before they can be quoted as
control-volume results — and the Re 100 one needs the clean-p protocol below,
because a 40 000-step run at `niter = 6` pollutes the stored pressure that the
budget reads. The `empty` and determinism gates were re-run and pass as stated.

Cross-validation of the new statistic, on the converged Re 40 state
(`cvpz_20301.h5`, step 20310, box 4-8 x 6.5-9.5): runtime C_D = 1.70387475,
C_L = -1.2e-3, reproduced to NINE significant digits per border by an
independent reassembled-global-plane evaluation, and to 0.011 % by the
offline `tutorials/naca/rans/postProcess/cv_forces.py --boxes 1.5 --nose 5.5
8.0 --span-z` (its residual is the first-order one-sided gradients it must
use at block edges, having no halos). Box scatter on that same field:
C_D = 1.698 (margin 1.5 D) / 1.704 / 1.726 / 1.765 (6 D) — use a TIGHT box.

Gates (python3 check_cylinder.py ...):

- `steady forces_re40.txt` — PASS (2026-07-13, PENALIZATION series):
  C_D = 1.6924 +- 2.7e-5 over the last 20% (t = 80..100), |C_L| = 4.0e-4.
  The value sits ABOVE the unbounded-flow 1.5-1.6 band as expected for this
  setup: ~6% Dirichlet-far-field blockage at 16D plus the first-order
  penalization's effective diameter (~D + h = 1.03D); the hard gate band is
  1.4-1.7. The control-volume budget on the converged state reads 1.704 for
  the committed box, i.e. inside the same band and 0.7 % above the
  penalization value.
- `strouhal forces_re100.txt` — PASS (2026-07-13, PENALIZATION series):
  St = 0.168 (spectral FUNDAMENTAL of C_L over t = 100..200), mean
  C_D = 1.448, mean C_L = 2.4e-4 (symmetric), C_L amplitude 0.51. Shedding
  self-starts at t ~ 40 from the 0.02 D offset. NOTE: the confined C_L
  carries a 3rd harmonic of comparable power (period exactly T0/3), so
  zero-crossing counting reports 3x St — the metric takes the lowest
  spectral peak within 35% of the largest.
  Control-volume cross-check over t = 201.5..205 on a clean-p restart of the
  same state: mean C_D 1.4469 vs the penalization 1.4486 (0.12 %), C_L range
  [-0.572, +0.486] vs [-0.536, +0.483]. This is the gate that exercises the
  budget's unsteady term — WITHOUT it the C_L mean flips sign (+0.503 vs
  -0.031), since d/dt of the box momentum is 1.15 in C_L units here.
- `empty forces_empty.txt` — PASS (exact): C_L = C_D = 0.0, and the
  aoa = 5 deg freestream preserved exactly (case-composed twin of
  validation/freestream/oblique.ini). The committed `cv_box` is deliberately
  one cell INSIDE the block boundaries so every border face is
  block-interior and the closed-box sum cancels EXACTLY; on a
  block-aligned box the same run gives 8.7e-16 (the one-sided pressure
  branch at a block's low edge, see below).
- force determinism — PASS: forces file for 1 vs 4 CPU ranks BYTE-IDENTICAL
  (the ordered global-id reduction); CPU vs GPU sampled C_L/C_D difference
  0.0 on this case. Re-verified 2026-08-03 for the control-volume budget on
  both `empty.ini` and the Re 40 case.

## Clean-p protocol (REQUIRED for control-volume forces on a long run)

Long IBM runs at production niter = 6 accumulate a large VELOCITY-NEUTRAL
oscillating mode in the stored pressure (std ~4e2 here; the same family as
the known channel pn drift, memory: pressure-volume-average-drift). u, the
dynamics and the old u-only penalization force were untouched by it, but the
control-volume budget reads that pressure directly: on the committed
`cyl_re40_20001.h5` it returns C_D = 261. The mode is SPATIALLY VARYING, so
subtracting a reference pressure does not rescue it.

Getting a clean-p state: copy the converged restart, ZERO its pn (h5py),
restart with niter = 60 for ~300 steps (the converged projection rebuilds the
physical p in a few substages — stagnation +0.56, outlet column pinned to
2e-5). Do NOT restart with the polluted p at niter = 60 (the accumulated
spurious grad-p loses its self-consistent sloppy-projection compensation ->
violent transient, C_L ~ O(100), dt collapse), and do NOT run the whole case
at niter = 60 (~15x production cost here). `cvpz_20301.h5` is such a state
for Re 40.

From a clean pressure the budget then HOLDS at production niter: 2000 steps
of Re 40 at niter = 6 kept C_D = 1.692 +- 0.001 with no growth. It is the
accumulation over a full multi-10 000-step run that eventually poisons it.

Note: the control-volume force is the TOTAL force (pressure + friction
combined; modeled turbulent stress enters through the eddy-viscosity part of
tau) — the split needs surface integration, which is not planned.
