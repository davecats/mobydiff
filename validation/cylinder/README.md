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

    ./setup.sh          # needs build_cpu/moby_prepare; the ibmc venv
                        # (trimesh + shapely + mapbox_earcut) only to
                        # regenerate the committed cylinder.stl

(Before 2026-09-26 setup.sh used the deleted `mobygrid` + the retired
`mobygeom.py stl-ibm-coeff`; the recorded results below come from those
legacy-format files. The moby_prepare files keep all 4096 blocks too.)

Since 2026-10-01 the two `cyl_re*.ini` name the STL and `remove_solid =
false` themselves and `setup.sh` prepares FROM them: the solver compares the
case file's input echo with its own ini, and a case file prepared from a
sed-derived variant was refused as stale (this gate had not been run since
the step-7 contract).

Runs (GPU recommended; one job at a time):

    mpirun -n 1 ../../build_cpu/moby_prepare empty.ini       # its case file (the cyl_* inis name theirs)
    mpirun -n 1 ../../build_gpu/moby_solve cyl_re40.ini     # steady drag
    mpirun -n 1 ../../build_gpu/moby_solve cyl_re100.ini    # vortex shedding
    mpirun -n 1 ../../build_cpu/moby_solve empty.ini        # empty-domain zero force

## The force statistic changed (2026-08-03)

These gates were originally built on the PENALIZATION integral
`F = int coef*u dV`, which has been REMOVED from the solver: it is exact
bookkeeping only while the solid interior is present, and production runs
remove the buried core (`[blocks] remove_solid`). The runtime now reports the
control-volume budget, and every ini here carries a `cv_box`.

CONSEQUENCE: the original `forces_re40.txt` and `forces_re100.txt` are the
OLD penalization series (generated files, not in git). The Re 40 gate below
is RE-MEASURED on the control-volume series (2026-10-01,
`run_f7_gates.sh steady` writes `f7/forces_st_new.txt`). The Re 100
`strouhal` gate still quotes penalization numbers: regenerating it needs the
clean-p protocol below, because a 40 000-step run at `niter = 6` pollutes the
stored pressure that the budget reads — and the outlet pressure mode of the
last section limits how long a clean-p leg stays clean. The `empty` and
determinism gates were re-run and pass as stated.

Cross-validation of the new statistic, on the converged Re 40 state
(`cvpz_20301.h5`, step 20310, box 4-8 x 6.5-9.5): runtime C_D = 1.70387475,
C_L = -1.2e-3, reproduced to NINE significant digits per border by an
independent reassembled-global-plane evaluation, and to 0.011 % by the
offline `tutorials/naca/rans/postProcess/cv_forces.py --boxes 1.5 --nose 5.5
8.0 --span-z` (its residual is the first-order one-sided gradients it must
use at block edges, having no halos). Box scatter on that same field:
C_D = 1.698 (margin 1.5 D) / 1.704 / 1.726 / 1.765 (6 D) — use a TIGHT box.

Gates (python3 check_cylinder.py ...):

- `steady f7/forces_st_new.txt` — PASS (2026-10-01, CONTROL-VOLUME series,
  2000 steps from the clean-p converged state, t = 101.5..111.5):
  C_D = 1.6924 +- 6.7e-4 over the last 20%, |C_L| = 1.6e-4 (over the last
  1000 steps: 1.69234525 +- 8.2e-4, the number the penalization-factor
  section below tracks). The value sits ABOVE the unbounded-flow 1.5-1.6
  band as expected for this setup: ~6% Dirichlet-far-field blockage at 16D
  plus the first-order penalization's effective diameter
  (~D + h = 1.03D); the hard gate band is 1.4-1.7. History: the 2026-07-13
  PENALIZATION series read 1.6924 +- 2.7e-5 over t = 80..100.
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

## Steady-state stop (`[case.airfoil] steady_tol`)

The budget's unsteady term doubles as a convergence monitor: `steady_tol`
stops the run, and writes the final field, once `|2 dmom/dt| / qref` (the
unsteady term in coefficient units, larger of the two components) stays below
the tolerance for `steady_samples` consecutive samples.

Exercised here 2026-08-03, both directions:

- Re 40 (steady), clean-p restart, `steady_tol = 1e-3`: stopped itself after
  779 steps at t = 105.4 with the measure at 9.3e-4, against a configured
  `t_final = 120` — about 2900 steps saved — and wrote the final snapshot.
- Re 100 (shedding) must NOT trigger, and does not come close: over 99
  samples spanning ~0.84 shedding periods the measure ranges 3.1e-2 to 5.7
  (median 0.89), never dipping below even `1e-2`. Each component's d/dt does
  cross zero twice per period, but the two are out of phase, so taking the
  larger of them never vanishes; the consecutive-sample requirement is a
  second line of defence.

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

## The penalization factor (numerics review step 9, F7) — measured 2026-10-01

The factor now shipped is the RATIONAL third-order form
(`ibm.f90 penal_*`, `../penalization/README.md`); the `dt` study below was
run with the exact exponential that preceded it by a few hours, and its
conclusion (the factor is not what limits this flow) carries over: on a
time-varying right-hand side the two are indistinguishable.

`run_f7_gates.sh` runs two legs with TWO binaries from the same state (the
restarts live in a sibling checkout, `SIB`; `REF` = `~/step8_ref_binaries`,
the implicit-Euler factor; `NEW` = the exact factor); `dt_study.py` analyses
the second.

**steady — the drag does not move.** Re 40 from the clean-p converged state
(`cvpz_20301.h5`), 2000 steps at the production settings:

| | C_D (last 1000 steps) | C_L |
|---|---|---|
| implicit Euler | 1.69234524 ± 8.2e-4 | −3.9e-4 |
| exact exponential (step 9) | 1.69234525 ± 8.2e-4 | −3.9e-4 |
| rational P3 (2026-10-01, the shipped factor) | 1.69234525 ± 8.2e-4 | −3.9e-4 |

As it must: all three factors have the same steady fixed point,
`lambda q = R`. Exponential against rational, the same leg with
`REF = ~/step9_ref_binaries`: mean C_D 1.6923452515 / 1.6923452519, and the
two sampled series differ by at most 1.0e-8 in C_D and 2.2e-9 in C_L (the
file's print precision is 1e-8).

**dt — the factor is not what limits the time accuracy of this flow.** Re
100 shedding state, pressure cleaned (zero `pn`, 300 steps at niter 60), then
ONE time unit at fixed dt = 5e-3 / 2^k, k = 0..3, Chebyshev niter 60.
2208 cut DOFs with `lambda dt_gamma` from 1e-4 to 44 at dt = 5e-3.

| dt | old − new, cut rms | cut max | fluid rms | body-box rms vs the finest run | C_D(T) old / new |
|---|---|---|---|---|---|
| 5.0e-3 | 2.5e-6 | 1.1e-5 | 6.3e-8 | 9.7e-3 | 1.50084002 / 1.50084006 |
| 2.5e-3 | 2.2e-6 | 8.6e-6 | 4.5e-8 | 2.2e-3 | 1.52468122 / 1.52467630 |
| 1.25e-3 | 1.8e-6 | 6.7e-6 | 3.0e-8 | 6.4e-4 | 1.45727014 / 1.45727675 |
| 6.25e-4 | 2.2e-6 | 8.0e-6 | 2.3e-8 | (reference) | 1.46381330 / 1.46381342 |

The two factors differ by ~1e-5 in velocity at every step size, three orders
of magnitude below the dt-to-dt differences near the body. A body at rest
keeps its cut cells quasi-steady (u ~ 0, slaved to the neighbours), so the
first-order error of the old factor multiplies a time derivative that is
nearly zero. **The expectation "the cut-cell time error drops faster" cannot
be observed here because that error was never visible**; the place where the
factor's order IS measured is `validation/penalization/` (first order →
round-off). A moving or impulsively started body would be the flow where it
matters.

**FOUND, pre-existing and identical for both binaries: this study is NOT a
clean time-convergence test, because the stored pressure pollutes within one
time unit.** From the cleaned state (p rms 0.30, max 1.4):

| dt | p rms at T | max |p| | C_L range over the run | where the largest velocity difference to the finest run sits |
|---|---|---|---|---|
| 5.0e-3 | 0.89 | 3.0 | −0.07 .. +1.16 | outlet corner (x = 15.9, y = 15.9), 0.15 |
| 2.5e-3 | 0.58 | 2.3 | −0.19 .. +0.73 | outlet corner, 2.4e-2 |
| 1.25e-3 | 0.79 | 3.1 | −0.75 .. +1.33 | outlet corner, 1.7e-2 |
| 6.25e-4 | 1.65 | 6.0 | −3.09 .. +3.90 | — |

The physical C_L amplitude is 0.5. The dt-to-dt velocity differences peak in
the outlet corners and the lateral far field, not in the wake, and the fluid
rms is non-monotone in dt (2.7e-2, 1.2e-3, 2.2e-3): this is the
velocity-neutral-then-active pressure mode at the Dirichlet-p outlet under a
Chebyshev projection that does not converge the long waves in 60 iterations
(`lmin` ~ 2.5e-5 on N = 512), the family of the clean-p section above and of
the B0 outlet finding. Only the body-box column falls like a time error
(ratios 4.4 and 3.4). Open; it belongs with the outlet / projection work, not
with F7.
