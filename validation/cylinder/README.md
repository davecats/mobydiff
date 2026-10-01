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

## The outlet pressure mode, run down (2026-10-01; handout item 3 after step 9)

`run_outlet_mode.sh` repeats the one-time-unit leg of the table above with
different PROJECTIONS, everything else fixed (the same clean-p state, fixed
dt); `outlet_mode.py` analyses the legs. Inputs as for `run_f7_gates.sh`.

    NEW=/abs/build_gpu/moby_solve [OUT=outlet] ./run_outlet_mode.sh [<name>:<k>[:<T>] ...]

**1. It is not the projection.** dt = 5e-3, 200 steps, final state:

| projection | p rms | max \|p\| | p rms: outlet band (x > 15) / lateral bands / body box / rest | div rms left | velocity vs red-black 240: max |
|---|---|---|---|---|---|
| Chebyshev, niter 60 | 0.891 | 3.09 | 2.013 / 1.046 / 0.267 / 0.727 | 1.4e-5 | 1.6e-4 |
| Chebyshev, niter 240 | 0.889 | 3.05 | 2.012 / 1.044 / 0.268 / 0.725 | 1.4e-6 | 1.9e-4 |
| plain Jacobi, niter 60 | 0.901 | 3.11 | 2.027 / 1.060 / 0.268 / 0.736 | 3.5e-4 | 6.0e-3 |
| red-black SOR, niter 60 | 0.894 | 3.04 | 2.014 / 1.053 / 0.270 / 0.730 | 1.7e-4 | 1.1e-3 |
| red-black SOR, niter 240 | 0.890 | 3.05 | 2.013 / 1.045 / 0.270 / 0.725 | 3.4e-5 | (reference) |

Three projections at two iteration counts, residuals a factor 250 apart, give
the same pressure to 1 % and the same velocity to 1e-3. The state is the
converged solution of the discrete problem: not the Chebyshev bounds, not
unconverged long waves. Its 2-dx content is 9.6e-4 in every leg: it is a
smooth field.

**2. What it is.** The stored pressure is ANTISYMMETRIC in y and largest in
the LAST CELL COLUMN: the y-rms of p grows from 0.19 at the inlet to 2.2 at
x = 15.98, half a cell from the face where the pressure is held at zero. Over
ten time units (Chebyshev 60, dt = 5e-3) it spikes at every vortex passage
through the outlet, with alternating sign and growing peaks:

| t | 202 | 203 | 204 | 205 | 206 | 207 | 208 | 209 | 210 | 211 |
|---|---|---|---|---|---|---|---|---|---|---|
| p rms | 0.63 | 0.22 | 0.15 | 0.79 | 0.38 | 0.28 | 0.91 | 0.49 | 0.36 | 0.99 |
| antisymmetric amplitude | +0.53 | −0.19 | +0.12 | −0.68 | +0.33 | −0.24 | +0.79 | −0.44 | +0.32 | −0.86 |

(amplitude: half the difference of the mean pressure of the upper and the
lower quarter of the domain, x < 15). So the one-time-unit numbers of the
table above are a snapshot of an oscillation, taken near a peak.

**3. The mechanism (from the code).** Every substage `apply_bc(...,
outflow_copy = .true.)` RESETS the outlet face to its interior neighbour,
`u_out := u_nx`. The cell on the outlet then needs
`u_out − u_nx = −dx (dv/dy + dw/dz)` for continuity, and the projection
supplies exactly that through the outlet face's correction `2 phi_nx/dx`,
again in every substage, because the next reset throws it away. The stored
pressure is incremental (`p += phi/dt_gamma`), so in the last column it
accumulates `−(dx²/2 dt_gamma) dv/dy` per substage: a time integral with no
`p_face = 0` to anchor it (the outlet face's "predictor" carries no pressure
gradient). Its y-gradient acts on v in that column until `dv/dy` vanishes
there. The discrete outlet is therefore not "do nothing" but "the outflow
must be parallel", enforced by a constraint force that is stiffer the
smaller dt. Measured at the same instant (t = 202.5) from the same state:

| dt | rms dv/dy at x = 15.0 | in the last column | ratio | p y-rms, last column |
|---|---|---|---|---|
| 5.0e-3 | 0.0528 | 0.0120 | 0.23 | 2.2 |
| 2.5e-3 | 0.0532 | 0.0047 | 0.09 | 1.4 |
| 1.25e-3 | 0.0534 | 0.0055 | 0.10 | 2.0 |
| 6.25e-4 | 0.0535 | 0.0042 | 0.08 | 4.1 |

That is why the pollution is not monotone in dt, why the largest
dt-to-dt velocity differences sit at the outlet, and why the steady Re 40
flow is immune (it settles with dv/dy = 0 at the outlet, where phi vanishes).
The same reading covers the Lamb-Oseen exit of `../freestream/` (energy
regrowth during the exit, reflected fraction 2.2e-2) and the factor 2 in its
Poiseuille last-cell pressure.

**4. A prototype that removes it** (branch `proto/outflow-incremental`,
commit `4908e7c`, NOT on main). The outlet face gets a predictor of its own,
built from its neighbour's: the neighbour's increment over the substage, with
the neighbour's pressure gradient taken out and the face's own put in,
against the held outlet pressure:

    u_out* = u_out + (u_nx* − u_nx) + dt_gamma [ dp/dx|face nx − dp/dx|outlet face ]

It is one small kernel over the boundary point list (`save_outflow_gap`)
that fills the per-point constant of the existing `dst = src + C` row before
`momentum()`; cases without an outlet are untouched by construction. A first
version without the pressure terms removed the reset but left the last
column's pressure unanchored (a frozen relic and a drifting outflow
imbalance): the outlet face needs its own pressure gradient.

| | main | prototype |
|---|---|---|
| p rms after one time unit, dt 5e-3 / 6.25e-4 | 0.89 / 1.65 | **0.080 / 0.073** |
| max \|p\| | 3.1 / 6.0 | 0.72 / 0.73 |
| p rms, outlet band | 2.0 / 3.7 | 0.040 / 0.011 |
| p y-rms, last column | 2.2 / 4.1 | 0.000 / 0.000 |
| rms dv/dy, last column over x = 15.0 | 0.23 / 0.08 | 1.12 / 1.12 |
| velocity, dt 5e-3 against 6.25e-4 at T: max | 0.15 | **6.0e-5** |
| ten time units: p rms | 0.15 … 0.99, spiking | 0.072 … 0.075 |
| ten time units: control-volume C_L (t > 203.5) | −0.65 … +0.62, sample-to-sample rms 0.046 | −0.40 … +0.42, 0.015 |
| mean C_D (t > 203.5) | 1.4478 | 1.4464 |
| **100 time units at the production settings** (Chebyshev niter 12, dt 5e-3), last 50: St | 0.1670 | 0.1744 |
| … mean C_D | 1.4474 | 1.4383 |
| … control-volume C_L amplitude | **0.655** | **0.393** |
| … stored p rms over the run | 0.20 … 1.13, periodic | 0.073 … 0.078 |
| `../freestream` oblique | exact | exact |
| `../freestream` Poiseuille: profile dev at x/lx 0.5 / 0.9 | 1.639e-3 / 1.539e-3 | 1.637e-3 / 1.521e-3 |
| Poiseuille: last-cell p (exact: G dx/2 = 1.25e-3) | 2.41e-3 | **1.247e-3** |
| Poiseuille: drift | 1.3e-15 | 5.6e-16 |
| Lamb-Oseen exit, E/E0 at t = 2.5 | 2.2e-2, regrowing (0.049 → 0.228) | **9.8e-5**, monotone |
| `../blasius` (outlets at x_max and at the top), worst theta / H | 1.33 % / 0.45 % | 1.65 % / 0.82 % (gate 2 %) |
| `../blasius` worst du/Ue / dv/v_edge | 2.40e-3 / 0.124 | 1.48e-3 / **0.013** |
| `../blasius` top entrainment `v_top` at x/lx 0.15 … 0.7 | +0.013 +0.002 −0.041 −0.202 | +0.025 +0.030 +0.028 +0.033 |

The body-box pressure (0.27) is the same in both: what the prototype removes
is the part of the stored pressure that was never physical. With it, the
time-convergence statement the F7 study could not make becomes possible: the
two step sizes agree to 6e-5 in velocity after one time unit.

**Not done, and why it is not on main.** It changes the numerics of every
outlet case. Before it ships: the outlet face at RESTART (it is not in the
field file; today it is re-copied, which with the prototype is a one-time
kick — reconstruct it from continuity instead), the penalization factor of a
face next to a body, 2:1 interfaces on an outlet face, and a re-measurement
of every outlet case (this one incl. the Strouhal gate, `../blasius`, the
boundary-layer tutorial, `tutorials/naca/rans`, the sailplane).

**The Re 100 `strouhal` gate on the control-volume series** (the part of
handout item 4b that was open): `run_outlet_mode.sh cheb12:0:100` then
`check_cylinder.py strouhal <OUT>/forces_om_cheb12T100_0.txt`.
MAIN: St 0.1670, mean C_D 1.4474, mean C_L −1.0e-2, C_L amplitude **0.655**
— the old penalization series read 0.168 / 1.448 / 0.51 on the same flow, so
the control-volume lift of main is inflated by the outlet mode (its pressure
sits on the box borders), by 28 % in amplitude. PROTOTYPE: St 0.1744, mean
C_D 1.4383, mean C_L +1.2e-2, C_L amplitude **0.393**. The prototype's flow
is a different one, not only a cleaner statistic: the Strouhal number moves
+4 % and the lift amplitude falls to 0.39 (literature, unbounded: St ~0.165,
C_L' ~0.33; this domain has ~6 % blockage). Which of 0.167 and 0.174 is the
right number for THIS confined domain is not settled here; it needs the
lateral far field varied.
