# penalization — the penalization factor, isolated (numerics review F7)

A triply periodic 8^3 box with UNIFORM velocity, a uniform body force `f`
and a UNIFORM, moderate penalization coefficient `lambda` on every velocity
DOF, plus one passive scalar with a uniform `source`, `ibm_wall = dirichlet`
and the same rate (`coef_p = lambda Pr`). Convection, diffusion and the
pressure correction vanish identically, so the solver integrates the one
thing the factors are for:

    du/dt = f - lambda u,              u(t) = f/lambda + (u0 - f/lambda) e^{-lambda t}
    ds/dt = S - lambda (s - s_body),   s(t) = s_inf + (s0 - s_inf) e^{-lambda t},  s_inf = s_body + S/lambda

`run_decay.sh` prepares the case file from `decay.ini` (its analytic wall
exists only so that `moby_prepare` writes coefficient tiles), overwrites the
tiles with the uniform `lambda`, and runs T = 0.2 at dt = 0.04 / 2^k,
k = 0..3, for `lambda = 20` (lambda dt_gamma up to 0.43: the range of a
graded cut cell) and `lambda = 200` (up to 4.3).

    [SOLVER=../../build_cpu/moby_solve] [REF=/abs/older/moby_solve] ./run_decay.sh

**The factors** (`ibm.f90 penal_incr_factor` / `penal_state_factor` /
`penal_state_minus_incr`, x = lambda dt_gamma): the AMPHIBIOUS rational form,
`e^x` replaced by its third-order Taylor polynomial `P3`,

    state = 1/P3,    incr = (1 + x/2 + x^2/6)/P3,    P3 = 1 + x + x^2/2 + x^3/6,

for the velocity (`ibm%mu = incr`, the state completed by
`step.f90 add_penalization_state_correction`) and for the Dirichlet scalar
(`scalar.f90`, inline).

**Gate**, for the velocity AND the scalar: at `lambda = 20` every observed
order >= 2.7 and `|err| <= 2e-4` at dt = 0.04; at `lambda = 200` (the steady
fixed point, lambda T = 40) `|err| <= 1e-13` at every dt; the fields uniform
to round-off, v = w = 0.

Recorded 2026-10-01, `lambda = 20`, `|value(T) - exact|` (CPU and GPU builds
give the same digits):

| velocity | dt 0.04 | 0.02 | 0.01 | 0.005 | observed order |
|---|---|---|---|---|---|
| rational P3 (this tree) | 1.01e-4 | 1.48e-5 | 2.00e-6 | 2.61e-7 | 2.77, 2.88, 2.94 |
| exact exponential (step 9, `~/step9_ref_binaries`, `878868f`) | 1.4e-17 | 4.2e-17 | 5.6e-17 | 5.6e-17 | — (round-off) |
| implicit Euler (`~/step8_ref_binaries`, `0811811`) | 1.23e-2 | 5.97e-3 | 2.93e-3 | 1.45e-3 | 1.04, 1.03, 1.01 |

| scalar | dt 0.04 | 0.02 | 0.01 | 0.005 | observed order |
|---|---|---|---|---|---|
| rational P3 (this tree) | 7.46e-5 | 1.09e-5 | 1.48e-6 | 1.92e-7 | 2.77, 2.88, 2.94 |
| implicit Euler (step 8 AND step 9: the scalar was not touched by step 9) | 9.06e-3 | 4.40e-3 | 2.16e-3 | 1.07e-3 | 1.04, 1.03, 1.01 |

At `lambda = 200` the solution has reached its steady value by T = 0.2: the
rational factors read 2.9e-15 (velocity) and 2.2e-15 (scalar) at dt = 0.04
and round-off below, implicit Euler 9.8e-9 → 2.5e-15. All three factor
families share the steady fixed point `lambda u = R` (`state + x incr = 1`),
which is why the steady cylinder drag and the developed IBM channel do not
move between them (`../cylinder/README.md`,
`../channel_interface/les_ibm/README.md`).

**Why the rational form and not the exact exponential.** This gate is the
one problem on which the exponential is exact: a right-hand side frozen over
the substage. As soon as the forcing varies in time, the RK3 coupling sets
the error and the two are indistinguishable (single-DOF model of our RK3
substages, `lambda = 20`, `f = cos 2 pi t`: rational 1.1e-4 / 2.5e-5 /
5.9e-6 / 1.4e-6, exponential 1.2e-4 / 2.7e-5 / 6.1e-6 / 1.5e-6, implicit
Euler 1.1e-3 → 1.7e-4, order 1). The rational form costs one divide where
the exponential cost an `exp` per DOF and substage (measured on `les_ibm`,
see the numerics review, step 9), has no small-x series branch, and calls no
libm on the device, so CPU == GPU by construction.

The factor functions themselves are unit-tested against exact rational
reference values by `build_cpu/penalization_test` (`src/test_penalization.f90`).
