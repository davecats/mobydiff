# penalization — the exact penalization factor, isolated (numerics review F7)

A triply periodic 8^3 box with UNIFORM velocity, a uniform body force `f`
and a UNIFORM, moderate penalization coefficient `lambda` on every velocity
DOF. Convection, diffusion and the pressure correction vanish identically,
so the solver integrates the one thing the factor is for:

    du/dt = f - lambda u,     u(t) = u0 e^{-lambda t} + (f/lambda)(1 - e^{-lambda t}).

`run_decay.sh` prepares the case file from `decay.ini` (its analytic wall
exists only so that `moby_prepare` writes coefficient tiles), overwrites the
tiles with the uniform `lambda`, and runs T = 0.2 at dt = 0.04 / 2^k,
k = 0..3, for `lambda = 20` (lambda dt_gamma up to 0.43: the range of a
graded cut cell) and `lambda = 200` (up to 4.3).

    [SOLVER=../../build_cpu/moby_solve] [REF=/abs/older/moby_solve] ./run_decay.sh

**Gate:** `|u(T) - exact| <= 1e-13` at every dt and both `lambda`, field
uniform, v = w = 0. It is an identity, not a convergence test: over a
substage the frozen right-hand side is exactly what
`penal_incr_factor`/`penal_state_factor` (ibm.f90) integrate, and
`alpha + beta = gamma` in every RK substage.

Recorded 2026-10-01 (CPU and GPU builds give the same digits):

| `lambda = 20` | dt 0.04 | 0.02 | 0.01 | 0.005 | observed order |
|---|---|---|---|---|---|
| exact factor (step 9) | 1.4e-17 | 4.2e-17 | 5.6e-17 | 5.6e-17 | — (round-off) |
| implicit Euler (`~/step8_ref_binaries`, `0811811`) | 1.23e-2 | 5.97e-3 | 2.93e-3 | 1.45e-3 | 1.04, 1.03, 1.01 |

At `lambda = 200` the solution has reached its steady value `f/lambda` by
T = 0.2 (lambda T = 40): the exact factor reads 8.7e-19, and the old one
9.8e-9 → 2.5e-15 — both factors share the steady fixed point, which is why
the steady cylinder drag and the developed IBM channel do not move with F7
(`../cylinder/README.md`, `../channel_interface/les_ibm/README.md`).

The factor functions themselves are unit-tested against `expm1` reference
values by `build_cpu/penalization_test` (`src/test_penalization.f90`).
