# The convective limit of the time step

`./run_gate.sh` (CPU, about a minute). A uniform flow in a triply periodic
box (`box.ini`) at Re = 1e6, seeded with a random perturbation of 1e-6
(`seed_ic.py`, written into a solver-minted template): convection is the only
thing that can go unstable, and the perturbation is small enough to stay
linear.

## What is measured

Central (skew-symmetric) convection linearised about a uniform flow has the
purely imaginary spectrum `i Σ_d u_d sin(k_d h_d)/h_d`, of largest modulus
`Σ_d |u_d|/h_d` (the 4h wave in every direction the flow has a component in),
and the RK3 step is stable on the imaginary axis up to √3. So

    dt × Σ_d |u_d|/h_d  ≤  1.732 .

`[time] cflmax` bounds exactly this sum (since 2026-10-02). It used to bound
the largest SINGLE component `dt |u_d|/h_d`, whose limit depends on the
direction of the flow:

| flow | Σ_d \|u_d\|/h_d | limit on the component number | below (sum 1.65) | above (sum 1.80) |
|---|---|---|---|---|
| (1, 0, 0) | 1 \|u\|/h | 1.732 | perturbation never above its seed | NaN |
| (1, 1, 0) | 2 \|u\|/h | 0.866 | 1.3 × seed | NaN |
| (1, 1, 1) | 3 \|u\|/h | 0.577 | 1.5 × seed | NaN |

(fixed step, 600 steps; growth per step of the worst mode at 1.80 is 1.034.)
Under the old meaning `cflmax = 0.8` was above the limit for a flow along the
diagonal (sum 2.4) and for most of the range in between; it held in practice
because the ratio sum / largest component of a turbulent flow is 1.2–1.5
(1.26 and 1.44 printed by the boundary-layer cases), which made the usual 0.8
an effective 1.0–1.2.

## What changed in the solver (2026-10-02)

`get_timestep_rates` reduces `max over cells of Σ_d |u_d|/h_d` (the three
staggered faces of the cell's index) and `cflmax` bounds `dt` times it. The
per-component report of the init log (`cfl:` line, review finding F8) is gone
with the thing it reported; a `cflmax` above √3 is announced at init.

**What moves.** A run bound by `dtmax` or by the diffusion limit is
bit-identical. A Courant-bound run takes a smaller step by the ratio sum /
largest component of the cell that binds. **The inis were not changed**:
`cflmax = 0.8` is now 46 % of the limit, where it used to sit at 60–70 % of
it in a turbulent flow. `cflmax = 1.2` restores the old step with a 30 %
margin. The two limits are bounded separately; RK3's stability region is not
a rectangle, and near `pecletmax = 0.5` the convective margin is smaller
(z = −2 + 1.2 i has |R| = 0.99).

## Gate

Six fixed-step legs (the table) and six adaptive legs, the same two values
for the three flows (`cflmax` 1.60: the perturbation stays below twice its
seed; 1.85: it grows by more than 1e3 and the solver announces the value at
init). Recorded 2026-10-02: 12/12 PASS.
