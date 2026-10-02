# The explicit-diffusion limit of the time step

`./run_gate.sh` (CPU, about a minute). A decaying Beltrami flow at Re = 1 in a
triply periodic box (`box.ini`): the Courant number stays below 0.15 and the
flow decays to round-off, so diffusion is the only thing that can go unstable.

## What is measured

The RK3 step is stable on the negative real axis down to z = −2.5127 and the
extreme eigenvalue of the discrete Laplacian is −4 Σ_d ν/h_d², so explicit
diffusion needs

    dt × Σ_d ν/h_d²  ≤  0.628 .

`[time] pecletmax` bounds exactly this sum (since 2026-10-02), so one value
below 0.628 is safe on any grid. It used to bound `dt × ν/h²` of the finest
SINGLE direction, and how much of the 0.628 that left depended on how many
directions of a cell are that fine:

| grid | Σ_d ν/h_d² | limit on `dt ν/h²` | stable leg | unstable leg |
|---|---|---|---|---|
| 16 × 16 × 16 | 3 ν/h² | 0.2094 | 0.205: decays to 1.5e-8 | 0.215: NaN by step 500 |
| 16 × 16 × 4 | 2.06 ν/h² | 0.3046 | 0.295: decays to 1.8e-10 | 0.315: NaN by step 400 |
| 16 × 4 × 4 | 1.13 ν/h² | 0.5584 | 0.545: decays to 2.3e-16 | 0.570: NaN by step 600 |

(fixed step, 600 steps, Chebyshev 12; the unstable mode starts from round-off.)
The limit is the sum, bracketed to ±2.5 % on three grids. Under the old
meaning `pecletmax = 0.5` was safe only where one direction dominates the sum
(a wall-normal line much finer than the other two: every channel and
boundary-layer case) or where the Courant limit binds first; on a grid with
two or three equally fine directions and a step set by diffusion it was
beyond the limit. There the adaptive step does not necessarily end in NaN:
an isotropic run at a sum number of 1.5 settled into a bounded state of
amplitude ~20 where the Courant limiter, fed by the growing velocity, holds
the step at the stability boundary. That is the "sustained O(0.65)
disturbance inside the patch" of `docs/next_session_outlet.md` (side
findings), reproduced without a body, an outlet or a level jump.

## What changed in the solver (2026-10-02)

- `pecletmax` bounds `dt × max over cells of Σ_d ν_eff/h_d²` (molecular,
  eddy viscosity and the most diffusive scalar alike).
- The conjugate body's rate was already half the Gershgorin diagonal at cut
  cells; it is now that everywhere, which IS the sum for a regular cell. The
  cut-cell special case is gone (its step is unchanged).
- A `pecletmax` above 0.628 is announced at init. A step beyond the limit,
  fixed or adaptive, is warned about once:

      WARNING: dt x sum_d(nu/h_d^2) =  0.660 exceeds the RK3 limit of explicit diffusion 0.628 at step 1
               unstable unless the field stays exactly uniform in enough directions: lower dt, or set pecletmax below the limit

  It is the limit of the GRID: a field exactly uniform in a direction never
  excites the modes that vary along it, so a 1D or 2D test problem with a
  fixed step can run on beyond it.

**What moves.** A run whose step is set by `dtmax` or by the Courant limit
is bit-identical. A diffusion-bound run takes a smaller step by the ratio
sum / largest direction of the cell that binds: a few per cent on a
wall-clustered grid, up to 3 on an isotropic one (where the old step was
beyond the limit). The suites say which is which (recorded in
`docs/next_session_outlet.md`, follow-up).

## Gate

Six fixed-step legs (the table: decay and no warning below, blow-up and the
warning above) and six adaptive legs, the same two values on all three grids
(`pecletmax` 0.60: no warning, decay; 0.66: the warning, and the run leaves
the decaying solution). Recorded 2026-10-02: 12/12 PASS.
