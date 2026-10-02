# The explicit-diffusion limit of the time step

`./run_gate.sh` (CPU, about a minute). A decaying Beltrami flow at Re = 1 in a
triply periodic box (`box.ini`): the Courant number stays below 0.15 and the
flow decays to round-off, so diffusion is the only thing that can go unstable.

## What is measured

The RK3 step is stable on the negative real axis down to z = −2.5127 and the
extreme eigenvalue of the discrete Laplacian is −4 Σ_d ν/h_d², so explicit
diffusion needs

    dt × Σ_d ν/h_d²  ≤  0.628 .

`[time] pecletmax` bounds `dt × ν/h²` of the finest SINGLE direction (the
largest over all cells and directions). How much of the 0.628 that leaves
depends on how many directions of a cell are that fine:

| grid | Σ_d ν/h_d² | limit on `dt ν/h²` | stable leg | unstable leg |
|---|---|---|---|---|
| 16 × 16 × 16 | 3 ν/h² | 0.2094 | 0.205: decays to 1.5e-8 | 0.215: NaN by step 500 |
| 16 × 16 × 4 | 2.06 ν/h² | 0.3046 | 0.295: decays to 1.8e-10 | 0.315: NaN by step 400 |
| 16 × 4 × 4 | 1.13 ν/h² | 0.5584 | 0.545: decays to 2.3e-16 | 0.570: NaN by step 600 |

(fixed step, 600 steps, Chebyshev 12; the unstable mode starts from round-off.)
The limit is the sum, bracketed to ±2.5 % on three grids. **`pecletmax = 0.5`
is safe only where one direction dominates the sum** — a wall-normal line much
finer than the other two, which is every channel and boundary-layer case — or
where the Courant limit binds first. On a grid with two or three equally fine
directions (an isotropic refined patch, a body-fitted finest level) and a
step set by diffusion it is beyond the limit.

With the adaptive step the instability does not necessarily end in NaN: the
16³ leg at `pecletmax = 0.5` (effective sum number 1.5) settles into a bounded
state of amplitude ~20 where the Courant limiter, fed by the growing velocity,
holds the step at the stability boundary. That is the "sustained O(0.65)
disturbance inside the patch" of `docs/next_session_outlet.md` (side
findings), reproduced without a body, an outlet or a level jump.

## What the solver says

`pecletmax` keeps its meaning (the limiter is unchanged, fields are
bit-identical). The solver reports, in the manner of the `cfl:` line:

     peclet: max single-direction diffusion rate  6.485E+00 (pecletmax bounds dt x this); worst-case SUM over directions  1.945E+01 =  3.00 x the direction max
          RK3 with explicit diffusion is stable for dt x SUM < 0.628: pecletmax x ratio =  1.50 is the effective sum bound (molecular viscosity only)

and, once, at the first step taken beyond the limit (fixed or adaptive step):

     WARNING: dt x sum_d(nu/h_d^2) =  1.500 exceeds the RK3 limit of explicit diffusion 0.628 at step 1
              unstable unless the field stays exactly uniform in enough directions: a 3D field on this grid needs pecletmax <= 0.209 (or a smaller dt / dtmax)

The sum is the molecular one (with scalars, of the most diffusive scalar): an
eddy viscosity or a conjugate body tightens the true limit further. It is the
limit of the GRID. A field that is exactly uniform in a direction never
excites the modes that vary along it, so a 1D or 2D test problem runs on
beyond it: `validation/scalar/conduction.ini` (4 × 16 × 4, `pecletmax = 0.4`,
sum number 1.19) warns and passes its gate, because its field varies in y
only and the per-column arithmetic keeps it so to the last bit. That is the
one case of the 7- and 9-case suites and the outlet suite that warns.

## Gate

Six fixed-step legs (the table: decay and no warning below, blow-up and the
warning above) and two adaptive legs on 16³ (`pecletmax` 0.5: effective sum
bound 1.50 reported, warning, blow-up; 0.2: 0.60, no warning, decay).
Recorded 2026-10-02: 8/8 PASS.
