# Boundary-layer track (branch `boundaryLayer`)

STATUS 2026-07-19: laminar Blasius precursor DONE and gated
(`tutorials/turbulentBoundaryLayer/`); one solver-level OPEN item (the
Chebyshev finding below) and the 3D/turbulent extension pending.

## What exists

- `[boundary] <x face>_<u|v>_profile = blasius` (boundary.f90): shooting-
  solved similarity table at init (secant on f''(0), RK4, uniform-eta
  lookup with the outer asymptote f = eta - 1.7208 beyond the table);
  u = U_inf f'(eta), v = the entrainment component; eta = beta*y/theta with
  `[boundary] blasius_theta` and beta = 2 f''(0). The `_value` key of every
  blasius row is U_inf; Re_theta = U_inf*theta*[flow] re. x faces, u/v
  only, wall at y = 0 (validated in prepare_blasius_profile).
- `[grid.<d>] one_sided` key for the natural distribution (the
  `natural_one_sided` flag existed but had no config key). NOT persisted in
  snapshot metadata: restart inis must keep the [grid.y] section.
- Case `tutorials/turbulentBoundaryLayer/`: Re_theta,in = 100, 400x100x4
  theta units, 384x160x4 (natural one-sided y: dy_wall 0.18, 31 pts in
  delta99_in), inlet/outlet/wall/outlet faces, plain damped Jacobi.
  Flow: template.ini (mint) -> make_blasius_ic.py (analytic Blasius field
  IC via block_scatter) -> blasius2d.ini (hold the solution, t = 2000) ->
  compare_blasius.py (independent ODE, gates at x/lx <= 0.7: theta < 2%,
  H < 2%, du/Ue < 1e-2, dv/v_edge < 0.15). Validated: theta <= 1.4%,
  H <= 0.8%, du <= 2.4e-3, steady to 4e-4; last 15% of the domain =
  outlet influence zone, entrainment aloft ~9% below Blasius (p = 0 top,
  finite height).

## OPEN: Chebyshev-Jacobi instability with Dirichlet-p outlets

`accel = chebyshev` + niter = 6 in THIS configuration grows a 2-dx
pressure mode from round-off: e-fold ~36 t.u., IC-independent (uniform or
analytic Blasius), saturates at O(0.2) u-wiggles with O(1) interior
divergence (projection collapse). Peak near y ~ 23 (above the BL edge);
pressure leads, velocity follows. Isolation table (analytic IC, t = 1500,
full table in the tutorial README): chebyshev/6/dt 0.5 RINGS with outlet
OR free-slip top; STABLE = plain Jacobi/6 (shipped default), chebyshev/60,
and chebyshev/6/dtmax 0.25 — the pump needs chebyshev AND the sloppy
niter = 6 projection AND a large-enough dt: a per-STEP resonance. The
6-term Chebyshev residual polynomial oscillates in SIGN across
[lmin, lmax] (mode-by-mode over/undershoot) and the sign-flipping
leftover resonates with the per-step momentum/BC map; halving dt detunes
it; damped Jacobi's (1 - 0.8 lambda)^6 is monotone-positive; at niter=60
the residual is too small. The pn family is velocity-neutral without
pressure-Dirichlet outlets (the A2 cylinder caveat is that neutral mode);
the x_max outlet alone makes it velocity-active, and the ~36 t.u. e-fold
was invisible over the few-hundred-t.u. freestream/cylinder gate
horizons. Next: dump the projection residual spectrum per iteration
here; map the stability boundary in (niter, dt); check long-horizon
chebyshev+outlet production runs (cylinder/naca) for latent growth.

## Next steps (the turbulent case)

1. 3D: widen z (lz ~ 30-40 theta, nz accordingly), trip or inflow
   perturbation to transition; or recycle-rescale inflow (new machinery).
2. Higher Re_theta,in (>= 300) so the TS/bypass route is available;
   revisit domain height (delta99 grows like x^0.8 turbulent).
3. Consider [blocks] refinement near the wall (the validated 2:1
   machinery) once the laminar base is 3D.
