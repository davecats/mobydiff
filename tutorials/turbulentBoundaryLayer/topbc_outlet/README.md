# topbc_outlet — outlet top boundary

The top (`y_max`) is an **outlet**: zero-gradient velocity + Dirichlet
p = 0, so the displacement-driven entrainment leaves through the top. See
the parent `../README.md` for the shared setup and the run recipe.

The domain is 100 theta_in tall (~13 delta99 at the inlet), so the induced
acceleration is small; the residual is the entrainment running ~9% below
Blasius aloft (the p = 0 top under-drives it). The `compare_blasius.py`
v gate is taken inside the layer (eta <= 8); the far-aloft deficit shows
up as `v_top` (info only, up to ~-0.19 near the outlet zone) and in the
lower entrainment near the top of `blasius2d.png`.

## Results

- `blasius.png` — u/Ue, v (similarity form), theta growth vs Blasius.
  Gate PASS (t = 2000, steady to ~6e-4): theta <= 1.13%, H <= 0.34%,
  du/Ue <= 2.4e-3, in-layer dv <= 6.6e-2.
- `blasius2d.png` — z-averaged u, v, p over the x-y plane.
- `convergence.png` — steady-state check. The field drift rate falls ~10x
  (3e-5 -> 2.5e-6) over the first ~500 t.u. and plateaus; the theta error
  oscillates in a fixed +/-1% band with no trend (identical at t = 1000 and
  t = 3000); interior divergence flat at ~5e-5. Clincher: the single-step
  drift rate (2.9e-5/t.u.) is ~10x the 100-t.u.-window rate (2.9e-6/t.u.) —
  the per-step changes cancel over long windows, a stationary residual, not
  development.
