# topbc_dirichletuv — full-Dirichlet-velocity (pin both) top

The top (`y_max`) prescribes the **full Blasius velocity**: u = U_inf
Dirichlet (pins the freestream/edge velocity) AND v Dirichlet from the
x-varying Blasius entrainment, with Neumann pressure. This is the "pin
both" answer to the `topbc_displacement` tradeoff: prescribing v alone
(with Neumann u) imposed the entrainment but let the freestream drift
because nothing held u = U_inf; pinning u fixes that. See the parent
`../README.md` for the shared setup and run recipe.

## The top BC

```
y_max_patch     = patch        ; generic face -> explicit per-variable rows
y_max_u_type    = dirichlet    ; u = U_inf (edge velocity pinned)
y_max_u_value   = 1.0
y_max_v_type    = dirichlet    ; v = the x-varying Blasius entrainment ...
y_max_v_profile = blasius
y_max_v_value   = 1.0          ; U_inf (the profile reference velocity)
y_max_w_type    = dirichlet    ; w = 0 (quasi-2D)
y_max_w_value   = 0.0
y_max_p_type    = neumann      ; zero-gradient pressure
```

At y = 100 the similarity variable eta ~ 44, where f'(eta) = 1 to machine
precision, so the constant u = U_inf IS the exact Blasius u there. The full
velocity is prescribed, so the projection leaves every top face pinned
(`face_grad` returns 0); the `x_max` outlet remains the compliant boundary
and pressure pin. Because the imposed (u, v) is the exact divergence-free
Blasius trace, the boundary is globally mass-consistent.

## Result: pinning u does NOT fix the drift -- the pressure is the culprit

Well-posed, stable and converged (drift 2.3e-6/t.u.). Pinning u = U_inf
does hold the top-FACE velocity (U_e at the top cell only 1.001-1.004 vs
the displacement top's 1.011), so the "free edge velocity" of the
displacement case is gone. But the Blasius match is the WORST of the three
tops:

| station x/lx | U_e(top) | theta err | H     |
|--------------|----------|-----------|-------|
| 0.15         | 1.0011   | -3.8%     | 2.67  |
| 0.50         | 1.0030   | -18.0%    | 3.03  |
| 0.70         | 1.0036   | -26.4%    | 3.28  |

**Why (diagnosed):** the pressure is still **Neumann** on the top, so the
same weak favorable dp/dx ~ -4e-5 as the displacement case develops
(measured; the p field floats near -1, pinned only at the x_max outlet
edge). Pinning u = U_inf at the top face does not remove that gradient --
it only prevents the TOP from accelerating, so the acceleration is forced
into an interior **freestream bump**: u overshoots to ~1.010 mid-domain
while the top face is held at 1.000 (overshoot grows +0.001 -> +0.007
downstream). That non-uniform freestream corrupts the integral
diagnostics -- in the overshoot region u/U_e > 1 so the deficit integrand
(1 - u/U_e) goes negative, which is what drives the spurious theta
thinning (-26%) and the inflated H (3.28). The near-wall u profile itself
is still within ~2-3% of Blasius; the damage is all in the outer flow.

**Takeaway across the three tops:** the freestream is governed by the
PRESSURE, not by imposing velocity at the top. Only `../topbc_outlet`
(Dirichlet p = 0 along the top) pins the freestream to ZPG and reproduces
Blasius (theta 1.1%). Both Neumann-p tops leak a favorable dp/dx: the
displacement top lets it drift the edge velocity (+1%), and pinning u on
top of that (this case) merely converts the drift into an interior bump
(worse). The right way to pin BOTH the entrainment and ZPG would be a
Dirichlet-p top with a prescribed v -- i.e. keep the outlet's pressure
pin and add the entrainment -- not a velocity-only top.

## Figures

- `blasius.png` — u/Ue, v (similarity form), theta growth vs Blasius.
- `blasius2d.png` — z-averaged u, v, p over the x-y plane.
