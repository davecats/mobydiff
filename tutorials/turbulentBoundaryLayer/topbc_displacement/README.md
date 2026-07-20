# topbc_displacement — displacement (prescribed-entrainment) top

The top (`y_max`) is a **displacement boundary**: Neumann u and p (zero
gradient) with the wall-normal v **prescribed** from the Blasius
entrainment. Instead of letting an outlet guess the top flux, the exact
far-field entrainment velocity is imposed. See the parent `../README.md`
for the shared setup and the run recipe.

## The top BC

```
y_max_patch     = patch        ; generic face -> explicit per-variable rows
y_max_u_type    = neumann      ; slip (zero-gradient tangential u)
y_max_v_type    = dirichlet    ; prescribed normal velocity ...
y_max_v_profile = blasius      ; ... = the x-varying Blasius entrainment
y_max_v_value   = 1.0          ; U_inf (the profile's reference velocity)
y_max_w_type    = neumann      ; quasi-2D, w == 0
y_max_p_type    = neumann      ; zero-gradient pressure
```

The `blasius` value profile, on a y face, imposes v(x, y=leng(2)) =
0.5 U_inf sqrt(theta/(Re_theta (x+x_v))) (eta f' - f) with
eta = leng(2) sqrt(Re_theta/(theta (x+x_v))), x_v = Re_theta theta/beta^2 —
i.e. the exact Blasius wall-normal velocity at the domain top, decaying
like x^-1/2 downstream (~5.7e-3 at the inlet to ~3.4e-3 at the outlet).

**Why the projection stays solvable:** a prescribed-normal-velocity /
Neumann-pressure face is treated by the projection exactly like the inlet
or a wall — `face_grad` returns 0, so the face contributes nothing to the
pressure denominator and its velocity is never corrected (pinned to the BC
value). The `x_max` outlet remains the single compliant boundary and the
pressure pin. Because the imposed top v is the exact (divergence-free)
Blasius entrainment, it is globally mass-consistent with the inlet, and the
outlet simply carries the Blasius outlet profile.

## Result: entrainment imposed, but the freestream is no longer pinned

The BC is well-posed, stable and converged (drift rate 2.3e-6/t.u. at
t = 2000, same as the outlet case), and the imposed top v tracks the
Blasius entrainment to < 0.5% (5.71e-3 vs 5.714e-3 at the inlet, decaying
like x^-1/2 to 3.4e-3 at the outlet). BUT the Blasius match is WORSE than
the outlet top:

| station x/lx | U_e     | theta err | H err  | du/Ue   |
|--------------|---------|-----------|--------|---------|
| 0.15         | 1.0016  | +0.3%     | +0.6%  | 0.7%    |
| 0.50         | 1.0066  | +6.9%     | -0.9%  | 2.8%    |
| 0.70         | 1.0093  | +9.1%     | -1.4%  | 3.5%    |

**Why (diagnosed, not a bug):** the outlet top carried Dirichlet p = 0,
which pinned the freestream pressure to ~0 all along the top -> true ZPG
-> U_e = 1.0000 exactly. This displacement top is **Neumann p**, so the
freestream pressure floats (pinned only at the x_max outlet EDGE). A weak
FAVORABLE gradient dp/dx ~ -4e-5 develops over the domain (measured), which
accelerates the edge velocity ~1% (Bernoulli: dU_e ~ -dx*dp/dx ~ 0.016)
and fills the profile (H just below 2.591, interior v below Blasius in
`blasius.png`), pushing theta ~9% above Blasius downstream. It is a
converged steady state, not a transient.

So there is a real tradeoff: prescribing the entrainment v imposes the
exact wall-normal velocity but, with Neumann pressure and Neumann u, does
NOT pin the freestream velocity to U_inf; the outlet top instead pins the
freestream (ZPG, U_e = const, theta within 1.1%) at the cost of a
self-selected entrainment ~9% under Blasius aloft. Pinning the edge
velocity as well (u = U_inf Dirichlet, the `../topbc_dirichletuv` variant)
does NOT recover ZPG: the favorable dp/dx comes from the Neumann-p top, so
pinning u only converts the drift into an interior freestream bump (worse).
Reproducing Blasius cleanly needs the top PRESSURE pinned -- i.e. the
outlet top -- so the true "pin both" would be a Dirichlet-p top plus a
prescribed v, not a velocity-only top.

## Figures

- `blasius.png` — u/Ue, v (similarity form), theta growth vs Blasius.
- `blasius2d.png` — z-averaged u, v, p over the x-y plane (note the weak
  favorable pressure ramp along the top, absent in the outlet case).
