# Turbulent boundary layer (precursor: laminar Blasius)

Zero-pressure-gradient flat-plate boundary layer. This first stage is the
quasi-2D **laminar Blasius** case that validates the inflow/outflow setup;
the turbulent 3D extension builds on it.

Three variants of the **top** boundary condition live in subdirectories,
sharing all of the physics, grid, and the Python tooling in this parent
directory:

- **`topbc_outlet/`** — outlet: zero-gradient velocity + **Dirichlet p = 0**.
  The entrainment simply leaves; the Dirichlet p pins the freestream to ZPG.
- **`topbc_displacement/`** — Neumann u and p, v **prescribed** from the
  x-varying Blasius entrainment.
- **`topbc_dirichletuv/`** — full Blasius velocity Dirichlet (u = U_inf AND
  v = entrainment), Neumann p.

Key result (details + mechanism in each subdir's README): **the freestream
is governed by the top PRESSURE, not by imposing velocity there.**

| top BC              | p on top   | U_e     | theta match | note |
|---------------------|------------|---------|-------------|------|
| outlet              | Dirichlet 0| 1.0000  | **1.1%**    | ZPG pinned; entrainment ~9% low aloft |
| displacement        | Neumann    | 1.011   | 9%          | favorable dp/dx drifts the edge velocity |
| dirichletuv         | Neumann    | 1.004*  | 26%         | dp/dx persists -> interior freestream bump |

*top-face value; the interior freestream overshoots to ~1.010. Only the
outlet top (Dirichlet p) reproduces Blasius; both Neumann-p tops leak a
weak favorable dp/dx (~-4e-5) that the velocity conditions cannot cancel.

Each subdirectory has its own `blasius2d.ini`, `template.ini`, README and
result figures; the executables and scripts are referenced from the parent.

## Nondimensionalization

Lengths with the **inlet momentum thickness** theta_in, velocities with
U_inf, so `[flow] re` is the inlet momentum-thickness Reynolds number
Re_theta,in. Both cases use Re_theta,in = 100: the virtual plate origin
sits at x_v = Re_theta/beta^2 = 226.7 upstream of the inlet
(beta = 2 f''(0) = 0.664115), delta99_in ~ 7.5, and Re_theta at the outlet
(x = 400) is ~166 — comfortably below the Tollmien-Schlichting threshold,
so the flow stays steady laminar.

## Boundary conditions (shared)

- `x_min`: **inlet**, Blasius similarity profile for u AND v (the
  entrainment component) via `_profile = blasius`. The `_value` key of every
  blasius row is U_inf; `[boundary] blasius_theta` is theta_in. The solver
  shooting-solves the Blasius ODE at init and evaluates
  u = U_inf f'(eta), v = U_inf beta (eta f' - f)/(2 Re_theta) at
  eta = beta y/theta.
- `x_max`: **outlet** (zero-gradient velocity, Dirichlet p = 0) — the single
  compliant boundary + pressure pin in both cases.
- `y_min`: no-slip **wall** (the plate).
- `y_max`: the variant (see the subdirectories and the table above).
- `z`: periodic (quasi-2D, 4 cells).

The `blasius` value profile also covers the **top (y) face** for v: it
imposes the x-varying entrainment at the domain top y = leng(2), station
x + x_v (used by `topbc_displacement`). A y-face blasius row is v-only.

The y grid is the one-sided `natural` stretching (new `[grid.y] one_sided`
key): dy_wall ~ 0.18 theta, 31 points inside the inlet delta99. NOTE: the
`one_sided` flag is not stored in snapshot metadata — a restart ini must
keep the full `[grid.y]` section.

## Run and verify

From a case subdirectory (paths reach the parent scripts and the repo build):

```bash
cd topbc_outlet          # or topbc_displacement
mpirun -n 4 ../../../build_cpu/moby_solve template.ini   # mint template_1.h5
python3 ../make_blasius_ic.py                             # -> IC_blasius.h5
mpirun -n 4 ../../../build_cpu/moby_solve blasius2d.ini   # restarts from the IC
python3 ../compare_blasius.py blasius2d_4001.h5 --plot blasius.png
python3 ../plot_fields.py     blasius2d_4001.h5 --out blasius2d.png
```

Starting from the analytic Blasius field means the run only has to HOLD
the solution — the gate then measures pure discretization + BC error, and
the impulsive-start transient (many flow-throughs of washout) is skipped.

`plot_fields.py` draws the z-averaged u, v, p fields over the x-y plane
(`blasius2d.png`; y clipped to `--ymax`, default 40).

`compare_blasius.py` solves the Blasius ODE independently (RK4 + secant
shooting; beware np.interp clamping — beyond the table it must follow the
asymptote f = eta - 1.7208, which is what `blasius_eval` does) and, at
stations x/lx <= 0.7 (the last ~15% is the outlet influence zone), gates
the measured momentum thickness (< 2% vs theta(x) = beta sqrt(nu (x +
x_v)/U)), shape factor H (< 2% vs 2.5911), u profile (< 1% of U_e) and
the v profile inside the layer, eta <= 8 (< 15% of the entrainment
scale); the far-aloft deficit is reported per station as `v_top`.

`check_convergence.py conv_*.h5 --plot convergence.png` verifies steady
state: run from the IC with frequent output (`field_interval = 200`) and it
reports, per consecutive snapshot pair, the field drift RATE
max|u(t2)-u(t1)|/(t2-t1) with the mid-domain theta error + interior
max|div|. Converged = the drift rate falls to a floor and the gate metrics
are stationary (details in `topbc_outlet/README.md`).

## FINDING: Chebyshev-Jacobi is unstable in this configuration

With `[pressure] accel = chebyshev`, this case grows a **2-dx pressure
mode** (peak near y ~ 23, above the BL edge): exponential with e-folding
~36 t.u. from round-off, independent of the IC (uniform or analytic
Blasius), until it saturates at O(0.2) u-wiggles and O(1) interior
divergence — the projection collapses. **Plain damped Jacobi (`accel`
unset, sor = 0.8) is dead stable** over the same horizon (p-mode flat at
1.6e-7 for 1500 t.u.), which is why these inis do not set `accel`.

Isolation runs, all restarted from the analytic IC, t = 1500 horizon
(p-mode rms: stable = flat at ~1.5e-7, ringing = saturated ~0.3):

| accel     | niter | dtmax | top       | result |
|-----------|-------|-------|-----------|--------|
| chebyshev | 6     | 0.5   | outlet    | RINGS  |
| chebyshev | 6     | 0.5   | free-slip | RINGS  |
| chebyshev | 60    | 0.5   | outlet    | stable |
| chebyshev | 6     | 0.25  | outlet    | stable |
| (none)    | 6     | 0.5   | outlet    | stable |

So the pump needs chebyshev AND the sloppy niter = 6 projection AND a
large-enough dt -- a per-STEP resonance: the 6-iteration Chebyshev
residual polynomial OSCILLATES IN SIGN across the eigenvalue interval
(the projection over/under-shoots mode by mode), and the sign-flipping
residual left in the field resonates with the step-to-step momentum/BC
update; halving dt detunes it, damped Jacobi's residual (1 - 0.8
lambda)^6 is monotone-positive, and at niter = 60 the residual is too
small to pump. The free-slip-top run shows the x_max outlet alone
suffices (the pn family is velocity-neutral without pressure-Dirichlet
faces -- the A2 cylinder caveat documents exactly that neutral mode).
Why not seen in the validated chebyshev cases (channels, freestream
gates, cylinder): e-fold 36 t.u. is invisible over their few-hundred-
t.u. horizons. Solver-level root-cause (which eigenmodes, interval
tuning) is OPEN -- until then, on long steady runs with Dirichlet-p
outlets use plain damped Jacobi, or crank niter, or shrink dt.
