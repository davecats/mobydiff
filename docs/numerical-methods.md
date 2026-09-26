# Numerical methods

This page describes the discretization and algorithms `mobydiff` uses. The authoritative
design records for the block-refinement machinery are in the
[design notes](index.md#design-notes--internal-history).

## Governing equations

`mobydiff` integrates the incompressible Navier–Stokes equations, non-dimensionalized by a
reference length and velocity so that viscosity enters as $1/\mathrm{Re}$:

```math
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\,\mathbf{u} = -\nabla p + \frac{1}{\mathrm{Re}}\,\nabla^2\mathbf{u} + \mathbf{f}, \qquad \nabla\cdot\mathbf{u} = 0 .
```

Here $p$ is the kinematic pressure and $\mathbf{f}$ collects the constant forcing
`[flow] forcing_*` and the optional spatially varying `[force]` field.

## Spatial discretization

The domain is a Cartesian box discretized with **second-order finite differences on a
staggered (MAC) grid**: pressure lives at cell centres and each velocity component lives on
the corresponding cell face. Staggering couples pressure and velocity on the tightest
stencil and avoids odd–even (checkerboard) pressure decoupling.

The momentum convection term is discretized in **skew-symmetric** form (the average of the
divergence and advective forms), which conserves kinetic energy discretely and is what
keeps the 2:1 interfaces stable; it is hardwired.

Each direction can be stretched independently (`[grid.x/y/z] distribution`):

- `uniform` — constant spacing.
- `cosine` / `tanh` — symmetric clustering toward the boundaries.
- `natural` — the Pirozzoli–Orlandi near-wall stretching used for wall-bounded turbulence,
  with the first off-wall spacing set in wall units (`natural_dyw_plus`); `one_sided`
  clusters at the low wall only.
- `blayer` — natural clustering up to `outer_height`, geometric coarsening above (boundary
  layers); `geometric`; or an explicit node line from `nodes_file`.

Metric terms from the stretching are carried per cell so the difference operators stay
second order on the non-uniform mesh.

## Time integration

Time advancement uses an **explicit low-storage three-stage Runge–Kutta (RK3)** scheme. The
convective and viscous terms are advanced explicitly; incompressibility is enforced at each
stage by a pressure projection. The immersed-boundary penalization is treated implicitly
(see [below](#immersed-boundary-method)), so the geometry adds no time-step restriction.

The step size adapts to stability limits set in `[time]`: it is the largest `dt` satisfying
the convective **CFL** limit (`cflmax`) and the viscous **Péclet** limit (`pecletmax`),
capped by `dtmax`.

## Pressure projection

Each RK stage produces an intermediate velocity $\mathbf{u}^{\ast}$ that is not divergence
free. A projection removes its divergence: solve a variable-coefficient Poisson equation for
a pressure correction $\phi$ and subtract its gradient,

```math
\nabla\cdot\left(\frac{1}{\rho}\nabla\phi\right) = \frac{1}{\Delta t}\,\nabla\cdot\mathbf{u}^{\ast}, \qquad \mathbf{u}^{n+1} = \mathbf{u}^{\ast} - \Delta t\,\nabla\phi .
```

The Poisson problem is solved **iteratively with a damped-Jacobi smoother**, optionally
accelerated by a **Chebyshev–Jacobi** polynomial iteration (`[pressure] accel = chebyshev`).
The Jacobi diagonal absorbs the inhomogeneous coefficients from grid stretching and from
the refinement interfaces, which is what lets a single global iteration stay robust across a
non-uniform, block-refined mesh. The discrete operator is kept **symmetric and positive
definite (SPD)** everywhere, including at 2:1 interfaces (composite face-gradient stencil
with conservative copy reconciliation) — SPD is what keeps the Chebyshev iteration stable
across the interface.

- `niter` sets the smoother iterations per projection.
- `sor` is the damping factor; plain damped Jacobi diverges above ≈ 0.8.
- `cheb_lmin` / `cheb_lmax` bound the operator spectrum for the Chebyshev polynomial and are
  auto-derived when left at their `-1` defaults.

> The projection replaced an earlier coupled red–black SOR scheme as the default. Red–black
> SOR remains selectable (`[pressure] solver = redblack`, with over-relaxation `sor ≈ 1.5`)
> and also runs across 2:1 interfaces; it cannot be combined with Chebyshev acceleration.

## Block-structured grid and 2:1 refinement

The grid is organized into **equal-size blocks** (BCM-style, after Nakahashi & Kim 2004 and
Jansson et al. 2019). A block is an `nb × nb × nb` box of cells; `[blocks] nb` sets the edge
length (even, ≥ 4, and must divide the global grid). Blocks are numbered along a **Z-order
(Morton) space-filling curve** and split linearly over the MPI ranks. Because each block
redundantly sweeps its open halo layer with its owner, the results are **exactly independent
of the block count and the number of ranks**.

Two capabilities build on the block layout:

- **Solid-block removal** (`remove_solid`, default on) — a block fully buried inside an
  immersed body (solid at cell centres and all staggered locations, dilated by one halo
  cell) is dropped from the computation. Its faces become exact zero-flux walls.
- **2:1 local refinement** — refinement is requested by box (`refine`, up to 16 boxes, with
  `refine_levels` levels) or by geometry (`refine_body`, refining blocks that touch the
  immersed surface plus a one-block buffer, optionally capped by `refine_body_levels` and
  raised locally by `refine_body_box`). Neighbouring blocks differ by at most one level
  (2:1 balancing). Refined node lines come from midpoint subdivision, so a fully refined
  dyadic region is bitwise identical to running at the doubled resolution. With
  `refine_dims = xz` the blocks refine in x and z only (quadtree) and y keeps one global
  node line — the natural choice for quasi-2D wall-bounded and airfoil cases.

### The 2:1 interface

At an interface between a coarse and a fine block, halo data is transferred conservatively:
restriction averages the fine-side face values feeding a coarse cell, and prolongation
injects the covering coarse value into the fine halo. The transfer is arranged so the wire
always carries destination-point values, and coarse faces are fed by up to four fine
sub-faces in a fixed child order.

The production interface uses an **energy-conserving constant-½ velocity transfer**. This is
a deliberate order-for-stability trade: the truncation-optimal higher-order reconstruction
breaks interface energy conservation and destabilizes turbulent runs, whereas constant-½
keeps the interface stable and conservative. The residual cost is a mild, localized
dissipation of the smallest fine-scale fluctuations crossing into the coarse mesh — a loss,
not a spurious band. This treatment is validated in developed turbulent channel flow,
including across edges and corners and in combination with LES; details and the accepted
residuals are in the [design notes](index.md#design-notes--internal-history).

## Immersed boundary method

Solid geometry is imposed by **volume penalization**. A per-cell coefficient field forces the
velocity toward zero inside the body through a source term that is treated **implicitly** in
the momentum update (roughly $\mu = 1/(1 + \Delta t\,\text{coef})$), so an arbitrarily strong
penalization introduces **no time-step restriction**.

Analytic walls (`[ibm] wall_shape`) are classified inline. STL geometry is classified by the
MPI-parallel `moby_prepare`, which runs the solver's own grid, block and coefficient code
and writes the case file the solver reads: per-leaf coefficient tiles (graded
sharp-interface coefficients in cut cells), the per-level touch/buried block masks, and the
wall distance used by the RANS models. See [Running](running.md#preprocessing-immersed-bodies-moby_prepare)
and the [`sailplane` tutorial](tutorials.md).

## Large-eddy simulation

An optional subgrid-scale model (`[les] model`) adds a turbulent eddy viscosity to the
diffusion term. Both a classical **Smagorinsky** model and the **WALE** model (which gives
the correct near-wall $\nu_t \to 0$ scaling without dynamic procedures) are available. The
model is **IBM-aware** by default (`ibm_aware`): the subgrid viscosity is zeroed inside solid
cells so the model does not read the penalized velocity drop as resolved strain. The eddy
viscosity steps by the physical filter-width ratio across 2:1 interfaces with no spurious
band, and this has been validated across refinement interfaces and IBM walls.

## RANS and IDDES

`[turbulence] model = rans` advances the **k-ω SST** model (`[rans]`): resolved walls or
wall functions, optional **γ–Re_θt transition** (Langtry–Menter, with separation-induced
transition), wall distance from the case file or computed from the analytic indicator.
The turbulence scalars use first-order upwind convection. `model = iddes` blends SST near
walls with the `[les]` SGS model away from them (Gritskevich SST-IDDES shielding and
elevating functions).

## Passive scalars and conjugate heat transfer

Any number of passive scalars (`[scalar]`, `[scalar.N]`) are transported with the same
staggered second-order scheme (central convection, no upwind option; eddy diffusivity
ν_t/Pr_t from whichever turbulence model is active). At an immersed body a scalar is held
at a value (`dirichlet`), insulated (`adiabatic`), or — with `ibm_wall = conjugate` —
solved inside the solid as well, with its own conductivity and capacity, coupled through a
harmonic-mean face coefficient at the cut faces. The derivation is in
`docs/conjugate/conjugate_ibm.tex`.
