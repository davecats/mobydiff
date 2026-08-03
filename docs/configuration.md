# Configuration reference

`mobydiff` reads a single INI-style configuration file, passed as the sole command-line
argument. This page documents every section and key.

## File syntax

- Sections are written `[name]`; keys are `key = value`.
- Section and key names are **case-insensitive** (lowercased on read).
- Comments start with `;` or `#` and run to end of line.
- Blank lines are ignored. **Unknown sections and keys are silently ignored** — check
  spelling carefully.
- String values may be single- or double-quoted; the quotes are stripped.
- Booleans accept `true` / `.true.` / `1` / `yes` and `false` / `.false.` / `0` / `no`.

## Required vs. optional

For a fresh run (no restart), these keys are **required** and must be positive:

- `[grid] nx, ny, nz` and `lx, ly, lz`
- `[flow] re`
- `[time] dt`, plus at least one of `nsteps` / `t_final`, and `dtmax > 0`

Everything else is optional and falls back to the defaults listed below. If
`[restart] file` is set, full validation is skipped and the run continues from that field.

---

## `[case]` — case selector

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | `generic` | Flow case: `channel` (honours `[case.channel]`), `airfoil` (honours `[case.airfoil]`) or `generic`. |

## `[case.channel]` — channel-case parameters

Only read when `[case] name = channel`.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `n_walls` | int (1 or 2) | 2 | One- or two-wall channel. |
| `natural_blend_index` (alias `jb`) | real | 40.0 | Pirozzoli–Orlandi wall-normal stretch blend index. |
| `mean_profile_sine_amplitude` | real | 0.0 | Amplitude of a `sin(2πy/Ly)` term added to the mean streamwise profile. |
| `large_disturbance_amplitude` | real | 1.0e-2 | Large-scale initial disturbance amplitude. |
| `small_noise_amplitude` | real | 1.0e-3 | Small random-noise amplitude. |
| `stats_sample_interval` | int | -1 (off) | Steps between statistics samples. |
| `stats_write_interval` | int | -1 (off) | Steps between statistics writes. |
| `stats_file` | string | `channel_stats.h5` | HDF5 statistics output. |
| `runtime_file` | string | `runtimedata.txt` | Runtime log file. |

## `[case.airfoil]` — quasi-2D immersed-body case (airfoil, cylinder, …)

Only read when `[case] name = airfoil`. The case composes its own far field
from patch types: `x_min` and both lift-direction faces become Dirichlet
inlets carrying the freestream, `x_max` a pressure outlet, and the span
direction is periodic. Explicit `[boundary]` keys still win.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `aoa` | real | 0.0 | Angle of attack [deg]. Sets the inlet freestream `(U cos α, U sin α)` in the chord–lift plane and the wind axes the reported coefficients use. |
| `u_inf` | real | 1.0 | Freestream speed. |
| `chord` | real | 1.0 | Reference chord for the coefficients. |
| `span` | `z` or `y` | `z` | Periodic (extrusion) direction. `z`: chord x, lift y. `y`: chord x, lift z — the orientation for `[blocks] refine_dims = xz`, which then never refines the span. |
| `force_sample_interval` | int | 10 | Steps between force samples. `0` (or negative) disables forces entirely — the only way to run without `cv_box`. |
| `cv_box` | 4 reals | — | **Required** unless forces are disabled: `c0 c1 l0 l1`, the control volume for the force budget, in physical coordinates along the chord and lift axes (the span is always the full periodic extent). See below. |
| `steady_tol` | real | 0.0 (off) | Stop the run once the flow is steady — see below. |
| `steady_samples` | int | 3 | Consecutive qualifying samples `steady_tol` needs before stopping. |
| `runtime_file` | string | `forces.txt` | Where `iteration time cl cd` is appended. |

### Forces: the control-volume momentum budget

`C_L`/`C_D` come from a momentum balance over the box, never from the body:

```
F = - d/dt ∫_V u dV - ∮_S [ u (u·n) + (p - p_∞) n - τ·n ] dS
τ = (ν + ν_t)(∇u + ∇uᵀ)
```

Only the four lateral borders are integrated — the box spans the full
periodic extent, so the two span faces cancel. Choosing the box:

- **It must enclose the body**, with the borders in reasonably clean flow.
  A margin of ~1.5 chords is a good default and matches
  `tutorials/naca/rans/postProcess/cv_forces.py --boxes 1.5`, the offline
  twin of this statistic.
- **Tighter is more accurate.** Longer borders accumulate more collocation
  and gradient error: on the Re 40 cylinder the same converged field reads
  C_D = 1.698 / 1.704 / 1.726 / 1.765 for margins from 1.5 D out to 6 D.
- Borders are **snapped at setup** to a face of the coarsest level they
  cross (a coarse node is a node at every finer level), so a border never
  lands on a fine-only face and the box stays closed. The snapped box is
  printed at startup — check it.
- Each physical face is integrated exactly once, by the block on its east
  side at that block's own level, so a border may cross a 2:1 interface
  freely. A uniform flow through the box integrates to zero to round-off
  (exactly zero when no border coincides with a block boundary).
- The reported force is **independent of the rank count and of CPU vs GPU
  by construction** (per-block partials scattered into the global block
  table, exact allreduce, ordered final sum).

Two things to know before trusting the numbers:

- **The budget reads the stored pressure**, which accumulates a large
  velocity-neutral drift over long runs at production `[pressure] niter`
  (the same mode as the known channel `pn` drift). Dynamics are unaffected,
  but a polluted `pn` destroys this statistic — on a 20 000-step Re 40
  cylinder run it returned `C_D = 261`. From a *converged* pressure the
  budget holds: 2000 steps at `niter = 6` kept C_D to ±0.001. If a long run's
  forces look wrong, re-converge the pressure (zero `pn`, restart at high
  `niter` for a few hundred steps) before believing them.
- **The unsteady term is differenced between consecutive samples.** At the
  first sample after a start or restart it is unknown and the budget is
  reported without it; it is also the term that limits accuracy on rapidly
  varying flows (on the shedding Re 100 cylinder the sampled C_L tracks
  the reference to ~20 % of its amplitude, with the mean C_D within 0.12 %).

### Stopping at steady state

`steady_tol > 0` ends the run as soon as the flow stops changing, writes the
final field (the usual end-of-run snapshot) and returns. The criterion is the
budget's own unsteady term, expressed in the units of the reported
coefficients:

```
|2 · d/dt ∫_V u dV| / (U∞² c L_span)   componentwise, take the larger
```

so `steady_tol` reads like a C_L/C_D increment — `1e-4` means the unsteady
term no longer moves the fourth digit of the coefficients. It is only
testable once that term exists, i.e. from the SECOND force sample onward, and
it must hold for `steady_samples` consecutive samples (default 3). The run of
samples matters: in an oscillating flow each component's d/dt passes through
zero twice per period, and a single-sample test could stop a perfectly
unsteady run at a turning point. Taking the larger of the two components
already guards against this — they are out of phase, so on the shedding Re 100
cylinder the measure never falls below 3e-2 — but the consecutive-sample
requirement is the cheap insurance.

The decision is taken from an allreduced quantity, so every rank stops at the
same step. Choosing the tolerance: watch the coefficient trace first, then set
it near the residual wobble you are willing to accept. On the steady Re 40
cylinder `steady_tol = 1e-3` stopped the run at t = 105.4 instead of the
configured t_final = 120.

There is no alternative force statistic. The penalization integral
`F = ∫ coef·u dV` this case used to report was removed: it is exact
bookkeeping only while the solid interior is present, and production runs
remove the buried core (`[blocks] remove_solid`, the default), which puts
its share of the pressure-dominated loading outside that bookkeeping.

## `[grid]` — global grid size and domain

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `nx`, `ny`, `nz` | int | — | Global cell counts. **Required.** |
| `lx`, `ly`, `lz` | real | — | Domain lengths. **Required.** |

## `[grid.x]`, `[grid.y]`, `[grid.z]` — per-axis grid controls

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `distribution` | enum | `uniform` | `uniform`, `cosine`, `tanh`, or `natural` (aliases `pirozzoli_orlandi` / `po`). |
| `stretch` | real | 0.0 | Stretching parameter for `cosine` / `tanh` / `natural`. |
| `natural_dyw_plus` (aliases `dyw_plus`, `dy_wall_plus`, `dyw+`) | real | 0.05 | First off-wall spacing in wall units (natural grid). |
| `n` | int | — | Per-axis alternative to `[grid] nx/ny/nz`. |
| `length` | real | — | Per-axis alternative to `[grid] lx/ly/lz`. |
| `subdivided` | bool | false | Mark the axis as subdivided. |

## `[mpi]` — domain decomposition

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `dims` | 3 ints | `0 0 0` | MPI process grid; `0` lets MPI choose the factor for that axis. |

## `[flow]` — physics and initial condition

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `re` | real | — | Reynolds number. **Required.** |
| `forcing_x/y/z` | real | 0.0 | Constant body force / mean pressure gradient per direction. |
| `initial_u/v/w` | real | 0.0 | Uniform initial velocity components. |
| `initial_noise` | real | 0.0 | Initial random-noise amplitude. |
| `initial` | string | `uniform` | Initial-field mode. |

## `[time]` — time stepping

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `dt` | real | — | Nominal time step. **Required.** |
| `nsteps` | int | 0 | Number of steps (one of `nsteps` / `t_final` must be > 0). |
| `t_final` | real | 0.0 | Final time. |
| `cflmax` | real | 0.0 | Maximum CFL for the adaptive step (≥ 0). |
| `pecletmax` | real | 0.0 | Maximum cell Péclet number (≥ 0). |
| `dtmax` | real | — | Hard cap on `dt`. **Must be > 0.** |

## `[pressure]` — pressure projection / Poisson solver

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `niter` | int | 3 | Damped-Jacobi iterations per projection (≥ 0). |
| `sor` | real | 0.8 | Damping/over-relaxation factor (plain Jacobi diverges above ~0.8). |
| `accel` | enum | `jacobi` (none) | `chebyshev` / `cheb` enables Chebyshev–Jacobi acceleration. |
| `cheb_lmin` | real | -1.0 (auto) | Chebyshev lower eigenvalue bound. |
| `cheb_lmax` | real | -1.0 (auto ≈ 2.0) | Chebyshev upper eigenvalue bound. |

## `[boundary]` — boundary conditions

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `periodic_x/y/z` | bool | false | Periodicity per direction. |
| `{x,y,z}_{min,max}_{u,v,w,p}_type` | enum | `dirichlet` | Face BC type: `dirichlet` / `0` or `neumann` / `1`. |
| `{x,y,z}_{min,max}_{u,v,w,p}_value` | real | 0.0 | The Dirichlet / Neumann value for that face and variable. |
| `{x,y,z}_{min,max}_patch` | enum | (inferred) | Face patch type for the RANS layer: `wall` or `patch`. Absent = infer a wall from Dirichlet tangential velocities. Non-periodic faces only (config error otherwise). |

Every non-periodic face defaults to a homogeneous Dirichlet condition. See the
[`sailplane`](tutorials.md#sailplane-external-aerodynamics-with-ibm) tutorial for a full
inflow/outflow/symmetry set. The patch type only affects the turbulence
models (wall distance, ω pinning, scalar wall ghosts) — declare a
Dirichlet velocity inlet `patch` so RANS does not read it as a no-slip
wall.

## `[ibm]` — immersed boundary method

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | true | Enable the volume-penalization IBM. |
| `coeff_file` | string | (empty) | HDF5 coefficient file produced by `mobygeom`. |

## `[blocks]` — block decomposition and 2:1 refinement

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `nb` | int | 0 (auto) | Block edge in cells; if set, must be even and ≥ 4, and divide the grid. |
| `remove_solid` | bool | true | Drop blocks fully buried inside a solid body. |
| `refine` | 6 reals | — | Refinement box `xmin xmax ymin ymax zmin zmax`; repeatable up to 4 boxes. |
| `refine_levels` | int | 1 | Number of refinement levels. |
| `refine_body` | bool | false | Refine blocks touching the immersed body (geometry-driven). |

See [Numerical methods](numerical-methods.md#block-structured-grid-and-21-refinement) for
what these do.

## `[force]` — spatially varying volumetric body force

Optional force `f(x)` added **on top of** the constant `[flow] forcing_*`.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | false | Enable the `f(x)` term. |
| `type` | string | `profile` | `profile` (analytic), `file` (read a field), or `custom` (code hook). |
| `profile` | string | `constant` | `constant` or `sine`. |
| `amp_x/y/z` | real | 0.0 | Force amplitude per direction. |
| `wavenumber_x/y/z` (aliases `k_x/k_y/k_z`) | real | 0.0 | Sine-profile wavenumbers. |
| `dir` | int | 1 | Force direction index. |
| `file` | string | (empty) | Force-field HDF5 file (velocity layout) for `type = file`. |

## `[turbulence]` — model family

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `model` | enum | `none` | Family: `none`, `les`, `rans`, or `iddes`. Absent + a configured `[les] model` implies `les` (an explicit `none` wins). |
| `fd_force` | real | -1 (off) | IDDES validation hook: force the DDES shielding function to a constant (0 = pure-SGS limit, 1 = pure-RANS limit). |

`les` needs an SGS model in `[les]`; `rans` needs `[rans] model = sst`;
`iddes` needs BOTH (SST transport near walls, the SGS model where the
DDES shielding releases the flow to LES). In RANS/IDDES runs the output
`p` is a modified pressure (the −2/3 k δij part of the Boussinesq stress
is absorbed into it).

## `[les]` — subgrid-scale model (used by `les` and `iddes`)

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `model` | enum | `none` | `none`, `smagorinsky` (`smag`), or `wale`. |
| `cs` | real | 0.10 | Smagorinsky constant. |
| `cw` | real | 0.325 | WALE constant. |
| `delta_scale` | real | 1.0 | Filter-width scaling (> 0). |
| `ibm_aware` | bool | true | Zero the subgrid viscosity inside solid (IBM) cells. |

## `[rans]` — k-ω SST (used by `rans` and `iddes`)

The section's presence alone builds the SST geometry state (wall
distance, IBM wall cells) so it can be inspected before any transport
runs; `[turbulence] model = rans|iddes` additionally advances k/ω.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `model` | enum | `none` | `sst` enables the transport equations. |
| `wall_treatment` | enum | `resolved` | `resolved` (y⁺₁ ≲ 1) or `wall_function` (Weber/OpenFOAM ω+νt wall functions). Rejected under `iddes` (not validated). |
| `transition` | bool | false | γ–Re_θt transition model (Langtry–Menter 2009); resolved walls only; rejected under `iddes`. |
| `tu` | real | 5.0 | Initial/freestream turbulence intensity in percent. |
| `nut_ratio` | real | 10.0 | Initial ν_t/ν, sets ω = k/(nut_ratio ν). |
| `dump_geometry` | bool | false | Write `<prefix>_ransgeom.h5` (dwall/yeff/wallcell + coordinates). |
| `dwall_tol` | real | 1e-10 | Analytic wall-distance polish tolerance. |

## `[output]` — field output

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `field_interval` | int | 0 | Steps between field dumps (≥ 0; 0 disables). |
| `field_prefix` | string | (empty) | Output filename prefix (`<prefix>_<n>.h5`). |

## `[restart]` — restart

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `file` | string | (empty) | Restart HDF5 field file. When set, the run continues from it and full validation is skipped. |

---

**Notes**

- There is no `[run]` section — run control lives in `[time]` and `[output]`.
- `[case.channel]` keys are consumed only when `[case] name = channel`; the `generic` case
  ignores them.
