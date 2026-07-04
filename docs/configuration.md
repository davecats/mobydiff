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
| `name` | string | `generic` | Flow case: `channel` (honours `[case.channel]`) or `generic`. |

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

Every non-periodic face defaults to a homogeneous Dirichlet condition. See the
[`sailplane`](tutorials.md#sailplane-external-aerodynamics-with-ibm) tutorial for a full
inflow/outflow/symmetry set.

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

## `[les]` — large-eddy simulation

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `model` | enum | `none` | `none`, `smagorinsky` (`smag`), or `wale`. |
| `cs` | real | 0.10 | Smagorinsky constant. |
| `cw` | real | 0.325 | WALE constant. |
| `delta_scale` | real | 1.0 | Filter-width scaling (> 0). |
| `ibm_aware` | bool | true | Zero the subgrid viscosity inside solid (IBM) cells. |

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
