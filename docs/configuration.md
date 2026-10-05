# Configuration reference

`moby_solve` reads a single INI-style configuration file, passed as the sole command-line
argument (`moby_prepare` reads the same file, plus the geometry keys marked below). This
page documents every section and key.

## File syntax

- Sections are written `[name]`; keys are `key = value`.
- Section and key names are **case-insensitive** (lowercased on read).
- Comments start with `;` or `#` and run to end of line.
- Blank lines are ignored. **Unknown sections, and unknown keys in most sections, are
  silently ignored** — check spelling carefully (`[boundary]`, `[rans]`, `[scalar]` and the
  `[case.*]` sections print a warning for an unknown key).
- String values may be single- or double-quoted; the quotes are stripped.
- Booleans accept `true` / `.true.` / `1` / `yes` and `false` / `.false.` / `0` / `no`.

## Required vs. optional

For a fresh run (no restart), these keys are **required** and must be positive:

- `[grid] nx, ny, nz` and `lx, ly, lz` (or `[grid.x/y/z] n` / `length`)
- `[flow] re`
- `[time] dt`, plus at least one of `nsteps` / `t_final`, and `dtmax > 0`

Everything else is optional and falls back to the defaults listed below. If
`[restart] file` is set, full validation is skipped and the run continues from that field.

---

## `[case]` — case selector

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | `generic` | Flow case: `generic`, `channel` (`[case.channel]`), `boundarylayer` (`[case.boundarylayer]`) or `airfoil` (`[case.airfoil]`). An unknown name warns and falls back to `generic`. |
| `file` | string | `<field_prefix>.case.h5` | The **case file**: grid, block size, leaf table, face kinds, IBM coefficients and wall distance. `mpirun -n N moby_prepare input.ini` writes it, `moby_solve` only reads it (a missing file stops with the prepare command printed, a stale one with the offending `key = value` named). The pre-2026-09-30 alias `[ibm] coeff_file` is gone. |

A case sets its own defaults (boundary faces, grid distribution, forcing, …) before the
rest of the file is applied, so explicit keys elsewhere in the `.ini` still win.

## `[case.channel]` — channel-case parameters

Only read when `[case] name = channel`.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `n_walls` | int (1 or 2) | 2 | One- or two-wall channel. |
| `natural_blend_index` (alias `jb`) | real | 40.0 | Pirozzoli–Orlandi wall-normal stretch blend index. |
| `mean_profile_sine_amplitude` (alias `mean_sine_amplitude`) | real | 0.0 | Amplitude of a `sin(2πy/Ly)` term added to the mean streamwise profile. |
| `large_disturbance_amplitude` | real | 1.0e-2 | Large-scale initial disturbance amplitude. |
| `small_noise_amplitude` | real | 1.0e-3 | Small random-noise amplitude. |
| `stats_sample_interval` | int | -1 (off) | Steps between statistics samples. |
| `stats_write_interval` | int | -1 (off) | Steps between statistics writes. |
| `stats_file` | string | `channel_stats.h5` | HDF5 statistics output. |
| `runtime_file` | string | `runtimedata.txt` | Runtime log file. |

## `[case.boundarylayer]` — spatially developing boundary layer

Only read when `[case] name = boundarylayer`. The case declares a Blasius inlet at
`x_min`, outlets at `x_max` and `y_max`, a wall at `y_min`, periodic `z`, the `blayer`
wall-normal grid, and enables the `[force] type = trip` forcing with the defaults below
(an explicit `[force]` section still overrides them).

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `u_inf` | real | 1.0 | Free-stream velocity. |
| `theta_in` (alias `blasius_theta`) | real | 1/H | Inlet momentum thickness of the Blasius profile. |
| `natural_blend_index` (alias `jb`) | real | 45.0 | Wall-normal stretch parameter. |
| `dyw_plus` | real | 0.15 | First off-wall spacing in wall units. |
| `resolved_height` (alias `outer_height`) | real | 30.0 | Height of the wall-resolved region; the grid coarsens geometrically above it. |
| `trip_enabled` | bool | true | Enable the trip forcing. |
| `trip_x0`, `trip_lx`, `trip_ly`, `trip_amp`, `trip_amp_s`, `trip_ts`, `trip_nmodes`, `trip_seed` | | 10, 4, 1, 0.18854, 0, 4, 16, 1 | Trip parameters (see `[force]`). |
| `stats_sample_interval` / `stats_write_interval` / `runtime_interval` | int | -1 (off) | Statistics sampling / writing, runtime-log interval. |
| `stats_file` | string | `bl_stats.h5` | Statistics output. |
| `runtime_file` | string | `runtimedata.txt` | Runtime log file. |

## `[case.airfoil]` — airfoil / bluff body in a free stream

Only read when `[case] name = airfoil`. The span axis is periodic; the case declares an
outlet at `x_max` and inlets (free-stream velocity at angle α in the chord–lift plane) on
every other face, and samples lift and drag from a control-volume momentum budget.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `aoa` | real | 0.0 | Angle of attack in degrees. |
| `u_inf` | real | 1.0 | Free-stream speed. |
| `chord` | real | 1.0 | Reference length for C_L / C_D. |
| `span` | enum | `z` | Periodic span axis: `z` (lift along y) or `y` (lift along z, for `refine_dims = xz`). |
| `cv_box` | 4 reals | — | `c0 c1 l0 l1`: control-volume extent along the chord (x) and lift axes (full span). **Without it force sampling is disabled** (with a warning). Borders snap to cell faces. |
| `force_sample_interval` | int | 10 | Steps between force samples (0 = off). |
| `runtime_file` | string | `forces.txt` | Force history (`iteration time cl cd dmomdt`). |
| `steady_tol` | real | 0.0 (off) | Stop once the unsteady term of the budget stays below this (coefficient units)… |
| `steady_samples` | int | 3 | …for this many consecutive samples. |

## `[grid]` — global grid size and domain

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `nx`, `ny`, `nz` | int | — | Global cell counts. **Required.** |
| `lx`, `ly`, `lz` | real | — | Domain lengths. **Required.** |

## `[grid.x]`, `[grid.y]`, `[grid.z]` — per-axis grid controls

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `distribution` | enum | `uniform` | `uniform`, `cosine`, `tanh`, `natural` (aliases `pirozzoli_orlandi` / `po`), `blayer` (natural clustering up to `outer_height`, geometric coarsening above), `geometric`, or `cans` (aliases `tanh_wall`, `tanh_lower`: the one-sided tanh clustered at the LOW wall, `y = L [1 + tanh(a (s − 1))/tanh(a)]`, `a = stretch` — CaNS `gtype = 2`, `gr`). |
| `stretch` | real | 0.0 | Stretching parameter for `cosine` / `tanh` / `natural`. |
| `natural_dyw_plus` (aliases `dyw_plus`, `dy_wall_plus`, `dyw+`) | real | 0.05 | First off-wall spacing in wall units (natural grid). |
| `one_sided` (alias `natural_one_sided`) | bool | false | Natural grid: cluster at the low end only (boundary layer) instead of both ends. |
| `outer_height` (aliases `natural_outer_height`, `resolved_height`) | real | — | `blayer`: height of the wall-resolved region. |
| `nodes_file` (aliases `node_file`, `nodes`) | string | — | Read this axis's node line from a file instead of generating it. |
| `n` | int | — | Per-axis alternative to `[grid] nx/ny/nz`. |
| `length` | real | — | Per-axis alternative to `[grid] lx/ly/lz`. |
| `subdivided` | bool | false | Build the line as the midpoint subdivision of the half-resolution line (bitwise-identical to a refined level of that grid). |

## `[mpi]` — domain decomposition

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `dims` | 3 ints | `0 0 0` | MPI process grid; `0` lets MPI choose the factor for that axis. |

## `[flow]` — physics and initial condition

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `re` | real | — | Reynolds number. **Required.** |
| `convection` | enum | `skew` | Momentum convection form: `skew` (skew-symmetric, the production form, energy-neutral for any advecting field) or `divergence` (alias `conservative`, `cons`: the CaNS/AMPHIBIOUS form, re-admitted 2026-10-05 for the `turbulentBoundaryLayer` code comparison; NOT energy-neutral under the incremental projection, so watch the divergence residual on a long run). Independent of `[scalar] convection`. |
| `forcing_x/y/z` | real | 0.0 | Constant body force / mean pressure gradient per direction. |
| `initial_u/v/w` | real | 0.0 | Uniform initial velocity components. |
| `initial_noise` | real | 0.0 | Initial random-noise amplitude. |
| `initial` | string | `uniform` | Initial-field mode (generic case): `uniform`, or the analytic test fields `beltrami` (3D ABC flow), `tgv` (2D Taylor–Green), `tgv3d` (manufactured 3D field, not an NS solution). |

## `[time]` — time stepping

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `dt` | real | — | Nominal time step. **Required.** |
| `nsteps` | int | 0 | Number of steps (one of `nsteps` / `t_final` must be > 0). |
| `t_final` | real | 0.0 | Final time. |
| `cflmax` | real | 0.0 | Maximum Courant number for the adaptive step (≥ 0): `dt` is bounded by `cflmax / max_cells Σ_d |u_d|/Δx_d`, the SUM over the three directions of a cell. RK3 with central convection is stable for `dt × Σ_d |u_d|/Δx_d < √3 = 1.73` whatever the direction of the flow (`validation/courant_limit/`); a value above it is announced at init. Until 2026-10-02 this key bounded the largest single component `|u_d|/Δx_d`, which in a turbulent flow is 1.2–1.5 times smaller than the sum (up to 3 for a flow along the diagonal): the usual 0.8 was effectively 1.0–1.2 on the sum and is now 46 % of the limit, so a Courant-bound run takes that many more steps at the same key; 1.2 restores the old step with a 30 % margin. |
| `pecletmax` | real | 0.0 | Maximum diffusion number for the adaptive step (≥ 0): `dt` is bounded by `pecletmax / max_cells Σ_d ν_eff/Δx_d²`, the SUM over the three directions of a cell (ν_eff includes the eddy viscosity and the most diffusive scalar). The RK3 limit of explicit diffusion is `dt × Σ_d ν/Δx_d² ≤ 0.628` on any grid (`validation/diffusion_limit/`), so the usual 0.5 is 80 % of it; a value above 0.628 is announced at init, and a step beyond the limit (a fixed `dt` too) is warned about once. Until 2026-10-02 this key bounded the single-direction number `dt ν/Δx²`: on a grid with one dominant direction (a fine wall-normal line) the two differ by a few per cent, on an isotropic grid by a factor 3. |
| `dtmax` | real | — | Hard cap on `dt`. **Must be > 0.** |

## `[pressure]` — pressure projection / Poisson solver

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `solver` | enum | `jacobi` | Smoother: `jacobi` (damped Jacobi) or `redblack` (aliases `red-black`, `sor`; coloured Gauss–Seidel / SOR). |
| `niter` | int | 3 | Iterations per projection (≥ 0). One red-black iteration is two colour sweeps. |
| `sor` | real | 0.8 | Damping / over-relaxation factor for plain damped Jacobi (diverges above ~0.8) and red-black SOR (over-relaxation, 1.5 is the usual setting; **0.8 is the wrong value for it**). **Ignored under `accel = chebyshev`**: the Chebyshev semi-iteration takes its damping from `cheb_lmin`/`cheb_lmax`, which depend only on the grid and are auto-derived. |
| `accel` | enum | `jacobi` (none) | `chebyshev` / `cheb` enables Chebyshev–Jacobi acceleration. Mutually exclusive with `solver = redblack` (Chebyshev needs a stationary linear operator). |
| `cheb_lmin` | real | -1.0 (auto) | Chebyshev lower eigenvalue bound. Auto = (2/3) sin²(π/N_max), N_max the largest base global grid size (the domain-scale Jacobi–Poisson mode; refinement does not change it). Printed at init. Use `niter = 12` with Chebyshev. |
| `cheb_lmax` | real | -1.0 (auto = 2.0) | Chebyshev upper eigenvalue bound: 2 by Gershgorin for any Jacobi-preconditioned Poisson operator, stretched grids and 2:1 interfaces included. |

Both smoothers run on 2:1-refined grids (red-black since R1 —
`validation/redblack_interface/`). They need **different `niter` for the same
quality**: measured on the production boundary layer, damped Jacobi needs
`niter ≥ 12` to keep the pressure zero-mode away where red-black stays clean
below 6, so comparing the two at a fixed `niter` compares nothing useful.
`solver = redblack` additionally requires even global sizes in periodic
directions, for a consistent colouring.

## `[boundary]` — boundary conditions

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `periodic_x/y/z` | bool | false | Periodicity per direction. |
| `{x,y,z}_{min,max}_{u,v,w,p}_type` | enum | `dirichlet` | Face BC type: `dirichlet` / `0` or `neumann` / `1`. |
| `{x,y,z}_{min,max}_{u,v,w,p}_value` | real | 0.0 | The Dirichlet / Neumann value for that face and variable. |
| `{x,y,z}_{min,max}_{u,v,w,p}_profile` | enum | `constant` | Per-point boundary values: `constant`, `parabola` (Poiseuille inlet) or `blasius` (Blasius similarity profile; the `_value` is U∞). |
| `blasius_theta` | real | — | Inlet momentum thickness for the `blasius` profile. |
| `{x,y,z}_{min,max}_patch` | enum | (inferred) | Face patch type: `wall`, `patch` (alias `generic`), `inlet` or `outlet`. Absent = infer a wall from Dirichlet tangential velocities. Non-periodic faces only (config error otherwise). |

Every non-periodic face defaults to a homogeneous Dirichlet condition. The patch type is
the one user-facing face axis: it derives the per-variable rows that are not set
explicitly (`outlet` = the face-normal velocity is an unknown the predictor advances, with zero-gradient tangential velocities and a Dirichlet pressure; `inlet` =
Dirichlet velocity, with RANS / scalar free-stream values), and tells the turbulence models
where the walls are (wall distance, ω pinning, scalar wall ghosts). An explicit `_type` key
that contradicts the declared patch is a hard error; explicit `_type` / `_value` keys also
override those stored in a restart file.

An immersed body may cross an inlet or an outlet plane (a rough or immersed wall
in a boundary layer, a plate leaving through the outlet). On an outlet the
face carries its own penalization (`step.f90 predict_outlet_faces`); on a
Dirichlet face the velocity datum is zero wherever the face's own staggered
location is inside the solid (the inlet does not blow through the body; the
solver reports how many non-zero data it zeroed). The ghost rows half a cell
OUTSIDE the domain are sampled from the geometry too, so an STL body that is
meant to cross a plane must extend beyond it (`validation/body_outlet/`).

## `[ibm]` — immersed boundary method

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | true | Enable the volume-penalization IBM. |
| `wall_shape` | enum | `wavy` | Analytic wall: `wavy` or `eggcarton` (alias `rough`; 3D sinusoidal roughness). |
| `amp_x`, `amp_z` | real | 0.025 | Analytic wall amplitudes. |
| `n_wave_x`, `n_wave_z` | int | 1 | Wavelengths per domain length. |
| `phase_x`, `phase_z` | real | 0.0 | Phase offsets. |
| `wall_offset` | real | 0.01 | Mean plane of the analytic wall (`wavy`: the trough height; `eggcarton`: the mid-plane, so `wall_offset = amp_x` puts the troughs on the domain wall and the surface is MacDonald et al.'s `k (1 + cos cos)`). |
| `wall_x_start`, `wall_ramp`, `wall_x_end` | real | — | `eggcarton` only: a smooth wall upstream of `x_start`, the full egg-carton beyond `x_start + ramp`, a smoothstep of the whole surface (mean plane included) in between; `wall_x_end` (optional) ramps it out again over `[x_end, x_end + ramp]`. Upstream the surface is held one crest-to-trough height below the domain wall so that the graded coefficient band does not reach the first fluid row; it crosses the wall about 2/3 of the way through the ramp. Unset = the wall everywhere. Both are recorded in the case file's input echo only when set, so older case files stay valid. |
| `stl_file` | string | — | **`moby_prepare` only.** STL body (binary or ASCII); repeatable, up to 8. The solver accepts it only beside a case file that carries the same geometry. |
| `stl_scale`, `stl_translate` | real, 3 reals | 1.0, `0 0 0` | **`moby_prepare` only.** Transform `v·scale + translate`. |
| `band_filter` | bool | false | 3-point low-pass on the predicted velocity in a thin band around the body (damps the staircase cell-Reynolds fan). |
| `band_width` | int | 3 | Band width in cells. |
| `band_theta` | real | 0.5 | Filter strength, in [0, 0.6]. |

## `[blocks]` — block decomposition and 2:1 refinement

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `nb` | 1 or 3 ints | 0 (auto) | Block edge in cells. One value broadcasts (`nb = 16`), three set each direction (`nb = 64 44 48`). Each must be even, ≥ 4, and divide the grid; set all three or none. |
| `remove_solid` | bool | true | Drop blocks fully buried inside a solid body. |
| `refine` | 6 or 7 reals | — | Refinement box `xmin xmax ymin ymax zmin zmax [level]`; the optional 7th value caps the box's level (default `refine_levels`). Repeatable, up to 16 boxes. |
| `refine_levels` | int | 1 | Number of refinement levels. |
| `refine_body_levels` | int | -1 (= `refine_levels`) | Cap on body-driven refinement. |
| `refine_body_box` | 7 reals | — | `xmin xmax ymin ymax zmin zmax level`: raises the body cap to `level` inside the box (only the surface band refines, not the volume). Repeatable, up to 16. |
| `refine_dims` | enum | `xyz` | `xyz` (octree) or `xz` (quadtree: blocks refine in x and z only, y keeps one global node line at every level). |
| `refine_body` | bool | false | Refine blocks touching the immersed body (geometry-driven). |
| `keep_buried` | bool | false | Keep blocks buried inside the body. Required by conjugate heat transfer with `refine_body` (the solid is an unknown there). |

**Choosing `nb`.** Each block carries a one-cell halo, so the cost of a layout is
its *allocated* volume, `(1+2/nb_x)(1+2/nb_y)(1+2/nb_z)` — measured to predict
the step time to better than 1 % across a 2× range of block size. Almost all of
that tax is device-local halo-copy traffic, not arithmetic. Prefer the largest
`nb` the refinement placement tolerates, per direction: under
`refine_dims = xz`, y is never subdivided, so `nb_y` buys nothing and is pure
overhead. On the 4096×176×192 boundary-layer grid, `nb = 16` costs 1.424 while
`nb = 64 44 48` costs 1.123.

The trade-off is granularity: a 2:1 interface can only sit on a block boundary,
so a larger `nb_y` leaves fewer legal interface heights (at `nb_y = 44`, y-index
44/88/132 instead of every 16 cells).

See [Numerical methods](numerical-methods.md#block-structured-grid-and-21-refinement) for
what these do.

## `[force]` — spatially varying volumetric body force

Optional force `f(x)` added **on top of** the constant `[flow] forcing_*`.

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `enabled` | bool | false | Enable the `f(x)` term. |
| `type` | string | `profile` | `profile` (analytic), `file` (read a field), `custom` (code hook `update_bodyforce`), or `trip` (random wall-normal trip, CaNS/SIMSON form). |
| `profile` | string | `constant` | `constant` or `sine`. |
| `amp_x/y/z` | real | 0.0 | Force amplitude per direction. |
| `wavenumber_x/y/z` (aliases `k_x/k_y/k_z`) | real | 0.0 | Sine-profile wavenumbers. |
| `dir` | int | 1 | Force direction index. |
| `file` | string | (empty) | Force-field HDF5 file (velocity layout) for `type = file`. |
| `trip_x0`, `trip_lx`, `trip_ly` | real | 0, 4, 1 | `type = trip`: streamwise centre and Gaussian widths of the envelope. |
| `trip_amp`, `trip_amp_s` | real | 0, 0 | Unsteady and steady trip amplitudes. |
| `trip_ts` | real | 4.0 | Time between random-walk redraws. |
| `trip_nmodes`, `trip_seed` | int | 16, 1 | Spanwise modes and RNG seed. |

## `[turbulence]` — model family

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `model` | enum | `none` | Family: `none`, `les`, `rans`, or `iddes`. Absent + a configured `[les] model` implies `les` (an explicit `none` wins). |
| `fd_force` | real | -1 (off) | IDDES validation hook: force the DDES shielding function to a constant (0 = pure-SGS limit, 1 = pure-RANS limit). |
| `iddes_cdt1` | real | 20.0 | IDDES shielding constant C_dt1 (8 = Spalart's DDES value). |
| `iddes_clip` | bool | false | Spalart's `max(0, l_RANS − l_LES)` clipping in the hybrid length. |
| `iddes_delta` | enum | `iddes` | LES length scale: `iddes` (wall-aware formula) or `cbrt` (`(ΔxΔyΔz)^⅓`). |

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
| `ambient_sustain` | bool | false | Rumsey free-stream sustaining sources (keeps the ambient k/ω from decaying). |
| `kpin_box` | 6 reals | — | `x0 x1 y0 y1 z0 z1`: k pinned to 0 inside (forced-laminar zone). Repeatable. |
| `ktrip_box` | 7 reals | — | `x0 x1 y0 y1 z0 z1 rate`: volumetric k source in fluid cells inside (trip strip). Repeatable. |
| `dump_geometry` | bool | false | Write `<prefix>_ransgeom.h5` (dwall/yeff/wallcell + coordinates). |
| `dwall_tol` | real | 1e-10 | Analytic wall-distance polish tolerance. |

## `[scalar]`, `[scalar.N]` — passive scalars

`[scalar]` holds the global keys; each scalar has a numbered section `[scalar.1]`,
`[scalar.2]`, … Scalars are written to the snapshots under their `name`.

| `[scalar]` key | Type | Default | Meaning |
|-----|------|---------|---------|
| `count` | int | number of `[scalar.N]` sections | Number of scalars (must match the sections if given). |
| `convection` | enum | `divergence` (`advective` when any scalar is `conjugate`) | `divergence` (conservative), `skew`, or `advective` (preserves a uniform scalar for any velocity field). With an `ibm_wall = conjugate` scalar an unset key resolves to `advective` and an explicit `divergence` is a hard error: the masked cut-face flux makes a uniform scalar drift in fluid cut cells under the divergence form (numerics review F2, 2026-09-27). |
| `stats_sample_interval`, `stats_write_interval` | int | off | In-solver scalar statistics. |
| `stats_file` | string | — | Statistics output (read with `tools/scalar_stats.py`). |
| `stats_layout` | enum | `profile` | `profile` (wall-normal rows, per level) or `plane` (x-y plane, z-averaged). |
| `heat_interval`, `heat_file` | int, string | off | Body heat-release diagnostic (Nusselt number for conjugate scalars). |
| `indicator_interval` | int | off | Conjugate cut-face tangential-error indicator. |

| `[scalar.N]` key | Type | Default | Meaning |
|-----|------|---------|---------|
| `name` | string | `sN` | Dataset name in the snapshots. |
| `pr` (aliases `prandtl`, `sc`, `schmidt`) | real | 1.0 | Molecular Prandtl / Schmidt number. |
| `prt`, `prt_model` | real, enum | 0.85, `constant` | Turbulent Prandtl number (`constant` or `kays`). |
| `initial` (alias `init`), `init_profile` | real, enum | 0.0, `uniform` | Initial value / profile (`uniform` or `linear_y`). |
| `source`, `source_type`, `source_dir` | real, enum, x/y/z | 0, `uniform` | Volumetric source; `velocity` (alias `kasagi`) scales it with the velocity along `source_dir`. |
| `inlet` | real | 0.0 | Value at `inlet` patch faces. |
| `{x,y,z}_{min,max}_type` / `_value` | enum / real | Neumann 0 | Domain-face BC (`dirichlet` / `neumann`); defaults follow the patch types. |
| `ibm_wall` | enum | `dirichlet` | Immersed-wall treatment: `dirichlet`, `adiabatic` (alias `neumann`) or `conjugate` (solid conduction). |
| `ibm_value` | real | 0.0 | Body value for `dirichlet`. |
| `solid_k`, `solid_rhocp`, `solid_init`, `solid_source`, `contact_resistance` | real | 1, 1, 0, 0, 0 | `conjugate` only: solid conductivity, capacity, initial value, source, interface contact resistance. |
| `solid_thickness`, `solid_outer_k`, `solid_outer_rhocp` | real | 0, 0, 1 | `conjugate` only: an optional second (outer) solid material band. |
| `tangential_correction` | bool | false | `conjugate` only: experimental tangential-flux term (measured to be worse; leave off). |

Conjugate scalars are hard errors with block removal (set `[blocks] remove_solid = false`),
with `refine_body` unless `keep_buried = true`, with `ibm_value`, and with RANS wall
functions.

## `[output]` — field output

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `field_interval` | int | 0 | Steps between field dumps (≥ 0; 0 disables). |
| `field_prefix` | string | (empty) | Output filename prefix (`<prefix>_<step>.h5`). |
| `profile` | bool | false | Per-phase step timing (`profiling.f90`): three nested profilers (`step_timing` / `proj_timing` / `exch_timing`) printed after the loop, plus a coverage line against the loop timer. Diagnostic only — it reads clocks, so fields are bit-identical either way. |
| `exchange_barrier` | bool | false | Diagnostic. Puts an `MPI_Barrier` immediately before every `MPI_Waitall`, so the exchange's arrival **skew** lands in its own `skew_barrier` bucket and whatever remains in `mpi_wait` is **transfer**. Aggregate `mpi_wait` contains both and cannot separate them. Needs `profile = true`. Bit-exact (verified: fields identical with it on and off), but it serialises the exchange — never quote a step time from a barrier run. |

## `[restart]` — restart

| Key | Type | Default | Meaning |
|-----|------|---------|---------|
| `file` | string | (empty) | Restart HDF5 field file. When set, the run continues from it and full validation is skipped. |

---

**Notes**

- There is no `[run]` section — run control lives in `[time]` and `[output]`.
- `[case.*]` keys are consumed only by the matching `[case] name`; other cases ignore them.
- Removed keys: `[blocks] interface_constant_half`,
  `momentum_reflux` and `interface_skew` no longer exist (and, being unknown, are silently
  ignored).
