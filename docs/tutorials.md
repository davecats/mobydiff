# Tutorials

This guide walks through two representative cases shipped in `tutorials/`: a wall-bounded
turbulent flow (`channel_kmm180`) and an external-aerodynamics flow around an immersed body
(`sailplane`). Each tutorial directory contains a ready-to-run `input.ini`.

Run every case through `mpirun`, even on a single rank.

---

## `channel_kmm180` — turbulent plane channel

A canonical incompressible turbulent channel at friction Reynolds number
$\mathrm{Re}_\tau = 180$ (the classic Kim–Moin–Moser setup). It exercises the channel flow
case, near-wall grid stretching, constant streamwise forcing, and turbulence-statistics
accumulation — no immersed boundary.

### Setup

The domain is $L_x \times L_y \times L_z = 4\pi h \times 2h \times 2\pi h$ on a
$288 \times 136 \times 288$ grid. Streamwise ($x$) and spanwise ($z$) are uniform and
periodic; the wall-normal direction ($y$) uses the **natural** near-wall stretching with the
first off-wall spacing at $\Delta y_w^+ \approx 0.05$:

```ini
[case]
name = channel

[case.channel]
n_walls = 2
natural_blend_index = 16
large_disturbance_amplitude = 1.0e-2
small_noise_amplitude       = 1.0e-3
stats_sample_interval = 500
stats_write_interval  = 5000
stats_file = channel_kmm180_stats.h5

[grid]
nx = 288
ny = 136
nz = 288
lx = 12.566370614359172   ; 4*pi
ly = 2.0
lz = 6.283185307179586    ; 2*pi

[grid.y]
distribution = natural
stretch = 16
natural_dyw_plus = 0.05

[flow]
re = 180.0        ; Re = Re_tau
forcing_x = 1.0   ; unit mean pressure gradient balances wall friction

[boundary]
periodic_x = true
periodic_y = false
periodic_z = true
```

At $\mathrm{Re}_\tau = 180$ this resolution gives $\Delta x^+ \approx 7.9$ and
$\Delta z^+ \approx 3.9$. The flow is driven by a unit streamwise body force
(`forcing_x = 1.0`), which under the friction scaling equals the mean pressure gradient that
balances the wall shear. The channel initializer seeds a laminar mean profile with a
large-scale disturbance plus small noise to trip transition to turbulence.

### Run

```bash
# 8-rank CPU run (build_cpu is the reference build)
mpirun -n 8 ./build_cpu/main tutorials/channel_kmm180/input.ini

# or single GPU
mpirun -n 1 ./build_gpu/main tutorials/channel_kmm180/input.ini
```

A pre-computed restart (`restart.h5` / `channel_kmm180_restart.h5`) is provided so you can
start from a developed field instead of waiting through transition — point `[restart] file`
at it.

### Inspect

Turbulence statistics accumulate into `channel_kmm180_stats.h5` at the configured intervals.
Post-process with the channel tools (see the [tools reference](tools.md)):

```bash
python3 tools/plot_channel_stats.py stats.png channel_kmm180_stats.h5:kmm180
python3 tools/channel_loglaw.py loglaw.png channel_kmm180_field_50000.h5:kmm180
```

The mean profile should collapse onto the law of the wall ($U^+ = y^+$ in the viscous
sublayer, $U^+ \approx 2.44\ln y^+ + 5$ in the log layer), and the rms fluctuation profiles
should match the reference DNS.

---

## `sailplane` — external aerodynamics with IBM

Flow around a sailplane geometry supplied as an STL mesh, imposed with the volume-penalization
immersed boundary method. It exercises the `generic` flow case, inflow/outflow/symmetry
boundary conditions, and the full STL → IBM-coefficient preprocessing workflow.

The source STL is in millimetres and symmetric about its `y = 0` plane. The tutorial works in
metres, keeps the symmetry plane at computational `y = 0`, and simulates only the
`y ≥ 0` half-domain.

### Setup

```ini
[case]
name = generic

[grid]
nx = 400
ny = 450
nz = 100
lx = 43.22162499838677   ; 5 * length
ly = 48.72317            ; 2.5 * span (half-domain in y)
lz = 9.750915            ; 5 * height

[flow]
re = 1.0e5

[ibm]
enabled = true
coeff_file = sailplane_ibm_coeff.h5

[boundary]
periodic_x = false
periodic_y = false
periodic_z = false

; Inlet: prescribed uniform velocity (1,0,0), pressure Neumann
x_min_u_value = 1.0
x_min_p_type  = neumann
; Outlet: velocity Neumann, pressure reference (Dirichlet 0)
x_max_u_type  = neumann
x_max_v_type  = neumann
x_max_w_type  = neumann
x_max_p_type  = dirichlet
; y = 0 symmetry: v = 0, Neumann for u, w, p  (similarly at the far y and z faces)
y_min_v_type  = dirichlet
y_min_v_value = 0.0
```

Uniform inflow enters at $x_{\min}$ with unit velocity; $x_{\max}$ is a pressure-reference
outflow; the $y$ and $z$ faces are symmetry/far-field. `[ibm]` points at the coefficient file
that encodes the solid geometry.

### Generate the IBM coefficients

The immersed body must be classified against the **exact** solver grid, so first export the
grid with `mobygrid`, then run the `mobygeom` STL-to-coefficient step:

```bash
cd tutorials/sailplane

# 1. export the exact node lines the solver will use
../../build_cpu/mobygrid input.ini sailplane_grid.h5

# 2. classify the STL against that grid and write the coefficient file
python3 ../../tools/mobygeom.py stl-ibm-coeff \
    --geometry "FRUE V0 ohneRundung.stl" \
    --output    sailplane_ibm_coeff.h5 \
    --grid-file sailplane_grid.h5 \
    --re 1.0e5 \
    --scale 0.001 \
    --translate 17.288649999032064 0.0 4.150549 \
    --check-fluid-points fluid_points.txt \
    --jobs 2 --tile-size 32 36 16
```

`--scale 0.001` converts the STL from millimetres to metres; `--translate` centres the
mirrored STL in the full domain (the solver then uses only the positive-`y` half).
`fluid_points.txt` is a set of known-fluid probe points used as a fail-fast sanity check on
the inside/outside classification. **Regenerate the coefficient file whenever the grid in
`input.ini` changes.**

### Run

```bash
mpirun -n 1 ./build_gpu/main tutorials/sailplane/input.ini
```

The provided `input.ini` is sized as a single-step smoke case tuned to fit a 6 GB GPU; raise
`[time] nsteps` (and `field_interval`) for an actual run. Field snapshots
(`sailplane_field_*.h5`, with companion `.xdmf`) can be opened in ParaView, and the
coefficient file itself has an `.xdmf` companion so you can visualize the classified geometry.

> See `tutorials/sailplane/README.md` for the exact bounding-box arithmetic and the
> preprocessing options (bounding-box culling, ambiguous-winding repair, ray-vote debugging).
