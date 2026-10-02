# Tutorials

The cases in `tutorials/` show how to set up and run real flows; each directory has its own
README with the exact commands and recorded results. The index of every runnable case
(tutorials and validation gates) is [`validation/README.md`](../validation/README.md).

| Tutorial | What it shows |
|----------|---------------|
| `channel_kmm180/` | Re_τ 180 turbulent channel (Kim–Moin–Moser), natural-stretched y, `[case] name = channel` — walked through below. |
| `min_channel/` | Minimal-span channel with a 2:1 wall-band refinement; also a regression-suite case (`input_gpu.ini`). |
| `turbulentBoundaryLayer/` | Spatially developing ZPG turbulent boundary layer (`[case] name = boundarylayer`, Blasius inlet, trip) validated against SIMSON/CaNS/AMPHIBIOUS; `overheadTest/` holds the performance campaign. |
| `naca/rans/` | NACA 0012, α = 5°, Re_c = 4e5, k-ω SST on an xz-refined grid vs OpenFOAM (control-volume C_L/C_D, Cp, Cf). |
| `cht/` | Conjugate heat transfer: Flageul turbulent channel (`channel/`) and Neuhauser pipe (`pipe/`) vs published DNS. |
| `sailplane/` | External flow around a CAD STL with the immersed boundary method, preprocessed by `moby_prepare` — walked through below. |

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
mpirun -n 8 ./build_cpu/moby_solve tutorials/channel_kmm180/input.ini

# or single GPU
mpirun -n 1 ./build_gpu/moby_solve tutorials/channel_kmm180/input.ini
```

No developed restart is shipped; to continue a run, point `[restart] file` at one of its
own snapshots (e.g. `channel_kmm180_field_50000.h5`).

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
boundary conditions, and the STL → case-file preprocessing with `moby_prepare`.

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
; [case] file = sailplane_case.h5 names the moby_prepare case file, see below

[boundary]
periodic_x = false
periodic_y = false
periodic_z = false

; Inlet: prescribed uniform velocity (1,0,0), pressure Neumann
x_min_u_value = 1.0
x_min_p_type  = neumann
; Outlet: the patch type sets every row (normal velocity predicted by the
; momentum step, zero-gradient tangential velocities, pressure Dirichlet 0).
; A Neumann row on the NORMAL velocity is a configuration error.
x_max_patch   = outlet
; y = 0 symmetry: v = 0, Neumann for u, w, p  (similarly at the far y and z faces)
y_min_v_type  = dirichlet
y_min_v_value = 0.0
```

Uniform inflow enters at $x_{\min}$ with unit velocity; $x_{\max}$ is a pressure-reference
outflow; the $y$ and $z$ faces are symmetry/far-field. `[ibm]` points at the file that
encodes the solid geometry.

### Prepare the case file

The immersed body is classified against the solver's exact grid by `moby_prepare`. Make a
copy of `input.ini` (say `prep_blocks.ini`) that adds a block layout and the STL:

```ini
[blocks]
nb = 10

[ibm]
stl_file = "FRUE V0 ohneRundung.stl"
stl_scale = 0.001
stl_translate = 17.288649999032064 0.0 4.150549
```

```bash
cd tutorials/sailplane
mpirun -n 2 ../../build_cpu/moby_prepare prep_blocks.ini sailplane_case.h5
```

then solve with `[blocks] nb = 10` and `[case] file = sailplane_case.h5`.
`stl_scale = 0.001` converts the STL from millimetres to metres; `stl_translate` centres the
mirrored STL in the full domain (the solver then uses only the positive-`y` half). The
committed `sailplane_ibm_coeff.h5` is an older single-level coefficient file, still usable
without `[blocks] nb`. **Regenerate the case file whenever the grid, the block layout or
`re` changes.**

### Run

```bash
mpirun -n 1 ./build_gpu/moby_solve tutorials/sailplane/input.ini
```

The provided `input.ini` is sized as a single-step smoke case tuned to fit a 6 GB GPU; raise
`[time] nsteps` (and `field_interval`) for an actual run. Field snapshots
(`sailplane_field_*.h5`) are in the block-table layout; reassemble them with
`tools/compare_fields.py --export-global` for visualization. The committed coefficient file
has an `.xdmf` companion for viewing the classified geometry.

> See `tutorials/sailplane/README.md` for the exact bounding-box arithmetic.
