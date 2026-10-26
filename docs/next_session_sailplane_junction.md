# Next session: the sailplane wing-fuselage junction (transition study)

STATUS: plan written 2026-10-25 on head `46b1b38`; step 1 re-based on
AeroSandbox 2026-10-26 (user's suggestion) and run once with NeuralFoil
section data (section 3.4: preliminary alpha_box); the XFOIL polars of
section 3.1 are still owed. The geometry measurements of section 1 and the
tools of section 1.3 are done. This is the study
decided 2026-10-24 (CLAUDE.md "Open, in order", item 3): the laminar-turbulent
transition at the junction of the sailplane wing and fuselage, DNS in the
junction region, LES / IDDES elsewhere, Re_c = 3e5 on the root chord, NO trip.
The strategy behind it was laid out in three session reports on 2026-10-24
(worklog); this file is the executable version. Section 2 is the ledger.

Read first: CLAUDE.md "State of the solver" (gate rules, hosts, landmines),
`tutorials/sailplane/README.md` (the present case: file IBM from the CAD STL,
one-step smoke case, 18 M cells), `tutorials/naca/rans/README.md` (the
airfoil protocol this study inherits: forces, `surface_cp_cf.py`,
`interp_restart.py`, the HoreKa chunk chain `submit_naca_stage.sh`),
`validation/body_outlet/` (a body through the outlet plane, needed in
section 4).

## 0. The question, and the decided configuration

**Question.** Where and how does the boundary layer become turbulent at the
wing-fuselage junction of this sailplane at Re_c = 3e5 with no trip: the
separation bubble of the root section, the corner flow at the fuselage, or
their interaction -- and what the transition line looks like across the span
near the root.

**Configuration (decided 2026-10-24, three reports).**

- Symmetry plane at y = 0 (zero sideslip: physical). Only y >= 0 is computed.
- The wing is TRUNCATED: it exits a free-slip face at y_cut = fuselage side +
  3 c (section 1: 2.95 m). The box therefore represents an infinite wing with
  a periodic row of fuselages every 2 y_cut; what that misses is handed in by
  the lifting line of section 3 (one number, the induced angle).
- The aft fuselage and fin leave the box: the outlet plane cuts the fuselage
  2.5 c behind the trailing edge (the body-on-outlet mechanism is in place and
  gated). The STL has NO horizontal tail (section 1), so nothing is lost.
- Zones by resolution, one model (`[turbulence] model = iddes`), plus a zonal
  model-off switch for the DNS core (section 5.2): DNS core = upper junction,
  wing root to 0.5 c span, the fuselage strip above the wing, wall cell
  Delta+ 3 (Delta+ 4 as the cheaper first run); wall-resolved LES ring =
  the rest of the truncated wing on both sides and the forebody, Delta+ 8-10;
  far field coarse. No RANS wall zone remains after the aft cut.
- No trip. Transition is selected by the bubble and the corner; the
  disturbance floor is the fixed-count projection residual and must be shown
  not to move the transition line (section 5.1).
- Cost (2026-10-24 estimate, Delta+ 4 core): 1.3-1.6e9 cells, ~30 A100,
  ~1 s/step, 3e5 steps, ~3 500 GPU-h per run, ~10 000 GPU-h for the study
  (Delta+ 3 core: 2.2-2.7e9 cells, 45-55 A100, 15 000-20 000 GPU-h). This is
  a project allocation on HoreKa, not the exasim share (D5).

## 1. The geometry, measured (2026-10-25)

`tutorials/sailplane/FRUE V0 ohneRundung.stl`, ASCII, millimetres, mirrored
about y = 0, x streamwise (nose at 0), z up. "ohne Rundung" = without fillet:
the junction is a sharp corner (D1). Measured with
`tutorials/sailplane/junction/stl_sections.py planform` (section 1.3):

| item | value |
|---|---|
| span b, half-span | 19.49 m, 9.745 m |
| wing area S (root chord carried to y = 0), aspect ratio | 12.6 m2, **30** |
| root chord c (at the fuselage side, y = 0.5 m) -- THE reference length | **0.831 m** |
| chord at y = 1 / 2 / 3 / 5 / 7 / 9 m, tip | 0.768 / 0.759 / 0.749 / 0.690 / 0.598 / 0.413 / 0.308 m |
| leading edge x at root / tip (forward sweep of the LE) | 2.630 / 2.321 m (about 1.8 deg) |
| incidence of the chord line to the fuselage axis | **+2.3 deg, constant** (no washout in this V0 model) |
| dihedral (LE z 0.445 m at y = 0.5, 0.908 m at y = 9) | 3.1 deg |
| section t/c | 0.129 (0.144-0.150 at y = 0.5-0.75: the root thickening) |
| section trailing-edge flap | hinge at about 0.87 c, deflected about **6 deg down** (both surfaces kink there) |
| fuselage length; max half-width; height at x = 2 m | 8.64 m; 0.373 m (x = 2.0); 0.97 m (z -0.25 .. 0.73) |
| fuselage at the wing (x 2.6-3.5 m) | z -0.20 .. 0.73; the wing LE at z 0.445: a **shoulder / upper-mid wing** |
| boom half-width at x = 5 / 7 m | 0.149 / 0.103 m |
| fin | x 7.7-8.6 m, z up to 1.70; **no horizontal tail in the STL** |

Derived, for the truncation argument (elliptic loading as the reference, AVL
gives the real one in section 3):

| error term of the truncated box | size at the junction |
|---|---|
| missing induced angle alpha_i = C_L / (pi AR) | 0.61 deg at C_L 1 (Delta C_l 0.07), uniform to first order |
| spanwise loading drop, flattened by the mirror | 0.4 % at y = 0.5 c, 4.7 % at the cut (y = 2.95 m) |
| image fuselages at 2 y_cut = 5.9 m = 7.1 c | (0.37 / 5.9)^2 = 0.4 % in velocity |
| top / bottom free-slip faces at 10 c | about 0.04 deg equivalent |

The fuselage's own effect on the root loading (cross-flow upwash, carry-over)
is inside the box and computed, not modelled.

### 1.1 Units and Reynolds number

Work in CHORD UNITS: `stl_scale = 0.001 / 0.831` so that the root chord is 1,
`re = 3.0e5`, U = 1; the Delta+ estimates and every length of this file are
in c. (The tutorial's `re = 1e5` is per metre with the STL in metres, i.e.
Re_c = 0.83e5 -- a different flow. The case file carries 1/Re: re-prepare.)

### 1.2 The angle of attack

There is no `[ibm] stl_rotate` key and none is wanted: the free-slip top and
bottom faces hold w = 0, so the inflow cannot carry an angle. The geometry
is rotated ONCE, offline, nose-up about the root quarter chord
(`stl_sections.py rotate FILE ALPHA OUT.stl`, pivot (x, z) = (2.84, 0.47) m
by default), then prepared. alpha in this file is the FUSELAGE-AXIS angle to
the free stream; the wing sits at alpha + 2.3 deg.

### 1.3 The tool

`tutorials/sailplane/junction/stl_sections.py` (numpy only; 2026-10-25):
`planform` prints the table above, `section FILE Y OUT.dat` chains the slice
segments at y = Y into one closed chord-normalised loop in the XFOIL
convention (TE, upper surface, LE, lower surface, TE), `rotate` writes the
pitched STL. The loop drops the tessellated trailing-edge stub (a 2.5e-4 c
three-point blunt TE) and coincident leading-edge points from shared STL
edges: a zero-length panel is a floating-point trap in the Debian XFOIL
build and a non-monotonic x is refused by AeroSandbox's `repanel`.
`tutorials/sailplane/junction/liftingline.py` is the step-1 model itself
(section 3).

## 2. Ledger

| step | what | status |
|---|---|---|
| 1 | Route A lifting line in AeroSandbox: NeuralFoil / XFOIL section polars + the nonlinear lifting line of the whole wing -> alpha_i(y), C_l(y), alpha_box (section 3) | PRELIMINARY 2026-10-26 (NeuralFoil; section 3.4); XFOIL polars and the Ncrit band owed |
| 2 | The calibration box: coarse truncated model at Re_c 3e5, alpha_box and +-0.5 deg, strip loads, cut at 3 c vs 4 c (section 4) | PENDING |
| 3 | Go / no-go physics tests: projection-residual floor on the SD7003 bubble, TS growth through the IBM (section 5.1) | PENDING |
| 4 | Solver increments: zonal model-off levels; strip-load post-processing (section 5.2) | PENDING |
| 5 | Production geometry and grid: fine STL (fillet per D1), the DNS core levels, prepare at 3-5 M leaves (section 5.3) | PENDING |
| 6 | Cross-node exchange performance on an 8-node slice of the real grid (section 5.4) | PENDING |
| 7 | Production run at Delta+ 4, then the Delta+ 3 core as the resolution check (section 5.5) | PENDING -- needs the allocation (D5) |

Decisions for the user, marked D1-D5, are collected in section 6.

## 3. Step 1: the lifting line in AeroSandbox (Route A)

Purpose: the loading of the WHOLE wing at the flight condition, from which
the box gets its one number: alpha_box = alpha_flight - alpha_i(root region).
The fuselage is NOT in the lifting-line model (the box supplies it): the
wing runs through the centreline with the root chord carried over, which
is the standard treatment and what the box's symmetry plane also implies.

### 3.0 The tool: AeroSandbox (decided 2026-10-26, replaces XFOIL + AVL by hand)

`pip install aerosandbox[full]` in a venv (4.2.10, NeuralFoil 0.3.3, Python
3.12, installed and tested on the workstation 2026-10-26; the venv is not in
git -- `python3 -m venv ~/asb && . ~/asb/bin/activate && pip install
aerosandbox[full]`). What it brings, all pure Python, no AVL binary:

- `asb.Airfoil(coordinates=...)` from the sliced sections; `.repanel()`.
- **NeuralFoil** (`get_aero_from_neuralfoil`, an XFOIL surrogate with Re,
  `n_crit`, `xtr`): C_l, C_d, C_m, transition locations in milliseconds --
  the sweep engine. Checked on the y = 1 m section at Re 3e5, Ncrit 12:
  C_l 0.79 / 0.91 at 2 / 3 deg against XFOIL's 0.83 / 0.95, transition
  0.68 / 0.67 c against 0.66 / 0.66: NeuralFoil reads this flapped section
  about 4-5 % low in C_l (a 0.4 deg bias in alpha at fixed C_L). Good enough
  for the design of experiments (sweeps over C_L, Ncrit, cut position), not
  for the final alpha_box: that takes the XFOIL polars (3.1).
- `asb.XFoil` wraps the installed `/usr/bin/xfoil` (polars with Ncrit, xtr,
  the Itr columns). On this workstation the Debian build (6.99.dfsg+1-3)
  needs an X display (`Cannot open display...aborting` without one) and
  with graphics off (`plop g`, which the wrapper sends) traps a floating-
  point exception inside `PACC` on the first converged point -- the wrapper
  reports every run as failed. Fix: run the wrapper from a session with a
  live `DISPLAY` (the 2026-10-25 polars came from one), or build XFOIL from
  source without `-ffpe-trap` and give the wrapper `xfoil_command=`. There is
  no `xvfb-run` on the workstation. (Replaces D4.)
- `asb.LiftingLine`: the NONLINEAR lifting line with viscous section data
  (NeuralFoil internally, Ncrit at its default 9 -- not settable through
  the class; the Ncrit 12 side of the band is done by polar inversion in
  `liftingline.py` and, finally, by the XFOIL polars), sweep and dihedral
  included; exposes the strip circulations, vortex centres and the induced
  velocity at any point, so alpha_i(y) comes out two ways (downwash at the
  bound vortex; polar inversion of the strip c_l).
- `asb.VortexLatticeMethod` (own VLM, inviscid) and `asb.AVL` (a wrapper
  that needs an `avl` binary) as cross-checks of the loading shape; the
  `asb.Opti` framework for the sweeps and for fixing C_L as a constraint.

### 3.1 The section polars

- Sections: `stl_sections.py section` at y = 0.5 (root, thickened), 1.0,
  2.0, 3.0, 5.0, 7.0, 9.0, 9.6 m (done; the scratch copies are not in git,
  regenerate them: seconds). Re(y) = 3e5 c(y) / 0.831, Mach 0.
- NeuralFoil polars at Ncrit 9 and 12 for every section (the band of D3),
  then XFOIL polars (`asb.XFoil`, `n_crit`, `max_iter 200`, alpha -4..10 by
  0.5) for the root region sections 0.5, 1.0, 2.0 at least -- these set the
  final alpha_box. Record C_l, C_d, C_m, Top_Xtr, Bot_Xtr; alpha_0 and the
  slope over the linear range. The XFOIL polar of the y = 1 m section at
  Ncrit 12 (2026-10-25): C_l 0.61 / 0.85 / 0.98 at 0 / 2 / 3 deg to the
  chord line, C_l,max 1.22 at 6 deg, upper-surface transition 0.66 c at low
  alpha (the bubble), lower surface laminar to the trailing edge. So at zero
  fuselage angle the wing already sits at C_l about 0.9: the flap is in its
  thermalling position (D2).
- Output: `tutorials/sailplane/junction/liftingline/` (polars as csv, the
  strip tables, one `README.md` with the tables).

### 3.2 The lifting line

`liftingline.py SECTION_DIR --cl C_L [--ncrit N] [--csv OUT]`: the planform
of section 1 is hard-wired (root chord carried to y = 0, incidence 2.3 deg,
the real sections), S_ref = the wing area, c_ref = 0.831 m. It sweeps alpha
(fuselage axis) -3..6, prints alpha for C_L 0.6 .. 1.1, and at the target
C_L the strip table with c_l(y), alpha_i(y) both ways and alpha_eff(y), the
root-region means over y = 0.55-1.05 m, the loading drop to the cut and
alpha_box. Spanwise resolution 8 strips per planform interval (every
0.06 m at the root). The AVL / VLM cross-check of the inviscid loading
shape is optional (3.0).

### 3.3 What step 1 settles

alpha_box for the chosen C_L; the size of the mirror's loading error at the
cut (whether 3 c suffices or 4 c is needed); the section's own transition
and bubble at the flight C_l (the XFOIL / NeuralFoil Xtr band, a first
reference for the DNS); and the confirmation that the fin-only STL loses
nothing. Route B (far field from a vortex lattice imposed on the box faces,
a new face kind; `asb.VortexLatticeMethod` would provide the field) is held
in reserve and taken only if step 2 shows the loading gradient matters.

### 3.4 First result (2026-10-26, NeuralFoil section data, Ncrit 9)

| quantity | C_L 1.0 | C_L 0.8 |
|---|---|---|
| alpha_fuselage (chord line) | 2.24 (4.54) deg | 0.08 (2.38) deg |
| root-region c_l (y 0.55-1.05 m) | 1.04 | 0.85 |
| alpha_i root region, downwash / polar | 0.45 / 0.66 deg | 0.46 / 0.41 deg |
| **alpha_box (fuselage axis)** | **1.6-1.8 deg** | **-0.35 deg** |
| loading drop y 0.6 m -> cut 2.95 m | 4.2 % | 7.5 % |
| C_L,max of the lifting line (NeuralFoil) | 1.11 at alpha 4 | |

Caveats: the NeuralFoil bias (3.0) moves alpha at fixed C_L by about
-0.4 deg and alpha_i by less than 0.1 deg; the strips over the thickened
root section (y 0.5-0.75) show a 0.3 deg wobble between the two alpha_i
estimates (the section polar changes there) -- the root-region mean is the
number to use; C_L 1.2 (D2) is NOT reachable: the wing stalls at C_L 1.1 in
this model, so the choice is between 0.8 (cruise-ish with this flap) and
1.0 (thermalling). The loading drop to the cut is 4-8 %, larger than the
elliptic estimate; at the DNS core (y < 0.9 m) it is under 1 %.

## 4. Step 2: the calibration box

A cheap truncated model at the production topology, run three ways. It
confirms the set-up (body through the outlet, the cut), measures the two
sensitivities that decide whether 0.1 deg in alpha_box is good enough, and
provides the initial field of the production run (`interp_restart.py`).

- Geometry: the STL rotated by alpha_box (section 1.2), `stl_scale` to
  chord units. Box in c: x from -4 (3.3 m ahead of the nose; the nose at 0,
  the wing LE at 3.2 c, the TE at 4.2 c) to the outlet at 6.7 c (TE + 2.5 c,
  through the boom: half-width 0.18 c there); y from 0 to y_cut = 3.55 c
  (2.95 m) -- and a second prepare at 4.55 c; z 20 c tall, the wing in the
  middle. Base grid Delta = c/16 (171 x 57 x 320 = 3.1 M), `refine_body` to
  c/512 at the wing and forebody (Delta+ about 26 at Cf 0.004: a wall-
  function RANS / IDDES wall, NOT resolved) -- about 120 M cells, 4 A100 at
  about 0.4 s/step, 10 c/U in 5 000 steps: under an hour per run.
- Model: `iddes` (RANS mode at this wall resolution, with the wall
  functions; transition under IDDES is an open item, so the calibration
  loading carries the fully-turbulent decambering error, section 4.1).
  Chebyshev `niter = 12` unless the outlet case needs red-black (the
  tutorial's note); `cflmax = 1.0`.
- Runs: alpha_box, alpha_box +- 0.5 deg at y_cut = 3.55 c; alpha_box at
  4.55 c. Measure the strip loads C_l(y) (section 5.2's script), the root
  Cp at y = 0.5 and 1.0 m, the junction separation line (the Cf sign on the
  wall), and the outlet's health with the boom through it (no backflow
  blow-up: the total-pressure datum).
- Read-out: dC_l/dalpha at the station (expected about 0.1 per deg);
  how much the junction Cp and separation line move per 0.5 deg (on a
  bubble section about 5 % chord per degree is the literature value, so
  0.1 deg is below the floor of section 5.1); the 3 c vs 4 c difference at
  the root (expected under 1 %: then 3 c stands).

### 4.1 Why alpha is TRANSFERRED, not fitted

The box's C_l at the station could be matched to the target by adjusting
alpha, but the calibration model's viscous decambering error (a wall-
function RANS misses the bubble: Delta C_l 0.03-0.05) is the same size as
the correction itself (2 pi alpha_i = 0.07) and would be fitted INTO alpha.
So alpha_box is taken from the lifting line (whose own C_l error is XFOIL's,
much smaller) and the box runs only verify consistency within the model's
error and measure sensitivities. The production run's own root C_l against
the step-1 target is the post-hoc check (section 5.5).

## 5. Steps 3-7

### 5.1 Step 3: the disturbance floor and the IBM (go / no-go)

Both from the 2026-10-24 reports, each a few A100-days, both BEFORE the
production grid exists:

- **Projection-residual floor.** The unforced SD7003 bubble at Re 6e4
  (`tutorials/naca/LES/` has the geometry and the strategy): the transition
  front and the bubble length at `niter = 12` vs 24 (Chebyshev) and under
  red-black with more sweeps. A moved front is a finding: the floor is then
  too high for a no-trip study and the projection in the DNS band needs
  more iterations (cost) -- decide then, never a compensating trip.
- **TS-wave growth through the IBM.** The Blasius case (`validation/blasius`)
  with a small harmonic inflow perturbation at a branch-I/II Reynolds number,
  wall as an immersed plate at Delta+ 3; the growth rate against Orr-
  Sommerfeld. The first-order cut-cell wall must not damp or seed the waves
  beyond a few percent of the inviscid growth. Fail = the DNS core needs a
  body-fitted-quality wall treatment, which this solver does not have:
  stop and report.

### 5.2 Step 4: solver increments (small)

- **Zonal model-off levels**: `[turbulence] dns_levels_from = L` sets the
  eddy viscosity to zero on blocks of level >= L, so the DNS core is
  literally model-free and the record says so (IDDES there would in any case
  be inert in laminar flow with k = 0 inflow and vanish at DNS resolution
  behind the breakdown, but the argument is not the record). Gated: bit-
  exact outside its levels (nofma, the 7/9/11 suites untouched), the IDDES
  case with the key at a level above the grid identical to without it.
- **Strip loads**: `tutorials/sailplane/junction/strip_loads.py` -- the
  surface pressure sampled off the wall along the STL normals (the
  `surface_cp_cf.py` way), integrated per spanwise strip over the triangles
  of the rotated STL -> C_l(y), C_m(y), plus the wall Cf sign for the
  separation line. Validated on the NACA (2D: every strip equal to the CV
  force's lift within the IBM sampling error).

### 5.3 Step 5: production geometry and grid

- D1 decides the fillet. Without one the corner is sharp and the corner
  transition is of the sharp corner; with one the STL must be re-tessellated
  (the V0 tessellation is 166 segments per wing section: fine for RANS, too
  coarse for a Delta+ 3 wall -- the level-set fraction sees the facets).
  Either way a fine STL is produced (the CAD source, not the STL, if it can
  be had) and checked by `stl_sections.py` (sections smooth at the DNS
  resolution).
- Levels: base c/16; `refine_body_box` nested: the DNS core box (x 3.0-4.4 c,
  y 0-0.5 c beyond the fuselage side, z around the upper surface) to
  Delta+ 4 = c/4096 (8 levels), the LES ring to c/2048, the forebody to
  c/1024, the far field at the base. `refine_dims = xyz`. The DNS / LES
  boundary sits OUTSIDE the laminar layer and away from the bubble, and no
  2:1 interface may cross the bubble's reattachment (the measured SD7003
  motivation for the dropped S5b): put the core's downstream face behind
  the reattachment line of the calibration run.
- Prepare on HoreKa CPU ranks at 3-5 M leaves (never done above 1 M: time
  it on the calibration grid first, report the scaling).
- `keep_buried` on (forces), `nb` dividing the grid, `report_partition`
  read on the first run (cross-node face cells).

### 5.4 Step 6: the exchange at 8-16 nodes

The measured cross-node stall (19-25 % of the step at 8-16 GPUs) is a
first-order cost at 30-55 A100: the red-black reduced exchange round and
the `phi_exchange` locality item (CLAUDE.md open list 1) move to the top
of the performance list and are measured on an 8-node slice of the real
grid (the calibration grid refined one level), before the allocation is
spent.

### 5.5 Step 7: production

- Run A: Delta+ 4 core, alpha_box, chunked on HoreKa (`submit_naca_stage.sh`
  generalised: it already chains, resumes the traces and primes the CV);
  start from the calibration field interpolated (`interp_restart.py`),
  10 c/U transient, then 10 c/U of statistics at the wall (Cf, Cp,
  transition line) and in the core (spectra of the wall-normal velocity in
  the bubble and in the corner).
- Check against step 1: the root C_l(y) from the strip loads vs the
  lifting-line target (within the XFOIL band); if it misses by more than
  the band, the remaining item is the mirror's loading gradient -> Route B.
- Run B: Delta+ 3 core, same everything: the resolution check of the
  transition line (within one bubble-length scatter = converged).
- Records: a `tutorials/sailplane/junction/README.md` in the RE-VALIDATED
  convention (binary hash, ini, settings, every number), the worklog bullet,
  `validation/README.md` row.

## 6. Decisions for the user

- **D1 fillet.** The STL has none ("ohne Rundung"). Study the sharp corner
  as built, or obtain / construct a fillet? (The real aircraft has one; the
  sharp corner makes the corner vortex and its transition a different
  problem.)
- **D2 flight C_L.** The flap is set for thermalling (C_l about 0.9 at zero
  fuselage angle; the lifting line stalls at C_L 1.1). Proposal: C_L = 1.0
  (the known root-separation regime of sailplanes in thermalling flight,
  alpha_box about 1.7 deg), with 0.8 (alpha_box about -0.35 deg) as the
  comparison if a second production run is ever affordable.
- **D3 XFOIL Ncrit.** Carry 9 and 12 as the band (proposal); the DNS floor
  measured in step 3 says which end applies.
- **D4 XFOIL binary.** The Debian build fails under the AeroSandbox wrapper
  without an X display (section 3.0): run the polars from a session with a
  live `DISPLAY` (proposal), or build XFOIL from source. AVL is not needed.
- **D5 allocation.** 10 000-20 000 A100-h on HoreKa: a project proposal has
  to be filed before step 7; steps 1-6 fit the exasim share and
  `hk-project-p0024333` (calibration under 100 A100-h, the physics tests a
  few hundred).
