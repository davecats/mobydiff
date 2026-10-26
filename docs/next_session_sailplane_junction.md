# Next session: the sailplane wing-fuselage junction (transition study)

STATUS (2026-10-26 evening, head `8aaca29` + this handout): step 1 DONE
(`tutorials/sailplane/junction/liftingline/README.md`: **alpha_box = 1.5 deg**
at C_L 1.0), decisions D1-D5 taken (section 6). **The next session starts
at section 4 (step 2, the calibration box), written as the executable
procedure: 4.0 the strip-load script first, then 4.1-4.4.** Nothing of
step 2 exists yet: no ini, no rotated STL, no case file, no script. The
HoreKa account for every job of this study is `hk-project-p0024333` (D5).
Five commits from `4d1759e` on are unpushed (push when the user says so).

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
| 1 | Route A lifting line in AeroSandbox: XFOIL section polars + the nonlinear lifting line of the whole wing -> alpha_i(y), C_l(y), alpha_box (section 3) | **DONE 2026-10-26**: alpha_box 1.5 deg (`liftingline/README.md`) |
| 2 | The calibration box (section 4): 2a `strip_loads.py` validated on the NACA; 2b four case files (alpha 1.0 / 1.5 / 2.0 at the 3 c cut, 1.5 at 4 c); 2c four SST wall-function runs on HoreKa; 2d read-out: C_l(y) vs the lifting line, dC_l/dalpha, the cut, the outlet through the boom | PENDING -- NEXT |
| 3 | Go / no-go physics tests: projection-residual floor on the SD7003 bubble, TS growth through the IBM (section 5.1) | PENDING |
| 4 | Solver increment: zonal model-off levels (section 5.2; the strip-load script moved to 2a) | PENDING |
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
  reports every run as failed. DONE 2026-10-26 (D4): XFOIL 6.99 is built
  from Drela's source by `tutorials/sailplane/junction/assets/xfoil/build.sh`
  (double precision, no FPE traps, `Xfoil/` and the tarball git-ignored),
  and `xfoil_polars.py` drives that binary DIRECTLY -- the wrapper also
  refuses the source build's polar files (the `CINC` column carries no
  header label, the Itr columns are absent). AVL is not needed.
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
- XFOIL polars (`xfoil_polars.py`: `PPAR N 240`, `ITER 300`, alpha -4..10
  by 0.5 as two `ASEQ` branches from 0) at Ncrit 9 and 12 for every section
  (the band of D3): DONE, `liftingline/polars/*.csv` with alpha, C_l, C_d,
  C_dp, C_m, Top_Xtr, Bot_Xtr; the summary table is in
  `liftingline/README.md`. The XFOIL polar of the y = 1 m section at
  Ncrit 12 (2026-10-25): C_l 0.61 / 0.85 / 0.98 at 0 / 2 / 3 deg to the
  chord line, C_l,max 1.22 at 6 deg, upper-surface transition 0.66 c at low
  alpha (the bubble), lower surface laminar to the trailing edge. So at zero
  fuselage angle the wing already sits at C_l about 0.9: the flap is in its
  thermalling position (D2).
- Output: `tutorials/sailplane/junction/liftingline/` (polars as csv, the
  strip tables, one `README.md` with the tables) -- DONE.

### 3.2 The lifting line

`liftingline.py SECTION_DIR --cl C_L [--ncrit N] [--polars DIR] [--csv OUT]`: the planform
of section 1 is hard-wired (root chord carried to y = 0, incidence 2.3 deg,
the real sections), S_ref = the wing area, c_ref = 0.831 m. It sweeps alpha
(fuselage axis) -3..6, prints alpha for C_L 0.6 .. 1.1, and at the target
C_L the strip table with c_l(y), alpha_i(y) both ways and alpha_eff(y), the
root-region means over y = 0.55-1.05 m, the loading drop to the cut and
alpha_box by two routes: the downwash route alpha_flight - alpha_i (inherits
the lifting line's NeuralFoil alpha at C_L) and the polar route
alpha_eff(c_l,root) - 2.3 with `--polars` on the XFOIL polars (independent
of it; THE number). Spanwise resolution 8 strips per planform interval
(every 0.06 m at the root). The AVL / VLM cross-check of the inviscid
loading shape is optional (3.0).

### 3.3 What step 1 settles

alpha_box for the chosen C_L; the size of the mirror's loading error at the
cut (whether 3 c suffices or 4 c is needed); the section's own transition
and bubble at the flight C_l (the XFOIL / NeuralFoil Xtr band, a first
reference for the DNS); and the confirmation that the fin-only STL loses
nothing. Route B (far field from a vortex lattice imposed on the box faces,
a new face kind; `asb.VortexLatticeMethod` would provide the field) is held
in reserve and taken only if step 2 shows the loading gradient matters.

### 3.4 Result (2026-10-26; `liftingline/README.md` has the full tables)

| quantity | C_L 1.0 (D2) |
|---|---|
| alpha_fuselage at C_L 1.0 (lifting line, NeuralFoil sections) | 2.24 deg |
| root-region c_l (y 0.55-1.05 m); at the cut (2.95 m) | 1.043; 1.022 (**-4.2 %**, < 1 % over the DNS core) |
| alpha_i root region (downwash at the bound vortex) | 0.45 deg |
| alpha_eff root region from the XFOIL polars, Ncrit 9 / 12 | 3.76 / 3.81 deg to the chord line |
| **alpha_box, polar route, Ncrit 9 / 12** | **1.46 / 1.51 deg -> 1.5 deg** |
| alpha_box, downwash route | 1.79 deg |
| C_L,max of the wing | 1.11 (C_L 1.2 unreachable) |

The polar route is taken (`alpha_eff(c_l,root) - 2.3`): at AR 30 the loading
is set by the planform, and the XFOIL polar does the rest; the downwash
route carries NeuralFoil's 0.3-0.4 deg bias on this flapped section. The
Ncrit band is worth 0.05 deg; the 0.3 deg between the routes is covered by
the +-0.5 deg runs of step 2, and if the junction turns out sensitive at
that level the refinement is a lifting line on the XFOIL polars themselves.
C_L 0.8 would be alpha_box about -0.35 deg (loading drop 7.5 %). First DNS
reference from the polars: the free-wing bubble transition at about 0.63 c,
the thickened root section at about 0.40 c, lower surface laminar.

## 4. Step 2: the calibration box -- the procedure for the next session

A cheap truncated model at the production topology (chord units, the STL
pitched by alpha_box, the 3 c cut, the boom through the outlet), run four
times with SST + wall functions. It confirms the set-up, measures the two
sensitivities that decide whether 0.3 deg in alpha_box matters (section
3.4), provides the initial field of the production run, and is the first
prepare above 20 000 leaves. Costs: four prepares on the workstation (CPU,
minutes each), four 1-hour `dev_accelerated` jobs on 4 A100.

Two facts found 2026-10-26 that shape it:

- **There is no runtime force output for this case.** `[case] name =
  airfoil` (the CV budget, `forces.txt`) requires a PERIODIC span and
  inlets on every non-outlet face; the sailplane is `name = generic` with
  a symmetry plane and free-slip faces. Every load of this study comes
  from the field file by post-processing -> 4.0 first.
- **Wall functions and the transition model are both REJECTED under
  `iddes`** (`docs/configuration.md` [rans]); the calibration wall is
  Delta+ about 27, so the model is `rans` + `sst` + `wall_function`, fully
  turbulent (no bubble: the decambering error is exactly why alpha is
  transferred, not fitted, 4.1). This also flags a production-run question
  for step 5: the LES ring at Delta+ 8-10 puts the first cell centre at
  y+ 4-5 under IDDES's `resolved` wall treatment -- acceptable for a
  wall-resolved LES, but it must be stated and checked on the SD7003 case
  of step 3, not assumed.

### 4.0 Step 2a: `tutorials/sailplane/junction/strip_loads.py`

Surface pressure and skin friction on the STL body from a block-table
field snapshot, integrated per spanwise strip:

    strip_loads.py FIELD.h5 --case CASE.h5 --stl ROTATED.stl --scale S --translate X Y Z
                   [--strips 0.1] [--chord 1.0] [--out strips.npz] [--cp-station Y ...]

- Geometry: the triangles of the rotated STL (ASCII reader of
  `stl_sections.py`), transformed exactly as `moby_prepare` does
  (`v * scale + translate`), centroids, outward normals, areas. The wing
  is the triangles with x < (x_TE + 0.3 c) and y > the fuselage side, the
  strips are bins in y of width `--strips` (chords) from the fuselage
  side to the cut; a "fuselage" bin takes the rest.
- Cp at each centroid: the solver's pressure `pn` sampled at centroid +
  1.5 Delta_local along the normal (the first fluid cell outside the
  penalization band; trilinear from the block holding the point -- the
  block table and the finest-lattice reassembly tools are in
  `tools/compare_fields.py` and `validation/prepare/compare_snapshots.py`;
  the NACA script `tutorials/naca/rans/postProcess/surface_cp_cf.py`
  samples the same way on an analytic 2D polyline and is the reference
  for the band logic). C_p = 2 (p - p_inf), p_inf from the inlet plane.
- Cf sign per centroid from the penalization-band estimator of
  `surface_cp_cf.py` (the band cells carry the wall gradient; a through-
  origin fit of the staggered velocity against the tangential distance):
  only the SIGN along the chordwise tangent is needed here (the
  separation line), the magnitude is a bonus.
- Per strip: c_l = sum(-C_p n_z dA)/(c_strip w) with the strip chord and
  width, c_m about the quarter chord, and the chordwise position where Cf
  changes sign on the upper surface (separation / reattachment).
- Cp stations (`--cp-station 0.6 1.0`, metres from the symmetry plane,
  converted with the scale): the upper and lower surface C_p(x/c) at that
  y, written as csv.
- VALIDATION before any sailplane use (the gate of 2a): the NACA nose case
  (`~/horeka_stage/naca/c11_nose_660000.h5` with
  `tutorials/naca/rans/assets/geometry/ibm_coeff_c11_nose.h5`; the STL is
  not there -- the NACA is an analytic section, so the script takes an
  `--analytic-naca` switch OR the 2D polyline is written as a one-strip
  STL by a 20-line helper): the single strip's c_l must match the runtime
  CV lift 0.5205 +- 0.0005 and `surface_cp_cf.py`'s C_p curve
  (`assets/reval_2026-10-24/cpcf_c11_nose_final_cp.dat`) to the IBM
  sampling error (a few 1e-3 in C_p, about 1 % in c_l). Record the
  numbers in the script's docstring.

### 4.1 Step 2b: geometry and case files (workstation, CPU build)

Chord units: c = 0.831 m, `stl_scale = 1.2034e-3` (= 0.001 / 0.831),
`re = 3.0e5`, U = 1. Four STLs / case files in
`tutorials/sailplane/junction/calib/` (the STLs and case files are NOT
committed -- 14 MB each / hundreds of MB; the inis and the rotation
commands are):

    cd tutorials/sailplane/junction
    for a in 1.0 1.5 2.0; do
      python3 stl_sections.py rotate "../FRUE V0 ohneRundung.stl" $a calib/frue_a$a.stl
    done

(`rotate` pitches about (x, z) = (2.84, 0.47) m in STL metres, nose-up
positive; the fuselage-axis angle of the run is the rotation.)

The ini `calib/calib_a1.5_y3.5.ini` (the other three differ in `stl_file`
and `ny`/`ly` only); every number below is in chords:

    [case]
    name = generic
    file = calib_a1.5_y3.5.case.h5
    [grid]
    nx = 176          ; 11.0 c: inlet 4 c ahead of the nose, outlet 2.8 c behind the root TE
    ny = 56           ; 3.5 c: symmetry plane to the cut (fuselage side at 0.54 c + 3 c) -- the 4 c variant: ny = 72, ly = 4.5
    nz = 320          ; 20 c, the wing at mid-height
    lx = 11.0
    ly = 3.5
    lz = 20.0
    [grid.x] / [grid.y] / [grid.z]  distribution = uniform
    [blocks]
    nb = 8                      ; 22 x 7 x 40 = 6160 base leaves of 512 cells
    refine_body = true
    refine_body_levels = 5      ; c/512 at the body: Delta+ about 27 (Cf 0.0045)
    refine_levels = 5
    refine_dims = xyz
    keep_buried = false         ; no penalization force is read at runtime; set TRUE if 4.0 ends up using the penalization integral
    [ibm]
    enabled = true
    stl_file = "frue_a1.5.stl"
    stl_scale = 1.2034e-3
    stl_translate = 4.0 0.0 9.434   ; nose at x = 4 c; y = 0 is the STL's symmetry plane; the root quarter chord at z = 10 c (0.47 m / 0.831)
    [flow]
    re = 3.0e5
    [turbulence]
    model = rans
    [rans]
    model = sst
    wall_treatment = wall_function
    tu = 1.0
    nut_ratio = 10.0
    ambient_sustain = false
    [time]
    nsteps = 100000000
    t_final = 15.0              ; 10 c/U transient + 5 c/U of samples (4.3)
    cflmax = 1.0
    pecletmax = 0.5
    dtmax = 2.0e-3
    [pressure]
    accel = chebyshev
    niter = 12                  ; fall back to solver = redblack, sor = 1.5, niter = 6 if the outlet case diverges (the tutorial's note)
    [output]
    field_interval = 2000
    field_prefix = calib_a1.5_y3.5
    [boundary]                  ; as tutorials/sailplane/input.ini: x_min inflow (1,0,0) Neumann p, x_max_patch = outlet, y and z faces free-slip

Prepare each on the workstation CPU build and keep the prepare logs
(leaf counts, the `report_partition` line, wall time -- the first data
point of the prepare scaling for step 5):

    mpirun -n 4 ../../../build_cpu/moby_prepare calib_a1.5_y3.5.ini 2>&1 | tee calib_a1.5_y3.5.prep.log

Expected: about 55 000 leaves (the c/512 band 41 000, c/256 10 000, the
rest coarser), 28 M cells, under 15 min on 4 ranks (the 18 000-leaf
tutorial prepare takes 2 min on 2). `moby_prepare --check-config` first.
Then a 20-step solve on the RTX 3060 on 1 and 2 ranks of the smallest
case (`max_abs 0` on `un vn wn pn`: rank identity, the smoke test of the
boom through the outlet).

### 4.2 Step 2c: the runs (HoreKa, 4 A100, one job each)

Ship the four case files and inis to `$W/calib_run/` (`$W` as in the
NACA re-validation; the worktree `moby-2to1-reval` is at `532ef54`, so
push the current head there first: `git push horeka:$W/moby-2to1-code
<sha>:refs/heads/junction` and `git -C $W/moby-2to1-code worktree add
--detach $W/moby-junction junction`). The job is one chunk of
`tutorials/naca/rans/horeka/submit_naca_stage.sh` with `NO_CHAIN=1`
(it builds the worktree if needed, restarts from the latest snapshot,
writes `provenance.txt`) -- copy it to `horeka/submit_calib.sh` with
`--account=hk-project-p0024333`, the `airfoil`-specific `done_already`
test replaced by "the log says main loop ended", and PREFIX = the ini's
`field_prefix`:

    export CODE_DIR=$W/moby-junction RUN_DIR=$W/calib_run HDF5_ROOT=$HOME/hdf5 \
           STAGE_INI=calib_a1.5_y3.5.ini PREFIX=calib_a1.5_y3.5 NSTEPS=12000 NO_CHAIN=1
    cd $RUN_DIR && sbatch --export=ALL --account=hk-project-p0024333 submit_calib.sh

Step time: the 17 M-cell NACA ran 0.062 s/step on 4 A100; 28 M cells
about 0.11 s/step; dt about 1.3e-3 at cflmax 1 (Delta = c/512) -> 15 c/U
in about 11 500 steps, 20-25 min: one `dev_accelerated` hour holds a run
with margin (if not, the chain mode resumes it). Four jobs, submitted in
sequence (one dev job per user at a time); the user's own HoreKa jobs are
not touched. Read `report_partition` and the first 100 steps' log before
leaving a job alone: the divergence residual, the Courant number, no
`WARNING ... exceeds the RK3 limit`.

### 4.3 Step 2d: the read-out

`strip_loads.py` on the final snapshot and on the 3 earlier ones of the
last 5 c/U (the mean and the scatter of each strip):

| measurement | from | expected / decides |
|---|---|---|
| c_l(y) over y = 0.54..3.5 c at alpha 1.5 | strips 0.1 c | shape vs the lifting line's (`liftingline/strips_cl1.0_ncrit12.csv`): the root-to-cut drop 4.2 % there; absolute c_l lower than 1.04 by the fully-turbulent decambering (expected 0.03-0.06) -- the SHAPE is the check, not the level |
| dc_l/dalpha at y = 0.6-1.0 m | alpha 1.0 / 1.5 / 2.0 | about 0.1 per deg (2 pi per rad x 0.9); linear over the 1 deg |
| junction C_p at y = 0.6 and 1.0 m, alpha 1.0 / 1.5 / 2.0 | `--cp-station` | how far the suction peak and the upper-surface pressure recovery move per 0.5 deg: this is the sensitivity the 0.3 deg uncertainty of alpha_box is measured against |
| separation line (Cf sign) on the upper surface near the root, alpha 1.0 / 1.5 / 2.0 | per strip | with SST wall functions a turbulent separation near the TE / corner at most; its motion per 0.5 deg |
| 3 c vs 4 c cut at alpha 1.5 | the two grids | root-region c_l within 1 % -> the 3 c cut stands for production; else 4 c |
| the outlet through the boom | the log, the snapshot | no backflow blow-up (the total-pressure datum), divergence residual as in the NACA (1e-5..1e-4) |
| prepare at 55 000 leaves | the prepare logs | wall time and memory on 4 workstation ranks: extrapolate to the production 3-5 M leaves (step 5) |

Record in `tutorials/sailplane/junction/calib/README.md` (RE-VALIDATED
convention: binary hash, ini, every number), the ledger row 2, the worklog
bullet. The alpha 1.5 / 3 c final snapshot is the initial field of step 7
(`interp_restart.py`), keep it on HoreKa and in `~/horeka_stage/junction/`.

### 4.4 Why alpha is TRANSFERRED, not fitted

The box's c_l at the station could be matched to the target by adjusting
alpha, but the calibration model's viscous decambering error (a wall-
function SST has no bubble: Delta c_l 0.03-0.06) is the size of the
correction itself (2 pi alpha_i = 0.05) and would be fitted INTO alpha.
So alpha_box comes from the lifting line (section 3.4) and the box runs
only verify consistency within the model's error and measure
sensitivities. The production run's own root c_l against the step-1
target is the post-hoc check (section 5.5).

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
- **Strip loads**: moved to step 2a (section 4.0) -- the calibration box
  cannot be read without it.

### 5.3 Step 5: production geometry and grid

- D1 decides the fillet. Without one the corner is sharp and the corner
  transition is of the sharp corner; with one the STL must be re-tessellated
  (the V0 tessellation is 166 segments per wing section: fine for RANS, too
  coarse for a Delta+ 3 wall -- the level-set fraction sees the facets).
  Either way a fine STL is produced (the CAD source, not the STL, if it can
  be had) and checked by `stl_sections.py` (sections smooth at the DNS
  resolution).
- The LES ring under IDDES: `wall_treatment = resolved` is the only one
  admitted under `iddes`; at Delta+ 8-10 (isotropic cells) the first cell
  centre sits at y+ 4-5. State it, and check it on the SD7003 case of step
  3 (the bubble's reattachment and the transition front at Delta+ 8 vs 4
  under IDDES) before the ring resolution is fixed.
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

## 6. Decisions -- TAKEN 2026-10-26 (user)

- **D1 fillet: the sharp corner as built, at first.** (A filleted variant
  is a later study; the fine STL of step 5 is cut from the same CAD.)
- **D2 flight C_L = 1.0** (thermalling; alpha_box 1.5 deg). 0.8 is the
  comparison if a second production run is ever affordable.
- **D3 Ncrit 9 and 12 carried as the band** (worth 0.05 deg in alpha_box).
- **D4 XFOIL built from source** under `tutorials/sailplane/junction/
  assets/xfoil/` (`build.sh`); AVL not needed (AeroSandbox's own lifting
  line), the user's offer of an AVL build stands if a cross-check is ever
  wanted (it would go to `assets/avl/`).
- **D5 the hours exist: project `hk-project-p0024333` on HoreKa** -- every
  job of this study, including the production runs, is submitted with
  `--account=hk-project-p0024333` (not the exasim share).
