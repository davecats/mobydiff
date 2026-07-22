# NACA 0012 Re 4e5 campaign — COMPLETE 2026-07-21 (C10 sweep done; results below)

STATUS (2026-07-21): the campaign RESUMED on the prepare/solve-split code
(moby_prepare built ibm_coeff_c10.h5 in 4m26s, 10615 leaves / 5.4M cells)
and the FULL C10 sweep (alpha = -2..5, Re 4e5, SST + `[rans]
ambient_sustain`, tu = 5 % / nut_ratio = 10, dt = 1e-4, 150k steps to
t = 15 per angle) COMPLETED across istmcorax (5/4/3/2), istmcetus GPU1
(1/-1) and the local 3060 (-2); aoa0 ran first as the verification.
KEY FIX on resume: the original nut_ratio = 1000 ambient DESTROYED the
first aoa0 (explicit eddy-diffusion dt bound ~1e-5 vs dt 5e-5; zombied at
dt 1e-8) — the Rumsey ambient-sustain sources (commit 94a9249) make
(k_inf, omega_inf) an exact fixed point so tu = 5 % arrives at the body
(measured EXACTLY 5.00 % at 2c upstream) with nut_inf = 10 nu only.
aoa0 gates ALL PASS: CV forces box-independent (C_D 0.01185/0.01187/
0.01240 at 1.5/2.5/4c), C_L(0) = -0.0003, NO interface artifacts (the R1
fan is GONE at the C10 interface distances; striping probe at distant
interfaces at/below the quiet-band floor; band_filter stays OFF), only
the surface-attached LE staircase jaggedness in Cf at x/c < 0.05 remains.
RESULTS (tutorials/naca: polars_c10.png, polar_mobydiff.dat,
cv_polar_raw.txt, cpcf_c10_aoa*_{cp,cf}.dat; XFOIL refs
xfoil_re4e5_n{1,9}.dat — the Debian xfoil FPE-crashes on any second
viscous point and on PACC: run ONE ALFA PER PROCESS and parse stdout):
C_L antisymmetric to 0.3 % (+-2: 0.1804/-0.1810); slope 0.0884/deg = 81 %
of 2pi (XFOIL n1: 0.108/deg); C_L(5) = 0.432+-0.022 vs XFOIL-n1 0.533
(-19 %, deficit grows with alpha); C_D(0) = 0.0128+-0.0018 vs n1 0.0104.
Suspected deficit drivers (unproven, next investigation): y+ 2.7-4.1
resolved-wall under-resolution (T3: first cells below y+ 30 carry the
10-19 % log-line error class), LE staircase suction-peak smearing
(measured peak -1.58 at alpha 5 vs XFOIL ~ -1.75), first-order scalar
upwind transition front.
OPENFOAM REFERENCE (compared 2026-07-21, tutorials/naca/
compare_openfoam.py + cpcf_vs_openfoam_aoa5.png; case
~/auswertung/20251027_MA_JannikWeber/run: simpleFoam kOmegaSST OF7,
alpha 5, Re 4e5, 48400-cell body-fixed 2D mesh, wall y+ avg 1.31,
decaying inlet turbulence tuned to Tu ~ 7 % at the LE): Cl = 0.5142,
Cd = 0.01339 — close to XFOIL n1 (0.533/0.0124). vs mobydiff C10
(0.432 +- 0.022 / 0.0157 +- 0.003): C_L -16 %, C_D +17 %. Surface
comparison: Cp SHAPES agree, deficit is the suction-side loading in
x/c < 0.4 (peak -1.59 vs -1.78, -11 %); pressure side matches; Cf
midchord levels agree (~0.0075) BUT OF shows a pseudo-laminar nose dip
(Cf ~ 0.002 at x/c 0.05-0.1, transition at ~0.12) while C10 is
turbulent from the staircase LE, and aft of x/c 0.4 the C10 suction Cf
sags below OF (thicker decelerated BL); TE separation onset 0.96 (C10)
vs 0.996 (OF). All consistent with the y+ 3-4 wall under-resolution +
LE staircase smearing hypothesis — the OF case resolves y+ 1.3.
L11 RESOLUTION STUDY (2026-07-22, aoa 5 only): .prep_c11.ini ->
ibm_coeff_c11.h5 (refine_levels = 11, 16042 leaves / 8.2M cells,
Delta11 = c/6144, y+ ~ 1.4-2.1); IC INTERPOLATED from the converged L10
field (interp_restart.py, t = 15 carried, zero uncovered cells), run
t = 15..20 at dt = 5e-5 on cetus GPU1 (0.45 s/step). VERDICT: the wall
resolution was the DRAG driver, not the lift driver — C_D 0.0157 ->
0.0126 +- 0.002 (now matches OF 0.0134 / XFOIL-n1 0.0124 within CV
scatter); C_L 0.432 -> 0.445 +- 0.024 (deficit vs OF -16 % -> -13 %);
Cp_min -1.59 -> -1.61 vs OF -1.78 (the suction peak BARELY moved ->
the remaining lift gap is NOT wall resolution; LE staircase smearing is
the prime suspect); TE separation onset 0.959 -> 0.983 (OF 0.996).
NOTE OF's transition is FORCED (fvOptions: k = 0 for x < 0.09c both
sides + k-source 2e-5 trip strips at x = 0.10-0.12c hugging the
surface) — a matching [rans] k-pin/trip increment is the agreed next
step (dormant-off bit-exact gate) before attacking the staircase.
OOM LANDMINE FIXED (cv_forces.py): cv_force painted the WHOLE control
box on the finest lattice (~200 GB at L11 margin 4c; the machine has
62 GB) — repeatedly OOM-killed the post-processing. Now paints four
thin border strips (161 MB peak, C10 values reproduced EXACTLY).
plot_c10_turb_fields.py windows carry the same risk at L11+ depths:
estimate cells x 6 x 8 B before choosing a window.

The pre-resume pause notes (fan analysis, R1 history) follow unchanged.

## The fan question: is it a 2:1-interface problem?

Short answer: the fan's ORIGIN is (with strong control-run evidence) NOT
the 2:1 interface — it is the staircase-IBM seed amplified by the
CENTRED convection scheme at cell-Reynolds >> 2. But the interfaces play
a REAL, distinct, secondary role: they partially reflect and alias the
parasite (near-Nyquist content CANNOT be transmitted to a coarser level
by any consistent transfer), and in snug nests — where the first
interface sits inside the fan core — the two mechanisms compound and are
hard to separate. The evidence, in order:

1. ORIGIN is seed + scheme, not interface (R1 controls,
   validation/naca0012 README):
   - the fan's dominant ripple wavelength is 2.44 fine cells with 68 %
     of its energy near the grid Nyquist — the stationary parasite of
     centred convection at high cell-Re;
   - it correlates purely with cell-Reynolds (Re 1e5 at cell-Re 146:
     fan; Re 1e3 at cell-Re 1.5: none), independent of the turbulence
     model and of niter;
   - the "fan-box" control (fine MEDIUM around an unchanged coarse
     ring) did NOT reduce the near fan: the RING SEED controls it;
     refining the ring (L5) collapsed it 15x;
   - the near strips (<= 0.023c) persist inside a single-level window —
     no interface needed to sustain them.

2. The INTERFACE is responsible for the far-field pattern and scatter:
   - R1 verdict (3): "the FAR strips were largely LEVEL-INTERFACE
     artifacts, not parasite: with the whole analysis window at one
     level they drop to exactly 0.0";
   - the Re 4e5 / L6-xz interface-overlay figure
     (tutorials/naca/flow_aoa0_interfaces.png) shows ripple STREAKS
     locked onto the interface lines below the LE (z ~ 5.90-5.93):
     coarse-side striping where the fan hits the L6->L5->L4 jumps;
   - mechanism: a 2*Delta_fine wave is beyond the coarse level's
     Nyquist. NO consistent 2:1 transfer can transmit it; the
     const-1/2 restriction necessarily returns it as part reflection,
     part aliased coarse content. The interface does not CREATE the
     parasite, but it scatters it — and a snug nest places that
     scattering surface right at the fan core (0.005-0.01c from the
     nose at L6-xz/Re4e5).

3. The interface MACHINERY itself is validated clean for RESOLVED
   content: uniform flow exact (0.0) through every orientation; Beltrami
   interface error BELOW the octree reference at equal base; Re_tau 180
   xz wall-band channel with NO interface band (u'/v'/w'/-u'v' jump
   ratios 0.995-1.039 vs the reflux-era 1.56/0.31); SD7003 L4-xz ==
   L4-3D to every checker digit. The artifacts under discussion arise
   only when GRID-SCALE (unresolvable) content reaches an interface.

4. Consequences observed at Re 4e5 / L6-xz (the snug-nest campaign,
   stopped): C_L(0) = -0.008 offset (fan-locked asymmetric mean state;
   the fan-free Re 1e5 L5-xz run had -0.0017), C_L rms 0.005-0.025
   (20x the fan-free class; f ~ 48 U/c = staircase-TE shedding at
   St_h ~ 0.2 on ~0.004c effective bluntness), Cf contaminated over the
   first ~1/3 chord. Cp remained essentially clean (stagnation +1.00,
   Cp_min -0.41 = XFOIL at aoa 0).

DISENTANGLING TEST (ready to run on resume): the B11/C10 protected-box
designs place the first interface 0.2-0.5c from the body. If the near
fan persists there unchanged -> pure seed+scheme (as R1 says); any
fluctuation bands at the DISTANT interfaces are then cleanly
attributable to interface scattering of whatever reaches them. That is
exactly the "interface-artifact test bed" geometry.

## Campaign state (all committed)

- Solver/tooling infrastructure (commit 6a14132 + earlier):
  - [blocks] refine = x0 x1 y0 y1 z0 z1 [level] — per-level refine
    boxes (16 max), mobygeom mirror;
  - WINDOWED per-level lattices end-to-end (solver lidOf, mask file
    format with per-level window attrs, mobygeom builder/masks): deep
    quadtrees (finest lattice 7.2e9 blocks) now build in ~1 GB; legacy
    full-raster files unchanged (suite bit-exact incl. les_ibm_ref);
  - block-table tile generation: far-field bbox shortcut (coef rows
    exact zeros; dwall exact igl batched 400 leaves/call) + --jobs
    worker pool. SD7003 regeneration BITWISE identical in 36 s.
  - [case.airfoil] span = y; make_airfoil_stl --span y --n (nose facet
    sagitta: n=2880 -> 0.25 fine cells at c/6144);
  - surface_cp_cf.py (constrained-LS Cf using only fluid points + the
    exact no-slip BC; LS wall-extrapolated Cp — validated: stagnation
    Cp = +1.00), plot_cp_cf.py, plot_polars.py, cv_forces.py
    (control-volume forces, momentum balance with full effective-stress
    terms, nested-box independence check — REQUIRED because the interior
    is removed), plan_decomposition.py, show_decomposition.py,
    make_ic_tests.py.

- B11 case (built, verified, NOT run): 128c x 96c, nose (50,48), 11
  levels, y+(D/2) ~ 2.06, 31088 leaves = 15.9 M cells (~6.7 GB);
  ibm_coeff_b11.h5 on disk (untracked, regenerate via setup_b11.sh);
  1-step 4-rank dry run PASSED (row cross-check + windowed lidOf).
  ~16-20 h/angle on the RTX 5090 (dt 5e-5, t_final 15).

- C10 proposal (template c10_base.ini committed; table NOT generated —
  generation was stopped at pause): 10 levels, y+(D/2) ~ 4.1 LE / ~2.7
  midchord, 0.25c protected box at L6 (c/192) + doubling wake-skewed
  margins; estimated ~11.5k blocks ~ 6 M cells, dt 1e-4 ->
  ~3-4 h/angle. Regenerate with:
      cd tutorials/naca && bash setup_b11.sh   # grid + STL if absent
      then the mobygeom block-table call from c10_base.ini's boxes
      (see gen_c10.log for the exact command line).

- IC findings (2026-07-16): homogeneous (1,0,0) impulsive start is the
  MILDEST of the three tested; zero-IC explodes (inviscid inlet front at
  cell-Re 146); body-masked IC (zero in body + 2-cell layer) rings
  WORSE (the shell discontinuity is a stronger seed than the interior
  impulse).

- Physics settings for the rerun (agreed): fully-turbulent SST with
  SUSTAINED ambient (tu = 5, nut_ratio = 1000 -> tu ~ 1.5-2 % arriving
  at the body; record it and match the XFOIL/OpenFOAM transition
  setting); niter = 18 (CV-pressure cleanliness); buried interior
  REMOVED; forces ONLY from cv_forces.py on snapshots; aoa sweep
  -2..5.

## Resume checklist

1. Regenerate/decide the decomposition (B11 vs C10 or the user's new
   geometric preparation), show it, dry-run 1 step.
2. Single aoa0 verification: check (a) fan strips vs the R1 tables at
   the new interface distances, (b) interface-locked striping present/
   absent at the distant interfaces (the artifact test), (c) C_L(0)
   offset, (d) CV forces box-independence, (e) ambient tu arriving at
   the body, (f) [ibm] band_filter decision (validated option if the
   seed-fan needs damping — R1: 3-5x at ~2 % cost).
3. Then the sweep (one host per angle) + polars + Cp/Cf vs XFOIL and
   OpenFOAM.
