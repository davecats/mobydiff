# NACA 0012 Re 4e5 campaign — PAUSED 2026-07-16 (fan/interface observations + state)

Paused by user decision pending a simplification/improvement of the
geometric preparation. Everything below is committed; no simulations or
generations are running on any host.

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
