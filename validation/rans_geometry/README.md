# RANS T1/T1b geometry gates — wall distance (dwall/y_eff) + IBM wall cells

Gates for IDDES phase T1 (docs/next_session_iddes.md): the `[rans]` section
builds the SST geometry state at init (rans.f90) and `dump_geometry = true`
writes `<prefix>_ransgeom.h5` (blocks table, interior dwall / yeff / wallcell,
per-block cell-centre coordinates). Nothing consumes the state yet, so these
runs advance the flow exactly like their `[rans]`-less twins.

Cases (all cold-start, 1 step — only the init matters):

- `flat_l1.ini`     file-based IBM plane-wall channel (the les_ibm walls at
                    y = 0.259375 / 2.259375 on the 64x80x64 grid), single
                    level, nb = 8. Reference: exact point-to-box distance to
                    the two wall slabs + the y domain walls + the half-cell
                    floor.
- `flat_refine.ini` same geometry with `refine_body` (levels 2): per-leaf
                    dwall at each leaf's level.
- `wavy.ini`        the analytic wavy-wall IBM (single bottom wall,
                    amp 2.5e-2, offset 1e-2) on a 64x64x16 section.
                    Reference: independent high-precision minimization of the
                    distance to the sine curve (scipy) + domain walls + floor.
- `wavy_refine.ini` T1b gate c: the same analytic case with `refine_body`
                    (levels 2), per-level dwall on the leaves' own centres.

Since phase T1b the analytic path is GEOMETRY-AGNOSTIC (walldist.f90):
dwall is computed from the isInBody indicator alone (finest-level surface
point cloud + kd-tree nearest + local polish to `[rans] dwall_tol`,
default 1e-10), so the wavy runs exercise the generic machinery — the
wavy-specific `body_surface_distance` minimization was deleted and the
scipy minimization in `check_rans_geometry.py` is the surviving
specialized reference. T1b results (2026-07-07):

- wavy.ini generic vs scipy: max|dwall-ref| = 2.28e-11 (default tol);
  sweeping dwall_tol 1e-2/1e-4/1e-6/1e-8 gives 9.8e-4/2.5e-5/2.6e-7/2.6e-9
  — monotone convergence with polish depth.
- wavy_refine.ini: level 0 = 2.35e-11, level 1 = 2.44e-11.
- second geometry through the SAME machinery: `mpirun -n 1
  build_cpu/walldist_test` (src/test_walldist.f90) — a sphere straddling
  the periodic x boundary vs the closed form |r - R|; errors track the
  tolerance (1.5e-4/1.6e-7/1.6e-10 at tol 1e-4/1e-7/1e-10).

Workflow:

```bash
./setup.sh                                   # block-table coeff files with dwall_blocks
mpirun -n 1 ../../build_cpu/main flat_l1.ini
mpirun -n 1 ../../build_cpu/main flat_refine.ini
mpirun -n 1 ../../build_cpu/main wavy.ini
mpirun -n 1 ../../build_cpu/main wavy_refine.ini
mpirun -n 1 ../../build_cpu/walldist_test    # T1b gate b (sphere)
python3 check_rans_geometry.py --mode flat  flat_l1_ransgeom.h5
python3 check_rans_geometry.py --mode flat  flat_refine_ransgeom.h5
python3 check_rans_geometry.py --mode wavy  wavy_ransgeom.h5
python3 check_rans_geometry.py --mode wavy  wavy_refine_ransgeom.h5
```

The mobygeom `block-table` files are generated from the les_ibm STLs/grid
(`../channel_interface/les_ibm/`); `setup.sh` also verifies that the
regenerated levels-2 file carries byte-identical `coef_blocks`/masks/blocks
to the committed `les_ibm/ibm_coeff_blocks.h5` (the dwall_blocks dataset is
a pure addition).
