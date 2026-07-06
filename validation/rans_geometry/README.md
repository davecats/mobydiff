# RANS T1 geometry gates — wall distance (dwall/y_eff) + IBM wall cells

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

Workflow:

```bash
./setup.sh                                   # block-table coeff files with dwall_blocks
mpirun -n 1 ../../build_cpu/main flat_l1.ini
mpirun -n 1 ../../build_cpu/main flat_refine.ini
mpirun -n 1 ../../build_cpu/main wavy.ini
python3 check_rans_geometry.py --mode flat  flat_l1_ransgeom.h5
python3 check_rans_geometry.py --mode flat  flat_refine_ransgeom.h5
python3 check_rans_geometry.py --mode wavy  wavy_ransgeom.h5
```

The mobygeom `block-table` files are generated from the les_ibm STLs/grid
(`../channel_interface/les_ibm/`); `setup.sh` also verifies that the
regenerated levels-2 file carries byte-identical `coef_blocks`/masks/blocks
to the committed `les_ibm/ibm_coeff_blocks.h5` (the dwall_blocks dataset is
a pure addition).
