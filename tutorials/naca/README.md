# NACA 0012 polar sweep at Re_c = 4e5 (fully-turbulent SST, L6-xz IBM)

Production polar for comparison against XFOIL and OpenFOAM RANS: aoa =
-2, -1, 0, 1, 2, 3, 4, 5 at Re = 4e5, on the R2D 2D-refinement stack
(docs/next_session_refine2d.md): span along y (`[case.airfoil] span =
y`: chord x, LIFT z), `[blocks] refine_dims = xz`, `refine_levels = 6`.

Physics/resolution choices (see naca_base.ini comments):

- **L6 surface Delta = 3.66e-4 c -> y+_1 ~ 2-3** at the turbulent
  midchord (cf ~ 4e-3 at Re 4e5): resolved-wall SST, credible Cf. The
  same grid in 3D octree mode would be ~121 M cells (it OOM'd a 49 GB
  A6000 at Re 1e5); the xz quadtree makes it ~15 k blocks.
- **Fully turbulent** (`[rans] transition` off): matches XFOIL run
  fully-turbulent (or with forced transition at the LE) and avoids the
  gamma-Re_thetat fine-grid unreliability documented in
  validation/sd7003 ("R2D-3 follow-up").
- **Keep-buried coefficients**: load-bearing for the penalization
  C_L/C_D (validation/naca0012 README).
- **12c box, Dirichlet far field**: costs ~10-15 % of the 2pi lift
  slope (blockage); compare slopes with that in mind, or compare
  against XFOIL shifted by the measured ratio.

## Run

```bash
PY=$HOME/ibmc/bin/python ./setup.sh   # STL + grid.h5 + L7-xz block table
# one host per queue, e.g. (cetus NEEDS the login shell -l):
ssh istmcorax 'setsid bash -l tutorials/naca/run_sweep.sh corax -2 -1 0 1 &'
ssh istmcetus 'setsid bash -l tutorials/naca/run_sweep.sh cetus 2 3 &'
setsid bash run_sweep.sh local 4 5 &
```

~50 min/angle on the RTX 5090 (t_final = 10 at dt = 1e-4), ~2x that on
the A6000 / RTX 3060.

## Polars

```bash
python3 plot_polars.py --xfoil xfoil_polar.txt --openfoam of_polar.csv
# -> polars.png (C_L-alpha + C_L-C_D), polar_mobydiff.dat
```

XFOIL reference (fully turbulent, matching this setup):
`xfoil -> naca 0012 -> oper -> visc 4e5 -> vpar (xtr 0.01 0.01) ->
pacc -> aseq -2 5 1` and save the polar dump.

## Surface Cp / Cf

```bash
python3 surface_cp_cf.py naca_aoa4_100000.h5 --plot cpcf_aoa4.png
python3 plot_cp_cf.py cpcf_naca_aoa4_100000.npz \
    --xfoil-cp cpx_a4.txt --of-cp of_cp_a4.csv --of-cf of_cf_a4.csv
```

- **Cf** is a least-squares wall gradient over NEAR-WALL FLUID CELLS
  ONLY, constrained through the origin by the EXACT immersed no-slip
  condition u_t(0) = 0 (no solid/forced value enters): span-averaged
  tangential velocities of cells within 2.5 fine cells of the analytic
  section are fitted as u_t = g d (weights 1/d), then re-restricted to
  d+ <= 5 with u_tau = sqrt(nu g) so only the viscous sublayer feeds
  the fit. Cf = 2 nu g, signed TE-ward (negative = reversed flow).
- **Cp** is a least-squares linear WALL EXTRAPOLATION of the same
  cells' pressure to d = 0 (keeps dp/dn ~ 0 where the BL approximation
  holds, still captures the finite normal gradient at the curved LE);
  p_inf from a level-0 far-upstream box. A pure zero-gradient
  (nearest-cell) Cp is the degenerate < 3-cell fallback.
- Outputs: `cpcf_*.npz` + XFOIL-comparable `*_cp.dat` / `*_cf.dat`
  (x/c value; upper block then lower).

External reference formats accepted: XFOIL CPWR (x [y] Cp) and CF
dumps; OpenFOAM as any '#'-commented CSV/whitespace table with x/c in
column 1 and the value in the last column.
