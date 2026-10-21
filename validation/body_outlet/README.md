# Body at an outlet face

Gates for an immersed body that crosses an inlet and an outlet plane
(`docs/next_session_body_at_outlet.md`, increments B1 and B2). Poiseuille
flow between two IMMERSED plane walls (`make_slabs.py`: solid below
`y = 0.259375` and above `1.259375`, both 0.3 dy off a cell face, the gap
exactly 1.0, both slabs padded beyond every domain face) in an
inflow/outflow channel 4 x 1.5 x 0.125 on 64 x 48 x 4 cells, `nb = 4` with
`remove_solid` (32 of 192 blocks buried, 16 of them along each slab, some at
the outlet), `re = 10`, uniform inlet `u = 1` at `x_min`, Dirichlet-pressure
outlet at `x_max`.

    SOLVER=/abs/build_cpu/moby_solve [RANKS=4] ./run_gates.sh [profile|mirror|permute|limits|ranks|restart|twin|pocket|all]

Every solve is prepared beside the solver by `tools/prepare_if_missing.sh`.
`SHORT=1` cuts the runs to 2000 steps (the exactness groups do not need the
developed state; `profile`/`mirror`/`permute` do).

What each group measures, and what it read (2026-10-05, `build_cpu`,
numbers to be found in the section below):

- **profile** -- `pois_body.ini` against its periodic twin
  `pois_body_periodic.ini` (same walls, `forcing_x = 12 nu U/h^2 = 1.2`) at
  `t = 40` (viscous time 10): the velocity ON the outlet face equals the
  twin's developed profile; the pressure along the gap centre is linear with
  the twin's gradient and the last fluid cell holds `G dx/2` above the
  outlet value; every solid-centred velocity (interior AND the 32 solid
  faces on the outlet plane, read from `un_xmax`) is at round-off; the flux
  through the outlet equals the inlet flux; the inlet rows inside the slabs
  are reported zeroed by the solver (B2: `boundary: non-zero Dirichlet
  velocity data zeroed inside the body: 64`).
- **mirror** -- `pois_body_low.ini`, the outlet on the LOW face (inlet at
  `x_max` with `u = -1`): the `x -> lx - x` image of the high case.
- **permute** -- `pois_body_y.ini` (stream y, walls normal x) and
  `pois_body_z.ini` (stream z, walls normal y): the x case under the
  coordinate permutation, including the `vn_ymax` / `wn_zmax` planes.
- **limits** -- `pois_body_ends.ini` (the slabs END at `x = 3.96`, between
  the last interior u face and the outlet face: `n` solid, `f` fluid,
  `r = incr_f/incr_n ~ 1e27`) and `pois_body_starts.ini` (the slabs exist
  only for `x >= 3.97`: `f` solid, `n` fluid, `r ~ 1e-27`): stable, solid
  faces at round-off, the mass balance closed.
- **ranks** -- `RANKS` ranks == 1 rank at `max_abs 0`.
- **restart** -- 1000 + 1000 steps == 2000, `max_abs 0` (the outlet face
  is state; a low face sits in `un`, a high one in `un_xmax`).
- **twin** -- a case file prepared WITH the slabs and then every
  coefficient zeroed (the body crossing the outlet with `coef = 0`) equals
  the body-free run at `max_abs 0`: the one-statement form of
  `predict_outlet_faces` is the unpenalized statement when the factors are
  exactly 1, measured rather than argued.

- **pocket** -- `pocket.ini` + `check_pocket.py`: the trough of a rough
  wall ON the outlet plane (`docs/next_session_rough_outlet.md`, section
  4). A laminar boundary layer (Blasius inlet, no trip) over the MacDonald
  egg-carton of the production rough case on its own wall line and dx
  (`cans`, dy_1 0.0105, dx 0.205; 106 x 192 x 22 on 21.7 x 50 x 4.33, one
  spanwise wave, five streamwise), the roughness grown in over x = 5..10
  and running through the outlet, 2500 steps = 50 t.u., `cflmax 1.0`,
  red-black 6. The last column at the spanwise trough is a dimple open to
  the outlet: on the plane the coefficients read 134 / 0.072 / 0.012 in
  rows 1..3, one cell behind solid / solid / 12.8, two cells behind solid.
  With the 2026-10-01 face equation alone the row-3 face read +0.35 at
  t = 6 (physical: 0.02), drifted into inflow, grew at an accelerating rate
  and blew up at **t = 34** -- the laminar form of what killed the
  production rough case and its 4 M-cell cut-down (t = 52). The gate: every
  field finite at t = 50, `max|u|` on the plane in the dimple rows below
  0.5, interior rms divergence below 1e-3. 4 CPU ranks, ~5 minutes.
  MEASURED 2026-10-06 with the total-pressure datum: the dimple face
  −0.077 at t = 50 (rows 1..6: −1e-5, −0.013, **−0.077**, −0.001, −0.001,
  −0.003; the inflow saturates at −0.094 around t = 38 and decays, −0.032
  by t = 150), rms divergence 1.0e-4 at t = 50 (the cold-start transient
  at red-black 6; 7.9e-8 by t = 150), CPU == GPU `max_abs 0` at 1500
  steps (nofma pair, same case file, the datum live from t ≈ 25).

The uniform inlet meets the cut row of each slab at the inlet plane, where
the flux it injects into a half-solid cell is deflected: `max|v| = 0.41` at
`i = 0`, `4e-3` by `x = 1`, `9e-6` in the last eight cells. That is the
entry flow of a uniform profile against a staircase wall, not an outlet
matter, and it is why the domain is four gap widths long.

## Measurements

RE-VALIDATED 2026-10-21, commit 532ef54 (banner), production flags (the
`gpu` group with the nofma pair), workstation CPU 4 ranks + RTX 3060;
`run_gates.sh all` 43 min, `pocket` 4 min, `gpu` + `twin` 2 min: ALL PASS.
profile: bulk-flux ratio 0.999102, `max|u_outlet - u_twin|` 1.347e-3,
developed-ness 2.6e-5, G 1.199866, last-cell p 3.7466e-2 / ratio 0.99921,
solid u 9.4e-28, plane solid faces 1.6e-28, flux imbalance -5.6e-16 -- the
2026-10-05 numbers below to the digit; mirror u 2.9e-15 (3.3e-15 then), p
5.0e-13, outlet face 3.3e-15; permute y 0.0, permute z 2.4e-15 / p 6.0e-13;
limits as recorded; ranks 1 == 4 EXACT (incl. `un_xmax`); restart ==
continuous EXACT; GPU == CPU EXACT at nofma; the zeroed-coefficient twin ==
body-free EXACT; pocket max|u| in the dimple rows 7.69e-2 (bound 0.5),
interior rms divergence 1.012e-4 (bound 1e-3; 1.0e-4 recorded). Settings
changed for the rerun: none.

All numbers 2026-10-05, `build_cpu` (production flags) unless said otherwise,
head of the session's branch; `run_gates.sh all` plus the `gpu` group with the
nofma pair.

**profile** (20000 steps, t = 40): bulk flux through the fluid rows of the
outlet face 1.000000 against the twin's 1.000898 (the inlet sets exactly 1);
`max|u_outlet - u_twin|` **1.347e-3** on a peak of 1.5007 (0.09 %, of which
the bulk ratio accounts for 0.09 %); the outlet face against the same run at
x/lx = 0.75: 2.6e-5 (developed). Pressure along the gap centre: gradient
**1.19987** against the twin's forcing 1.2, last-cell p **3.7466e-2** against
G dx/2 = 3.7496e-2 (**ratio 0.99921**), the last four cells 0.2623, 0.1873,
0.1124, 0.0375 (steps of G dx = 0.075 down to the outlet value 0). Solid
locations: `max|u|` **9.4e-28** (2048 u faces), `|v|` 5.4e-28, `|w|` 0.0;
the 32 solid faces ON the outlet plane **1.6e-28**; `max|p|` in the solid
8.3 (the decoupled Brinkman pressure, bounded; fluid 11.1). Flux in = out to
**-2.2e-16**. The solver reports 32 non-zero inlet data zeroed (the 64
u faces of the slabs' rows at the inlet, half of them in removed blocks).

**mirror**: u **3.3e-15**, v 1.7e-15, w 0.0, p 5.4e-13, outlet face 2.9e-15
against the x_max image; its own solid faces on the x_min plane 1.6e-28,
flux in = out to 0.0.

**permute**: stream along y -- **0.0 in every field and on the `vn_ymax`
plane**; stream along z -- 2.4e-15 (stream), 1.5e-15 (normal), 0.0 (span),
1.1e-13 (p), 2.9e-15 (`wn_zmax`).

**limits**: `ends` (n solid, f fluid): 20000 steps stable, the 64 solid
neighbour faces at **3.4e-28**, the fluid outlet faces behind the body's end
carry an outflow of 0.233, flux imbalance -4.4e-16. `starts` (f solid, n
fluid): stable, the 64 solid faces on the plane at **9.8e-28**, flux
imbalance 0.0 (the inlet is now 1.5 high: flux 0.1875, exiting through the
fluid rows at up to 2.0).

**ranks**: 4 == 1 at **max_abs 0** incl. `un_xmax` (2000 steps).
**restart**: 1000 + 1000 == 2000 at **max_abs 0** incl. `un_xmax`, 4 ranks.
**gpu**: `build_gpu_nofma` == `build_cpu_nofma` at **max_abs 0** (2000
steps, same case file); isolated: the body-free outlet case and the
periodic body case are each CPU == GPU at 0.0 too. (A production CPU run
against a nofma GPU run reads 2e-15 / p 6e-14: the FMA class, not a
defect -- the group now runs its own CPU leg.)
**twin**: zeroed coefficients == body-free at **max_abs 0** incl. `un_xmax`.

**The rough-wall boundary layer** (not part of `run_gates.sh`; inputs and
results in `~/outlet_runs/body_outlet/{rough,smooth}_tbl`, run on istmcetus):
a production-shaped cut-down of `overheadTest/horeka/configs/rough_jacobi.ini`
-- the production y line (176 cells, dy_1 = 0.0103, so the egg-carton
roughness +-0.05 spans rows 1..6 at the crests and rows 0..1 are solid on
BOTH x planes), 320 x 176 x 24 on 60 x 100 x 8, wavelengths 15 / 8, nb 64 44
24 -- and its smooth twin on the same grid: 200000 steps from rest, then
100000 with statistics (t = 200..300), 13 ms/step on one A6000. Both run
(the head of 2026-10-02 refused the rough one at init); the solver reports 34
non-zero inlet data zeroed. The momentum-thickness growth dtheta/dx
(`assets/postpro/momentum_integral.py`'s loader; c_f from the first cell is
meaningless on the rough wall): the SMOOTH twin holds **1.20e-3 (x 55..57)
-> 1.19 -> 1.17 -> 1.16 -> 1.15e-3** into the last half unit (the old outlet
collapsed it to 5 % of its value in the last cell on the production case);
the ROUGH one oscillates with the roughness phase and in the last unit reads
1.66e-3 against **1.31e-3 / 1.42e-3 at the same phase one and two
wavelengths upstream** (+26 % / +17 %), where the smooth twin reads +12 % /
+4 % in the same band -- the layer's own slow streamwise modulation, present
on both, makes a finer separation impossible on this case; theta grows
uniformly through the last six cells (0.45503 -> 0.45691, no kink). What
this does NOT show: both layers are still LAMINAR at t = 300 (H = 2.58 /
2.62; the trip at x0 = 10 with amp 0.03 does not transition them within 60
units), so the turbulent balance of the handout's gate is not measured here;
it needs the production domain on the cluster.
