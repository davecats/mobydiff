# A3 INCREMENT 1 — RANS scalar inlet values (gate b)

The T4 gap closed: rans.f90's per-scalar ghost-mode tables now key on the
DECLARED patch type — `SCALAR_BC_VALUE` at `PATCH_INLET` faces (freestream
values computed once at transport init: k_inf = 1.5 (tu/100 U_inf)^2 with
U_inf the face's Dirichlet velocity magnitude, omega_inf = k_inf/(nut_ratio
nu), gamma_inf = 1, Re_thetat~_inf = the T4 lambda = 0 tu correlation) and
`SCALAR_BC_COPY` at `PATCH_OUTLET` — pure functions of the patch type, no
re-inference from the velocity BC rows. nut wall-ghost handling untouched.

`inlet_channel.ini`: SST plane channel (Re = 100, 64x32x8), Poiseuille
inlet at x_min (parabola, peak U_inf = 1), outlet at x_max, declared walls
at y, tu = 5 %, nut_ratio = 10 => k_inf = 3.75e-3, omega_inf = 0.0675.

Gates (run_gates.sh):

- report_patch_types shows x_min inlet / x_max outlet / y walls, all
  "(declared)";
- after 500 steps the first interior column at mid-channel holds the
  freestream: k and omega within 10 % of k_inf/omega_inf (the ghost pins
  the face value by the midpoint identity; one cell of convective
  adjustment);
- 1 == 4 ranks EXACT (compare_fields --tolerance 0 on un vn wn pn k omega
  nut).

Gate (a) — the standard 7-case suite (min_channel 4-rank CPU + GPU,
les_ibm ± refine_body, Beltrami y-slab, turb180, wf180_y30, lam30t)
bit-exact vs the pre-change binary (nofma, max_abs 0 incl. all RANS
scalars, CPU AND GPU) — runs from the session gate driver; the new modes
are inlet-face-gated and dormant in every channel.

## Results (2026-07-13, all PASS)

- gate (a): all 7 cases BIT-EXACT (max_abs 0) vs the dd37937 binaries,
  CPU AND GPU, incl. k/omega/nut (+gamma/rethetat on lam30t).
- gate (b): report_patch_types prints inlet/outlet/wall "(declared)";
  first-column mid-channel k = k_inf to 0.21 %, omega = omega_inf to
  2.3 % after 500 steps; 1 == 4 ranks EXACT (max_abs 0 on un vn wn pn
  k omega nut). CPU vs GPU on the ACTIVE inlet path <= 1.8e-10 over 500
  steps (omega; velocities 1.7e-12) — the established A0 in/outflow
  CPU/GPU ulp class (the freestream gates gate GPU at 1e-12/200 steps),
  not a scalar-inlet artifact: the dormant-path suite is exactly 0.
