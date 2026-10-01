# blasius — laminar Blasius boundary layer (the laminar stretched-line gate)

Quasi-2D laminar zero-pressure-gradient flat-plate boundary layer,
Re_theta,in = 100 (lengths in inlet momentum thicknesses, velocities in
U_inf). It is the B0 precursor of `tutorials/turbulentBoundaryLayer/`, kept
here because it is a GATE: a one-sided natural (stretched) y line, the
Blasius u + v inlet, a Dirichlet-p outlet at `x_max` AND at the top, and a
run that starts FROM the analytic similarity field and must simply hold it.
It is the laminar stretched-line case the F1 flux-form viscous stencil was
validated on (numerics review step 8).

Domain 400 x 100 x 4, grid 384 x 160 x 4 (dy_wall ~ 0.18 theta, 3.8 % max
growth), t = 2000 at the dt cap 0.5 (4000 steps). One rank on the local GPU
takes ~41 s; 4 CPU ranks are fine too (`RANKS=4 SOLVER=../../build_cpu/moby_solve`).

## Run

```bash
./run_blasius.sh          # SOLVER (default ../../build_gpu/moby_solve), RANKS (default 1)
```

The driver does, in order: derive `template.ini` from `blasius2d.ini`
(nsteps 1, no restart, prefix `template`); prepare the ONE case file both
inis name (`[case] file = blasius.case.h5`); mint `template_1.h5`;
`make_blasius_ic.py` fills it with the Blasius field (`IC_blasius.h5`,
staggered u/v rows, same virtual origin as the inlet table); run
`blasius2d.ini`; `compare_blasius.py` measures theta, H and the u/v
profiles at x/lx = 0.15..0.7 against its own RK4 + shooting Blasius solve.
Gates: theta and H within 2 %, u within 1 % of U_e, v within 15 % of the
local entrainment scale inside the layer. The last 15 % of the domain is
the outlet zone and is not sampled; the p = 0 top under-drives the
entrainment aloft (`v_top`, informational).

**Landmine, the reason the template is derived and never hand-edited:** a
restart file carries its boundary rows (`bc_type`/`bc_value`) and the ini
overrides only the rows it sets explicitly. A template minted from a sibling
variant (a Dirichlet top, say) silently transplants that variant's faces —
the outlet-top run restarted that way diverged by t ~ 14 on every binary
(`docs/next_session_after_step8.md`).

## Recorded results (2026-09-30, branch head after F1, GPU build, 1 rank)

`worst` over the four stations, t = 2000:

| projection | theta | H | du/Ue | dv/v_edge | gate |
|---|---|---|---|---|---|
| damped Jacobi, niter 12 (**the shipped ini**) | 1.33 % | 0.45 % | 2.40e-3 | 0.124 | PASS |
| damped Jacobi, niter 6 | 1.74 % | 0.72 % | 2.23e-3 | 0.096 | PASS |
| red-black SOR 1.5, niter 6 | 1.55 % | 0.58 % | 2.41e-3 | 0.115 | PASS |
| Chebyshev-Jacobi, niter 12 | 1.55 % | 0.58 % | 2.41e-3 | 0.115 | PASS |

Red-black 6 and Chebyshev 12 agree to every printed digit: that is the
converged-projection answer, and the ODE comparison shows the ~1.5 %
theta / ~0.6 % H discretisation error of this grid at the x/lx = 0.7
station (the first three stations sit within 0.4 % / 0.2 %). Jacobi 12
reads closer to Blasius (1.33 / 0.45) by an under-converged projection's
luck, not by accuracy; it stays the shipped configuration because it is the
one B0 documented. B0 itself (2026-07-19, pre-review code, niter 6) recorded
theta 1.13 %, H 0.34 %, du/Ue 2.4e-3, dv 0.066.

**Chebyshev is stable here.** The B0-era finding that Chebyshev at niter 6
pumps a 2-dx pressure mode (e-fold ~36 t.u.) is not reproducible on main
(review F3, retracted 2026-09-27): the Chebyshev row above ran to t = 2000
at dt 0.5 with `|p|` bounded and `dt` at its cap throughout. The stale
warning was removed from the ini on 2026-09-30.

Figure: `blasius.png` (u/Ue and v profiles in similarity form, theta growth)
from the shipped-ini run.

## 2026-10-01: the outflow-predictor prototype on this gate

`../cylinder/README.md` (last section) traces the outlet pressure mode to the
predictor's reset of the outlet face; a prototype (local branch
`proto/outflow-incremental`, not on main) gives that face a predictor of its
own. This case has two outlets (x_max and the top), so it was run with both
binaries, the shipped ini (damped Jacobi, niter 12), GPU, 1 rank:

| | main (`~/step9b_ref_binaries`) | prototype |
|---|---|---|
| theta error at x/lx 0.15 / 0.30 / 0.50 / 0.70 | 0.24 / 0.39 / 0.08 / −1.33 % | 0.26 / 0.74 / 1.27 / 1.65 % |
| H error | −0.12 / −0.19 / −0.09 / 0.45 % | −0.12 / −0.35 / −0.63 / −0.82 % |
| du/Ue, worst | 2.40e-3 | 1.48e-3 |
| dv/v_edge, worst | 0.124 | **0.013** |
| `v_top` (relative entrainment error aloft) | +0.013 +0.002 −0.041 −0.202 | +0.025 +0.030 +0.028 +0.033 |
| gate | PASS | PASS |

The wall-normal velocity is ten times closer to Blasius and the top
entrainment deficit ("the p = 0 top under-drives the entrainment aloft",
above) is gone: it was the outlet treatment, not the p = 0 condition. The
momentum thickness now grows 0.3 → 1.7 % above Blasius along the plate where
main reads within 0.4 % up to x/lx = 0.5 and −1.3 % at 0.7; both are inside
the 2 % gate, and which of the two is the grid's own error is not settled by
one run (main's converged-projection value is 1.55 %).
