# Minimal-span channel: cheap full-step reproduction of the 2:1 interface blow-up

> **HISTORICAL.** This case and the numbers below date from the original
> diagnosis, run with the old `sor = 1.5 / niter = 3`. That `sor = 1.5` ALONE
> diverges (simple Jacobi needs `sor <= 0.8`), so it confounded the interface
> study. The inis have since been updated to the current defaults
> (`sor = 0.8, niter = 6, accel = chebyshev`) and the interface treatment has been
> reworked, so the runs below no longer reproduce as written -- re-measure if used.

This is the **Step-1 confirmation gate** for the 2:1 interface projection
instability (`docs/interface_review.md` §vi, `docs/next_step_prompt.md`). It
reproduces the `validation/channel_interface` turbulent blow-up at ~1/200th the
cost by keeping the *exact* Reynolds number, wall-normal/streamwise grid and
velocity scale of the Re_tau=180 channel but shrinking the spanwise extent to a
minimal-flow-unit slab.

Why the cheap laminar reproductions failed first: the zero-coarse-average
interface mode is high-k and viscously damped, so on a low-velocity laminar base
it always decays. It only blows up when the **channel velocity (~18)** makes the
nonlinear advection of the amplified interface-normal velocity outrun viscosity.
The minimal-span channel keeps that velocity; the laminar Poiseuille (Umax=1) did
not, which is why it stayed stable at every Re.

## The two cases (run on GPU; CPU is ~80x slower)

    module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
    mpirun -n 1 build_gpu_nofma/main input_gpu.ini        # refined: 2:1 interface
    mpirun -n 1 build_gpu_nofma/main input_gpu_noref.ini  # control: uniform, no interface

Both: 128 x 64 x 16, Re_tau=180, natural-y grid, forcing_x=1, the channel
initialiser's mean profile + transition disturbance. `input_gpu.ini` refines
both near-wall bands to level 1 (flat y-interfaces at y+=112, as in the
validation); `input_gpu_noref.ini` is identical with `[blocks] refine` removed.

Inspect with (max|u|, max|v| vs step):

    python3 tools/check_interface_shear_mode.py <rundir> --prefix channel_field

## Result (current code, the defect)

| step | refined (2:1 interface) | control (uniform) |
|------|-------------------------|-------------------|
| 250  | max\|u\|=18.4, max\|v\|=0.38 | 18.4, 0.14 |
| 500  | 21.5, 3.36              | 18.4, 0.11 |
| 750  | **245, 441**            | 18.5, 0.09 |
| 1250 | **5e5**                 | 18.6, 0.06 |
| 3500 | (10^80, dead)           | **19.2, 0.015** (bounded) |

Identical inputs; the *only* difference is the 2:1 interface. The refined case
blows up at the interface around step 750 (t~0.23); the uniform control stays
bounded. This isolates the blow-up to the interface, full-step, at the real
channel velocity.

## Step-2 acceptance

The projection-coupling fix must make `input_gpu.ini` stay **bounded** (turbulent
O(20), like the control) instead of diverging — while every exact gate
(uniform-flow 0.0, channel nb=4 bit-exact, interface_decay contraction, mass
round-off, MOBY_HALO_AUDIT clean) still passes.
