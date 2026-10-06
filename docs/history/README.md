# History: closed handouts and superseded design notes

Moved here on 2026-10-06 (second cleaning pass, `docs/cleanup_review_2026-10-06.md`
section 6.3) with every citation redirected (`git mv`, so `git log --follow`
works). Each file is the record of a finished piece of work: the plan, every
deviation, and the gate numbers. Nothing here is a current instruction; the
live handouts stay in `docs/`. One line each, in rough chronological order.

| file | what it recorded | closed |
|---|---|---|
| `momentum_interface_handout.md` | the 2:1 momentum-interface investigation (reflux band) | superseded 2026-06 |
| `jacobi_interface_handout.md` | the damped-Jacobi / Chebyshev projection across the 2:1 interface | 2026-06 |
| `interface_review.md` | review of the first interface treatment | superseded |
| `corner_reconstruction_strategy.md` | cubic edge/corner halo reconstruction (code deleted 2026-07-01) | superseded |
| `next_session_les.md`, `next_session_edges_les.md`, `next_session_les_ibm.md` | LES across refinement, edge/corner validation, LES-IBM coupling | 2026-06-30 (item 3 of edges_les, a less-dissipative transfer, was never started: an accepted residual) |
| `next_session_interface_normal.md` | the coarse-owns-face v' excess, accepted as the const-1/2 price | 2026-06-30 |
| `next_session_bodyforce.md` | the spatially varying body force (`[force]`) | 2026-07-01 |
| `next_session_iddes.md` | k-omega SST, wall functions, transition, IDDES (T0-T5) | 2026-07-10 |
| `next_session_airfoil.md`, `next_session_cv_forces.md`, `next_session_naca_rerun.md`, `next_session_naca_re4e5.md` | the airfoil case, inlets/outlets, control-volume forces, the NACA campaigns | 2026-07/08 |
| `next_session_refine2d.md` | quadtree (`refine_dims = xz`) refinement | 2026-07-15 |
| `next_session_prepare_everything.md`, `next_session_after_step7.md`, `next_session_after_step8.md` | prepare does all preprocessing (step 7) and its loose ends; F1 | 2026-09-30 / 10-01 |
| `next_session_boundary_layer.md` | the laminar Blasius precursor of the boundary-layer tutorial | superseded by `tutorials/turbulentBoundaryLayer/` |
| `next_session_scalar.md`, `next_session_verification.md` | passive scalars S0-S5a and the host/device staleness audit | 2026-08-07 |
| `next_session_conjugate.md`, `next_session_pipe_cht.md` | conjugate heat transfer C1-C3 and the Neuhauser pipe | 2026-08-28 / 09-28 |
| `next_session_block_overhead.md`, `next_session_2to1_penalty.md`, `next_session_apply_registers.md`, `next_session_divergence_halos.md` | the performance campaign (block tax, the 2:1 penalty, registers, divergence halos); results in `tutorials/turbulentBoundaryLayer/overheadTest/results_*.md` | 2026-09 |
| `next_session_merge_to_main.md`, `next_session_port_finish.md`, `next_session_skew_convection.md` | branch consolidation, the jacobi-interface feature port, the skew lockdown | 2026-09-25/26 (port_finish section 3 still lists the 14 airfoil inis without a `cv_box`, an undecided follow-up) |
| `next_session_numerics_steps.md` | execution of the 2026-09-26 review, steps 0-6 | 2026-09-27 |
