# Next session: code simplification (the review's section 9, continued)

STATUS: handout written 2026-10-20 on head `68c6360`. Nothing below is
started. It continues `docs/cleanup_review_2026-10-06.md` section 9 (the
ranked findings) after the items DONE on 2026-10-08: the config hardening
(F11/F12), the body-force support restriction, the apply k1+k2 and
`cheb_combine` fusions, the boundary face-index helpers, the shared
projection setup/tail, the log-law constants in `init.f90`. Every item here
is a performance-neutral simplification unless marked PERF; every one is
gated the same way (section 6). The user's standing order: simplification
first, then the one-by-one numerics re-validation of `validation/README.md`
(a separate handout when this one is done).

Read first: the review's sections 2 (duplication), 3 (simplifications, with
the "moves nothing / moves an expression" verdict per item), 4 (risks) and
8 (performance hypotheses). Line numbers below are of `68c6360`.

## 1. Order of work

Ranked by (what a future edit would trip over) / (risk). Items 1-4 are a
day each; 5-8 an afternoon each; 9-11 an hour each.

### 1. `scalar_transport_kernel` split by mode  (review §3 item 13)

`scalar.f90:2769-3373`, ~600 lines, ~100 private variables; the per-run
flags `useIbm / anyConj / useNut / nondiv / wallfn` and the per-scalar flags
`adiab / conjug / tang / banded / prtModel / srcType` are re-evaluated per
CELL inside the collapsed loop. Target: a PLAIN kernel (dirichlet /
adiabatic; nut, wall function, Kays; ~57 privates) and a CONJUGATE kernel
(~70; ~50 if the shipped-disabled C2 tangential block becomes a third
optional kernel), each looping over its own scalar list built once in
`validate_conjugate_config` (`sc%plainList` / `sc%conjList`).

What moves VERBATIM (bit-exact): the face velocities and `divu`
(`:2886-2895`), `cl*` (`:2903-2908`), `sol*` (`:2913-2920`), `nt*`
(`:2954-2963`), `m*` (`:2959-2964`), the masked and unmasked `conv`
(`:2972-2986`), `dm`, the `useNut.and.wallfn` / `useNut` ladders
(`:3253-3315`), `fw/fe/diff` (`:3317-3328`), `srcVal` (`:3335-3343`),
`rhs = -conv + diff + srcVal`, `ss` (`:3368`), the penalization
(`:3380-3387`). Lines that change: the `is` loop header, `if (adiab .or.
conjug)` -> `if (adiab)` in the plain kernel, `cm* = m*` (a logical
identity there), `divuse = divu`, the conjugate `rhs` branch deleted from
the plain kernel. Moving the C2 tangential `diff + corr` out of the fused
`rhs` is NOT bit-exact for `tangential_correction = true` runs (only for
the default-off ones) -- keep it inside the conjugate kernel unless a
measurement says otherwise.

GATE FIRST: there is no case with one plain AND one conjugate scalar.
Write it (a `[scalar.1] ibm_wall = dirichlet` + `[scalar.2] ibm_wall =
conjugate` twin of `validation/conjugate/slab.ini`'s geometry) and add it to
`run_bitexact_s3.sh`'s list before touching the kernel; then the 9-case
suite + `validation/conjugate/run_gates_c{1,2,3}.sh` at nofma, CPU and GPU.
Then ncu on the two kernels (registers, occupancy) against the single one --
the `jacobi_apply` lesson: the dormant branches cost occupancy; this is the
lever, and it is unmeasured.

The same cut applies to `rans_transport_kernel` (`rans.f90:1606-2129`, 103
privates, the transition block alone 190 lines): a gamma/Re_theta kernel
FIRST writing `gams/rets` and a per-cell `gameff` scratch, the k/omega
kernel reading `dkfac = min(max(gameff,0.1),1)` and `pk = gameff*pk` from it
-- the pure-RANS branch stays verbatim (the T2-T5 bit-exactness argument).
Second, after the scalar one.

### 2. One face-diffusivity helper  (review §3 item 14)

The ladder `if (conjug) [banded|plain] + nut-not-at-cut / else if (useNut
.and. wallfn) / else if (useNut)` is written three times --
`scalar.f90:3035-3315` (transport), `scalar_stats.f90:414-446` (sample),
`:833-860` (heat) -- and the convective-mask ladder twice (`scalar.f90:
2951-2969`, `scalar_stats.f90:402-410`). S4's design claim "the statistics
use the TRANSPORT kernel's own face diffusivity" is an identity by COPY
today; one `declare target` `scalar_face_diffusivity(...)` makes it one by
construction. The C2 tangential stencil exists three ways (six direction
blocks `scalar.f90:3108-3230`; `indicator_pass` `:2431-2544` in offset form;
`face_corr_rate`): the kernel blocks take the `indicator_pass` form
verbatim. ~250 lines less. Identity is per-cell operand order: nofma gate,
the 9-case suite + conjugate C1-C3 + the S4 heat gates
(`validation/scalar/run_gates_s4.sh`).

### 3. `comm.f90`: one enumeration  (review §3 item 7)

`init_block_exchange` (`comm.f90:547-876`, 330 lines): the local, receive
and send enumerations are the same ~40-line body three times (resolve ->
round filter -> piece loop -> `candidate_boxes` -> fill), and the receive
loop re-runs `resolve_neighbors` over ALL local blocks once per peer.
Enumerate candidates ONCE per (block, dir) for my blocks and once per peer
block, bucket by owner. The canonical wire order is preserved if the
bucketing keeps (peer, round, block, dir, cand, piece) order -- that order
is what both ends derive independently, so the gate is the load-bearing one:
`min_channel` 1 == 2 == 3 == 4 ranks EXACT (four peer topologies, the gate
the divergence-list work was validated with), `validation/redblack_interface`
1 == 4, the 7-case suite at 4 ranks. A copy is a copy: production flags
suffice, run nofma anyway. Also here: `interface_normal_dim` /
`div_entry_var` / `iface_restrict_normal` (`:1226-1320`) are one "pure +axis
face: which dim" loop with three op tests -> one `plus_face_dim(off)`.

### 4. The lattice probe and the node-line subdivision, once each  (review §2)

The same-level -> parent -> children occupancy probe exists four times:
`face_kind` (`blocks.f90:1167-1208`), `occupied_any_level` (`:1213-1233`),
`partition_surface` (`:1525-1586`), `resolve_neighbors` (`comm.f90:948-
1011`); `face_kind` is literally `resolve_neighbors`' op (COPY -> OPEN,
PROLONG -> COARSE, RESTRICT -> FINE, none -> CLOSED). Move the probe to
blocks.f90 beside `occupied_any_level` and derive the others from it.
Midpoint subdivision of a node line exists three times: `subdivide_node_line`
(`blocks.f90`), the `subdivided` branch of `build_node_line` (`init.f90`),
and `ibm.f90`'s re-implementation of `build_level_lines` -- move
`subdivide_node_line` to `init.f90` (it needs nothing from blocks) and use it
from the other two. Gates: P0 29/29 (`validation/prepare/run_gates.sh`,
incl. the egg-carton classification), P1 15/15, `run_gates_order.sh`
43/43, the refined suite cases (`min_channel`, `les_ibm_refine`,
`beltrami_yslab`). Integer and copies: production flags suffice.

### 5. `interface_correct`: one launch, static table  (review §3 item 5; PERF)

Three near-identical kernels (`pressure_solver.f90:823-913`), 3 launches x
nIter per substage, `face_grad(...)` evaluated per face point although it
depends on `(d, side, b)` only -- a static table like `cfLow` (the
2026-09-14 hoist's own bit-exactness argument). One kernel over a face
list. Worth ~1 ms/step at 16 ranks, 0.7 % of the refined 1-rank device time
(profile section 7.3). Gate: the refined suite cases + `redblack_interface`
at nofma; ncu before/after on `refined_yp82_rect_jacobi`
(`horeka/exchange/run_ncu.sh`, `KERNEL_RE=interface_correct`).

### 6. The per-iteration phi mirror on outlet cases  (review §3 item 16; PERF)

`pressure_solver.f90:319`: `apply_scalar_bc(blk, bc, phi, phiMode)` after
every phi exchange, 36 surface-list launches per step at niter 12.
`jacobi_apply` reads `phi(i) - phi(ip)` with `phi(ip) = -phi(i)`, which is
`2*phi(i)` EXACTLY in IEEE, so the mirror is removable IF nothing else
reads the mirrored ghost. The comment at `:316-318` says the exchange's
tangential extension can write physical halos; the proof owed is a halo
READ-SET proof (which kernels read `phi` at a physical ghost index, incl.
edge halos `phi(nx+1, ny+1, k)`) -- write it down before the change, the
way `sync_divergence_halos`'s was. Gate: the outlet suite
(`validation/freestream/run_bitexact_outlet.sh`, 1 and 4 ranks, CPU+GPU),
`run_gates.sh` groups `low refined restart`, `validation/body_outlet` all
groups incl. `pocket`.

### 7. The red-black sweep's diagonal from the tables  (review §3 item 17; PERF)

`redblack_sweep` (`pressure_solver.f90:1102-1222`): 26 privates, per-cell
`iLo/jLo/kLo/colorOffset`, six branchy `face_grad_denom` per cell and one
fp64 divide -- the whole lever set that gave the Jacobi kernels -14 to -30 %
was never applied; its register count is unmeasured. The sweep visits index
0 (the redundant halo layer) and the tables `dnLow/dnHigh/d1P` start at 1:
they need a `0:nb` extension -- new code, same expression. Red-black is the
production smoother of the boundary-layer tutorial (single level), so this
one lands on a production case. `zero_phi` (9.3 % of the refined red-black
device time) runs only with a 2:1 interface and is NOT on this list.
Gate: `validation/redblack_interface` 1 == 4, the 7-case suite with
`solver = redblack` variants of min_channel/les_ibm (make them: the suite
has no red-black case), `tutorials/turbulentBoundaryLayer` 20 steps vs
`~/pass2_ref_binaries` at nofma; ncu before/after.

### 8. The two "which value wins" decisions  (review F13, F14; NEEDS THE USER)

F13, restart precedence (`io.f90:730-811`): the ini wins for `cflmax
pecletmax dtmax t_final re forcing_* niter sor ibm_enabled` and the BC rows;
the FILE always wins, undocumented, for `globalSize`, `leng`, `isPeriodic`,
`distribution / stretch / natural_dyw_plus` (while `one_sided / subdivided
/ nodes_file / outer_height` come from the ini: a mixed-source node line),
`step`, `t_current`, `cfl`; `dt` deliberately. Grid size, lengths and
periodicity belong to the CASE FILE now, so the ini should have to AGREE
(error on mismatch) rather than be overridden; a `[grid] nx` mismatch today
surfaces as "case file grid differs from [grid] nx", blaming the ini for the
restart file's value. Propose: config-authority for everything except
`step`, `t_current`, `cfl`, `dt`; a mismatch is an error naming both; one
table in `configuration.md`. Gate: `validation/freestream/run_restart_outlet.sh`
(restart == continuous), the S2 restart legs, a deliberate mismatch stops.

F14, the Courant reduction over ghosts (`step.f90:977-1051`): the Courant
loop runs `0..n+1`, the Péclet loop `1..n`; a tangential Dirichlet ghost is
`2v - q`, so at an inlet whose datum exceeds the interior velocity the GHOST
sets dt (conservative, never unstable; `(nb+2)^3/nb^3` wasted cells). Interior
only is NOT bit-exact where a ghost currently binds: gate with `dtmax`-bound
cases (max_abs 0) and expect the inlet cases (blasius2d, the freestream
Poiseuille, the TBL tutorial) to move by the step alone -- show it with
`dtmax` binding on both sides.

### 9. Latent robustness  (review F15, F16)

- `comm.f90:228` `MPI_Cart_create(..., reorder=.true.)` while the block
  split and `select_target_device` assume cart rank == world rank and
  node-contiguous ids (`node_id` is world-rank, `report_partition`'s `nodeOf`
  cart-rank). `MPI_Comm_dup(MPI_COMM_WORLD)` removes the premise; the
  Cartesian topology is a rank-box leftover (`c%coords` is gone, `c%dims`
  survives for the nb rule only).
- `blocks.f90:1318-1336` xz-mode `leaf_key`: the y tile is not masked to
  `LEAF_KEY_YBITS` and `morton2_key` drops bits >= 21; `key_covers` guards
  xyz only -> extend it.
- `blocks.f90:382, 771` `levelOffset` product in default integer: a guard
  (error stop above 2^31), not a kind change (`lidOf` is int32-indexed).
- `comm.f90` component identity by literal (`var /= 1/2/3`, `var <= 3`):
  `VAR_U/V/W`, `NVEL`.
- One OpenMP guard convention: `pressure_solver.f90` wraps every `target`
  construct in `#ifdef USE_OPENMP_OFFLOAD` (26), `step.f90` / `boundary.f90`
  / most of `scalar.f90` none; `turbulence.f90` and `bodyforce.f90` issue
  `enter/exit data` UNGUARDED; `bodyforce.f90` guards its trip kernels so
  the CPU build runs them serially. Equivalent today only because the CPU
  build links no OpenMP (`CMakeLists.txt`). Pick: guard `enter/exit data`
  and `update` only, let `target` constructs fall back. Production-flag
  gate (no arithmetic moves).

### 10. Small leftovers (an hour, one commit)

- `set_constrained_cells_host` (`rans.f90:1221-1255`) duplicates the device
  twin without its `kfix` branch; `rans_prepare` runs the device twin before
  step 1 -- delete after verifying no host read of `k`/`omg` between
  `init_rans_transport` and the first `rans_prepare`.
- `read_case_layout`'s `has_layout` flag (`io.f90:823-944`): the only caller
  stops on false -- fold the stop in, drop the flag and the `moby_solve`
  branch.
- `ibm_body_blocks` runs twice per run on identical data (`update_ibm_mu` ->
  `bodyBlocks`; `narrow_rdenom_blocks` -> `rdenomBlocks`): pressure_solver
  reuses the public `bodyBlocks`.
- `tools/prepare_if_missing.sh:29-33` re-implements `case_file_name` in
  shell and defaults the prefix to `field` where the Fortran defaults to ""
  (a hidden `.case.h5` vs `field.case.h5`): make the Fortran print the
  resolved name (`moby_prepare --case-file-name input.ini`) and let the
  script call it.
- The duplicated ini readers (`scalar.f90:1512-1557` vs `config.f90`, plus
  `flow/case_config_helpers.f90`); the ini is parsed three times.
- `apply_bc` rewrites the BCK_FACE rows (`w = 0`, `0*q + C = C`) on every
  call (review §6 item 2 of 2026-09-26, still alive): skip after init.
- The per-block face-start table replacing `momentum_face_start(...)` x3
  per cell in five kernels (`step.f90:245, 466, 681, 806` + one): integer,
  unmeasured -- measure on the trip case before keeping it.
- `les.f90:55-125`: `update_generic_les_viscosity` has one live case;
  collapse it and drop the mapped `model`/`cs`.
- `scalar_grow` (`scalar.f90:1198-1304`) copies 30 arrays twice through a
  temporary: parse `count` first.
- `check_conjugate_refinement` (`scalar.f90:2007-2088`): six plane loops ->
  one over `(d, side)`.

### 11. Modules and documents

- `scalar.f90` (3428 lines) into `scalar_types` (type + constants) ->
  `scalar_conjugate` (face coefficient, bands, `plane_box_fraction`,
  tangential, indicator, Peclet rate, `init_scalar_conjugate`) +
  `scalar_wallfn` (Kays, Jayatilleke) -> `scalar` (init, maps, transport,
  finish) + `scalar_config` (parse + validate); `scalar_stats` imports the
  face helpers from `scalar_conjugate`. Pure move: production-flag gate.
- CLAUDE.md (190 kB, loaded every session): move the "Active work" log to
  `docs/history/worklog_2026.md` and leave conventions, build/run,
  verification and a two-page state-of-the-solver with pointers. THE USER'S
  CALL (review 6.5); ask before doing it.
- The fourteen airfoil inis without a `cv_box` (`docs/history/
  next_session_port_finish.md` section 3): decide per ini (a box, or the
  warning stays). THE USER'S CALL.

## 2. What is NOT on this list, and why

- `zero_phi` of the red-black path: runs only with a 2:1 interface; the
  production red-black case (the TBL tutorial) never has one.
- The fused predictor split into per-component kernels: the register lever
  is exhausted on it; needs its own ncu-first session (review 7.5 item 6).
- The fp64 divide in `compute_rdenom` (body blocks only) and the reduced
  divergence round for red-black: performance items with their own read-set
  proofs, after the list above.
- The numerics re-validation of `validation/README.md` at niter 12 (the long
  tier) and the review's feature tail (steps 10-12): separate handouts.

## 3. Reference set and gates

Reference: `~/pass2_ref_binaries` (commit `d821341`; banners print
`ac722ec`, identical `src/`), the four builds + `moby_prepare` beside each +
`tools/prepare_if_missing.sh`. **Gate at NOFMA against it** (the fused
projection statements can contract differently at production flags).
Cut a new set from a COMMIT (never a mid-edit tree) after any item that
moves an expression, and record the hash in its PROVENANCE.

The standard round, two streams (~1.5 h): CPU -- `validation/scalar/
run_bitexact.sh` (RANKS=4), `run_bitexact_s3.sh`, `validation/freestream/
run_bitexact_outlet.sh` (1 and 4 ranks), `validation/body_outlet/run_gates.sh
pocket`, `validation/bodyforce/run_gates.sh`, `validation/prepare/run_gates.sh`
(29); GPU -- the 7-case, 9-case, outlet and bodyforce gates with the
`build_gpu_nofma` binaries. Stage the binaries in a scratch `bin/<label>/`
after checking every banner; absolute `NEW=`/`REF=`; `module load` in the
same shell; `OMPI_MCA_hwloc_base_binding_policy=none`; never edit a running
driver; a gate that reuses another group's output compares whatever binary
that group ran. Memory notes: `gate-batch-workflow`, `cleanup-pass2-state`.

Per item, the additional gate named above (the mixed scalar case, the rank
topologies, the restart pair, the dtmax-bound pair, ncu before/after).

## 4. Landmines, this tree

- `config.f90` STOPS on an unknown key/section/value since 2026-10-08: an
  older binary run on a current ini, or a current binary on a `@KEY@`
  template, stops with the key named -- `moby_prepare --check-config
  input.ini` parses an ini without building anything.
- A `[case.<x>]` section of another case than `[case] name` stops the run
  (the body-force gate once ran the GENERIC case silently).
- The case-file echo records one `wall_amp` value; files written with two
  are accepted (`io.f90` compatibility rule) -- do not "simplify" that rule
  away while any pre-2026-10-08 analytic case file is in use.
- `compile.sh` takes ONE mode; `./compile.sh cpu gpu` builds cpu only.
- The predictor's `if (skew)` guard (`skew = dns%conv_skew`) is a REGISTER
  guard (128 vs 100 registers, -1.5 % step): do not fold it.
- `d1P` must stay a separate factor in the diagonal (`(dnLow*mu +
  dnHigh*mu)*d1P`); `x/6.0d0` in a per-DOF kernel is a second divide: write
  reciprocals of non-dyadic constants as multiplies.
- The nut consumer chain (`add_eddy_viscosity_correction`) is the
  bit-exactness argument of T0-T5: a helper there is a nofma gate on every
  LES/IDDES case, or leave it.
