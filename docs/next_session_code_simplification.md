# Next session: code simplification (the review's section 9, continued)

STATUS: handout written 2026-10-20 on head `68c6360`. **Items 1 (the scalar
kernel split), 2 (one face-diffusivity helper), 3 (comm.f90's one
enumeration) and 4 (the lattice probe and the node-line subdivision, once
each) are DONE 2026-10-20, item 5 (`interface_correct` as one launch over
a static face list) DONE 2026-10-21** -- see "Item N, done" at the end of
this file for the measurements and what each one cost; items 6-11 are not
started.
It continues `docs/cleanup_review_2026-10-06.md` section 9 (the
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

## Item 1, done (2026-10-20)

`scalar_transport_kernel` is two kernels, `scalar_transport_plain`
(dirichlet / adiabatic, the wall-function and Kays ladders, the
penalization) and `scalar_transport_conjugate` (the conjugate scheme with
its band and C2 tangential branches), each over its own scalar list
(`sc%plainList` / `sc%conjList`, built in `validate_conjugate_config`,
1-cell dummies when empty, a kernel launches only for a non-empty list).
Every per-cell expression moved verbatim; the lines that changed are the
ones the handout listed (`if (adiab .or. conjug)` -> `if (adiab)`, the
masks, `divuse`, the `rhs` branch), plus three per-run flags the config
guarantees for a conjugate scalar and that the conjugate kernel therefore
no longer tests (`useIbm`, `nondiv`, `wallfn`). The C2 tangential block
stayed inside the conjugate kernel (not a third kernel).

GATE FIRST, as asked: `validation/scalar/mixed.ini` (dirichlet + Kays,
adiabatic + Kasagi source, conjugate; the analytic wavy wall, flow on) and
`mixedles.ini` (the same under WALE) are the first cases with one plain AND
one conjugate scalar; both are in `run_bitexact_s3.sh`. Gates vs
`~/pass2_ref_binaries` (`d821341`), all `max_abs 0`: the 11-case S3 suite
(CPU 1 rank, GPU), the 7-case suite (CPU 4 ranks, GPU), the mixed pair at
CPU 1 and 4 ranks and GPU, 4 == 1 ranks, CPU == GPU; at PRODUCTION flags
the CPU is `max_abs 0` too on `mixedles` and `turbsst` (every dataset, the
split moves no contraction), so **the reference set is unchanged**;
conjugate C1 / C2 / C3 with the new nofma binaries: see the commit message
for the counts.

MEASURED (RTX 3060, production flags, 128x128x32 = 524k cells, nb 32, one
substage launch, nsys over 36 launches; `cuobjdump -res-usage` for the
static resources -- ncu needs performance-counter permission this
workstation does not grant):

| kernel | registers | stack | ms/launch before | after | change |
|---|---|---|---|---|---|
| fused (one dirichlet scalar, constant Pr_t, WALE) | 255 | 88 B | 1.849 | 1.086 (plain) | **-41 %** |
| fused (one dirichlet scalar, KAYS Pr_t, WALE) | 255 | 88 B | 3.965 | 3.198 (plain) | -19 % |
| fused (one conjugate scalar, WALE) | 255 | 88 B | 1.932 | 1.881 (conjugate) | -3 % |
| fused (dirichlet + adiabatic + conjugate, WALE) | 255 | 88 B | 5.824 | 3.974 + 1.882 | +0.6 % |

Registers did NOT move: all three kernels sit at the 255 cap (the
`jacobi_apply` lever is not available here -- occupancy is register-capped
identically before and after), what moved is the stack frame (88 -> 40 B in
both new kernels, i.e. the spills) and the dormant work per scalar. The
production shapes gain (one plain scalar under LES/RANS, the S2 campaigns:
-41 %; one conjugate scalar, the CHT tutorials: -3 %); a case with scalars
of BOTH families pays the hoisted per-cell loads (face velocities, masks,
nu_t averages) once per kernel and is at parity. TWO FINDINGS the
measurement produced: (1) **`prt_model = kays` TRIPLES the plain kernel**
(3.20 vs 1.09 ms per launch at one scalar): `prt_kays` is evaluated six
times per cell per scalar, once per face from the face Peclet number; a
cell-centred Pr_t (one evaluation, faces averaging it) is a numerics change,
not bit-exact -- the review section 8 item, now sized. (2) The dormant
wall-function ladder in the plain kernel costs NOTHING: a build with it
compiled out runs 1.095 / 3.204 ms on the two plain cases (vs 1.086 /
3.198), so a third "wall-function family" kernel is not warranted.
Not done: the same cut on `rans_transport_kernel` (the handout's "second,
after the scalar one") -- its own session.

## Item 2, done (2026-10-20) -- and what it costs

THE HELPERS. `plain_face_diffusivity(dm, useNut, wallfn, nutL, nutR, ypL,
ypR, pr, prt, prtModel, re, P, y+_T)` is the plain family's ladder
(molecular; + the thermal wall function's face value; + nu_t,face/Pr_t),
`conjugate_face_total(dm, phL, phR, ks, rc, invd, addNut, nutFace, ...)` and
its F2 twin `band_face_total(...)` the conjugate family's (the C1 / band
coefficient + the eddy part where `addNut` = a fluid cell's uncut face
under a turbulence model). The transport kernels AND the two statistics
kernels of scalar_stats.f90 (sample, heat) call them, so the S4 statement
"the statistics report the flux the kernel applied" is one definition, not
three copies. The band dispatch is the CALLER's, once per scalar around its
six faces -- see the cost table for why. The kernels hoist the nu_t and
(under wall functions) y+ stencils once per cell, GUARDED, because both
arrays are 1-cell dummies when their model is off and a helper argument
is evaluated at the call site; the helper forms the 0.5 (L + R) average in
the order the kernels used.

NOT DONE, measured and dropped: the C2 tangential stencil as ONE
offset-form loop in the conjugate kernel (the handout's "six blocks ->
indicator_pass form"). It is bit-exact (gated with `tangential_correction
= true`) but costs the conjugate kernel a further +3.6 % (stack frame 40 ->
160 B: the loop's index arithmetic on a kernel already at the 255-register
cap). The six blocks stay; what is unified is the CONVENTION: the two
tangential axes in increasing order everywhere (`indicator_pass` and
`face_corr_rate` used (z, x) for a y face, the kernel (x, z); the dot
products in `conjugate_tangential` round differently) -- the printed
indicator moves at round-off on y faces, invisible at its 4 printed digits.

COST (RTX 3060, production flags, 524k cells, one substage launch, nsys,
36 launches; item 1's kernels are the baseline):

| variant | plain kernel (1 dirichlet, WALE) | conjugate kernel (1 conjugate, WALE) |
|---|---|---|
| item 1 (inline ladders) | 1.085 ms, stack 40 B | 1.883 ms, stack 40 B |
| helpers, band dispatch INSIDE the helper, + the tangential loop | 1.109 (+2.2 %) | 2.043 (+8.5 %), stack 160 |
| helpers, band dispatch inside, six blocks | 1.111 (+2.4 %) | 1.973 (+4.8 %), stack 80 |
| **shipped: helpers, band dispatch by the caller, six blocks** | **1.110 (+2.0 %), stack 64** | **1.929 (+2.4 %), stack 96** |

So the one definition costs about 2 % of EACH transport kernel: on a
scalar case the kernel is ~10 % of the step, i.e. ~0.2 % of the step, at
the edge of what the campaign matrix resolves. Every form tried keeps the
255-register cap and raises the stack frame; the inlined helper changes the
schedule, not the arithmetic. **This is the user's call: the commit stands
alone, a one-commit revert returns the inline ladders.** The handout
labelled the item performance-neutral; it is not quite.

GATES (all at nofma vs `~/pass2_ref_binaries`): fields `max_abs 0` on CPU
and GPU for the mixed pair, `wavy.ini` with `tangential_correction = true`,
`wfs180_y30` (the wall-function ladder), a banded slab (`band_slab.ini`,
F2) and `ibmwf180`; the STATISTICS FILE and the HEAT FILE byte-identical on
the CPU for `ibmwf180` (wall functions + heat) and `wavy` + stats (the
conjugate Nusselt columns); the 11-case and 7-case suites CPU + GPU; S4
gates; conjugate C1-C3 -- counts in the commit message; the S4 statistics
groups (stats accum plane levels restart det noeffect tools, 25 checks)
PASS once `turbles_67600.h5` -- the S2 campaign snapshot they restart from,
which lives in the sibling checkout `mobydiff.scalar`, not in this tree --
is linked into `validation/scalar/`. FOUND, PRE-EXISTING (identical on the
reference binary): the S4 `heat` group's check (u), `check_scalar_stats.py
heat ... --balance`, reports the graded column as 3.8e+25 against the
solver's 0.25 -- its Python form of the penalization term is stale (the
2026-10-01 rational factor, or the case file's coefficient column), while
the `--balance` closure at 1 % passes. One for the one-by-one
re-validation pass; not touched here. LANDMINE FOUND: **the
GPU statistics and heat files are not bit-reproducible run to run with ONE
binary** (4e-16, the `atomic update` accumulation order; the fields are
exactly reproducible), so a statistics identity is a CPU statement and a GPU
one is round-off.

## Item 3, done (2026-10-20)

`init_block_exchange` enumerates ONCE. `collect_pieces` is the body the
three families used to repeat (resolve -> round filter -> piece loop ->
`candidate_boxes`), appending `piece_type` records -- what a transfer piece
is before any family stores its part -- in the canonical (direction,
candidate, piece) order for one destination block, one source owner (or
any) and one round. My blocks are enumerated once for EVERY owner, in
(round, block, direction, candidate, piece) order: the local entries are
that list filtered to my own rank, the recv entries of peer p the same list
filtered to p -- so the per-peer re-run of `resolve_neighbors` over all my
blocks is gone and the order both ends rely on is visibly one list. The
peers' blocks are enumerated once per peer with my rank as the owner (the
send list, peer-major by construction). Three small `fill_*` routines hold
what differs (local: both slots, gather map, blend, normal, phi row; recv:
destination slot and box, blend, normal; send: source slot, destination
box, gather map), `allocate_entries` the allocation block, and the pass
loop (count, then fill) is gone with it: the records are counted. Also
`axis_dim(off)` (the signed axis of a pure axis direction) is the one loop
behind `interface_normal_dim` / `div_entry_var` / `iface_restrict_normal`.
Lines: +315/-229 in comm.f90 -- not shorter, since the fills and the
record type are new text, but the 40-line enumeration body exists once and
the three families are its three uses.

Integer bookkeeping and copies, so a copy is a copy; gated at nofma anyway
vs `~/pass2_ref_binaries`: `min_channel` at 1, 2, 3 AND 4 ranks each
`max_abs 0` vs the reference and 2 == 3 == 4 == 1 ranks (the load-bearing
gate: four peer topologies, the one the divergence-list work was validated
with); `validation/redblack_interface` 1 == 4 ranks EXACT and the uniform
3-level patch 0.0, CPU and GPU; the 7-case suite CPU 4 ranks and GPU 7/7;
the 11-case suite CPU and GPU 11/11.

## Item 4, done (2026-10-20)

`lattice_candidates(blk, level, cl, off, n, ids, ops, subs)` in blocks.f90
is THE 2:1-aware probe: the same-level leaf (OP_COPY), else the coarser
leaf covering the cell (OP_PROLONG), else the finer leaves inside it that
are adjacent across `off` (OP_RESTRICT, fixed child order, their parity),
`off = 0` keeping every child. The four copies are its four uses:
`face_kind` (off = 0: ANY finer leaf makes the face FINE -- the one place
the former copies differed, so the probe has the mode rather than the
callers being made "equivalent"), `occupied_any_level` (off = 0, n > 0),
`partition_surface` (off = +e_d, the children on the near side) and
comm.f90's `resolve_neighbors` (off, plus owner/slot and the PROLONG
parity). The OP_* constants moved to blocks.f90 with it: what the probe
finds IS the exchange operation. `subdivide_node_line` moved to init.f90
(it needs nothing from blocks) and is called by `build_node_line`'s
`subdivided` branch, by blocks' `subdivide_level_line` -- now public, and
what ibm.f90's level-line loop calls instead of re-implementing
`build_level_lines` -- i.e. one midpoint formula for the three. Integer
bookkeeping and the same formula: gated at nofma anyway vs
`~/pass2_ref_binaries`: P0 29/29 (`validation/prepare/run_gates.sh`, the
egg-carton classification included), P1 15/15 (`run_gates_stl.sh`),
`run_gates_order.sh` 42/43 -- the one failing leg compares today's
PRODUCTION binary with the 2026-10-01 one at production flags and differs
by 1e-13, the FMA class the 2026-10-08 projection fusions introduced
(pre-existing, the same on the previous commits) -- the 7-case suite (CPU
4 ranks, GPU) 7/7, the 11-case suite (CPU, GPU) 11/11,
`validation/redblack_interface` CPU and GPU all PASS.
CAUGHT BY THE GPU GATE, FIXED BEFORE THE COMMIT: the probe's first form
sized its result arrays 4, the exchange's maximum of ADJACENT children, but
`face_kind`'s "every child" mode returns up to 8 -- a stack overflow inside
`face_kind` on any refined cell that clobbered the caller's saved registers.
The CPU builds passed every gate by stack layout; the GPU builds (same host
code, `-mp=gpu`) segfaulted on every refined case, in `build_block_metadata`
right after the call, with a bounds-checked -O0 build running clean. The
arrays are `MAX_CANDIDATES = 8` now and `resolve_neighbors` copies the
at-most-4 adjacent ones into the entry arrays with a guard. LESSON: a gate
that passes on one build is not a gate; and gdb on a deterministic
segfault (callee-saved registers holding small integers = a stack
overflow in the callee) beats reasoning about the source.


## Item 5, done (2026-10-21)

`interface_correct` is ONE kernel over a face list built once
(`init_interface_faces`, called from both projections beside the
`hasIface` test): per 2:1 interface face of the rank its block, normal
direction, face index (1 or nb+1) and `ifCoef = face_grad(kind, .true.,
d1f, refd(d))` -- the cfLow argument: the same expression on the same
inputs, evaluated once on the host instead of per point on the device.
The collapsed loop runs (face, t2, t1) over the largest tangential extents
of the three orientations and skips a face's excess (nothing on cubic
blocks); low and high faces are one statement, `phi(f-1) - phi(f)` with f
the face index. `refd` left the kernel's signature; a rank without an
interface face (nIfaces = 0) does not launch at all, where it used to
launch three kernels over every block. 52 registers, no stack (was 3 x 54).
MEASURED on the RTX 3060, before = `0882695`, production GPU builds, 100
steps, two repetitions each (spread < 0.7 %):

| case | step before -> new | `proj_timing apply` | kernel device time |
|---|---|---|---|
| `min_channel` GPU (refined bands, 1600 leaves, Chebyshev 12) | 36.90 -> 35.63 ms (**-3.4 %**) | 14.36 -> 13.00 ms/step (-9.5 %) | 0.829 -> 0.304 ms/step (-63 %) |
| `redblack_interface/refined_channel` (red-black 3) | 14.37 -> 13.80 ms (**-4.0 %**) | 1.001 -> 0.419 ms/step (-58 %; the bucket IS this kernel there) | 0.254 -> 0.108 ms/step |
| `block_nb/beltrami_rect` (nb 8 4 8, 3D patch, 1600 leaves) | 16.01 -> 14.75 ms (**-7.8 %**) | 6.37 -> 5.16 ms/step (-19 %) | 0.657 -> 0.339 ms/step (-48 %) |

The step moves by MORE than the device time removed (min_channel: 1.3
ms/step off the step, 0.5 off the device): these steps are launch-bound
(the review's 55 % device busy), so the two launches per call saved on the
host are the larger half -- red-black's per-call `apply` wall time 55.6 ->
23.3 us says the same, one launch's host cost instead of three. The
rectangular case pays the skip (its one kernel 9.4 us against the three
at 4.4/6.9/7.0) and still halves the device time. ncu on
`refined_yp82_rect_jacobi` was NOT run (HoreKa): the per-launch floor is
the mechanism and the workstation has no counter permission; the handout's
"~1 ms/step at 16 ranks" stands as an estimate, this one measured at 1 rank
with no MPI. NEW GATE CASE, because none existed: `validation/block_nb/
beltrami_rect.ini` -- a real flow (phi /= 0) through a 3D patch on
NON-CUBIC blocks, i.e. the only case whose interface faces have unequal
tangential extents; the uniform-flow cases of that directory are blind to
the correction (phi = 0). It is the 8th case of
`validation/scalar/run_bitexact.sh` (the "7-case suite" of the records).
GATES (nofma vs `~/pass2_ref_binaries`): the 8-case suite CPU
4 ranks AND GPU 8/8 `max_abs 0` each; the 11-case suite CPU AND GPU 11/11;
`min_channel` and `beltrami_rect` at 1 and 2 CPU ranks `max_abs 0` (with the
suite's 4 and the GPU's 1: four rank counts, the face list is per rank);
`validation/redblack_interface` CPU (incl. 1 == 4 ranks EXACT) and GPU ALL
PASS; `validation/block_nb` CPU and GPU ALL PASS. ALSO at PRODUCTION flags:
`min_channel`, `beltrami_rect`, `les_ibm_refine` on the GPU and the first two
on the CPU `max_abs 0` -- the coefficient is the same value from the same
expression, nothing contracts differently -- so **the reference set is
UNCHANGED (`~/pass2_ref_binaries`)**. Landmine met, not caused:
`validation/block_nb/gate_unset.case.h5` was a ONE-rank file from
2026-09-29 and stopped the `nb unset` leg at `RANKS=4` with "rank owns no
blocks" (the prepare log beside it says `rank layout 1 1 1`); deleted and
rerun, PASS. An nb-less case file is pinned to the rank count that prepared
it -- a gate directory keeps the last run's.
