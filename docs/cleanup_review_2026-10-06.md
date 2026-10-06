# Second cleaning and profiling pass (2026-10-06, head `b1f41be`)

Scope: the whole solver source as it stands after the outlet work (O0-O5,
B1-B2, the backflow datum), step 7 (prepare does all preprocessing), F1/F7,
the sum limiters, the trip port and the rough-wall keys -- everything that
landed since `docs/numerics_review_2026-09-26.md`, whose sections 1-5 and 9-10
are closed and whose section 6 is the open list this pass re-ranks. Read in
full: `step.f90`, `pressure_solver.f90`, `boundary.f90`, `comm.f90`,
`blocks.f90`, `init.f90`, `config.f90`, `ibm.f90`, `io.f90`, `prepare.f90`,
`moby_prepare.f90`, `moby_solve.f90`, `scalar.f90`, `scalar_stats.f90`,
`turbulence.f90`, `les.f90`, `bodyforce.f90`; the kernels of `rans.f90`;
every driver under `validation/`; the `docs/next_session_*.md` set. Nothing
in `src/` was modified by this pass; every item below is a finding with a
file and line, ranked in section 9, and the profile (section 7) is the fresh
HoreKa baseline of the current head (job 5184720), not a comparison.

Rule carried in from the previous session: a gate's bound is taken from the
state the gate READS (the pocket divergence bound was set from t = 150 and
failed at t = 50).

---

## 1. Dead code and leftovers

Small and mechanical in the numerics core; the outlet work left the files
clean (no `outflow_copy`, `BCK_FACE_OUTFLOW` or outlet argument of
`interface_correct` survives anywhere in `src/`, and no comment still claims
`apply_bc` runs per iteration). What is left:

**Unused symbols**

| where | what |
|---|---|
| `step.f90:16-19` | `NVEL`, `FACE_COARSE`, `FACE_FINE` imported, never used |
| `pressure_solver.f90:169` | `character(len=32) :: env`, the removed `CHEB*` env toggles |
| `pressure_solver.f90:214` | `pressure_projection(ps, blk, dns, ...)`: `dns` never referenced |
| `pressure_solver.f90:656` | `jacobi_apply(ps, ...)`: `ps` never referenced |
| `boundary.f90:11` | `SIDE_MAX` unused |
| `boundary.f90:162` | local `DIRICHLET = 0` shadows `BC_DIRICHLET` |
| `boundary.f90:485-488` | `case (BC_NEUMANN, BC_OUTFLOW)` on a tangential row: unreachable, the comment says so |
| `boundary.f90:513-527` | `dn` formed for every row, read by the Neumann branch only |
| `comm.f90:49, 232` | `c%coords` filled by `MPI_Cart_coords`, read nowhere (the Cartesian topology is a rank-box leftover; `c%dims` survives for the nb rule only) |
| `comm.f90:561, 769` | `ierr` unused; `if (allocated(divVarL)) deallocate` on a fresh local |
| `comm.f90:2114, 2225, 2325, 2404` | `c%nLocal/nSend/nRecv/nPeers` in `map(to:)` of kernels that never read them |
| `comm.f90:2226` | `copy_local_scalar_entries` maps `blk%x/y/z` the body never reads (the deleted cubic reconstruction); with them the `blk` argument of that routine and of `exchange_scalar_halos` (`:1700`) |
| `gpu_runtime.f90:34`, `moby_solve.f90:136` | `enter_grid_data` maps `g%xNode/yNode/zNode` to the device; no target region in `src/` reads `g%` (verified: the readers in scalar_stats/rans are host code). A dead device mapping since Phase 0 |

**Stale comments** (these mislead and should be fixed first)

| where | says | is |
|---|---|---|
| `pressure_solver.f90:784-785` | `interface_correct` writes "the high face at index nb" | the kernel writes `nx+1` (`:828, 855, 882`) |
| `pressure_solver.f90:64-66, 550` | `compute_rdenom` "formed once per substage" | body blocks only, once per run body-free (`:137-146` of the SAME comment block says so) |
| `pressure_solver.f90:653-655` | `jacobi_apply` writes "only the cell's own LOW faces" | also the outlet high faces (`:730-744`) |
| `comm.f90:1623-1624` | phase order justified by "the prolong tangential interpolation" | that branch is gone (`:1991-1998` says so) |
| `comm.f90:453`, `pressure_solver.f90:347-350` | "the DIVERGENCE prefix" | an index list since 2026-09-17 |
| `comm.f90:1677` | `exchange_halos(..., interp=)` | the flag means "include cross-level entries"; only red-black passes `.false.` |
| `blocks.f90:13-17, 561` | "ONE block per MPI rank ... `field_type`", "mirroring `enter_field_data`" | both gone |
| `blocks.f90:181-183`, `init.f90:94-95` | `active` "only meaningful with [blocks] nb"; "0 = one block per rank box" | nb is always derived now |

**Elsewhere** (scalar/turbulence/rans/les/bodyforce/ibm/io):

| where | what |
|---|---|
| `rans.f90:1617, 1665` | `nutc` private in `rans_transport_kernel`, never used |
| `rans.f90:17-20, 1879-1880` | "gamma_sep not yet implemented" -- `lm_gamma_sep` exists (`:1049`) and is applied (`:2070`) |
| `scalar.f90:46-49, 1727-1731, 1736-1738` | pre-step-7 comments (two dwall producers; "analytic coefficients computed on the device") |
| `scalar.f90:2570-2572` | `enter_scalar_data` maps `initValue/inlet/initProfile`, host-only arrays |
| `scalar.f90:304, 277` | `scalar_count` public with no external caller; `conjugate_face_local_k` public "for scalar_stats", never called |
| `les.f90:85, 111-115` | unreachable `return`; a one-case `select` with `model`/`cs` mapped for the dead dispatch |
| `les.f90:167` | literal `1.0d20` where `SOLID_FACE_THRESHOLD` is "the ONE definition" (`init.f90:22`) |
| `ibm.f90:295` | a second, Re-dependent definition of "solid": `thresh = 0.5*SOLID/re` |
| `bodyforce.f90:29-31` | header still names the Schlatter & Orlu trip; the CaNS port is below it |
| `ibm.f90:1310-1311, 273-275` | "must read the DEVICE copy of coef: the analytic path fills it there" -- false since 7-4, coefficients come from the file on the host |
| `ibm.f90:449` | `if (len_trim(dns%case_file) == 0) return` unreachable (`moby_solve.f90:93` always sets it) |
| `io.f90:838-863, 821, 831, 1110` | `read_case_layout`'s `has_layout` and three legacy-file notes for paths the solver refuses; "[ibm] case file" half-edit |
| `moby_solve.f90:162-168` | "keeps the parameters for the diagnostics that print them" -- nothing prints them (the call stays for key validation) |
| `walldist.f90:286-288` | `bisect_indicator` justified as a host twin of a device routine that is now host code (`ibm.f90:1145`); not mergeable bit-exactly (different exit convention) -- fix the reason |
| `prepare.f90:185` | "exactly what init_rans_geometry computed inline" (gone) |
| `ibm.f90` | the module has no `private` statement: 15 helpers are public by default |

**Two predicates for one property**: `io.f90:604, 647` test
`faceBcType(dir, f) == BC_OUTFLOW` where `bc%faceOutlet(f)` encodes the same
fact; `BC_OUTFLOW` as a row TYPE is otherwise consistent (`boundary.f90:215,
327, 472` derive it, flag the face, resolve it to no row).

## 2. Duplication

The arithmetic is written once; what is duplicated is INDEX BOOKKEEPING and
ENUMERATION, which is where a future edit goes wrong by touching one copy.

**Index bookkeeping (step/boundary)**

- The outlet face/neighbour pair `(1, 2)` / `(nb+1, nb)` is spelled five times
  (`step.f90:605-611`; `boundary.f90:846-856, 913-919, 507-511, 964-968`) and
  the ghost pair `(0, 1)/(nb+1, nb)` once more (`:1022-1028`).
- `dir = (face_id + 1)/2; side = modulo(face_id - 1, 2)` is inlined in four
  kernels (`step.f90:598-599`; `boundary.f90:838-839, 909-910, 1018-1019`)
  while `boundary_face_dir/side` exist (`:780-790`) but are not `pure` +
  `declare target`.
- `uStartX = momentum_face_start(blk%physLow(1,b))` (x3) is recomputed per
  CELL in five kernels (`step.f90:244-246, 465-467, 686-688, 737-739,
  811-813`) for a per-block quantity.
- The per-direction `d1` pick exists three times (`pressure_solver.f90:515-539`
  host; `step.f90:617-626`).
- The bilinear edge-`nut` interpolation is spelled 12 times in
  `add_eddy_viscosity_correction` (`step.f90:821-871 ...`). CAVEAT: the nut
  consumer chain is the bit-exactness argument of the turbulence phases; a
  helper may move FMA contraction -- nofma gate if touched at all.

**Two projections, one tail**: `pressure_projection` and `redblack_projection`
share verbatim setup (`refd/hasIface` `:244-245` vs `:1002-1003`,
`outLow/outHigh` `:254-261` vs `:1007-1010`) and tail (final `apply_bc` + full
exchange + profiling, `:325-353` vs `:1036-1048`). The projection diagonal is
written twice: from the static tables in `compute_rdenom` (`:586-589`) and as
six branchy `face_grad_denom` calls per cell in `redblack_sweep`
(`:1154-1163`) -- red-black never received the 2026-09-14 table hoist.
`interface_correct` is three near-identical kernels (`:812-891`);
`append_boundary_face_points` three copies of an 8-line body
(`boundary.f90:717-756`).

**Enumeration (comm/blocks)**

- `init_block_exchange`: the local, recv and send enumerations are the same
  ~40-line body three times (`comm.f90:575-630, 646-690, 696-742`); the recv
  loop re-runs `resolve_neighbors` over ALL local blocks once per peer.
- `interface_normal_dim` / `div_entry_var` / `iface_restrict_normal`
  (`comm.f90:1228-1252, 1316-1326`): one "pure +axis face: which dim" loop
  with three op tests.
- The lattice occupancy probe (same level -> parent -> children) exists four
  times: `face_kind` (`blocks.f90:1167-1208`), `occupied_any_level`
  (`:1213-1233`), `partition_surface` (`:1557-1576`), `resolve_neighbors`
  (`comm.f90:950-1013`); `face_kind` is literally `resolve_neighbors`' op
  (COPY -> OPEN, PROLONG -> COARSE, RESTRICT -> FINE, none -> CLOSED).
- `wrap_lattice` (`blocks.f90:1136-1155`) vs `neighbor_origin`
  (`comm.f90:1385-1406`).
- Midpoint subdivision of a node line three times: `subdivide_node_line`
  (`blocks.f90:1623-1637`), the `subdivided` branch of `build_node_line`
  (`init.f90:396-400`), and `ibm.f90:779-801`, which re-implements
  `build_level_lines` (`blocks.f90:655-671`) wholesale.
- `init_block_set` / `init_block_set_from_table` share a verbatim tail
  (`blocks.f90:222-233` vs `:401-417`); the box-intersection predicate twice
  (`:819-823`, `:838-842`); `av = c%activeVars` + its 5-line comment four
  times (`comm.f90:2040, 2110, 2343, 2413`).

**Physics written twice (scalar/rans/turbulence)**

- The face-diffusivity ladder (`conjug ? banded|plain, nut not at a cut :
  useNut.and.wallfn : useNut`) is written three times: the transport kernel
  (`scalar.f90:3035-3315`), the statistics sample (`scalar_stats.f90:414-446`)
  and the heat diagnostic (`:833-860`); the convective-mask ladder twice
  (`scalar.f90:2951-2969`, `scalar_stats.f90:402-410`). The S4 design ("the
  stats use the TRANSPORT kernel's own face diffusivity") is an identity by
  copy today; one `declare target` helper makes it an identity by construction.
- The C2 tangential stencil three ways (`scalar.f90:3108-3230` six direction
  blocks; `indicator_pass :2483-2530` offset form; `face_corr_rate :2316-2370`).
- SST: `arg2/f2/nutloc` (`rans.f90:1409-1412` == `:1744-1747`), the strain
  `s2` (`:1398-1401` == `:1703-1707`), the k-based y+ (`:1403, 1498, 1766`);
  the flux-divergence line 12 times (`:1812-1875, 1941-1996`), the first-order
  upwind block 4 times (`:1806-1830, 1920-1938`);
  `set_constrained_cells_host` (`:1221-1255`) vs `rans_set_constrained_cells`
  (`:1260-1303`), the host twin lacking `kfix`.
- The velocity-gradient tensor is CALLED from five kernels (`les.f90:115,
  199`, `turbulence.f90:353`, `rans.f90:1403, 1702`): under IDDES up to four
  evaluations per cell per substage of the same `q`.
- Constants: `IDDES_KAPPA` (`turbulence.f90:33`) duplicates `WF_KAPPA`
  (`rans.f90:168`) because turbulence sits below rans; `WF_CMU25` is
  `SST_BETA_STAR**0.25`; `pi` five times; `solid_threshold` copied into a local
  and mapped per launch (`rans.f90:1640`, `les.f90:79`) where `scalar.f90:2903`
  uses the parameter directly (the pattern to keep).
- Config readers: `scalar.f90:1512-1557` duplicates config.f90's
  (`:1538` says "read_bool verbatim"); `case_config_helpers.f90` has a third
  tokenizer; the ini is parsed three times (`read_case_name`,
  `flow%read_config`, `read_runtime_config`).
- `ibm_body_blocks` runs TWICE per run on identical data (`ibm.f90:1258` ->
  `bodyBlocks`; `pressure_solver.f90:390` -> `rdenomBlocks`).
- `case_file_name` has a shell re-implementation in
  `tools/prepare_if_missing.sh:29-33` that defaults the prefix to `field`
  where the Fortran defaults to `""` (a hidden `.case.h5` vs `field.case.h5`
  for an ini without `field_prefix`).

**Fortran <-> Python lockstep pairs** (name both ends when either moves):
`blocks.f90` key functions <-> `tools/partition_analysis.py:149-196, 240`;
`init.f90` natural line <-> `tools/make_channel_restart.py:31ff` (whose
`morton_key` at `:118` is the LEGACY interleave -- harmless only because
`block_row_map` matches restart rows on (origin, level)); `build_leaf_table`
<-> the retired `tools/mobygeom.py:1652-1744`.

## 3. Simplifications (performance-neutral)

Each item says whether it moves an expression. "Moves nothing" = integer or
control-flow only, gateable at production flags; "moves" = the nofma gate.

**Bit-exact, moves nothing**

1. Remove the dead symbols and fix the stale comments of section 1.
2. `boundary_face_dir/side` made `pure` + `declare target` and used from the
   four kernels; one `face_pair(side, n, f, nbr)` helper for the five
   `(1,2)/(nb+1,nb)` spellings.
3. A static per-block `faceStart(3, b)` table (from `physLow` at init)
   replacing the per-cell `momentum_face_start` triple in five kernels.
4. `outlet_flags(bc, outLow, outHigh)` + a shared `projection_finish(blk, bc,
   c)` for the two projections' identical setup and tail.
5. `interface_correct`: `face_grad(...)` depends on `(d, side, b)` only yet is
   evaluated per face point (`pressure_solver.f90:825, 830, ...`); a static
   table like `cfLow` (the 2026-09-14 hoist's own argument). Pairs with the
   open "three kernels -> one launch over a face list".
6. `apply_bc` writes the BCK_FACE rows (`w = 0`, `boundary.f90:470-471`) on
   every call; `0*q + C` is `C` exactly and nothing modifies a pinned face
   between calls -- the 2026-09-26 §6 item 2 is still alive in this form.
7. `comm.f90`: enumerate candidates once per (block, dir) and bucket by owner
   (order preserved as (peer, round, block, dir, cand, piece)); one
   `plus_face_dim(off)` helper for the three `+axis` functions; `face_kind`
   derived from `resolve_neighbors`' op; `finish_block_set` for the two
   builders' shared tail; `subdivide_node_line` moved to `init.f90` and used
   by `build_node_line` and `ibm.f90`.
8. `blk` declared `intent(inout)` in read-only kernels (`scalar_stats`
   sample/heat, both LES kernels, `rans_assemble_nut`, `compute_iddes_fd`) ->
   `intent(in)`; `turbulence_is_enabled(turb) .and. allocated(turb%nut)` --
   the second test is implied by `init_turbulence` (four sites).
9. `les.f90:55-125`: `update_generic_les_viscosity` has one live case
   (Smagorinsky); collapse it and drop the mapped `model`/`cs`.
10. Delete `set_constrained_cells_host` (`rans.f90:1221-1255`): `rans_prepare`
    runs the device twin before step 1 anyway (verify no host read of `k`/`omg`
    between `:931` and the first `rans_prepare`).
11. `scalar.f90`: `check_conjugate_refinement`'s six plane loops -> one over
    `(d, side)`; `scalar_grow` copies 30 arrays twice through a temporary --
    parse `count` first and allocate once.
12. Constants in one place: `WF_KAPPA`/`WF_E` to `init.f90` beside
    `SOLID_FACE_THRESHOLD` (removes `turbulence.f90:33`'s `IDDES_KAPPA` copy and
    `scalar.f90`'s only dependency on `rans`); `les.f90:167`'s literal `1.0d20`;
    `pi` (five definitions: `bodyforce.f90:267, 475`, `init.f90:675`,
    `ibm.f90:525`, `pressure_solver.f90:188`).

**Bit-exact at nofma (expressions move or control flow changes contraction)**

13. **Split `scalar_transport_kernel` by mode** (review §6 item 6, now sized):
    the kernel has ~100 private variables (not the ~60 the review guessed)
    and re-evaluates per CELL the per-run flags `useIbm/anyConj/useNut/
    nondiv/wallfn` and the per-SCALAR flags `adiab/conjug/tang/banded/prtModel/
    srcType`. A plain kernel (dirichlet/adiabatic, nut, wall function, Kays;
    ~57 privates) and a conjugate kernel (~70, ~50 if the shipped-disabled C2
    tangential block becomes a third optional kernel) over per-mode scalar
    lists built in `validate_conjugate_config`. The per-cell expressions move
    verbatim (list in the agent notes: `uw..wt/divu`, `cl*`, `sol*`, `nt*`,
    `m*`, the two `conv` forms, `dm`, the `useNut` ladders, `fw/fe/diff`,
    `srcVal`, `rhs`, `ss`, the penalization); the lines touched are the `is`
    loop header, `if (adiab .or. conjug)` -> `if (adiab)`, `cm* = m*`,
    `divuse = divu`. Gate: 9-case scalar suite + conjugate C1-C3 at nofma, plus
    a case with one plain AND one conjugate scalar, which does not exist yet
    and must be written first. The same cut applies to `rans_transport_kernel`
    (103 privates, the transition block alone 190 lines): a gamma/Re_theta
    kernel first writing `gams/rets` and a per-cell `gameff`, the k/omega
    kernel reading it -- the pure-RANS branch stays verbatim.
14. One `declare target` `scalar_face_diffusivity(...)` for the ladder written
    three times (`scalar.f90:3035-3315`, `scalar_stats.f90:414-446, 833-860`)
    and the six C2 tangential blocks (`scalar.f90:3108-3230`) in the
    offset-vector form `indicator_pass` already uses (`:2483-2530`). ~250
    lines less; the identity is per-cell operand order, so nofma.
15. The bilinear edge-`nut` helper in `add_eddy_viscosity_correction` -- the
    nut consumer chain is the bit-exactness argument of T0-T5; nofma gate on
    every LES/IDDES suite case, or leave it.
16. **Per-iteration phi mirror on outlet cases**: `jacobi_apply` reads
    `phi(i) - phi(ip)` with `phi(ip) = -phi(i)`, which is `2*phi(i)` EXACTLY in
    IEEE, so the per-iteration `apply_scalar_bc(blk, bc, phi, phiMode)`
    (`pressure_solver.f90:319`; 36 launches per step at niter 12) is removable
    IF nothing else reads the mirrored ghost. The comment at `:316-318` says
    the tangential extension can write physical halos; a halo read-set proof
    comes first (no kernel found reading the edge halos `phi(nx+1, ny+1, k)`).
    Candidate, not asserted.
17. Red-black diagonal from the `dnLow/dnHigh/d1P` tables instead of six
    `face_grad_denom` per cell: the sweep visits index 0 (the redundant halo
    layer, `:1128-1130`) and the tables start at 1 (`:475-477`), so they need
    a `0:nb` extension -- new code, same expression; production-flag gate on
    `validation/redblack_interface` + the red-black suite cases.

## 4. Correctness risks found by reading

None is a defect in a shipped configuration; each is a premise the code
relies on without checking, or an inconsistency that a later edit could turn
into one. Host/device staleness and map-bound-only variables were audited in
every file read: CLEAN (every host read of a mapped array sits before the
`enter_*_data` or behind a `target update from`; every section bound is also
an ordinary loop bound).

1. **The Courant reduction spans the ghost planes, the Péclet one does not**
   (`step.f90:1012-1014` `0..n+1` vs `:1043-1045` `1..n`, verified). A
   tangential Dirichlet ghost is `2v - q`, so at an inlet whose datum exceeds
   the interior velocity the GHOST sets `dt`; and the `(nb+2)^3/nb^3` extra
   cells are wasted. Conservative, never unstable. Restricting to `1..n` is
   not bit-exact where a ghost currently binds (the step changes): gate with
   `dtmax`-bound cases and expect the inlet cases to move.
2. **A floating-point weight used as a type tag** (`boundary.f90:962`:
   `bcW(var,f) /= -1.0d0` means "Neumann row"). Exact today (`w` is ±1), but
   the information is `faceBcType(var,f) == BC_NEUMANN`; the 2026-08-05 lesson
   ("never rely on a floating-point value as a classification") in spirit.
3. **The Cartesian communicator's `reorder = .true.`** (`comm.f90:225`): the
   block split and `select_target_device` assume cart rank == world rank and
   node-contiguous ids; `node_id` is a WORLD-rank id (`:235`) while
   `report_partition` indexes `nodeOf` by CART rank (`:538`). Open MPI ignores
   `reorder`, so latent; `MPI_Comm_dup(MPI_COMM_WORLD)` removes the premise
   (the topology itself is a rank-box leftover, `c%coords` is never read).
4. **xz-mode `leaf_key` is unchecked for width** (`blocks.f90:1333-1334`): the
   y tile is not masked to `LEAF_KEY_YBITS` and `morton2_key` drops bits >= 21
   of x/z; `key_covers` (`:1341`) guards xyz mode only. An xz table beyond
   2^21 cells per direction would alias keys silently at read.
5. **`levelOffset` is formed in default integer** (`blocks.f90:383, 772`:
   `int(product(int(winDims)))`): a per-level window above 2^31 cells wraps
   silently; the comment at `:98-99` cites 7e9-block finest lattices. `lidOf`
   is int32-indexed anyway, so the fix is a guard (error stop), not a kind
   change.
6. **Component identity by literal** (`comm.f90:2139-2149, 2179, 2373-2375,
   2446`: `var /= 1/2/3`, `var <= 3`) instead of `VAR_U/V/W`, `NVEL`, which
   comm.f90 does not import. Correct today, fragile if the layout moves.
7. **The OpenMP guard convention is split three ways**: `pressure_solver.f90`
   wraps every `target` construct in `#ifdef USE_OPENMP_OFFLOAD` (26 ifdefs);
   `step.f90` (9 kernels), `boundary.f90` (3) and most of `scalar.f90` do not;
   `turbulence.f90:155-186` and `bodyforce.f90:293-321` issue `enter/exit data`
   UNGUARDED while `scalar.f90:2562-2592` guards them; `bodyforce.f90:472, 531`
   guard the trip kernels so the CPU build runs those SERIALLY. Equivalent
   today only because the CPU build links no OpenMP at all
   (`CMakeLists.txt:211-245`). Pick one convention (guard only `enter/exit
   data` and `update`, let `target` constructs fall back) and apply it.
8. **Adiabatic mode keeps the unmasked `divuse = divu`** while its `conv` is
   masked (`scalar.f90:2987-2993`): the F2 inconsistency the conjugate mode
   fixed, kept deliberately for gated bit-exactness; a uniform scalar is not
   preserved at solid-face cells in adiabatic + advective. A documented trade,
   recorded here so it is not rediscovered.
9. `rans.f90:1221` host twin of the constrained-cell pinning lacks the `kfix`
   branch of the device twin (`:1286`); benign because the device twin runs
   before step 1 -- a reason to delete the host twin (section 3 item 10).
10. Unknown keys: `[scalar]`/`[scalar.N]` only WARN on an unknown key
    (`scalar.f90:1010-1012, 1152, 1176`); with the `config.f90` `case default`
    (section 5) both should error-stop for the same reason.

## 5. Configuration (`config.f90`)

**The `case default` landmine, sized.** Fourteen key blocks have no default
and silently ignore an unknown key: `[grid]` (`:126`), `[grid.x/y/z]`
(`:901`), `[flow]` (`:149`), `[time]` (`:190`), `[output]` (`:211`),
`[restart]` (`:228`), `[pressure]` (`:233`), `[ibm]` (`:275`), `[blocks]`
(`:338`), `[force]` (`:387`), `[turbulence]` (`:719`), `[les]` (`:875`),
`[mpi]` (`:703`), and `[case]` (`:427`, `if (key == "file")` only -- `name`
and the case's own keys are read by `flow_case.f90`, so a default here must
accept them). Three more only WARN: `[rans]` (`:864`), `[boundary]` (`:975,
992`, after trying the `<dir>_<side>_<var>_<field>` parse) and
`[scalar]`/`[scalar.N]` (`scalar.f90:1008`). The outer section `select`
(`:124`) has a default that routes `[scalar.N]` only, so an unknown SECTION
is dropped too; `[case.<name>]` sections ride that silent path and must be
whitelisted when it error-stops. The same class, one level down: seven VALUE
parsers warn and keep the previous value (`read_bool :1105`, `read_int :1360`,
`read_real :1431`, `read_bc_type :1185`, `read_grid_distribution :1217`,
`read_turbulence_model :1239`, `read_les_model :1263`), the `[pressure]`
`solver`/`accel` value selects warn (`accel = chebychev` runs plain Jacobi),
and `flow_case.f90:35` falls back to `generic` on an unknown `[case] name`.
The fix lands ONCE if `apply_config_value` (362 lines) is first split into one
routine per section like `apply_rans_value`; then every ini under
`tutorials/` and `validation/` goes through `moby_prepare` to catch the keys
that existed by accident (the handout's hour).

**A dead key that looks live**: `[ibm] amp_z` is parsed, documented
(`configuration.md:197`) and echoed into the case file (`config.f90:548`),
but `wavy_wall_height` (`ibm.f90:546-549`) uses `amp_x` for both sine
factors; a user varying it gets an unchanged wall and a "stale case file"
after re-preparing. Either wire it or remove it (MacDonald's roughness has
one amplitude; remove). `[rans] model` is read only by the validator (a
guard, fine). `config_seen_type%size/length/dt/nsteps` are set
(`config.f90:128-144, 926-929, 193, 196`) and never read.

**Restart precedence is half config-authority** (`io.f90:728-812`): the ini
wins for `cflmax pecletmax dtmax t_final re forcing_* niter sor ibm_enabled`
and the BC rows; the FILE always wins, undocumented, for `globalSize`,
`leng`, `isPeriodic` (`:805`), `distribution/stretch/natural_dyw_plus`
(`:779`) -- while `one_sided/subdivided/nodes_file/outer_height` come from the
ini, a mixed-source node line of the B0 class -- plus `step/t_current/cfl`;
`dt` deliberately (`:784`). A `[grid] nx` mismatch therefore surfaces as
"case file grid differs from [grid] nx/ny/nz", blaming the ini for the
restart file's value. Decide the list (grid size, lengths and periodicity
belong to the case file now, so the ini should have to AGREE, not be
overridden) and document it in one place.

**Input echo coverage** (`config.f90:493-560`): complete for everything that
shapes the case file except STL CONTENTS (a re-exported mesh at the same path
passes) -- acceptable, note it in `configuration.md`. `wall_x_end` without
`wall_x_start` is echoed (`:557`) but ignored by the geometry (`ramped`
false): echo and geometry should agree.

**Analytic-wall defaults exist three times** and must agree: `init.f90:166-178`
(dns), `ibm.f90:155-166` `set_ibm_geometry_defaults` (ibm), and the literal
`1.0d-2` in the echo condition `config.f90:554`. The middle one is also the
defect of section 4a below.

### 4a. One real defect: `moby_prepare` classifies analytic geometry against the DEFAULT wall

`prepare.f90:150/157` call `classify_refinement_masks` / `classify_active_mask`
BEFORE `init_ibm` + `set_ibm_geometry` (`:171-174`), and
`classify_active_blocks` (`ibm.f90:611`) and `classify_block_geometry`
(`:771`) each call `set_ibm_geometry_defaults(ibm)` themselves -- while the
indicator they are handed, `isInBody(x, ibm, dns)`, reads its geometry from
that `ibm`. So for `wall_shape = eggcarton`, or any wavy wall with non-default
amplitude/offset/phase/waves, the buried and touched masks describe the 2D
wavy wall of amplitude 0.025 at offset 0.01 while the coefficients describe
the configured wall. (Verified by reading all three routines; the
`set_ibm_geometry` comment at `ibm.f90:170-174` claims the opposite.) The
inline path before step 7 had the same order. Consequence today: NONE --
`rough_jacobi.ini` runs exactly this path (explicit nb, `remove_solid`
default) and buries no block under either wall, the TBL rough case has a
derived nb (no removal), and no egg-carton case uses `refine_body`. It was
never gated: the P0 classify gates run the default wall, the egg-carton gate
checked coefficients only. FIX = the simplification: delete
`set_ibm_geometry_defaults`, call `set_ibm_geometry(ibm, dns)` once in
`prepare_case` before the classification (and in the test drivers that call
`classify_*`). Bit-exact on every existing case (identical values for the
defaults); the gate is a `refine_body` egg-carton prepared before/after
(masks must DIFFER) plus P0 26/26 (masks must not).

## 6. Test drivers and documents

### 6.1 Drivers under `validation/` (76 files)

About 58 are current gates (named in `validation/README.md` or a per-case
README with results dated 2026-09-26 or later). The rest:

| driver | verdict | why |
|---|---|---|
| `prepare/run_gates_step7.sh` | RETIRE | its 7-0 inline reference predates F1 and F7: 18 pass / 26 fail with any current binary; what it covered is the suites + `run_gates_order.sh` |
| `cylinder/run_f7_gates.sh`, `cylinder/run_outlet_mode.sh` | one-off investigation drivers (step 9, the outlet mode); keep with a RETIRED header or move under the README's "how the numbers were produced" | need `~/step8_ref_binaries` and a sibling checkout |
| `channel_interface/core_patch/run_gate5.py` | unclear | untouched since 2026-06-30; points at `../run_reference.sh`, which lives in `developed/` |
| `prepare/run_gates_big.sh` | unclear | sailplane leg only, never re-run since P1b (2026-07-17); needs the retired mobygeom + venv |
| `scalar/compile_nofma.sh` | RETIRE | duplicates `./compile.sh cpu_nofma` / `gpu_nofma` (2026-09-14) |
| `scalar/{env,run}_cetus.sh` | host wrappers, not gates | `run_cetus.sh` hard-codes the `mobydiff.scalar` checkout |
| `conjugate/check_{coco,convergence,st}.py` | study scripts of the ON-HOLD tangential plan (`next_session_tangential.md`), not gates | keep with that plan |

Staleness inside current drivers (each a one-line fix):

- `/home/davide/ibmc/bin/python` is the default venv in `les_ibm/setup.sh:31`,
  `cylinder/setup.sh:18`, `rans_geometry/setup.sh:14`; it does not exist on
  this workstation (`~/ibmc` does); `block_nb/run_gates.sh` and
  `redblack_interface/run_gates.sh` hard-code `/home/ws/xt8786/ibmc`.
- Five drivers default to the `build_*/main` symlink (`freestream`, `iddes`,
  `multilevel_body`, `rans_sst`, `rans_inlet` `run_gates.sh`).
- `scalar/run_bitexact_s3.sh:21` cites `~/s2_ref`, which no longer exists
  (current sets: `~/backflow_ref_binaries`, `~/rough_ref_binaries`).
- Headers that describe removed behaviour: `channel_interface/run_validation.sh`
  ("niter = 6"), `developed/run_reference.sh` and `run_uniform128.sh`
  ("reflux_on/off"), `cylinder/check_cylinder.py` ("penalization C_D").
- `rans_sst/turb180.ini` keeps `niter = 6` under plain Jacobi; the T2 re-measure
  at niter 12 is recorded in `docs/history/next_session_after_step7.md` only, not in
  `rans_sst/README.md`; `iddes/README.md` has no entry after 2026-07-10.
- `prepare/run_gates_order.sh` writes a `cflmax = 0.8` literal into a generated
  ini.

Ini sweep over `tutorials/` and `validation/` (live files; the
`overheadTest/results_*` archive copies excluded): `coeff_file`,
`interface_constant_half`, `momentum_reflux`, `interface_skew`, `MOBY_*` hooks
-- **0 hits**. `[flow] convection` is a live key again (six boundary-layer
inis). **`sor` under a Chebyshev/Jacobi solver, where it is ignored: 60 live
inis.** `cflmax = 0.8`: 85 live inis (the rerun convention is 1.0; change it
when a case is rerun, not in a sweep). ~70 UNTRACKED generated inis
(`.s7_*`, `.bx_*`, `.cyl_*`, `.obl*`) with `coeff_file` litter `validation/`
and `tutorials/min_channel` -- gate leftovers to delete.

### 6.2 The C3 transient control rows

`validation/conjugate/README.md` ("(1) the transient two-material slab") and
`c3_transient.dat` carry three "C2, pointwise (the control)" rows (orders
1.71/1.09, 1.80/1.08, 1.28/0.54) produced by `run_gates_c3.sh transient`'s
`c2` leg with `CBIN=~/c2_ref_binaries/moby_solve_cpu` (commit `c243e56`,
archived 2026-08-28). They are NOT regenerable as written: the current binary
has no pointwise-capacity switch (`vfrac` is unconditional), and the archived
binary cannot read today's template -- at `c243e56` the config has
`[ibm] coeff_file` and no `[case] file`, and with no `case default` it would
silently run the ANALYTIC wavy indicator instead of the slab. Even with a
toggle, F2 (advective default), the sum limiter, F1 and F7 have moved since,
so the control would no longer differ "by one line". Decision for the pass:
**retire the rows and the `c2` leg** (keep the one-sentence finding that the
pointwise capacity costs a full order), or add a test-only pointwise toggle
and re-measure both legs with ONE binary. The first is recommended: the
finding is recorded, the toggle would be dead code in production.

### 6.3 The `docs/next_session_*.md` set (34 files + 4 superseded design notes)

Of the 38: **KEEP 7** -- `cleanup_pass2` (the live handout), `rough_outlet`
(section 2 owed, closed by this session), `redblack_interface` (its section-6
physics/performance gates were never run), `tangential` (ON HOLD, with its
three study scripts), `after_step9` (until its item 5 is folded here; 11
referrers incl. `test_keyorder.f90` and `blocks.f90`), `outlet` and
`body_at_outlet` (the design records of the current outlet; 21 and 10
referrers). **MOVE 29** to `docs/history/` (does not exist yet) and **FOLD 1**
(`boundary_layer`, superseded by the tutorial). Two have zero referrers
(`les_ibm`, `naca_rerun`). The superseded design notes
`jacobi_interface_handout`, `momentum_interface_handout`, `interface_review`,
`corner_reconstruction_strategy` move too (`docs/index.md` lists them).

The move needs REDIRECTS: about 15 source comments and 20+ driver headers
cite handout basenames (`scalar.f90`, `blocks.f90`, `config.f90`, `step.f90`,
`rans.f90`, `test_*.f90`, `partition_analysis.py`, the cht/pipe inis). A
one-line `docs/history/README.md` index plus a path fix in each citing
comment is the whole job; `git mv` keeps the history. CLAUDE.md pointers to
fix at the same time: "read `docs/history/next_session_edges_les.md` first" (a DONE
doc), NEXT -> `after_step8` / `after_step9` (both DONE); the only live NEXT
is `cleanup_pass2`. `next_session_port_finish.md`'s table of the 14 cv_box-less
airfoil inis is still the undecided follow-up it was on 2026-09-25 -- carry
it, do not lose it in the move.

### 6.4 The index docs

`docs/tutorials.md`: three stale statements (the sailplane
`sailplane_ibm_coeff.h5` "committed, still usable" -- retired at step 7-4;
the `prep_blocks.ini -> sailplane_case.h5` flow, contradicted by the
tutorial's own README; "`[ibm]` points at the file", now `[case] file`).
`docs/tools.md`: `prepare_if_missing.sh` and `set_case_file.sh`, both
load-bearing since step 7, are undocumented. `docs/index.md`: its last bullet
calls the `next_session_*` set "working handouts" and lists the superseded
design notes. `developer-guide.md`, `validation.md`: no broken path.

### 6.5 CLAUDE.md itself

187 KB, ~45k tokens, loaded into every session. Its "Active work" log is a
chronological record of 2026-06 to 2026-10 whose every bullet is also in a
handout or a results report; the conventions, the build/run recipe and the
verification rules are under 15 KB. Proposal (a documentation move, for the
user to decide since it changes what every session reads first): keep
conventions + build/run + verification + a 2-page "state of the solver" with
pointers, and move the log verbatim to `docs/history/worklog_2026.md`, where
the handouts moved in 6.3 sit beside it.

## 7. The profile of the current head

HoreKa job 5184720 (`horeka/exchange/submit_profile2.sh`, head `b1f41be`,
2 nodes, 8 x A100-40, nvhpc 25.3, 36 min): the six campaign configs at
1/2/4/8 ranks + `tutorials/min_channel` at 1 and 4, 200 steps each with
`[output] profile = true`; the nsys per-kernel timeline at 1 rank on five
configs (30 steps, `collect_kernels.py`); ncu on every kernel of one step
of `refined_yp82_rect_jacobi` (the `rough_jacobi` ncu leg profiled nothing:
its 3 steps hold fewer than the 450 launches the driver skips -- rerun with
`LAUNCH_SKIP=260`). Raw tables: `horeka/results_job5184720/{phases,kernels,
ncu}.md`. ONE column, no comparison: this is the baseline the pass works
from. The configs run Chebyshev niter 6 (red-black 3), the benchmark
settings, not the 2026-09-26 production 12.

### 7.1 Step and ratios (200 steps)

| config | cells | 1 rank | 2 | 4 | 8 | strong scaling 1->8 |
|---|---|---|---|---|---|---|
| `base_jacobi` | 138.4 M | 548.2 ms | 272.7 | 140.9 | 75.8 | 90 % |
| `rect_jacobi` | 138.4 M | 551.6 | 280.4 | 144.1 | 75.4 | 92 % |
| `rough_jacobi` | 138.4 M | 571.7 | 290.4 | 149.4 | 78.0 | 92 % |
| `refined_yp82_rect_jacobi` | 60.6 M | 246.2 | 127.9 | 68.2 | 38.0 | 81 % |
| `refined_yp82_rect_redblack` | 60.6 M | 216.7 | 113.7 | 61.3 | 35.5 | 76 % |
| `refined_big_rect_jacobi` | 242.2 M | -- | -- | 252.5 | 130.9 | 96 % (4->8) |
| `min_channel` | 0.82 M | 17.1 | -- | 15.1 | -- | 28 % (1->4) |

Block tax `rect/base` 1.006 / 1.028 / 1.023 / 0.994 at 1/2/4/8; the 2:1
machinery (`big/rect` per cell) 1.001 at 4 ranks, 0.992 at 8; red-black
0.880-0.935 of Jacobi at niter 3 vs 6 (not a like-for-like quality); the
roughness costs +3.6 % at 1 rank, +3.5 % at 8 over `rect`. All within the
2026-10-01 numbers after the predictor guard (`base` 75.93 ms at 8 ranks
then, 75.83 now).

### 7.2 Where the step goes (phase buckets)

| phase | base 1 rank | base 8 | rect 8 | refined 8 | red-black 8 | rough 8 | min_channel 1 |
|---|---|---|---|---|---|---|---|
| projection | 67.1 % | 68.7 | 70.8 | 73.0 | 71.2 | 69.2 | 91.7 |
| -- `apply` | 47.0 | 41.1 | 42.4 | 41.2 | 3.8 | 40.8 | 44.7 |
| -- `sweep` | 18.7 | 16.4 | 16.7 | 15.0 | 43.5 | 16.0 | 21.3 |
| -- `phi_exchange` | 0.4 | 5.3 | 6.3 | 9.8 | 10.5 | 6.2 | 14.0 |
| -- `vel_exchange` | 0.2 | 5.0 | 4.4 | 6.0 | 12.5 | 4.3 | 10.6 |
| momentum (incl. body force, copy) | 29.9 | 26.7 | 25.8 | 22.5 | 24.0 | 25.8 | 5.1 |
| exchange, all buckets | 0.4 | 11.3 | 11.4 | 16.5 | 25.2 | 11.1 | 27.1 |
| `mpi_wait` | -- | 4.6 | 2.9 | 3.7 | 6.1 | 2.9 | -- |

### 7.3 Device time per kernel (nsys, 1 rank, A100; share of device-busy time)

`base_jacobi`, 138 M cells, device busy 95.2 % of the 574 ms traced step:

| kernel | launches/step | ms/step | share | ms/launch |
|---|---|---|---|---|
| `jacobi_apply` k2 (`@707`, the face correction) | 18 | 204.1 | **37.3 %** | 11.34 |
| `jacobi_compute_phi` | 18 | 104.3 | 19.1 % | 5.79 |
| `step_momentum` (fused predictor) | 3 | 101.2 | 18.5 % | 33.7 |
| `jacobi_apply` k1 (`@671`, `p += phi`) | 18 | 54.8 | 10.0 % | 3.05 |
| **`step_add_bodyforce_correction`** | 3 | **46.3** | **8.5 %** | 15.4 |
| `step_momentum` qs -> q copy (`@457`) | 3 | 15.8 | 2.9 % | 5.26 |
| `bodyforce_fill_trip_kernel` | 3 | 11.0 | 2.0 % | 3.66 |
| `step_get_timestep_rates` | 1 | 5.3 | 1.0 % | 5.33 |
| everything else (BCs, outlet faces, copies, stats) | ~65 | 3.4 | 0.6 % | 0.03-0.2 |

`rough_jacobi` is the same table plus `update_ibm_mu` 9.7 ms (1.7 %),
`add_penalization_state_correction` 6.5 ms (1.1 %), `compute_rdenom` 5.6 ms
(1.0 %), and the local copies 2x longer per launch (15 M vs 6 M points).
`refined_yp82_rect_jacobi` (60.6 M): apply k2 35.6 %, phi 18.0, momentum
15.9, k1 9.8, body-force correction 8.5, copy 3.1, the same-level local copy
2.2, `interface_correct` (3 kernels x 18) 0.7 % in total -- at 1 rank the
2:1 interface's kernels are 3 % of the device time. `refined_yp82_rect_redblack`
(niter 3): `redblack_sweep` 97.7 ms = **45.7 %**, `zero_phi` **19.8 ms =
9.3 %** (36 full-array launches per step, zeroing what the sweep overwrites),
`copy_local_same_level` 7.9 %. `min_channel` (0.82 M cells, niter 12): the
device is busy **54.6 %** of the 21 ms step -- the rest is launch latency and
host time, ~200 launches per step at 20-100 us each; the exchange copies
alone are 27 % of the device time.

### 7.4 ncu, every kernel of one refined step (registers / occupancy / DRAM % of peak)

| kernel | regs | occupancy % | DRAM % | SM % | doubles/cell vs min |
|---|---|---|---|---|---|
| `jacobi_apply` k2 | 64 | 46.0 | 60.5 | 40.5 | 10.31 / 10 (1.03x) |
| `jacobi_apply` k1 | 59 | 45.5 | 62.5 | 66.2 | 3.08 / 3 |
| `jacobi_compute_phi` | 80 | 35.0 | 60.4 | 52.4 | 5.30 / 5 |
| `step_momentum` predictor | 128 | 23.8 | 34.5 | 40.4 | 16.84 / 16 |
| `step_add_bodyforce_correction` | 56 | 53.9 | **81.6** | 34.4 | 18.20 |
| `step_momentum` copy | 48 | 58.7 | 75.4 | 62.7 | 6.43 / 6 |
| `step_get_timestep_rates` | 74 | 35.5 | 35.0 | 62.9 | 3.38 |
| `copy_local_same_level` | 72 | 42.1 | 51.3 | 12.0 | -- |
| `copy_local_cross_level` | 124 | 23.8 | 9.7 | 37.2 | -- |
| `copy_local_scalar_entries` (cross-level) | 126 | 23.0 | 11.7 | 37.5 | -- |
| `predict_outlet_faces` | 130 | 16.6 | 31.9 | 18.1 | -- (surface) |
| `interface_correct` x3 | 54 | 48-51 | 0.1-43 | 46-59 | -- (surface) |
| `apply_bc` | 78 | 34.7 | 61.5 | 26.8 | -- |

Unchanged from 2026-09-14 where measured then (apply k2 64 regs / 46 %,
phi 80 / 35 %, momentum 128 / 24 %): the register work holds. New: the
body-force correction is a well-behaved streaming kernel (81.6 % of peak
DRAM) doing work that does not exist -- see 7.5.

### 7.5 What the profile says (ranking for section 8)

1. **The trip's body-force correction is 8.5 % of the step on every
   boundary-layer case, and almost all of it is zeros.**
   `add_bodyforce_correction` streams `f`, `qs`, `oldrhs` and `mu` for all
   three components over the whole domain every substage (18 doubles/cell,
   at 82 % of peak -- it is efficient at doing nothing), while the trip
   force is `v` only and non-zero in the active block list `fill_trip`
   already uses (40 of 256 blocks on `rect`; 4 of 16 on the tutorial). The
   2026-09-14 optimisation covered the FILL (now 2 %) and not the
   CORRECTION. Restricting the correction to the force's component and
   block list is bit-exact by construction (the skipped cells add exactly
   0.0) and removes ~7-8 % of the step on `base/rect/rough` and the
   production boundary layer. The general `[force] type = profile|file|
   custom` keeps the full kernel. **The cheapest, largest item of this
   pass.**
2. `jacobi_apply` k2 is 36-37 % of the device time at 1.03x its minimum
   traffic and 46 % occupancy: the remaining lever is fusing k1 into k2
   (one pass over `phi`, ~10 % of apply) and the Chebyshev `cheb_combine`
   into `compute_phi` (section 3, review §6 items 7-8). Nothing else about
   it is wasteful.
3. Red-black: `zero_phi` at 9.3 % of the device time is pure waste (zero
   only the halo planes the sweep does not write, or let the sweep write
   them); `redblack_sweep` at 46 % has never had the table/reciprocal
   hoists and runs six branchy `face_grad_denom` and a divide per cell
   (section 3 item 17): the two together are the red-black lever, worth
   measuring before the "red-black vs Jacobi" ratio is quoted again.
4. The 2:1 interface's own kernels are 3 % of the device time at 1 rank
   and `interface_correct`'s three launches 0.7 %: the "~1 ms/step at 16
   ranks" item is real but small; the per-launch floor matters on
   `min_channel`-sized cases, where the device idles 45 % of the step.
5. `get_timestep_rates` at 1 % with the ghost planes included (section 4
   item 1): ~20 % of that kernel for nothing.
6. `step_momentum` at 128 registers / 24 % occupancy is the known limit;
   the per-component split remains the untested idea.
7. Not yet profiled: the scalar and RANS kernels (no campaign config runs
   them; `min_channel` has none). The split of section 3 item 13 needs its
   own ncu on a scalar case before and after -- the LES/scalar suite cases
   on the workstation GPU suffice for registers and occupancy.

## 8. Performance items, re-ranked

Hypotheses from the read, ranked by the profile of section 7 (ncu reads
`launch__registers_per_thread` and `sm__warps_active`, never a step time
alone). THE PROFILE'S OWN TOP ITEM COMES FIRST: **the trip's body-force
correction kernel runs over the whole domain and all three components every
substage -- 8.5 % of the step on every boundary-layer config -- when the
force is `v` only on the trip's block list** (section 7.5 item 1; bit-exact
by construction). Then, grouped by what they cost.

**Launches that need not happen**

- The per-iteration `apply_scalar_bc(phi)` on outlet cases: 36 surface-list
  launches per step at niter 12 (section 3 item 16).
- `interface_correct`: 3 launches x nIter per substage with `face_grad`
  recomputed per point (CLAUDE.md's "~1 ms/step at 16 ranks"; section 3
  item 5).
- `zero_phi` per colour on red-black with interfaces (`pressure_solver.f90:
  1025`): two full-array launches per iteration zeroing halos the sweep
  overwrites.
- `scalar_finish`'s separate `qs -> q` copy per substage (`scalar.f90:
  3380-3404`): foldable into the velocity copy (`1, NVEL+nScalar`), a copy.
- `rans_substage` (`rans.f90:1552-1558`): 2-4 separate `exchange_scalar_halos`
  + `apply_scalar_bc` per substage where the passive scalars ride ONE batched
  exchange; `rans_wall_yplus` two launches + bc + exchange per substage.
- The exchange: four `av` automatic-array copies + `map(to: av)` per full
  exchange x 39 rounds (`comm.f90:2040, 2110, 2343, 2413`) -- three fixed
  variable sets resident on the device make it zero; the local same-level
  copy folded into the pack kernel (both before the Waitall; the CLAUDE.md
  "cheap concrete win", untaken); `apply_bc` maps an allocatable `vl` per
  call (small, 3 calls/step).

**Work done every substage on static data**

- Static conjugate face coefficients (`scalar.f90:3035-3061`): six
  `conjugate_face_diffusivity` calls (branches + a divide at cuts) per cell
  per conjugate scalar on `phi/ks/rc/invd`, all static -> three face arrays
  per conjugate scalar at init, bit-exact (a pure function of static data);
  the kernel then reads `kface` like the nut path and drops 7 `phi` loads and
  the `cut*` logic. NOT `invCap` (`*invCap` vs `/cap` is not bit-exact).
- `dm = ire/sc%pr(is)` (`scalar.f90:3022`): one fp64 divide per cell per
  scalar; the identical expression at init is bit-exact (the `x/6` lesson).
- `momentum_face_start` x3 per cell in five kernels (section 2).
- The velocity-gradient tensor up to four times per cell per substage under
  IDDES (section 2): store the nine components once, bit-exact.
- The Courant reduction over ghosts: `(nb+2)^3/nb^3` extra cells.

**Kernels never profiled**

- `redblack_sweep` (`pressure_solver.f90:1090-1210`): 26 privates, per-cell
  `iLo/jLo/kLo/colorOffset`, six branchy `face_grad_denom`, one fp64 divide
  per cell -- the whole lever set that gave the Jacobi kernels -14 to -30 %
  (tables, reciprocal) was never applied; its register count is unmeasured.
  Also its `copyOnly` pack/unpack still run the per-point `find_entry`
  binary search the divergence lists removed for Jacobi (`comm.f90:2355-2361,
  2423-2429`).
- `add_eddy_viscosity_correction` (`step.f90:789-971`): 25 privates, 12
  bilinear interpolations, 48 `nut` loads per cell.
- `scalar_transport_kernel` (~100 privates) and `rans_transport_kernel`
  (103): the `jacobi_apply` lesson says the dormant branches cost occupancy;
  the split (section 3 item 13) is the lever.
- `copy_local_cross_level` threads (point x var) with var innermost, the
  layout `copy_local_same_level` measured 25 % slower (`comm.f90:2048-2057`),
  and it is the ~128-register kernel.
- `prt_kays` evaluated 6x per cell per scalar (each face from both sides); a
  cell-centred Pr_t is a numerics change, not bit-exact -- measure first.

**Still open from the campaign** (unchanged): the fused predictor split into
per-component kernels (register lever exhausted; ncu before/after); the fp64
divide in `compute_rdenom` inside body blocks; a reduced divergence round for
red-black (needs its own read-set proof); the unexplained `phi_exchange`
locality shift; plan entry 7-5 only if a case is load-imbalanced.

## 9. Findings ranked, and the order of work

| # | severity | what | where | status |
|---|---|---|---|---|
| F10 | correctness (latent; no shipped case affected) | `moby_prepare` classified analytic geometry against the DEFAULT wavy wall (`classify_*` reset the geometry, and ran before `set_ibm_geometry`) | `prepare.f90:150-174`, `ibm.f90` `set_ibm_geometry_defaults` | **FIXED 2026-10-06** (section 4a; gate 4 of `validation/prepare/run_gates.sh`, `check_eggclass.py`: pre-fix 256 leaves = the default wall's, fixed 464, touch raster == the independent transcription; P0 29/29; every nofma suite `max_abs 0`) |
| F11 | configuration | 14 key blocks without `case default`, 3 that warn, 7 value parsers that keep the old value, unknown sections dropped | `config.f90` (section 5) | the handout's section 3 item 1: split `apply_config_value` per section, error-stop once, run every ini through `moby_prepare` |
| F12 | configuration | `[ibm] amp_z` parsed, documented, echoed, never used | `ibm.f90:546-549` | remove the key (one amplitude) |
| F13 | configuration | restart precedence undocumented and half config-authority (grid size, lengths, periodicity, distribution params from the FILE) | `io.f90:728-812` | decide: the ini must AGREE on what the case file owns; document in `configuration.md` |
| F14 | consistency | Courant reduction over ghosts, Péclet over the interior | `step.f90:1012, 1043` | interior only; not bit-exact where a ghost binds -- gate with `dtmax`-bound cases |
| F15 | robustness (latent) | `MPI_Cart_create(reorder=.true.)` + world-rank/cart-rank mix; unchecked xz key width; int32 window product | `comm.f90:225, 235, 538`; `blocks.f90:1333, 383` | `MPI_Comm_dup`; `key_covers` for xz; a guard on the product |
| F16 | hygiene | three OpenMP guard conventions | section 4 item 7 | one convention, applied once |
| -- | dead code / stale comments | section 1 | | **DONE 2026-10-06** for every item of section 1 that is a comment or an unused symbol (the `read_case_layout` flag, the `[rans]`/`[scalar]` warnings and the `vl` map are simplifications, left) |
| -- | drivers and docs | section 6 | | **DONE 2026-10-06**: `run_gates_step7.sh` and `compile_nofma.sh` retired, the C3 control leg retired (rows kept as the record), venv/reference/header fixes, 31 documents moved to `docs/history/` with 214 citations redirected, `tutorials.md`/`tools.md`/`index.md` fixed, CLAUDE.md's stale NEXT/read-first pointers |

**Order of work after this pass** (the handout's section 3, re-ranked by what
the read found; the profile of section 7 may reorder the performance items):

1. F11, the `case default` (one hour, the top item of the handout as well),
   taking F12 with it.
2. The bit-exact simplifications of section 3 that remove duplication the
   next edit would trip over: items 2-4, 7, 12 (index helpers, the two
   projections' shared tail, `comm.f90`'s one enumeration, constants in one
   place). Gate: production flags (nothing moves) + the nofma suites once at
   the end.
3. F13 and F14 together (both are "which value wins" questions with a
   measurable effect on `dt` or on a restart) -- one session, gated on the
   dtmax-bound cases and a restart pair.
4. The scalar kernel split (section 3 item 13) with its missing mixed-mode
   gate case, then the face-diffusivity helper (item 14): the two largest
   readability wins in the solver, and the lever for whatever the profile
   says about `scalar_transport_kernel`'s occupancy.
5. The performance items in the profile's order (section 7.5): the
   body-force correction on the trip's block list first (8.5 % of the
   boundary-layer step for an afternoon), red-black's `zero_phi` and sweep
   hoists next, then apply k1+k2 / `cheb_combine` fusion.
6. CLAUDE.md's log to `docs/history/worklog_2026.md` (section 6.5) -- the
   user's call.
