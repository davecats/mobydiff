# Handout — consolidate the branches and merge to main

> **STATUS 2026-09-25: BOTH MERGES DONE AND GATED. The `main` fast-forward and
> the branch deletions are NOT done — they are the user's call and are the only
> things left.**
>
> * `origin/boundaryLayer` merged at `965a92e`. §3 resolved as the handout
>   recommended: their CaNS/SIMSON trip maths kept whole, our block-list +
>   spanwise-table optimisation re-applied on top. The envelope re-check it
>   asked for came back usable — the CaNS port keeps the same separable
>   Gaussian and the same `ex < -50` cutoff, so `select_trip_blocks` is
>   unchanged and the printed fraction is 4/16 on the gate case (it differs
>   from the old 40/256 only because the default `trip_x0` moved 15 -> 10).
>   Gate: CPU nofma, trip case with `trip_ts` small enough to redraw ten times
>   inside the run, `max_abs 0` against a binary carrying their bodyforce.f90
>   verbatim, at 1 and 4 ranks; 1 rank == 4 ranks.
> * `origin/scalar` merged at `e3abcb1`. 11 conflicting files as measured. The
>   one non-mechanical resolution was NOT in §3's list and the handout did not
>   see it: **both branches had independently built a per-phase step timer on
>   the same `[output] profile` key.** Kept this branch's `profiling.f90`
>   (three nested profilers; the whole campaign is written in its buckets),
>   deleted the `scalar` branch's six `STEP_PROF_*` chron buckets and the
>   optional `prof` argument to `pressure_projection`, and preserved the one
>   thing theirs had that ours lacked as a new `PROF_SCALAR` bucket.
>
> **§5's premise was wrong in one place, and it matters.** The handout says the
> merge "changes no arithmetic on either side, so it should be bit-exact
> against BOTH parents". It is not, and must not be: `scalar` carries the
> 2026-08-05 fix for the cold-started RANS initial condition (`k` a factor 4
> low on the last plane of every block). So `turb180` / `wf180_y30` / `lam30t`
> are bit-exact against `origin/scalar` and differ from the pre-merge head by
> O(1e-2) — and a `max_abs 0` there would mean the fix had been LOST. The gate
> asserts the difference rather than tolerating it.
>
> **Gates run (GPU, job 5163314, `submit_merge_gate.sh`; nofma every side):**
>
> | leg | reference | result |
> |---|---|---|
> | min_channel 1 + 4 ranks, beltrami_slaby, les_ibm | pre-merge `f0a8fe0` | `max_abs 0` |
> | Pass G, `rect_jacobi` + `refined_yp82` (138 M / 60 M cells, 4 ranks, trip off) | pre-merge `f0a8fe0` | `max_abs 0` |
> | turb180, wf180_y30 | `origin/scalar` | `max_abs 0` |
> | lam30t | `origin/scalar` | `pn` 1.7e-28, `v`/`w` 1e-30, everything else 0 |
> | the same three | pre-merge `f0a8fe0` | DIFFER 0.18 / 0.88 / 2.09 — the RANS IC fix, as required |
>
> Plus, on the CPU: the scalar branch's own cases against `origin/scalar` (the
> scalar fields `max_abs 0` on conduction / prsweep / wave, 1.7e-18 on
> ibmwavy / ibmwavyr); `scalar_test` and `transition_test` ALL PASS; the
> prepare/solve P0 cases (wavy, wavy_refine, wavysolid) solve-from-case-file
> `max_abs 0` vs inline, which is what gates io.f90's 152 merged lines.
>
> The 1e-14 velocity / 1e-13 pressure differences against `origin/scalar` on
> Beltrami-flow cases are NOT the merge: they are this branch's 2026-09-02
> reciprocal-`rdenom` change, which deliberately gave up projection
> bit-exactness at exactly that magnitude. Controlled by comparing the two
> PARENTS directly, which reproduces them to the last digit.
>
> **§4 verdicts, checked rather than assumed.** `claude/jacobi-interface`'s six
> features are confirmed absent here (all seven keys: `kpin_box`, `ktrip_box`,
> `kpin_dwall`, `boostconv`, `steady_tol`, `refine_body_box`,
> `refine_body_levels`) — DO NOT DELETE IT. `claude/blocks` is now safe to
> delete: `validation/channel_interface_mfu` is salvaged into the tree (its
> whole tool chain survives here; carried over with a header saying it has not
> been re-run), and `validation/poiseuille` cannot be salvaged — it drives the
> case through `MOBY_POISEUILLE=1`, one of the 19 hooks deleted in the
> 2026-06-30 cleanup, and `tools/check_poiseuille.py` is gone too. `multiGPU`
> is a strict ancestor as stated. The stray local `origin` branch of §7 no
> longer exists.
>
> **STILL OPEN:** the `main` fast-forward; deleting `multiGPU`, `claude/blocks`
> and `bench/rdenom-always`; removing the `moby-2to1-{always,divlist,head4,
> premerge,scalarref}` worktrees; and the campaign-matrix re-run
> (`submit_matrix4.sh`, 4 nodes) that §5 asks for before any ratio in
> `results_horeka_2026-09-25.md` is quoted again.

Written 2026-09-25 at `c876c75`, from the session that finished the divergence
and `rdenom` optimisation campaign. **Everything below is measured, not
estimated: the conflict counts come from trial merges that were run and then
thrown away.**

## 0 — The branch map, measured

Divergence against `optimiseBlockRefinement_parentBoundaryLayer` (this branch),
after `git fetch --all --prune` on 2026-09-25:

| branch | tip | commits it has that we lack | verdict |
|---|---|---|---|
| `origin/main` | `7aa1c7b` 2026-05-18 | **0** — strict ancestor | merge target; the final step is a **fast-forward** |
| `origin/boundaryLayer` | `c7b4ef0` 2026-09-23 | **2** | merge — one conflict, and it is the interesting one (§3) |
| `origin/scalar` | `be1451e` 2026-09-21 | **59** | merge — 11 conflicting files, 230 clean (§2) |
| `origin/multiGPU` | `b51d398` 2026-06-09 | **0** — strict ancestor | **DELETE, fully contained** |
| `origin/claude/blocks` | `4962b10` 2026-06-22 | 51 | **almost certainly delete** — see §4 |
| `origin/claude/jacobi-interface` | `304765b` 2026-08-27 | 43 | **DO NOT DELETE YET** — carries unmerged features (§4) |

**Commit counts do not measure feature divergence, and this matters here.**
`convection = skew` is implemented on our branch, yet its originating commit
`57bd1e3` shows as "one we lack" — the content arrived by another route. Judge
by content, never by `git rev-list --count`.

## 1 — Recommended order

1. **`origin/boundaryLayer`** first. Two commits, one conflict, and it is a
   physics change to the trip that the performance work touched (§3). Settle it
   while it is the only thing in flight.
2. **`origin/scalar`** second. Much bigger, but the conflicts are small and
   mechanical (§2), and none of them is a correctness argument.
3. **`main`** last — a fast-forward, since `main` is a strict ancestor of this
   branch and stays one after both merges.

Do NOT merge into `main` first and then bring the features in: that turns the
fast-forward into a three-way merge for no gain.

## 2 — `scalar`: the measured conflict surface

A trial merge (`git merge --no-commit --no-ff origin/scalar`) gives **11
conflicting files against 230 that merge cleanly**:

| file | hunks | ~lines | note |
|---|---|---|---|
| `CLAUDE.md` | 1 | **821** | prose, not code — the biggest single item, and purely editorial |
| `src/moby_solve.f90` | 8 | 114 | the driver: `scalar` adds its calls, we added `set_ibm_geometry` |
| `src/modules/pressure_solver.f90` | 5 | 54 | |
| `src/modules/ibm.f90` | 1 | 35 | |
| `src/modules/comm.f90` | 1 | **30** | tiny, despite 905 lines of our churn — the exchange rewrite and the scalar halos barely overlap |
| `.gitignore` | 1 | 40 | |
| `config.f90`, `init.f90`, `chron.f90`, `moby_prepare.f90`, `docs/next_session_profiling.md` | 1 each | 5–20 | trivial |

**~270 lines of real code conflict across 8 source files.** `scalar.f90` (3376
lines), `scalar_stats.f90` (943), `test_scalar.f90` (398), `step.f90` and
`io.f90` all arrive clean.

`CLAUDE.md` deserves its own pass rather than a hurried resolution: both sides
grew long, carefully-worded sections, and the merged file is the project's
primary context. Resolve it last, deliberately, as an editing job.

## 3 — The one decision that is not mechanical

`boundaryLayer`'s `c7b4ef0` **"port the exact CaNS/SIMSON trip into
bodyforce.f90"** rewrites that file: −203/+91 lines. It is the only conflict in
that merge — and `bodyforce.f90` is exactly where this session's trip
optimisation lives (`results_stepwork_2026-09-14.md` §2: the envelope block list
plus the tabulated spanwise sums, worth −88.6 % of that bucket).

So the two sides are:

- **theirs**: a physics correction — the trip is now the exact CaNS/SIMSON form;
- **ours**: a performance optimisation of the trip that theirs replaces.

**Keep their maths and re-apply the optimisation on top.** Then re-derive its
gate, because the optimisation rests on an *envelope*: the old trip's `ex < -50`
cutoff confined it to 40 of 256 blocks on the rank that owns it, and the solver
prints that fraction. **The new trip may have a different envelope, so the
printed fraction must be re-checked, not assumed** — and if the CaNS/SIMSON form
has no equivalent cutoff, the block list does not apply at all and should be
dropped rather than forced.

## 4 — Deletion verdicts, with the evidence

**`multiGPU` — delete.** `git merge-base` with this branch *is* its tip
(`b51d398`), so it is a strict ancestor: every commit is already here.

**`claude/jacobi-interface` — keep for now.** Probing our tree for each feature
its 43 unmerged commits name:

| feature | in `src/`? | verdict |
|---|---|---|
| `convection = skew` | yes | present |
| `[ibm] band_filter` | yes | present |
| `cv_forces` | `tutorials/naca/cv_forces.py` | present (tool, not Fortran) |
| `[rans] kpin_box` | **no** (only named in a doc) | **missing** |
| `[rans] ktrip_box` | **no** | **missing** |
| `[rans] kpin_dwall` | **no** | **missing** |
| `[rans] boostconv` | **no** | **missing** |
| `[case.airfoil] steady_tol` | **no** | **missing** |
| `refine_body_box` / `refine_body_levels` | **no** | **missing** |

Six genuinely unmerged features, all RANS/airfoil-side: OpenFOAM-matching forced
transition, a surface-following laminar band, a steady-state residual
accelerator, a steady-state stop criterion, and refinement-box controls. Plus
the 2026-08 `naca/rans` tutorial state. **Decide deliberately whether those are
wanted before deleting the branch**; they are not reachable from anywhere else.

**`claude/blocks` — almost certainly delete, but check two cases first.** Its 51
unmerged commits are June work on the red-black-coupled projection that
`claude/jacobi-interface` was explicitly forked to *replace*, plus GPU exchange
optimisations that the `map(to: c)` / index-list rewrite has entirely
superseded. Only two artefacts are unique:
`validation/channel_interface_mfu` and `validation/poiseuille` (a Poiseuille
case with wall blowing/suction). Our `validation/freestream` already gates an
in/outflow Poiseuille, and we hold 15 validation directories it does not. **Look
at those two, salvage if they cover something `freestream` and
`channel_interface` do not, then delete.**

## 5 — Gates

The merge changes no arithmetic on either side, so it should be **bit-exact
against BOTH parents**, each on its own cases:

- vs this branch: the 7-case suite plus Pass G on both production cases,
  `max_abs 0` at production flags (`horeka/exchange/run_mapgate.sh`,
  `submit_rough.sh` for the body case).
- vs `origin/scalar`: its own scalar/CHT validation cases — `tutorials/cht/`,
  the conjugate pipe case, `src/test_scalar.f90`.
- vs `origin/boundaryLayer`: the turbulent boundary layer tutorial, **after**
  §3 is settled. This one is NOT bit-exact by construction if the trip maths
  changed; it needs a physics check against the CaNS/SIMSON reference instead.

Then re-run the campaign matrix (`horeka/exchange/submit_matrix4.sh`) on the
merged head. Nothing published in `results_horeka_2026-09-25.md` survives a
merge untested.

## 6 — A correction to carry forward

**Passive scalars DO exist — on `scalar`.** When the rough-wall benchmark was
designed (2026-09-19) I told the user the solver had no generic scalar
transport and that adding it would be its own multi-session track. That was true
of *this branch* and wrong as a statement about the project: `scalar` carries
`src/modules/scalar.f90` (3376 lines, "Passive scalars"), `scalar_stats.f90`,
conjugate heat transfer at the immersed interface (`SC_IBM_CONJUGATE`), a
`[scalar]` / `[scalar.N]` config surface, and validation against a body-fitted
NekRS DNS.

**Once `scalar` is merged, the MacDonald & Hutchins rough-wall
forced-convection benchmark needs no new physics** — only a case file combining
`[ibm] wall_shape = eggcarton` with a `[scalar.N]` section. That is the natural
follow-up to `configs/rough_jacobi.ini`, and it was wrongly written off as
expensive.

## 7 — Housekeeping

- A stray **local** branch named `origin` points at `7aa1c7b`. It is not a
  remote; delete it (`git branch -d origin`) so `origin/...` refs stop being
  ambiguous.
- `bench/rdenom-always` is a benchmark reference, not for merging. Delete it
  once `results_horeka_2026-09-25.md` is accepted.
- Worktrees currently live: `moby-2to1-{code,always,divlist,head4}`. Only `code`
  outlives this work.

## 8 — After the merge

The user's stated plan: new performance tests, then simplification sessions.
Both want the merged head, and the simplification pass in particular wants
`scalar` in, since `scalar.f90` is the largest module in the project and has
never been through one.
