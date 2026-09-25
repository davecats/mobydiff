# Handout — finish the jacobi-interface port, then merge it to main

Written 2026-09-25 on branch `port/jacobi-interface-features` (off `main` at
`68e8f16`). The previous session consolidated every branch into `main`
(`docs/next_session_merge_to_main.md`); this one ported four features out of
`claude/jacobi-interface` and dropped two. **Everything below is measured.**

## 0 — Where this stands

| commit | what |
|---|---|
| `b1798f2` | S3 skew lockdown + the `[scalar] convection` split |
| `1e97617` | `[blocks] refine_body_levels` / `refine_body_box` |
| `877e3ab` | `[rans] kpin_box` / `ktrip_box` |
| `e95ea37` | runtime CV forces + `[case.airfoil] steady_tol`, `tutorials/naca/rans`, docs |

CPU dormancy against the pre-port binary, completed after `e95ea37` was written
(that commit says these were still running): `min_channel` 4 ranks,
`beltrami_slaby`, `turb180`, `lam30t`, `conduction` -- all `max_abs 0`, and
`conduction`'s scalar field `s1` separately `max_abs 0`. That last check was
needed because `tools/h5maxdiff`'s default dataset list is velocity/pressure
plus the RANS scalars only: **a scalar case gated with no dataset arguments
silently compares four datasets and never looks at the scalar.** Name the
scalars explicitly. `les_ibm` remains unrun on CPU.

DROPPED on the user's instruction and the branch's own evidence: `boostconv`
(V1 negative on turb180 — best configuration 2.7x SLOWER than plain marching;
V2's win invalidated, the recombination suppressed the ktrip strip) and
`[rans] kpin_dwall` (patched an instability that skew convection fixed at its
root; no validated ini uses it). Both remain on `claude/jacobi-interface`.

Every gate below is CPU, `-Mnofma`, on the login node. **The GPU has not been
touched at all** — that is task 1.

## 1 — What is OWED, in order

1. **The GPU gate. Nothing here has run on a GPU.** The cc80 build compiles
   (that is all it proves). Run the equivalent of the four CPU gates plus the
   standard suite under `gpu_nofma`, against a binary built from `68e8f16`
   (`main` before this branch). Reuse `submit_merge_gate.sh` as the shape — it
   already knows how to build two references and run `run_mapgate.sh` — but its
   reference logic is merge-specific, so write a sibling rather than editing it.
2. **`les_ibm` dormancy**, the one case never attempted on CPU (too slow
   there). Against `68e8f16`, it must be `max_abs 0`. The other legs are DONE
   (see §0).
3. **`cv_box` for validation/naca0012, validation/sd7003, tutorials/naca.**
   See §3 — the largest piece of judgement left.
4. Then merge to `main` (fast-forward if nothing else has landed).

## 2 — The traps this port walked into, so the next one does not

**A shared config key can mean two different things.** `[flow] convection` was
read by BOTH `step.f90` and `scalar.f90`, and the two operators are NOT the
same: momentum subtracts HALF of `s (div u)` (skew-symmetric, energy-neutral),
the scalar subtracts the FULL amount (advective, uniform-preserving but no
longer exactly conservative — the property `validation/scalar/conserve.ini`
gates). Hardwiring the momentum form would have silently changed the scalar's.
The scalar now has `[scalar] convection`. **Before removing any key, grep for
every reader, not just the one the feature belongs to.**

**A passing gate can be vacuous.** The first scalar legs of the skew gate
(`conduction`, `prsweep`, `wave`) passed at `max_abs 0` — and would have passed
with the operator deleted, because those cases have no advection. Only
`ibmwavy` discriminates, and it was checked by confirming the new key moves
`theta` by 4.1e-10 between its two settings. **For every "no change" gate, show
the knob actually does something on that case.**

**The same applies to positive gates.** The first `refine_body_box` test used a
box spanning the whole domain and reproduced the uncapped leaf count exactly —
a correct result that discriminated nothing. The half-domain box (240 leaves,
128 refined, between the uncapped 256 and the capped 0) is the real test.

**A branch's own tree can be broken in a way its commit messages do not
mention.** `claude/jacobi-interface` makes a missing `cv_box` an `error stop`,
which means eighteen of its own inis cannot start. Nothing in its CLAUDE.md,
which documents the CV work at length, says so. **Run the other cases before
believing a feature is complete.**

## 3 — The one real decision left: cv_box for the legacy airfoil cases

The CV budget replaced the penalization integral, and a box has no sane
default — the case does not know where the body is. Eighteen inis
(`validation/naca0012/*`, `validation/sd7003/*`, `tutorials/naca/{naca,c10,
b11}_base.ini`, `validation/scalar/cylheat.ini`) were written for the old
statistic. Right now they start, print a loud warning and run with force
sampling off, which keeps their flow-field / Cp / transition gates working.

**Do not simply add boxes.** The budget is sensitive to two things the branch
records the hard way: `p_inf` must be subtracted PER FACE (a first attempt at a
signed-area form returned C_D = 261), and a border may legitimately cross a 2:1
interface but is snapped to a face of the coarsest level it crosses. A box
dropped in without re-validation gives a plausible wrong number.

The procedure that IS validated, from `validation/cylinder/`: pick the box,
then check it two ways — uniform flow through it must integrate to EXACTLY 0
when the borders are block-interior (`empty.ini` gates that; it reads 0 of 300
entries here), and the result must be reproduced by the offline
`postProcess/cv_forces.py` on a snapshot. `tutorials/naca/rans/c11_aoa5.ini`
shows the shape for a chord-1 body at (50, 48): `cv_box = 48.5 52.5 46.5 49.5`.

**And note the published-numbers problem**: those cases' READMEs quote forces
from the PENALIZATION integral. Re-enabling forces changes the statistic, so
the numbers must be re-measured, not just re-run. The branch flagged the same
thing for its own naca/sd7003/cylinder figures.

## 4 — The pn caveat that bites the CV budget

The statistic reads the STORED pressure. From a converged `p` it holds (2000
steps at `niter = 6` kept C_D to ±0.001), but a 20k-step `niter = 6` run
poisons it — the drift mode is spatially varying, so no reference subtraction
rescues it. Recipe in `validation/cylinder/README.md`: zero `pn` in a copy of
the converged restart and rerun ~300 steps at `niter = 60`. This is the same
family as the A2 velocity-neutral `pn` drift already in CLAUDE.md.

## 5 — What is still on claude/jacobi-interface

After this port: the `boostconv` module and its `[rans] boostconv*` keys, the
`[rans] kpin_dwall` band, the naca LES kickoff (`tutorials/naca/LES`,
`atzori.pdf`), the two BoostConv papers, and the C10/C11 analysis history
(`docs/next_session_naca_re4e5.md` came over; the per-commit dossiers did not).
**Do not delete that branch yet** — nothing here supersedes the LES kickoff.

## 5b — The scalar convection form is now selectable, and UNMEASURED

`[scalar] convection = divergence (default) | skew | advective` — all three are
`conv_div − f·s·(div u)|stencil` with f = 0, 1/2, 1:

| f | name | what it buys | what it gives up |
|---|---|---|---|
| 0 | divergence | exact global conservation (`conserve.ini` gates it to round-off) | a uniform scalar is not preserved when `div u` is only projection-small |
| 1/2 | skew | the momentum kernel's form: neutral in `∑s²` for any advecting field | leaves `½ s div u` on a uniform field |
| 1 | advective | `u·∇s`: a uniform scalar preserved exactly for any advecting field | conservation |

**WHICH IS BEST FOR A SCALAR IS OPEN.** Every gate in `validation/scalar` runs
the default divergence form — `validation/scalar/README.md` says so in as many
words ("the gates above all run the DEFAULT divergence form, so neither is
exercised by them"), and that was still true when f = 1/2 was added. The scalar
plan's item 8 RECOMMENDS a non-divergence form for production (`niter = 6`
leaves `div u` only projection-small) but by analogy with the momentum finding,
not from any scalar measurement.

The comparison to run, with the instrumentation that already exists:

1. `run_gates.sh conserve` under each form — it reports `∫s dV` drift directly
   (divergence gives −6.838e-19; the other two must degrade it, and BY HOW MUCH
   is the number that decides this).
2. `run_gates.sh det` under each form — it reports the 2:1-interface
   conservation residual (1.047e-05 under divergence) AND the rank/GPU identity.
3. A uniform-scalar case under each form: f = 1 must hold it EXACTLY, f = 1/2
   must leave `½ s div u`, f = 0 more. `run_gates.sh uniform` is that case.
4. Something advecting and non-trivial for accuracy — `ibmwavy`, or the
   turbulent channel scalar campaign.

Sanity check already done: on `ibmwavy` the three forms are distinct and the
f = 1/2 result sits at the MIDPOINT of the other two (divergence→advective
4.134e-10, divergence→skew 2.067e-10), which is what f = 0, 1/2, 1 predicts and
confirms the factor is wired correctly. Default is bit-exact vs the pre-change
binary on conduction / prsweep / ibmwavy (`max_abs 0` including the scalar
fields), so nothing published moves until someone sets the key.

## 6 — Unrelated but still owed from the previous session

The campaign matrix re-run (`submit_matrix4.sh`, 4 nodes). Nothing in
`results_horeka_2026-09-25.md` has been re-measured since the merges, and the
skew lockdown makes it worse: skew costs ~+3.4 % s/step on the branch's
measurement, and the campaign configs already set `convection = skew`
explicitly, so their numbers should be unaffected — **but that is an inference,
not a measurement, and it is exactly the kind that has been wrong before.**
