# The tangential-term escalation (route B) — implementation plan

STATUS: **stage 0 DONE (kill gate 1 PASSED), stage 1 DONE (gate 2 = tier 2).**
Written 2026-09-14 on branch `scalar`; stages 0-1 the same day, all in Python,
no solver and no Fortran. Measurements and tables in
`validation/conjugate/README.md` ("Stage 0 of the route-B escalation" and
"Stage 1"); `./run_gates_c2.sh stzero` re-runs the convergence gate.

**Where it stands.** The one-sided `s_t` converges on a curved interface
(order 0.83-1.08 against the shipped estimator's -0.05 to 0.37) and is
0.38-0.95 % across six decades of contrast against the shipped 1.65-2461 %.
With `k_area` it makes the scheme EXACT on a plane at every ratio and contrast
(6-9 orders below the baseline) and 7-63× better on a curved interface at
moderate contrast, improving with refinement. It is 1.4-11× WORSE at high
contrast on a curved interface, and `s_t` is no longer the cause -- the same
`s_t` with the `k_loc` multiplier sits at baseline parity. **The remaining
blocker is the MULTIPLIER on a curved high-contrast interface, and stages 2-5
below are on hold until it is resolved.** Do not start the Fortran work first:
its cost is identical either way and the gating quantity is the unknown one.

Prerequisite reading, in order: `validation/conjugate/README.md` sections
"C2 — measuring the tangential term", "the way out, measured" and "…and the
same fix on a CURVED interface, which is where it stops";
`docs/conjugate/conjugate_ibm_asbuilt.tex` §6 and §8; the DEVIATION comment in
`scalar.f90 conjugate_tangential`.

## What this is, in one line

The C1 baseline's cut-face flux is `F = a q_n + K s_t` where the exact flux is
`a q_n` alone. This plan removes the `K s_t` term by getting **both** factors
right. It is the only route that does not change the stencil, so it survives
the planned move to implicit direction-split diffusion unchanged.

## What is already known, and why this is a decide-first plan

Three measurements bound the outcome before any code is written:

| finding | number | source |
|---|---|---|
| `k_area` is the right multiplier | cut-cell \|div\| 1.3e2 → **1.33e-8** at r = ∞ | README "the way out" |
| …but alone it does not earn default-on | crossover r\* 0.027 → **0.065**, still above DNS 1e-2 | same |
| `s_t` on a curved interface does not converge | **43 / 54 / 46 %** at h = 1/64, 1/128, 1/256 | README table (b) |
| at high contrast there is no signal to estimate | exact `s_t` falls 90× from κ_s = 10 to 10³; the estimate's error does not fall at all | same |

So **`s_t` is the binding constraint and the plan must attack it first.** The
whole decision is reachable in pure Python on harnesses that already exist and
need no solver — `check_oblique.py` and `check_cylinder.py` rebuild `φ` from
the case file and evaluate the scheme's face quantities on the analytic field,
and `check_cylinder.py dipole` already scores `s_t` against the **exact**
`s_t_exact`. Stages 0–1 are therefore hours, and they carry kill gates. Do not
start stage 2 before they pass.

## The design idea: a one-sided `s_t` that fits in the EXISTING halo

The shipped estimator differences `T` across the interface and then de-biases
the straddle in closed form. That de-bias is derived for a PLANE, which is
exactly why it does not converge on a curved one. The fix is to not straddle:
estimate `s_t` from cells on ONE side only.

**The lever that makes this well posed is that `s_t` is continuous.** With
`∇T = ∇ₜT + (∂ₙT) n`,

```
s_t = e_d·∇T − n_d (n·∇T) = e_d·∇ₜT
```

and `T` is continuous across Γ, so its surface gradient is single-valued — only
`∂ₙT` jumps. Both sides estimate the *same* number, which gives a free
consistency indicator and, at high contrast, self-limiting behaviour: a
one-sided fit inside a nearly isothermal solid returns ≈ 0, which is the
correct answer the shipped estimator cannot produce.

**`conjugate_tangential`'s DEVIATION comment says same-side least squares needs
a two-deep halo. That is true of a full 3D fit and avoidable here.** With one
material layer available per side you can measure the two coordinate-tangential
derivatives `t1^σ, t2^σ` within that layer (no straddle, margin 1 — the stencil
`face_measure` already trims for), and close the missing normal component with
the face's own normal flux, `∂ₙT^σ = q_n/κ_σ`:

```
g_d^σ = (q_n/κ_σ − n_1 t1^σ − n_2 t2^σ)/n_d
s_t   = g_d^σ − n_d q_n/κ_σ
      = (q_n/κ_σ)(1 − n_d²)/n_d − (n_1 t1^σ + n_2 t2^σ)/n_d
```

`q_n` depends on `s_t` through the face balance, so this is linear in `s_t` and
solves in closed form — the same structure as the shipped
`st = (st − c·gtd)/(1 − c)`, and it reuses the existing `MIN_COSINE` grazing
guard for the `1/n_d`.

**Variant B, worth testing because it needs no face balance at all.** Demanding
that the two sides agree eliminates `q_n` directly:

```
q_n = [ n_1 (t1^f − t1^s) + n_2 (t2^f − t2^s) ] / [ (1 − n_d²)(1/κ_f − 1/κ_s) ]
```

i.e. the normal flux read off the JUMP in tangential derivatives. Elegant, and
degenerate exactly where it should be (κ_f = κ_s carries no information;
n_d → ±1 is grid-aligned, where `s_t ≡ 0` anyway) — but ill-conditioned NEAR
those degeneracies, so it is a candidate to measure, not to assume.

**Availability and fallback.** A side's tangential neighbours must be in the
same material. Use `φ > 0` / `φ < 0` masks: two available → centred difference
(O(h²)); one → one-sided (O(h)); none → fall back to the shipped estimator, or
to `s_t = 0`. Thin bodies and high curvature are where this bites, and the
fallback must be counted and reported, not silent.

## Stages

### Stage 0 — the estimator, in Python. KILL GATE 1.

`face_measure` computes on full arrays and masks by `cut` at the end, so a new
estimator is a pure function of its output dict alongside `debias()`; add the
per-side tangential differences and the availability masks to what it returns.

Score against the exact reference on both geometries, h = 1/64, 1/128, 1/256,
κ_s = 10 and 10³:

```
./run_gates_c2.sh cylinder     # already prints raw / shipped vs exact rms
```

**GATE: the cylinder `s_t` error must CONVERGE — observed order > 0, i.e. the
43/54/46 % must fall with h.** If it does not, STOP. Nothing downstream can
work, and the plan has cost an afternoon. Report which of variant A / B / both
converges, and the fallback fraction.

### Stage 1 — multiplier + face set, in Python. KILL GATE 2.

`k_area` is already implemented (`check_oblique.py face_area_fraction`, mode
`"area"`). Combine it with the stage-0 estimator and regenerate the two
decision tables: the plane crossover sweep, and the cylinder's ratio-banded
flux error plus cut-cell truncation at both contrasts.

**GATE, three tiers — the outcome decides how it ships:**
- **default-ON**: `r* < 1e-2` on plane AND cylinder, AND corrected ≤ baseline
  at κ_s = 10³ (the "never worse" property; the consistency indicator is the
  mechanism — blend the correction by confidence so it cannot lose).
- **default-OFF, documented regime**: `r*` between 1e-2 and 0.065. Same
  shipping decision C2 already made, but with a better term. Legitimate outcome.
- **abandon**: no improvement on 0.065, or still worse at high contrast.

### Stage 2 — Fortran: `k_area` and the face set. Only if gate 2 passed.

1. `face_area_fraction` in `scalar.f90` — the 2D plane-in-rectangle form,
   sibling of the existing `plane_box_fraction`; unit-test it in
   `test_scalar.f90` beside it, against the Python one.
2. **Extend the conjugate face set from marker-disagreeing to interface-CLIPPED
   faces.** This is the non-obvious half: the two sets differ, and by a lot —
   `k_area` at marker-cut faces only reads 3.92e+01 against 1.33e-08 at every
   clipped face. On a clipped face whose markers agree, `w` is 0 or 1 so
   `k_face = k_loc` and the C1 flux is unchanged; the correction is purely
   additive there. Verify that, do not assume it.
3. Keep the C2 discipline: everything stays behind
   `[scalar.N] tangential_correction`, inert BY CONSTRUCTION when off, not by
   a `+0.0` argument.

### Stage 3 — Fortran: the one-sided `s_t`.

New `declare target` function beside `conjugate_tangential` (keep the old one:
it is the control, and the README's numbers are quoted against it). Gather the
per-side tangential neighbours with `φ`-sign masks, apply the closed-form
closure, count fallbacks.

**HALO CHECKPOINT:** if stage 0 showed the 1-deep closure does not work and a
genuine two-deep stencil is needed, STOP HERE and re-scope. That is the
`comm.f90` halo-depth change shared with the S5b TVD increment — second upwind
cell, the per-dim affine gather maps, every 2:1 transfer — and it is its own
session, not a sub-task of this one.

### Stage 4 — re-gate. This is the expensive stage.

- C1 (`run_gates_c1.sh`): unchanged with the correction off — 1D slab sweep,
  conservation drift, 1 == 4 ranks, CPU == GPU max_abs 0, the 5 config guards.
- C2 (`run_gates_c2.sh`): re-run all, record the new numbers in the README
  tables beside the old ones. The old rows stay — they are the control.
- C3 (`run_gates_c3.sh`): fraction, transient ORDER, nusselt, budget. The
  transient order study (1.99/2.00) must not move.
- **`scalar_conjugate_peclet_rate` with the correction ON.** The correction is
  an added explicit flux divergence and C2 shipped disabled, so its stability
  cost has never been measured. C3 already found the cut-cell Gershgorin factor
  the hard way (two NaN campaigns); do not repeat that — measure the rate,
  sweep `pecletmax` on the worst case, adjust `share` if needed.
- **The missing order study**, which closes the accuracy story either way:
  global solution order for an oblique interface at r > 0, baseline vs
  corrected. `check_oblique.py field` is the harness. This is the one gap
  §6 of the as-built note currently declares openly.
- Bit-exactness: 7-case and 9-case suites vs `~/s5c_ref_binaries` (nofma, CPU
  AND GPU) with `[scalar] count = 0` and with every `ibm_wall /= conjugate`;
  plus every conjugate case with `tangential_correction = false`.
- Determinism: 1 == 4 ranks and CPU == GPU on both geometry paths.

### Stage 5 — decide and document.

Update `conjugate_ibm_asbuilt.tex` §6 (order of accuracy — the plateau either
goes away or gets a number), §8 (the disabled correction), §11 (validation);
`validation/conjugate/README.md`; the CLAUDE.md conjugate bullet.

## Cost, honestly

Stages 0–1 are an afternoon and are the whole decision. Stages 2–3 are a day
each. Stage 4 is the bulk and is unavoidable — it is the same re-gate C3 paid.
If gate 1 fails the plan costs an afternoon, which is the point of ordering it
this way.

## What this does NOT do

It does not make the oblique case second order unconditionally. That is route
C — carrying `T_Γ` as an interface unknown — which removes the tangential
expansion instead of correcting it, and which is hostile to the implicit
direction-split solve planned in [[solver-rewrite-plan]] (it breaks
tridiagonality at every cut cell). Route B is the one that stays compatible.
