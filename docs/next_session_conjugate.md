# Conjugate heat transfer at an immersed interface (branch `scalar`) — strategy

STATUS: **C1 AND C2 ARE DONE AND GATED** (C1 2026-08-27/28 on top of
`e1b5c2c`; C2 2026-08-28 on top of `735d360`, branch `scalar`). C3/C4 are not
started; §10 still describes them. Gates, the commands that produced every
number, and the three time-step findings are in
[`validation/conjugate/README.md`](../validation/conjugate/README.md).
This document is the plan; the derivation, the limit checks, the accuracy
argument and the sketches live in
[`docs/conjugate/conjugate_ibm.tex`](conjugate/conjugate_ibm.tex) (build with
`make` in that directory).

**THE C2 VERDICT — the tangential correction SHIPS DISABLED.** It is
implemented, gated and switched off by default (`[scalar.N]
tangential_correction`), and the numbers below say it should stay off for the
geometries this solver is for. The measurement, not the design, decided:
`e_face` and §3's closed form agree per face to every digit, and the term the
baseline drops costs ~3.6 % rms of the local interface flux at the DNS-like
ratio `|∇_tT|/|∂_nT| ~ 10⁻²` and moderate contrast — but putting it back is
better only above a measured crossover ratio (`r* ≈ 0.027` at κ_s = 10 on a
plane; a LOCAL ratio ≈ 0.1 on a body), it makes the SOLUTION worse even where
it makes the face flux exact, and at high contrast on a curved interface it is
**40–65× worse than dropping the term** — because `|∇_tT| = G|sin t|·2/(1+κ_s)`
vanishes in the isothermal limit, so the term the correction exists to add
disappears while the discrete estimate of it does not, and that noise is then
multiplied by `κ_s`. The INDICATOR, by contrast, ships as a permanent, free
diagnostic: any conjugate run can now report the error it is making.

## C2 — what landed

Two things, both at the SAME cut faces the C1 coefficient lives on, and
neither touching anything else:

1. **The indicator.** `[scalar] indicator_interval = N` reduces
   `e_face = h_d s_t/(T_R − T_L)` over every cut face and prints
   `max|e_face|`, its rms and the face count — at step 0 (on the initial
   field, which is how the manufactured gates use it) and every N steps
   after. It is §3's own expression, so it IS the relative error the
   baseline makes in the interface-normal flux at that face. Two passes, so
   that faces whose denominator vanishes (a cut face carrying no normal
   flux: `e_face = ∞`, meaning nothing) are excluded against a measured
   scale rather than silently dominating.
2. **The correction**, `F += s_t (k_loc − k_face)` behind `[scalar.N]
   tangential_correction` (**default off**), where `k_loc` is the
   conductivity of the material containing the face MIDPOINT — the sign of
   `φ_L + φ_R`, no division. It is added as its own flux divergence AFTER
   the C1 term rather than folded into it, so the fused expression is
   untouched and the feature is inert by CONSTRUCTION when off (the
   bodyforce lesson: not a `+ 0.0` argument).

`s_t` and its indicator go through ONE routine, `conjugate_tangential`, so
the number reported is the number applied — the invariant S4 exists for. The
explicit correction enters `scalar_conjugate_peclet_rate` as its Gershgorin
row sum (measured cost: 26 % of `dt` at κ_s = 10³). No new dataset, no
case-file change, no new stencil beyond the 3³ neighbourhood the 26-neighbour
halo already delivers.

## C2 — measured gate numbers

| §10 gate | measured |
|---|---|
| oblique plane, baseline error vs §3's closed form | agree **per face, to every printed digit** (`1.5080e+00` vs `1.5080e+00`, κ_s = 10, r = 0.1) |
| …vs the ratio `r` | baseline rms **strictly linear in r** and INDEPENDENT of h: `3.6 r` (κ_s = 10), `13 … 56 r` (κ_s = 10³) — the contrast amplifies it |
| …the corrected error | **flat in r**: rms `9.6e-2` (κ_s = 10) / `1.4e-1` (κ_s = 10³) ⇒ crossover `r* = 0.027` / `0.003 … 0.011` |
| …`r = ∞` (pure tangential) | corrected face flux **exact to 1e-11** (the STL floor / h); baseline 36 % rms |
| the solver's `e_face` vs the checker's | `1.0000E+00` max and rms, both, on 404 cut faces (the analytic value is exactly 1) |
| cylindrical shell, baseline | rms `1.61e-2 → 8.38e-3 → 4.09e-3`, **order 1.00** — §8's `O(κ_curv h/a)`, measured |
| …corrected | `4.1e-2 / 8.8e-2 / 5.8e-2`: **2.5–15× worse, and no longer converging** |
| cylinder in a UNIFORM GRADIENT (the curved case WITH a tangential gradient) | position-resolved the correction wins above a local ratio ≈ 0.1 (2.7–47× better); but `s_t` itself is **43–54 % wrong and does not converge**, and at κ_s = 10³ the correction is **40–65× worse** than the baseline (`|∇_tT| ∝ 2/(1+κ_s) → 0`) |
| the `r = ∞` BVP (field error) | baseline `L2 = 3.33e-3`; corrected `8.74e-3` — 2.6× worse, at the ratio most favourable to it and with the corrected face flux exact to 1e-11 |
| …and why (local truncation on the exact field) | cut-cell `|div|` rms **130** baseline vs **174** corrected; the area-weighted multiplier below gives **1.3e-8** |
| time-step penalty of the correction, κ_s = 10³ | `dt` × **0.743** |
| every C1 gate with the correction ON | **`ALL C1 GATES PASS`, 61 checks** — the grid-aligned ones digit for digit (`9.658940e-15`, `7.438494e-15`, `q = 1.332889036988`, rate order 1.000), the oblique one re-measured (conservation drift `6.939e-18`, guards 5/5, the 2:1 precondition both ways) |
| 1 == 4 ranks, CPU == GPU, correction ON | **max_abs 0** on all five datasets — and the indicator itself is bit-identical across both |
| the bit-exactness protocol | **32 case-runs, every one max_abs 0**, CPU and GPU, on the 7-case and 9-case suites |
| the S-suite physics groups | `run_gates{,_s2,_s3,_s4,_s5}.sh` — 91 checks, `rc = 0`, no failures |

## C2 — deviations from this plan

1. **The note's construction 1 is NOT what shipped, and gate 1 is why.**
   §7.3 argues that ordinary central differences suffice because the jump in
   `∇T` is purely normal and the projection removes it. That is right for the
   two TANGENTIAL differences, which weight the two sides `(½, ½)`, and WRONG
   for the arm difference, which weights them `(w, 1−w)`: the raw estimate
   carries a bias `μ n_d (½ − w)(1 − n_d²)` proportional to the jump `μ`, so
   it grows with the conductivity contrast. Measured, it makes the interface
   flux WORSE than dropping the term at every ratio ≤ 0.1 (at r = 0.01: rms
   0.359 against the baseline's 0.036 at κ_s = 10, and 1.74 against 0.129
   at κ_s = 10³). The note's own fallback — same-side least squares — needs the far
   cell's far neighbour, i.e. a TWO-DEEP halo, which is the comm-layer change
   the TVD increment is blocked on. What shipped instead removes the bias in
   CLOSED FORM: estimating `μ` from the face's own normal flux makes `n_d`
   cancel and leaves `s_t = (s_t^raw − c ∇_dT)/(1 − c)` with
   `c = κ_face(1/κ_R − 1/κ_L)(½ − w)(1 − n_d²)`. It needs no wider stencil,
   is inert exactly where it should be (equal materials, `w = ½`, a
   grid-aligned face, a purely tangential field), and is UNCONDITIONALLY well
   posed: whenever `c > 0` the geometry forces `κ_face ≤ 1/w`, hence
   `c ≤ ½`, hence `1 − c ≥ ½` for every material pair and every cut position.
2. **§10's "the baseline stalls near first order" is optimistic, by a full
   order.** §3's own bullet two paragraphs earlier has it right — "set by the
   flow, not by `h` … it does not converge away" — and that is what the gate
   measures. Over a 4× refinement at θ = 30°, κ_s = 10, the baseline's
   relative `q_n` error is `3.559e-02 → 3.774e-02 → 3.586e-02` (apparent order
   −0.08, +0.07) and its face-flux error `3.284e-02 → 3.331e-02 → 3.318e-02`
   (−0.02, +0.01). It cannot converge: `e_face = h_d s_t/(T_R − T_L)` is a
   RATIO OF TWO GRADIENTS, and both scale like `h`. In the usual accounting
   the cut-face flux is therefore **zeroth order** — inconsistent —
   `F_num − F_exact = s_t(k_face − k_loc)` being O(1) absolute; "first order"
   would mean O(h), which it is not. What makes the scheme usable anyway is
   not consistency but CONSERVATION plus sign alternation: `k_loc` is a
   staircase in `w` that flips as the interface sweeps past the face midpoint,
   so neighbouring faces carry errors of opposite sign and they largely cancel
   in the cell balances (the EJIIM Remark 24 regime the plan cites).
   A genuine first order does appear in the suite, but it is a DIFFERENT
   error: the cylinder gate, where `s_t ≡ 0` by symmetry and what is left is
   the curvature error in `w` — order `+0.94, +1.03` (κ_s = 10) and
   `+0.99, +1.00` (κ_s = 10³), i.e. §8's `O(κ_curv h/a)`. NOT MEASURED, and
   the natural C3 item: the order of the SOLUTION error, which needs two
   grids of the `r = ∞` BVP and was priced out here (`dt = 0.4 h²/κ_s` against
   a relaxation time `L²/α_f = 1` is 1e5 steps per leg at n = 64).
3. **§3's `≈ tanθ |∇_tT|/|∂_nT|` is the EQUAL-MATERIAL estimate.** The
   tangential term competes against `a q_n R`, and `R` is cut down by the
   high-conductivity leg, so the real error is amplified by the contrast: at
   `r = 0.01` it is 3.6 % rms at κ_s = 10 but 13–56 % at κ_s = 10³. The
   closed form in §3 is exact; only its `tanθ` reading is not.
4. **A finite-ratio BVP is not available, and that is a property of the
   problem.** Every `r < ∞` manufactured solution has a JUMP in `∇T`, so its
   exact boundary data differs BETWEEN THE TWO MATERIALS on the two domain
   faces the plane cuts, and the solver's per-face rows are constants. The
   finite ratios are therefore measured face by face on the analytic field
   (which needs no boundary condition at all and isolates the scheme), and
   `r = ∞` — the one ratio with no jump — carries the global statement.
   Adding per-point boundary data for a gate was rejected: it is solver
   surface no production case would use.
5. **A `45°` plane on a uniform grid is lattice-degenerate.** Every cut face
   of a direction gets the same `(w, a)`, so `max = rms` and the numbers
   track the offset `mod h`. Worse, with `y₀` ON the lattice a whole diagonal
   of cell centres lands EXACTLY on the plane, `w` collapses to 0 or 1 and
   the measured flux is off by the full contrast. The driver offsets `y₀`;
   the 30° numbers are the ones to read.
6. **Found by C2's geometry, but a C1 property: `pecletmax = 0.4` is
   marginal at an oblique high-contrast interface.** `scalar_conjugate_
   peclet_rate` grants `dt = pecletmax·3C/diag` at a cut cell while
   Gershgorin plus this RK3's real-axis limit allow `1.25 C/diag` — the
   default is at **96 % of the bound**. C1's gates never saw it because
   their interface is grid-aligned. Measured: the baseline at κ_s = 10³ on
   the oblique plane goes to NaN at 0.4 and is stable at 0.2. NOT changed
   here (moving `dt` invalidates every C1 number and is a C3-sized
   decision); recorded, with the recommendation to run such cases at
   `pecletmax ≤ 0.2`.

## C2 — the verdict, and what it rests on

**SHIP DISABLED.** Three independent measurements agree:

- **the crossover.** The correction's own error is flat in `r`, so it beats
  the baseline only above `r* ≈ 0.027` (κ_s = 10) or `0.003 … 0.011`
  (κ_s = 10³). §3's DNS-like ratio `10⁻²` sits AT or BELOW that.
- **curvature.** On a cylinder with a RADIAL field — where `s_t ≡ 0` — the
  correction injects 4–9 % where the baseline has 0.4–1.6 %, and destroys its
  first-order convergence. On a cylinder in a UNIFORM GRADIENT — the curved
  case that does carry a tangential gradient — the correction wins
  position-resolved above a local ratio ≈ 0.1, but its `s_t` is 43–54 % wrong
  and does NOT converge (the de-bias is derived for a plane), and at
  κ_s = 10³ it is **40–65× worse than the baseline**: `|∇_tT| = G|sin t|·
  2/(1 + κ_s)` vanishes in the isothermal limit, so the signal disappears
  while the estimate's error does not, and the difference is amplified by
  `(k_loc − k_face) ≈ κ_s`. **This is the strongest single argument in the
  increment**, and it is one no choice of multiplier can address.
- **the global test.** At `r = ∞`, where the corrected FACE flux is exact to
  1e-11 and the baseline is 36 % off, the corrected FIELD is 2.6× worse.
  That is not a contradiction, it is the lesson: the flux divergence of a
  finite-difference cell balance wants a face-AVERAGED flux, and `k_face`
  (smooth in `w`) approximates that better than the exact midpoint value
  `k_loc`, which is a staircase in `w`. **Making the pointwise flux exact
  makes the divergence worse.** §4 derives `F` as the flux AT THE FACE
  MIDPOINT; that is the quantity the escalation makes exact, and it is not
  the quantity the cell balance needs.

  MEASURED DIRECTLY, so it is a mechanism and not a story
  (`check_oblique.py residual`): the discrete flux divergence of the EXACT
  field at cut cells, whose true value is zero, has rms **130** with the
  baseline and **174** with the correction (`h = 1/64`, κ_s = 10; away from
  the interface both are `7e-12`). The correction makes the pointwise face
  flux exact and the local truncation error 34 % larger. And the same
  measurement points at the fix — see §13.

The last point is the one to carry into C4: a corner model anchored on the
same pointwise-flux premise will inherit the same problem.

## C2 — the way out, if the tangential term is ever wanted

Measured, not proposed on paper (`check_oblique.py residual`, and the table
in `validation/conjugate/README.md`). The three schemes differ in exactly ONE
number — the coefficient multiplying `s_t`, since `F = a q_n + K s_t` with
`a q_n` continuous:

| | `K` | what it is |
|---|---|---|
| C1 baseline | `k_face` | distance-weighted harmonic mean along the arm |
| C2 as shipped | `k_loc` | the material at the face MIDPOINT — a staircase in `w` |
| **the fix** | `k_area = f k_f + (1−f) k_s` | the exact face AVERAGE, `f` = fluid AREA fraction |

`k_area` is what the cell balance wants, in two lines: with exact face
averages the discrete divergence IS `∮ k∇T·n dS = ∫ ∇·(k∇T) dV = 0`.

**It needs no new data** — for a plane, `f` follows in closed form from `φ` at
the face centre (the mean of the two cell values, exact for a plane) and the
in-plane components of the face-centred `∇φ`, both of which `s_t` already
forms. And **one thing must change besides the multiplier**: it applies at
every face the interface CLIPS, not only where the two cell markers disagree
(different sets; the difference is `3.9e+01` vs `1.3e-08` in the residual).

Measured (rms cut-cell `|div|` on the exact field, `h = 1/64`, θ = 30°):

| κ_s, r | `k_face` | `k_loc` | `k_area` |
|---|---|---|---|
| 10, ∞ | 1.30e+02 | 1.74e+02 | **1.33e-08** |
| 10, 1 | 1.30e+02 | 1.72e+02 | **8.51e+00** |
| 10, 0.01 | **1.30e+00** | 7.27e+00 | 8.51e+00 |
| 10³, ∞ | 1.92e+04 | 1.93e+04 | **1.55e-06** |
| 10³, 0.01 | **1.92e+02** | 1.08e+03 | 1.21e+03 |

Three readings, and the third is the honest one:

1. `k_area` **removes the paradox** — exact wherever `s_t` is exact, ~10¹⁰
   down. "Making the flux exact makes the divergence worse" was the wrong
   MULTIPLIER, not a flaw in the escalation idea.
2. At finite `r` it hits a floor INDEPENDENT of `r` (8.5 / 1.2e3), which is
   the discrete `s_t`'s own residual error.
3. So the decision moves but does not flip: that floor still crosses the
   baseline at `r ≈ 0.065`, **above the DNS-like 10⁻²**. `k_area` alone would
   not earn the correction its default-on.

### …and on a CURVED interface the fix stops working, for a reason worth
### keeping

Measured on the **cylinder in a uniform far-field gradient** — harmonic on
both sides, no source (so the exact divergence is zero), and unlike the radial
log solution it carries a tangential gradient whose local ratio sweeps
`0 → ∞` around the body (`check_cylinder.py dipole`):

1. **Position-resolved, the correction pays over most of a body.** The
   crossover is at a LOCAL ratio ≈ 0.1 (κ_s = 10, h = 1/64): baseline wins
   only in the 0–0.1 band (32 of 512 faces); above it the correction is
   2.7–47× better. That is friendlier than the plane's global-ratio sweep.
2. **But `s_t` itself is 43–54 % wrong there, and does not converge**
   (h = 1/64, 1/128, 1/256 → 43 %, 54 %, 46 % of the signal's own rms). The
   de-bias is derived for a PLANE; curvature breaks its model.
3. **So `k_area` helps but no longer fixes**: cut-cell truncation rms
   `1.66e+1` (baseline) / `2.29e+1` (`k_loc`) / **`1.01e+1`** (`k_area`) at
   κ_s = 10, consistently ~1.6× better than the baseline across h.
4. **At κ_s = 10³ both corrections are 40–65× WORSE than the baseline**, and
   the reason is physics: `|∇_tT| = G|sin t|·2/(1 + κ_s) → 0` in the
   isothermal limit, so the term the correction exists to add VANISHES while
   the discrete estimate of it does not (measured: the exact `s_t` rms falls
   90× from κ_s = 10 to 10³, matching `2/(1+κ)`; the estimate's error falls
   not at all). That noise is then multiplied by `(k_loc − k_face) ≈ κ_s`.

**Conclusion: the fix is right and cannot rescue the feature.** `k_area` is
the correct multiplier — provably, and exactly so on a plane — but the binding
constraint is the `s_t` estimate on curved interfaces, and at high contrast
there is no signal left to estimate. Anyone reopening this should start from
`s_t`, not from the multiplier: the de-bias in `conjugate_tangential` corrects
the straddle of the ARM difference only, while the two TANGENTIAL differences
straddle too, and on a curved interface the plane model behind all of it needs
replacing (the same-side least squares of the note's construction 2, which
needs a two-deep halo). C4's corner model inherits the identical premise, so
this is the first thing to settle before attempting it.

**What ships enabled is the indicator.** It costs one reduction at a chosen
interval, it is the honest answer to "how much is this scheme's premise worth
on MY geometry", and its verdict on the wavy C1 case is worth recording:
`max|e_face| ~ 0.9`, `rms ~ 0.3`. That surface is resolved by 1.3–4.5 cells,
so it is a warning about the CASE, not about the scheme.

## C1 — what landed

`[scalar.N] ibm_wall = conjugate` (a THIRD branch, per §5), with `solid_k`
(κ_s), `solid_rhocp` (C_s), `solid_init`, `solid_source` and
`contact_resistance`. The scheme is the one face coefficient of §1:
`conjugate_face_diffusivity` (scalar.f90) is `κ·dm` within one material and
the distance-weighted harmonic mean `dm/(w/κ_L + (1−w)/κ_R + R_c dm/h)`
across a cut face, on `w = φ_L/(φ_L − φ_R)`. Around it: solid cells are real
unknowns (κ, C, pointwise — the fluid-fraction-weighted capacity is C3);
convection HARD-masked on solid and cut faces with the skew term's divergence
built from the SAME masked face velocities; `ν_t` in neither the solid nor a
cut face; no penalization. `sc%phi = ±dwall`, ghost-inclusive, signed by the
cell-centred IBM marker, built at init from the two EXISTING dwall producers.
**No new dataset, no case-file format change** — §8 held.

Everything is dormant without a conjugate scalar, and gated so: `[scalar]
count = 0` and every `ibm_wall /= conjugate` run is **max_abs 0**, CPU and
GPU, on both the 7-case standard suite and the 9-case scalar suite.

## C1 — measured gate numbers

| §10 gate | measured |
|---|---|
| 1D slab, `δ_L/h ∈ {0.05 … 0.95}` × `κ_s ∈ {10⁻², 1, 10, 10³}` | **max\|θ − exact\| ≤ 9.66e-15** over all 28 pairs (exactly 0 for every `κ_s = 1`) |
| …and `w` itself, from the case file alone | **max\|w − w_exact\| = 2.2e-14** (the STL distance's float floor, 72 cut arms per case) |
| cold start reaches the same profile | 2.39e-14 / 7.44e-15 / 4.77e-15 at `κ_s = 10⁻²/1/10`, residual between the last two writes **exactly 0** |
| Peclet limiter over materials (`α_s/α_f = 200`) | bounded with (`max\|θ\| = 0.869`), NaN without |
| capacity irrelevant at steady state (`C_s ∈ {0.5, 1, 8}`) | 8.85e-14 / 4.37e-14 / 4.88e-15 |
| contact resistance (`R_c ∈ {0, 0.25, 4}`) | 6.9e-17 / **0.0** / 4.3e-19 |
| `κ_s → ∞` == `dirichlet`, `κ_s → 0` == `adiabatic` | field difference **3.19e-4 / 3.19e-6 / 3.19e-8** and **2.87e-3 / 2.87e-5 / 2.87e-7**; interface-flux rate **order 1.000** both ways |
| `Σ C·T·ΔV` in an insulated composite box | drift **−6.94e-18**, relative **1.22e-16** (flow ON) |
| 1 == 4 ranks, CPU == GPU | **max_abs 0** on all five datasets, on BOTH geometry paths |
| the §12 config guards | 5/5 hard-rejected |
| the 2:1 precondition | checked at init, not assumed: a hand-placed box across the wall is rejected (224 bad faces); `refine_body` + `keep_buried` runs |

## C1 — deviations from this plan

1. **§11's "no `moby_prepare` change" needed one TRIGGER line.** `moby_prepare`
   computed `dwall` only when a `[rans]` section was present, so a conjugate
   case file came out without `dwall_blocks`. It now also builds them when
   `ibm_wall = conjugate` — the same arrangement §8 asks for in the solver,
   applied to the preprocessor. **No new dataset and no format change**: a
   file prepared with `[rans]` and one prepared with a conjugate scalar are
   the same file. §11's stronger claim ("no `moby_prepare` change") was too
   strong as literally written; the substance of it holds.
2. **The cut test is the MARKER DISAGREEMENT, not `φ_L·φ_R < 0`.** They are
   the same test — the sign of φ *is* the marker — but the marker form cannot
   be defeated by a cell centre lying exactly on the surface. `φ` is built so
   a solid cell is always STRICTLY negative, which makes `φ < 0` an exact
   material test and keeps `uncut ⇔ one material` a theorem rather than a
   hope.
3. **§7's escalation item 1 is wrong as written, and gate 1 proved it.**
   "Take the maximum over materials" of each material's own `α = κ/C` is NOT
   the explicit limit: a cut face carries `k_face` up to `max(κ_L, κ_R)` —
   which IS §7's own not-stiff bound — but it feeds the cell on the **other**
   side, whose capacity belongs to the other material. A fluid cell against a
   `κ_s = 1000` solid sees 1000× the fluid rate even when `α_s = α_f`
   exactly. `scalar_conjugate_peclet_rate` therefore builds the rate from the
   ACTUAL face coefficients after the interface exists (it cannot live in
   `precompute_peclet_rate`, which runs before `φ` does). That is also far
   *less* conservative than `max(κ)/min(C)`: at `w = ½`, `κ_s = 1000` the true
   penalty is 2×, not 1000×.
4. **A cut cell pays a Gershgorin factor the uniform interior never does.**
   ρ ≤ `2 A_ii/C_i`, and the shipped `pecletmax` convention is ~1.9× short of
   that — uniform runs survive only because their extreme modes are never
   excited. A cut cell's row is strongly asymmetric, so its worst mode is
   local and IS attained (this RK3's real-axis limit is 2.5). Measured:
   `(κ_s, C_s) = (0.01, 0.01)` at `w = 0.95` blows up at `pecletmax = 0.3`,
   stable at 0.2; `(1000, 1000)` at `w = 0.80` blows up at 0.4, stable at 0.2.
   The rate is doubled at cut cells only, which makes the nominal 0.4 behave
   as the measured-stable 0.2 in both.
5. **`solid_init` applies on a COLD START only.** On a restart the solid field
   is saved state and comes from the file like every other cell.
6. **`scalar_stats.f90` was extended but is only SMOKE-gated at C1.** Its
   y-face diffusivity and convective mask now take the same conjugate branch
   as the transport kernel, so the rows keep reporting the flux the kernel
   applied (the invariant S4 exists for). C1 checks only that the branch runs
   on CPU and GPU with finite output; the quantitative gate on the conjugate
   flux columns is C3's Nusselt increment, which is where §10 puts it.
7. **The grazing-arm guard is a fixed parameter, not a config key**
   (`CONJ_MIN_COSINE = 0.05` in scalar.f90). §8 asks for "a threshold on `a`";
   making it configurable would ship a knob no gate constrains.

## C1 — a pre-existing defect the gate run surfaced (fixed here)

Running the full §10 protocol turned up a regression in `21a7701` — the
commit this increment sits on — that has nothing to do with conjugate heat
transfer, and it blocked gate 5's `run_gates_s2.sh` / `run_gates_s4.sh`.

`trim_dt_for_final_time` signalled "the remaining time is round-off, not a
step" by setting `dns%dt = 0` and letting the main loop's `dt <= 0` exit
fire. But `dns%dt` is written into every snapshot's metadata and read back by
the restart, so the FINAL snapshot of a `t_final` run became unusable in a
NEW way — `config.f90: time step must be positive`. It traded one bad final
restart (the `1/dt`-amplified `pn`) for another. Measured:

```
relax/turbles_60000.h5   dt = 0.0   t = 29.999999999975433
RT_turbles.h5            dt = 0.0   t = 29.999999999975433
turbles.log              ERROR STOP time step must be positive
```

and because the failing leg deletes `turbles_*.h5` before it dies, S4 went
with it (`s4stats.ini` restarts from one of those snapshots).

`trim_dt_for_final_time` is now a LOGICAL FUNCTION that reports rather than
zeroing. The trajectory is unchanged — the loop exits at the same step and a
genuine final partial step still gets `dt = remaining` — only the recorded
metadata differs, and both bit-exactness drivers force `t_final = 0.0`, so
the routine returns immediately there.

**The commit's own claim was also wrong**, and that is the transferable part:
it says the fix is "inert on all 32 bit-exactness case-runs (every suite case
is nsteps-terminated)". `run_gates_s2.sh`'s `les_legs()` REWRITES `t_final`
in a generated variant, so its `les` and `band` legs are `t_final`-terminated.
A claim of the form "no case exercises this path" has to be checked against
the GENERATED inis, not only the committed ones.

## C1 — one result stronger than the plan predicted

§10's `κ_s → ∞` gate is written as a limit "to the discretisation's own
tolerance", on the expectation of an O(h) floor: the S3 modes put their
effective boundary on the STAIRCASE while the conjugate interface sits at its
true position. **There is no floor.** The difference from the `dirichlet`
twin falls like `1/κ_s` straight to `3.19e-8` — five decades below any `h²` in
that case — and the interface flux approaches `q_Dirichlet` at measured order
1.000. That is §5's algebraic claim confirmed rather than merely approached:
the S3 `dirichlet` mode, with its second-order graded `((d0−d)/d)/d0²`
coefficient, **is** Luchini's λ, and the `κ_s → ∞` limit of the cut-face
coefficient reproduces it. The `κ_s → 0` side is equally clean:
`|q|/κ_s → 4.000000`, the closed form.

This does NOT reopen §5's "one code path" question — re-expressing the two S3
modes through the cut-face arithmetic still could not be bit-exact, which is
what §5 actually turns on. It does mean the limits are a sharper regression
than expected, and worth keeping as one.

**REVISION (2026-08-05): the plan was restructured around a much simpler
baseline.** The scheme is unchanged in substance but the machinery collapsed:
the interface is now evaluated from the **existing `dwall` field**, so the
per-arm cut dataset, the interface normals, the least-squares weights and the
volume-fraction tiles — an entire preprocessing increment — are all gone. The
earlier, heavier formulation survives as the *escalation path* (§3), one
arithmetic term away. Read §0–§3 and you have the whole method.

Prerequisite: **S3** of `docs/next_session_scalar.md` (cell-centred IBM
coefficients, `[scalar.N] ibm_wall = dirichlet | adiabatic`). Conjugate is the
third mode of the same machinery — see §5. **The prerequisite is MET: S0–S5a
are all landed and gated** (S0–S2 2026-08-03; S3, S4 and S5a 2026-08-04, branch
`scalar`, HEAD `5756aa9`). Read the STATUS header of
`docs/next_session_scalar.md` before starting — S3 landed with two deviations
that matter here:

- **`ibm%mu` was NOT given the `VAR_P` extent** (only `ibm%coef` was); the
  scalar penalization factor `mu_s = 1/(1 + dt_gamma·coef_p/Pr)` is formed
  inline in the kernel because it is Prandtl-dependent. A conjugate cut-face
  path that wants a stored per-face coefficient must add its own array.
- **The body heat release cannot be measured as `∫coef_p (s_body − s) dV`** — a
  Dirichlet solid cell holds the body value to the last bit, so the product is
  `1e28 × 0 = 0` and ~37 % of the heat is invisible. The cancellation-free form
  (staircase solid/fluid face flux + graded-cell penalization) is implemented in
  `scalar_stats.f90` and in `validation/scalar/check_scalar_ibm.py surface`.
  A conjugate Nusselt diagnostic (§10, C3) must reuse it, not reinvent it.

S4 also means the statistics machinery C1–C3 need already exists: per-row
accumulators built with the TRANSPORT KERNEL's own face diffusivity, plus the
runtime heat file. A conjugate cut face must extend that same face-flux
expression, or the statistics stop reporting the flux the kernel applied.

Source papers, all in `literature/`:

- **Luchini et al.**, JCP 539 (2025) 114245 — `luchini-ibm-2025.pdf`. The IBM
  the solver already implements (`ibm.f90`), and the near-wall-1D argument (§3.1)
  that the baseline below rests on.
- **Wiegmann & Bube**, SIAM J. Numer. Anal. 37(3) (2000) 827–862 —
  `interface-method.pdf`. EJIIM: jump-corrected differences; supplies the exact
  cut-face flux and the "O(h) truncation on the interface band is enough"
  result (their Remark 24).
- **Cipelli et al.** (2025) — `cipelli-2025.pdf`. The corner correction (COCO);
  gives the scheme a second, independent derivation (§4) and the route to
  conducting sharp corners (increment C4).

---

## 0. Scope — two questions that may remove the problem

Ask both of every target case before implementing anything.

**Can the interface be made grid-aligned?** If the conducting wall is a plane —
a channel with conducting slabs, or the existing `les_ibm` geometry with the
wall moved onto a cell face — then put the interface *on* a face. The weight `w`
below is exactly ½, the face conductivity is the plain harmonic mean, and the
scheme is the textbook variable-coefficient finite-volume discretisation with
**no immersed-boundary involvement and no approximation**: second order, exactly
conservative, exact interface condition. This is how the canonical conjugate-DNS
channel studies are set up and it costs only a grid choice.

**Does the solid need to be solved at all?** If the solid is thin compared with
the thermal penetration depth, effectively isothermal, or quasi-steady on the
fluid time scale, replace it with a Robin condition
`k_f ∂_nT = h_eff(T_Γ − T_ref)` from a 1D solid model. No solid field, no
buried-block constraint, no solid time-step limit, no capacity fraction — and
the physics still enters through the same two dimensionless groups. This is a
modelling decision; make it deliberately rather than by default.

Everything below is for the case where the answer to both is no.

---

## 1. The scheme

At DNS resolution the near-wall field is locally one-dimensional in the
wall-normal direction — Luchini's own premise, transferred from momentum to
temperature. Under it, conjugate heat transfer reduces to **getting one
coefficient right per face**: the effective conductivity of a face whose two
cells are in different materials.

```
φ_L, φ_R = signed distance at the two cell centres     (= ±dwall; sign = solid marker)

if (φ_L·φ_R < 0):                                      ! THE only new branch
    w      = φ_L/(φ_L − φ_R)                           ! level-set fraction of the arm
    k_face = 1/( w/k_L + (1−w)/k_R )                   ! distance-weighted harmonic mean
else:
    k_face = k_L   (= k_R)                             ! today's arithmetic, untouched

F = k_face·(T_R − T_L)/Δ                               ! the existing flux line          (★)
```

**The obliquity lemma** (LaTeX note §6.2) is what makes this work at any
interface orientation. `dwall` stores the *perpendicular* distance, not the
distance along the arm; with `a = n·e_d`,

```
φ_L = a·δ_L ,   φ_R = −a·δ_R   ⇒   w = φ_L/(φ_L − φ_R) = δ_L/h_d
```

— the direction cosine **cancels**. The ordinary level-set zero-crossing
fraction is exactly the fraction the series resistance needs, with no normal
ever computed. Exact for a plane at any angle. Free by-product for §3:
`a = (φ_L − φ_R)/h_d`.

What (★) inherits, without further argument:

- **Exact in 1D** — it is two resistances in series, so exact for a grid-aligned
  interface at *any* cut position, and exact at any orientation whenever the
  local gradient is wall-normal.
- **Conservative** — it is a face flux; `Σ C·T·ΔV` conserved to round-off.
- **Correct limits** — `κ_s → ∞` reproduces Luchini's λ *algebraically*
  (including the `((d0−d)/d)/d0²` form `ibm.f90` stores); `κ_s → 0` gives zero
  flux. See §5.
- **Bit-exact when off** — an uncut face takes the else-branch, which is the
  current kernel line unchanged.

---

## 2. Governing equations

One `T` over fluid *and* solid, with `[T] = 0` and `[k ∂_nT] = 0` at the
interface. Both properties jump, in different combinations: the interface
condition involves the **conductivity** `k`, the transient involves the
**diffusivity** `α = k/(ρc)`. A formulation carrying only `α` cannot represent
the interface condition — get this right first.

Two pointwise ratios, both ≡ 1 in the fluid:

```
κ = k/k_f          C = (ρc)/(ρc)_f
∂T/∂t = (1/C)[ ∇·(κ ∇T)/(Re·Pr) − χ_f ∇·(u T) ]
```

In the fluid this is *exactly* what the current kernel solves — the structural
reason the feature is bit-exact when disabled and untouched away from the body.
Per scalar the solid is two numbers, `κ_s = k_s/k_f` and `C_s = (ρc)_s/(ρc)_f`,
whence `α_s/α_f = κ_s/C_s`.

---

## 3. What (★) costs, how to measure it, and the way back

(★) drops the tangential term of the exact cut-face flux. The error in the
physically meaningful quantity — the interface-normal flux `q_n` — is available
in closed form (LaTeX note §6.4):

```
(q_n^num − q_n)/q_n  =  h_d·s_t/(T_R − T_L − h_d·s_t)  ≈  tanθ·|∇_tT|/|∂_nT|
s_t = e_d·∇T − a(n·∇T)          cosθ = a = n·e_d
```

Properties, all of them load-bearing for the decision to drop it:

- Vanishes identically for a grid-aligned interface and for `κ_s = 1`.
- **Set by the flow, not by `h`** — it does not converge away. (★) is second
  order in the bulk and for normal-dominated interface transport, and first
  order in the local flux at obliquely cut faces carrying a tangential gradient.
  That is the honest statement and it is the whole trade.
- In a DNS thermal boundary layer the normal gradient acts over ~1 wall unit
  while `T_Γ` varies over tens of wall units along the surface, so
  `|∇_tT|/|∂_nT| ~ 10⁻²` and the local error is a few percent at strongly
  oblique faces.
- `∮ ∇_tT dS = 0` over a closed body or periodic wall ⇒ **the surface-mean heat
  flux (Nusselt) is unaffected at leading order**; the error is in the local,
  instantaneous distribution.

**Make it a measurement, not an assumption.** The middle expression is
computable at run time from quantities the solver can form:

```
e_face = h_d·s_t/(T_R − T_L)         evaluated at cut faces only
```

Report `max|e_face|` and its rms. This is a deliverable of C2, not an optional
diagnostic — it converts the premise of §1 into a number.

**Escalation.** If the indicator says the local flux matters, the route back to
the formally second-order scheme is one extra term at the *same* face:

```
F = k_face·(T_R − T_L)/Δ  +  s_t·(k_loc − k_face)                            (★★)
```

and `s_t` is already computed if the indicator is running. Nothing about the
data layout, the face loop, conservation or the limits changes. **Nothing new is
stored either** — see §8.

---

## 4. Why this is defensible (condensed; full argument in the LaTeX note)

**Luchini** gives the Dirichlet IBM: extrapolate the ghost through the first
fluid point and the *known* wall value, substitute back so only the Laplacian
centre weight changes (their λ), and integrate that term implicitly. It cannot
express CHT — their §3.4 says Neumann-type conditions are supported only where
the IBM is not involved, and a conjugate interface is a *transmission*
condition, with `T_Γ` an unknown fixed by the flux balance.

**EJIIM** supplies the missing piece: keep the standard stencil, correct it with
the jumps, closed by `[u] = 0`, `[u_ξ] = (ρ−1)u_ξ⁻` (their Eq. 42). Its cost —
the unknown jump as an extra variable, a Schur complement with GMRES, six-point
one-sided operators — is unaffordable here. **The bridge:** our scalars are
advanced *fully explicitly*, so at each substage every one-sided limit is known
data; the augmented system collapses to a local evaluation. That is what makes
the exact cut-face flux (★★) computable at all.

**COCO** (Cipelli et al.) reaches the same formula from the other side: apply
its "match the local analytic solution" recipe to conduction with two
complementary 180° wedges and the local model is
`T = T_Γ + a_t·s + q_n·d/k(side)` — the same `q_n` on both sides *is* the
matched-normal-derivative condition. Eliminating `T_Γ` returns (★★) with no
Taylor expansion. Corollary: **COCO ≡ EJIIM with the jump data supplied by an
analytic local solution** rather than by differentiating the jump conditions.

**The parameter count** settles why (★) cannot be formally second order at
oblique faces, and why the baseline is a deliberate trade rather than an
oversight:

| problem | local model | free parameters | equations from one arm |
|---|---|---|---|
| velocity, no-slip (Luchini/COCO) | `A·u_S(x)` with `u_S = 0` on Γ | **1** | 1 |
| conjugate interface | `T_Γ + a_t·s + q_n·d/k` | **3** | 2 |

One arm gives two equations for three parameters, so the third (`s_t`) must come
from transverse data — *always*. Two consequences: no arm-local scheme,
harmonic averaging included, can be formally second order on an oblique
conjugate interface; and **no single precomputed λ can exist for CHT**, because
COCO's ratio form requires the boundary condition to pin the local solution up
to one scale (the `κ_s → ∞` Dirichlet limit is the only case where it does).

**Two traps if COCO's form is copied literally** (they also explain two of the
landmines in §12):

- The conduction analogue of `u_S` is `T − T_Γ`, which is zero *on* Γ and
  **changes sign across it**, so a ratio-rescaled coefficient is unbounded and
  sign-indefinite. A negative λ is anti-diffusive and implicitness does not cure
  it — unlike the velocity case, where `u_S > 0` and `λ → ∞` is benign.
- A point-anchored correction gives the two sides of a face different local
  models, so it is a **cell source, not a flux**: the interface creates or
  destroys energy at O(h) per step. **Anchor the model to the face.**

**What COCO genuinely adds** is the case (★) cannot cover: a *conducting sharp
corner* (riblet tip, fin root, conducting trailing edge), where the locally flat
model fails at leading order. That is increment C4.

---

## 5. Three body thermal modes, one implementation

| `[scalar.N] ibm_wall` | limit | `k_face` | flux |
|---|---|---|---|
| `dirichlet` (S3) | `κ_s → ∞` | `k_L/w` | `k_L(T_body−T_L)/δ_L` = Luchini λ |
| `adiabatic` (S3) | `κ_s → 0` | `0` | zero (+ `k_L s_t` if escalated) |
| `conjugate` (this doc) | finite `κ_s` | Eq. (★) | (★) |

1. This document originally said S3's two modes **should be implemented as
   limits of the conjugate cut-face path**, one code path with three
   configurations. **That is now blocked, and deliberately so.** S3 shipped
   `dirichlet` as an inline penalization statement and `adiabatic` as six face
   masks, both gated (solid cell == `ibm_value` at 0.0 on 704–5984 cells;
   `∫s dV` drift 0.0), and §10's own bit-exactness rule requires every
   `ibm_wall /= conjugate` run to stay **max_abs 0** against the pre-increment
   binaries. Re-expressing those two modes through a new cut-face arithmetic
   cannot be bit-exact. So: **`conjugate` is a THIRD branch**, and the κ_s → ∞ /
   κ_s → 0 limits are gates that the third branch must reproduce *to the
   discretisation's tolerance* — which is what §10's C1 gate list already says.
   The "one code path" ideal is available later as a separate, explicitly
   non-bit-exact refactor, if the limits ever come out tight enough to justify
   it.
2. S3's `adiabatic` mode ("mask the six faces", landed as planned) is the
   `s_t = 0` truncation of the third row: it drops the tangential flux, so it is
   only first order for oblique interfaces — correct for grid-aligned walls.
   The masking is symmetric across a face, so the flux form still telescopes and
   the fluid conserves `∫s dV` exactly (measured: drift 0.000e+00 over 200 steps
   with a body present). Recorded in `validation/scalar/README.md`.

---

## 6. The discrete scheme

**Face loop** — a branch on face type inside the existing kernel:

| face type | flux |
|---|---|
| uncut, both cells fluid | unchanged: `D_face (T_R−T_L)/Δ`, `D = 1/(Re·Pr) + ν_t/Pr_t` |
| uncut, both cells solid | `κ_s (T_R−T_L)/(Re·Pr·Δ)`, no eddy part |
| cut (`φ_L·φ_R < 0`) | (★): `w` from the two signed distances, then `k_face` |
| `FACE_CLOSED` | masked, as today |

Fluid–fluid faces are byte-identical to today's arithmetic. The cell update
divides the flux divergence by the local capacity `C`. Structure the cut-face
branch as **one local model per face, always consumed as a face flux** — that is
what keeps conservation exact and what lets C4's wedge model drop in later.

**Convection** must be **hard-masked** on every face whose staggered velocity
node is solid and on every cut face. `u·n = 0` on Γ exactly, so masking is
consistent to second order. Do not rely on the penalized velocity being "small":
unlike a passive scalar inside a Dirichlet body, the solid now carries a real
temperature field that must not be advected.

**Capacity at cut cells**: `C_cell = φ + (1−φ)C_s` with `φ` the fluid volume
fraction, estimated from the same signed distance (`φ ≈ clip(½ + φ_c/h, 0, 1)`,
or the closed-form plane-in-cube expression). Irrelevant at steady state,
required for second-order transients.

**Eddy diffusivity**: at DNS there is none. With a model active it is fluid-side
only (`ibm_aware` zeroes `ν_t` in solid cells, and `ν_t → 0` at a resolved
wall), so cut faces are effectively molecular. `[rans] wall_treatment =
wall_function` with `conjugate` stays a hard config error.

---

## 7. Stability and the time step

**The conjugate interface is not stiff.** Since `δ_L + δ_R = h_d`,

```
R = δ_L/k_L + δ_R/k_R  ≥  h_d/max(k_L,k_R)  >  0     for EVERY cut position
```

The `δ → 0` singularity that forces Luchini's implicit `B(λΔt) =
λΔt/(e^{λΔt}−1)` in the Dirichlet case **does not exist** here: as the interface
approaches one cell centre that sub-segment's resistance vanishes but the other
takes over. Cells keep full volume (embedded boundary, not cut cell), so there
is no small-cell problem either. No exponential integrator at the interface.

**What does bite is the solid diffusivity**: `Δt ≤ h²/(6·max(α_f, α_s))`, and
`α_s/α_f ~ 10²` (metal against water) is ordinary. Escalation, in order:

1. extend `precompute_peclet_rate` / `get_timestep_rates` to take the maximum
   over materials, and accept the cost;
2. point-implicit the solid diagonal with **Luchini's** `B(λΔt)` verbatim — the
   same formula, reused for the solid interior rather than for the wall;
3. implicit solid diffusion on the damped-Jacobi / Chebyshev machinery
   `pressure_solver.f90` already owns.

Start at (1); escalate only on measured evidence. The escalation term (★★) is
an explicit spatial operator whose size grows with contrast, so measure its
time-step cost at `κ_s = 10³` if it is ever switched on. **Measured in C2**:
the correction costs `dt` × 0.743 at `κ_s = 10³` — small, and not what decided
against it. What C2 DID find here is that level (1) as shipped is marginal:
`scalar_conjugate_peclet_rate` grants `dt = pecletmax·3C/diag` at a cut cell
while Gershgorin plus this RK3's real-axis limit allow `1.25 C/diag`, so the
default `pecletmax = 0.4` is at **96 % of the bound** and an oblique interface
at `κ_s = 10³` goes to NaN there (0.2 is stable). C1's gates never saw it
because their interface is grid-aligned. That is C3's to settle.

---

## 8. Data: everything comes from `dwall`

`φ = ±dwall`, sign from the existing solid marker. **`dwall` is already computed
at every cell centre including ghosts, at the `VAR_P` position where the scalar
lives**, for both geometry paths:

| path | source |
|---|---|
| analytic geometry | `fill_body_distance_analytic` (rans.f90), from `isInBody` alone |
| file geometry | `dwall_blocks`, written by default by `moby_prepare`, read by `read_dwall_blocks` |

per leaf at that leaf's own level, exact to `[rans] dwall_tol`, and already
validated (the T1 gates measured it against closed-form references to 1e-11).

| quantity | baseline (★) | escalation (★★) |
|---|---|---|
| cut fraction `w` | `φ_L/(φ_L − φ_R)` | same |
| direction cosine `a` | not used | `(φ_L − φ_R)/h_d` |
| normal `n` | not used | `∇φ` (unit; central differences) |
| `∂_nT` | not used | `∇φ·∇T`, ordinary differences |
| `s_t` | not used | `e_d·∇T − a(∇φ·∇T)` |
| volume fraction `φ` | transients only | `clip(½ + φ_c/h)` |
| **new case-file datasets** | **none** | **none** |

The last row is the point. Because `φ` is a true distance function, `|∇φ| = 1`
and `∇φ` **is** the unit normal — the identity the RANS wall functions already
exploit for `sst%wnorm`. So the normal, the obliquity and the normal derivative
all follow from central differences of `φ` and `T`: **the entire method,
baseline and escalation, needs exactly one extra stored field, and that field
already exists.**

**Two arrangements are required, neither of them new data:**

1. **Availability without RANS.** Today the distance state is built only when a
   `[rans]` section is present (the T1 hook). A conjugate run is typically a DNS
   with no turbulence model, so the build must also be triggered by
   `ibm_wall = conjugate` — a dataset read for file geometry, an existing
   routine call for analytic geometry.
2. **A ghost-inclusive sign.** `φ` is needed on the halo layer so that cut faces
   on a block boundary see both signs. `dwall` is already ghost-inclusive; the
   solid marker must be too.

**Accuracy caveat on `w`** (LaTeX note §6.6): exact for a plane; for a curved
interface `φ_L` and `φ_R` are distances to *different* nearest points and `w`
carries a relative error `O(κh/a)`. At DNS resolution with `refine_body`,
`κh ≪ 1` — but note the `1/a`: the estimate degrades for **grazing** arms
(nearly tangent to the surface). Those carry little interface flux; the guard is
to fall back to the plain harmonic mean below a threshold on `a`. If exactness
in `w` is ever required, the per-arm bisection in `add_neighbor_coeff` supplies
it — at the price of the per-arm dataset this baseline exists to avoid.

---

## 9. Configuration surface

```ini
[scalar.1]
name        = theta
pr          = 0.71          ; FLUID Prandtl number (unchanged)
ibm_wall    = conjugate     ; dirichlet | adiabatic | conjugate
solid_k     = 100.0         ; k_s/k_f      — the interface condition
solid_rhocp = 3.0           ; (ρc)_s/(ρc)_f — the transient; α_s/α_f = solid_k/solid_rhocp
solid_init  = 0.0           ; initial solid temperature (default = `initial`)
solid_source = 0.0          ; volumetric source in the solid (e.g. Joule heating)
contact_resistance = 0.0    ; R_c; adds to R in (★), zero = perfect contact
tangential_correction = false ; the (★★) term; default off, see §3
```

`ibm_value` keeps its S3 meaning and is rejected with `conjugate` (the body
temperature is an outcome). `solid_k`/`solid_rhocp` must be positive;
`solid_k` with `ibm_wall /= conjugate` is a hard config error.

---

## 10. Increments and gates

Every increment ends with the standard bit-exactness gate: `[scalar] count = 0`
**and** a scalar run with `ibm_wall /= conjugate` must be max_abs 0 (nofma, CPU
AND GPU) versus the pre-increment binaries on the standard 7-case suite
(min_channel, les_ibm ± refine_body, Beltrami y-slab, turb180, wf180_y30,
lam30t). New cases and drivers go in `validation/conjugate/` with a README
recording commands and numbers.

**C1 — the baseline.** Solid cells become real unknowns (`κ`, `C`); `φ = ±dwall`
built without `[rans]` and ghost-inclusive; the cut-face branch (★); masked
convection; the Peclet limiter over materials.
*Gates*:
- **1D two-material slab**, cut position swept through a full cell
  (`δ_L/h ∈ {0.05 … 0.95}`), `κ_s ∈ {10⁻², 1, 10, 10³}`: **exact** — the steady
  solution is piecewise linear and the series resistance is then exact. This
  also validates `w` from the signed distances against the analytic cut
  position, which is the one genuinely new ingredient.
- **`κ_s → ∞` reproduces the `dirichlet` mode; `κ_s → 0` the `adiabatic` one**,
  to the discretisation's own tolerance.
- **Conservation**: `Σ C·T·ΔV` drift at round-off in an insulated composite box.
- 1 rank == 4 ranks EXACT; CPU == GPU EXACT.

**C2 — measure the trade, then decide.** Implement the indicator `e_face` and,
behind `tangential_correction`, the (★★) term.
*Gates*:
- **Oblique plane interface**, manufactured solution at 30°/45° with a
  *controlled* ratio `|∇_tT|/|∂_nT|`, `κ_s ∈ {10, 10³}`: measure the
  interface-flux error against that ratio and against `h`, and check it against
  the §3 prediction and against `e_face`. Expected: at DNS-like ratios (~10⁻²)
  the baseline is adequate; at O(1) ratios it stalls near first order and (★★)
  recovers second. **Record both curves** — this is the number that decides
  whether (★★) ever ships enabled.
- **Cylindrical shell**, exact log solution: isolates the `O(κh/a)` error in `w`
  in the interface flux (Nusselt).
- Grid-aligned cases from C1 unchanged to round-off (`s_t = 0` by construction).
- Time-step penalty of (★★) measured at `κ_s = 10³`.

**C3 — transients and production.** Fraction-weighted capacity; the `[f/β]`
second-derivative jump term only if C2 demanded it — **it did not**: C2's
verdict is that the tangential escalation ships disabled, so C3 carries no
second-derivative work. C3 also inherits C2's time-step debt (the shipped
`pecletmax = 0.4` sits at 96 % of the Gershgorin bound at a cut cell), which
belongs here because the capacity change moves `dt` at the same cells.
*Gates*:
- **Transient two-material slab/sphere** with a capacity jump vs the analytic
  solution, second order in the time-resolved interface flux.
- **Conducting channel wall**: reuse `validation/channel_interface/les_ibm/` —
  its off-grid plane walls sit mid-cell (`y = 0.259375`), so the cut fractions
  are non-trivial while the interface stays grid-aligned. Clean regression, a
  case where the baseline is *exact*, and it exercises the LES/IBM/2:1 stack.
- **Nusselt diagnostic** (`Σ` over cut faces of (★)) against the Gauss/CV
  border-flux cross-check, the way A2 validated `C_L`/`C_D`.

**C4 — conducting sharp corners. SEPARATE SESSION, optional.** The COCO wedge
model: two-material wedge eigensolution `γ(α, κ_s)` from `det M(γ) = 0` with
`[T] = 0`, `[k ∂_θT] = 0` on both wedge faces, tabulated at prepare time,
applied on flagged corner cells within a radius `≈ 2Δ`, anchored on
`T − T_corner` with `T_corner` reconstructed (the regular modes are not removed
here, so the parameter count of §4 grows to four).
*Gates*: **`wedge` → `flat` to round-off as `α → 180°`** (the analogue of COCO
reducing to Luchini's λ on a flat wall — cheap, and it validates the whole
eigensolution path); a conducting-wedge manufactured solution, second order in
the corner flux versus first with the flat model; conservation and the C1/C2
gates unchanged; a conducting riblet/fin as the production demonstration.

**Reference implementation (during C2, it is cheap):** a standalone 2D EJIIM
solve (augmented system + GMRES) for the manufactured cases, in
`validation/conjugate/reference/`. It separates "the scheme is wrong" from "the
implementation is wrong" — the role `mobygeom` plays as a cross-implementation
reference.

---

## 11. Files touched

| file | change |
|---|---|
| `src/modules/scalar.f90` | `κ`/`C` fields, the cut-face branch (one local model per face), config keys, solid init; C2 adds `s_t` + the indicator |
| `src/modules/rans.f90` (or a shared home) | build the distance state when `ibm_wall = conjugate` and no `[rans]` — no new computation, only a new trigger |
| `src/modules/init.f90` | Peclet limiter over materials |
| `docs/conjugate/` | the derivation note (already written) |
| `validation/conjugate/` | NEW — cases, drivers, reference implementation, README |

Note what is **absent**: no `moby_prepare` change, no `io.f90` / `field_hdf5.c`
change, no case-file format change, no new dataset, no "re-run `moby_prepare`"
error path. That is the whole benefit of the baseline.

---

## 12. Landmines (read before writing code)

- **`remove_solid` must be OFF and `keep_buried = true` for conjugate runs.**
  Buried blocks now carry the solid temperature field; removing them deletes the
  solid domain. Same class of trap as the A3 penalization-force finding
  (`validation/naca0012/README.md`), and it must be a **hard config error**.
- **`dwall` must be built without `[rans]`** and its sign marker must be
  ghost-inclusive (§8). A cut face on a block boundary needs both signs.
- **Never rescale by a nodal value** the way COCO does: the conduction analogue
  `T − T_Γ` vanishes on the interface and changes sign across it, so the
  coefficient is unbounded and sign-indefinite, and a negative λ is
  anti-diffusive (§4).
- **Conservation is the invariant to lean on**: write the interface term as a
  *face flux*, never as a cell source, and anchor the local model to the face.
  Then any error in `s_t` — including setting it to zero, which is the baseline
  — costs accuracy but never conservation, and the conservation gate stays a
  real test of the implementation.
- **`refine_body` is what keeps cut cells away from 2:1 interfaces**: it refines
  touched blocks plus a one-block 26-neighbour buffer to the finest level. Make
  this a **checked precondition** at init, not an assumption. It is also what
  keeps `κh ≪ 1` in §8.
- **Grazing arms**: guard `w` when `a = (φ_L − φ_R)/h_d` is small (§8).
- **Thin bodies** (thickness < h) are missed: both cell centres read fluid.
  Inherent to the cell-centred marker; document, do not fix.
- **`ν_t` must not enter the solid or the cut faces.**
- **Nusselt diagnostics need `keep_buried`** for the same reason the
  penalization forces do.

---
## 13. Next-session prompt (C3)

Written 2026-08-28 against branch `scalar`, C2 landed on top of `735d360`.
**Reference binaries: `~/s5c_ref_binaries`, commit `8f60944`** — read its
`PROVENANCE.txt`. It remains valid against C2: the correction and the
indicator are both dormant without a conjugate scalar, and C2 changed no
arithmetic on any other path (32/32 bit-exactness case-runs, `max_abs 0`).

> Implement increment **C3** of `docs/next_session_conjugate.md` — *the
> fraction-weighted capacity, the conjugate Nusselt diagnostic, and the
> time-step debt C2 left* — on branch `scalar`.
>
> **Read first, in this order:** the STATUS header of
> `docs/next_session_conjugate.md` (C1 and C2 are DONE; C2's six deviations
> and its verdict are the ground truth), then §6's capacity paragraph and
> §10's C3 entry, then `validation/conjugate/README.md` — the whole C2
> section, because three of its findings change what C3 has to do — then the
> "Active work" section of CLAUDE.md.
>
> **C3 IS THREE THINGS, AND THEY SHARE ONE COST.** The capacity change and
> the time-step fix both move `dt` at cut cells, so they must land together:
> the re-gate is paid once.
>
> 1. **The fluid-fraction-weighted capacity** (§6). A cut cell contains both
>    materials, so `C_cell = φ + (1 − φ)C_s` with `φ` the fluid VOLUME
>    fraction. It is irrelevant at steady state — C1 gated exactly that
>    (`8.85e-14 / 4.37e-14 / 4.88e-15` over `C_s ∈ {0.5, 1, 8}`) — and
>    required for second-order TRANSIENTS, which is the whole of this item.
>    **C2 leaves you the machinery**: `check_oblique.py face_area_fraction`
>    is the 2D plane-in-rectangle form, validated against brute force to the
>    quadrature error (9e-4 at 400² samples); the cell version is its 3D
>    sibling (plane-in-cube, same family), and both need only `φ` at the
>    centre and `∇φ` — data the solver already has. Do NOT reach for
>    `clip(½ + φ_c/h)` without measuring it against the closed form first.
> 2. **The Nusselt diagnostic** (§10). `Σ` over cut faces of (★), against the
>    Gauss/CV border-flux cross-check, the way A2 validated `C_L`/`C_D`.
>    **Reuse `scalar_stats.f90`'s cancellation-free form, do not reinvent
>    it** — the S3 finding (`docs/next_session_scalar.md`) is that
>    `∫coef_p(s_body − s)dV` is `1e28 × 0` in a penalized cell and loses
>    ~37 % of the heat, and the 2026-08-05 fix on top of it is that the
>    split must EXCLUDE solid cells explicitly rather than rely on a
>    floating-point cancellation. C1 left the conjugate flux columns
>    SMOKE-gated only (deviation 6); this is where they get their numbers.
> 3. **The time-step debt.** C2 measured that the shipped `pecletmax = 0.4`
>    sits at **96 % of the Gershgorin bound** at a cut cell
>    (`dt = pecletmax·3C/diag` against an allowed `1.25 C/diag`), and that an
>    OBLIQUE interface at `κ_s = 10³` therefore goes to NaN at 0.4 while 0.2
>    is stable. C1's gates never saw it because their interface is
>    grid-aligned. Decide it now, with the capacity change already in hand:
>    either tighten the cut-cell `share` (3 → 2 gives a 1.56× margin; → 1.5
>    reproduces the measured-stable 0.2) or leave it and make the guard a
>    documented config error. **Whatever you choose, C1's recorded numbers
>    move**, which is why it belongs here and not in a later increment.
>
> *Gates* (§10's C3 list, plus the two C3 inherits):
> 1. **Transient two-material slab** with a capacity jump against the
>    analytic solution, second order in the TIME-RESOLVED interface flux.
>    `validation/conjugate/slab.ini` already sweeps the cut position; the
>    transient reference is the standard two-material series solution.
> 2. **Conducting channel wall**: reuse
>    `validation/channel_interface/les_ibm/` — its off-grid plane walls sit
>    mid-cell (`y = 0.259375`), so the cut fractions are non-trivial while
>    the interface stays grid-aligned, i.e. a case where the baseline is
>    EXACT. It also exercises the LES/IBM/2:1 stack.
> 3. **Nusselt** against the Gauss/CV border flux. LANDMINE inherited from
>    A2: `niter = 6` IBM runs accumulate a velocity-neutral oscillating mode
>    in stored `pn`, so a border-flux cross-check needs a clean-`p` snapshot
>    (zero `pn` in a copy of the converged restart, rerun ~300 steps at
>    `niter = 60`; restarting the POLLUTED `p` at `niter = 60` transients
>    violently).
> 4. **Every C1 and C2 gate re-run**, since `dt` moves:
>    `./run_gates_c1.sh` and `./run_gates_c2.sh` must both come back — the
>    steady-state equalities to their recorded numbers, and the C2
>    measurements (which are evaluated on ANALYTIC fields and so must be
>    bit-stable) unchanged. If you tighten the rate, the `converge` group's
>    200 000 steps buy less physical time — check it still converges rather
>    than assuming it.
> 5. **The full bit-exactness protocol, unchanged**, against
>    `~/s5c_ref_binaries` (4 drivers) plus the five S-suite groups. C2's run:
>    32/32 case-runs `max_abs 0`, 91 S-suite checks, `rc = 0` throughout.
>    LANDMINE: `run_gates.sh`'s `det` group compares CPU vs GPU at TOLERANCE
>    0, so hand it the nofma pair.
>
> **Scope boundaries — do NOT:** touch the tangential correction. C2 measured
> it and the verdict is SHIP DISABLED; the analysis of what a working version
> would need is in the README's "the way out" and in the LaTeX note, and it
> starts with `s_t`, not with the correction. Do not attempt C4's wedge model
> — it rests on the same premise C2 falsified and is gated by the same
> measurement.
>
> **What C1 and C2 learned that will bite you:**
> - **the explicit limit at a conjugate interface is not a per-material
>   `α = κ/C`**, and a cut cell attains a Gershgorin factor the uniform
>   interior never excites. `scalar_conjugate_peclet_rate` encodes both, and
>   C2 showed the remaining margin is 4 %. The capacity change alters `C_i`
>   at exactly those cells, so re-derive the rate rather than assume it
>   survives;
> - **never rely on a floating-point cancellation as a classification** — the
>   2026-08-05 body-heat double count is the precedent, and the Nusselt
>   diagnostic is the same shape of computation;
> - **a claim of the form "no case exercises this path" must be checked
>   against the GENERATED inis**, not only the committed ones (C1's own
>   `t_final` regression);
> - **`check_oblique.py` and `check_cylinder.py` evaluate schemes on ANALYTIC
>   fields** with `φ` from the case file. That makes them cheap and
>   solver-independent — use the same style for the transient reference
>   rather than building another solver-in-the-loop gate;
> - REBUILD ALL FOUR BINARIES after any change that touches `dt`. The `det`
>   group compares three binaries against each other and passes happily when
>   all three are equally stale.
>
> **Conventions:** build both paths (`./compile.sh cpu && ./compile.sh gpu`,
> module `toolkits/nvhpc/25.9`) plus the nofma pair via
> `validation/scalar/compile_nofma.sh`; always launch through `mpirun`, even
> on one rank; new cases and drivers go in `validation/conjugate/` beside C1's
> and C2's; never declare the increment done with a failing build or an
> ungated result.
>
> Stop after C3's gates and report. Update this document's STATUS header with
> what landed, each gate's measured number, and any deviation from the plan —
> C1's and C2's entries are the model. **And state explicitly what the
> time-step decision was and what it cost**, because every conjugate number
> recorded before it is measured against a different `dt`.
