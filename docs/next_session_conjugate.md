# Conjugate heat transfer at an immersed interface (branch `scalar`) — strategy

STATUS: **C1 IS DONE AND GATED (2026-08-27/28, branch `scalar`, on top of
`e1b5c2c`).** C2/C3/C4 are not started; §10 still describes them. Gates, the
commands that produced every number, and the two time-step findings are in
[`validation/conjugate/README.md`](../validation/conjugate/README.md).
This document is the plan; the derivation, the limit checks, the accuracy
argument and the sketches live in
[`docs/conjugate/conjugate_ibm.tex`](conjugate/conjugate_ibm.tex) (build with
`make` in that directory).

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
time-step cost at `κ_s = 10³` if it is ever switched on.

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
second-derivative jump term only if C2 demanded it.
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

## 13. Next-session prompt (C1)

Refreshed 2026-08-27 against branch `scalar` HEAD `21a7701`. **The reference
binary set changed**: earlier revisions of this prompt named
`~/s5a_ref_binaries`, which its successor's `PROVENANCE.txt` documents as STALE
for cold-started RANS. Use `~/s5c_ref_binaries` (commit `8f60944`).

> Implement increment **C1** of `docs/next_session_conjugate.md` — conjugate
> heat transfer at an immersed interface, baseline scheme — on branch `scalar`
> (HEAD `21a7701`).
>
> **Read first, in this order:** `docs/next_session_conjugate.md` (§0–§3 are the
> whole method; §5, §8, §10 and §12 are the constraints), then
> `docs/conjugate/conjugate_ibm.tex` §6 for the derivation and the sketches
> (`make` in that directory), then the STATUS header of
> `docs/next_session_scalar.md` (S0–S5a is CONCLUDED and conjugate is the live
> scalar work), then the "Active work" section of CLAUDE.md, then
> `validation/scalar/README.md` for the gate machinery you inherit.
>
> **The scheme is one face coefficient.** At a face whose two cell centres
> straddle the interface (`φ_L·φ_R < 0`, with `φ = ±dwall` and the sign from the
> solid marker), replace the face diffusivity by the distance-weighted harmonic
> mean built on `w = φ_L/(φ_L − φ_R)`; every other face keeps today's kernel
> line verbatim. `dwall` already exists at every cell centre including ghosts,
> at the `VAR_P` position, for both the analytic and the file geometry path.
> **No new dataset, no `moby_prepare` change, no case-file format change** — if
> you find yourself adding one you have missed the obliquity lemma (§8, or §6.2
> of the LaTeX note).
>
> **Also in C1:** solid cells become real unknowns (`κ = k/k_f`,
> `C = ρc/(ρc)_f`, both ≡ 1 in the fluid); convection HARD-masked on cut faces
> and on faces whose staggered velocity node is solid; the Peclet limiter takes
> `max` over materials; the distance state is built when `ibm_wall = conjugate`
> even with no `[rans]` section, and its solid marker is ghost-inclusive.
>
> **Scope boundaries — do NOT:**
> - implement the tangential correction (★★) or the `e_face` indicator. That is
>   C2, and whether (★★) ever ships enabled is a question C2 answers with a
>   measurement, not a design decision to pre-empt.
> - re-express the existing `dirichlet` / `adiabatic` modes through the new
>   cut-face arithmetic. §5 records why: they shipped as an inline penalization
>   statement and six face masks, both gated, and re-expressing them cannot be
>   bit-exact. **`conjugate` is a THIRD branch**; the `κ_s → ∞` / `κ_s → 0`
>   limits are gates it must reproduce to discretisation tolerance, not a
>   refactor.
>
> **Two S3 deviations that bite here** (both in this document's header):
> - `ibm%mu` was never given the `VAR_P` extent — only `ibm%coef` was — because
>   the scalar penalization factor is Prandtl-dependent and formed inline in the
>   kernel. A stored per-face conjugate coefficient needs its own array.
> - The body heat release **cannot** be measured as `∫coef_p (s_body − s) dV`: a
>   Dirichlet solid cell holds the body value to the last bit, so the product is
>   `1e28 × 0` and ~37 % of the heat is invisible. The cancellation-free form
>   lives in `scalar_stats.f90` and `validation/scalar/check_scalar_ibm.py
>   surface`; a conjugate Nusselt diagnostic must reuse it. S4's per-row
>   accumulators are built from the TRANSPORT KERNEL's own face diffusivity, so
>   a conjugate cut face must extend that same expression or the statistics stop
>   reporting the flux the kernel actually applied.
>
> **Reference binaries: `~/s5c_ref_binaries`, commit `8f60944`.** Read its
> `PROVENANCE.txt` — it supersedes `~/s5a_ref_binaries` and `~/s5b_ref_binaries`
> and says exactly how each is stale. It was cross-checked max_abs 0 against the
> profile-timer commit, and the `t_final` fix at HEAD is provably inert on the
> suites (every suite case is nsteps-terminated), so it is valid against
> `21a7701`. If you ever archive a new set, cut it from a COMMIT, never from a
> working tree, and record the hash.
>
> **Bit-exactness protocol** — `[scalar] count = 0` AND any
> `ibm_wall /= conjugate` run must stay max_abs 0, CPU and GPU:
> ```
> cd validation/scalar
> REF=~/s5c_ref_binaries/moby_solve_cpu_nofma MODE=cpu ./run_bitexact.sh     # 7 cases
> REF=~/s5c_ref_binaries/moby_solve_gpu_nofma MODE=gpu ./run_bitexact.sh
> REF=~/s5c_ref_binaries/moby_solve_cpu_nofma MODE=cpu ./run_bitexact_s3.sh  # 9 cases
> REF=~/s5c_ref_binaries/moby_solve_gpu_nofma MODE=gpu ./run_bitexact_s3.sh
> ./run_gates.sh ; ./run_gates_s2.sh ; ./run_gates_s3.sh
> ./run_gates_s4.sh ; ./run_gates_s5.sh
> ```
> LANDMINE: `run_gates.sh`'s `det` group compares CPU vs GPU at TOLERANCE 0, so
> hand it the nofma pair:
> `BIN=$PWD/../../build_cpu_nofma/moby_solve GBIN=$PWD/../../build_gpu_nofma/moby_solve ./run_gates.sh`.
>
> **Landmines that will actually bite** (§12 has the rest):
> - `remove_solid = true` with `ibm_wall = conjugate` must be a HARD CONFIG
>   ERROR — buried blocks carry the solid temperature field.
> - `refine_body`'s one-block 26-neighbour buffer is what keeps cut cells away
>   from 2:1 interfaces and keeps `κh ≪ 1`: CHECK it at init, do not assume it.
> - Write the interface term as a **face flux, never a cell source**, and anchor
>   the local model to the face. That is what keeps conservation exact
>   regardless of the dropped tangential term, and what lets C2 and C4 drop in
>   without touching the surrounding machinery.
> - Guard `w` for grazing arms (small `a = (φ_L − φ_R)/h_d`): fall back to the
>   plain harmonic mean below a threshold.
>
> **Conventions:** build both paths (`./compile.sh cpu && ./compile.sh gpu`,
> module `toolkits/nvhpc/25.9`) plus the nofma pair via
> `validation/scalar/compile_nofma.sh`; always launch through `mpirun`, even on
> one rank; new cases and drivers go in `validation/conjugate/` (it does not
> exist yet) with a README recording what was run and measured; never declare
> the increment done with a failing build or an ungated result.
>
> **C1 gates — all must pass and be recorded:**
> 1. 1D two-material slab, cut position swept through a full cell
>    (`δ_L/h ∈ {0.05 … 0.95}`), `κ_s ∈ {10⁻², 1, 10, 10³}`: **EXACT**, and `w`
>    matches the analytic cut position — that weight is the one genuinely new
>    ingredient in this increment.
> 2. `κ_s → ∞` reproduces the `dirichlet` mode and `κ_s → 0` the `adiabatic`
>    one, to the discretisation's own tolerance.
> 3. `Σ C·T·ΔV` conserved to round-off in an insulated composite box.
> 4. 1 rank == 4 ranks EXACT; CPU == GPU EXACT.
> 5. The full bit-exactness protocol above.
>
> Stop after C1's gates and report. Update this document's STATUS header with
> what landed, each gate's measured number, and any deviation from the plan —
> the S3 deviations recorded in the header are the model for how that is done.
