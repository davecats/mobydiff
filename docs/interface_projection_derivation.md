# 2:1 interface pressure projection — correct SPD formulation and SOR recast

Status: **design note** (2026-06-22). Captures the correct formulation after three
failed bolt-on fixes (task #38). The diagnosis is settled; this note derives *why*
the bolt-ons fail and *what* the consistent scheme must be, so the next session can
implement it directly.

## 0. The observed failure (recap)

On the 2:1-refined Beltrami flow, fed the **exact divergence-free** velocity, the
projection still changes it (corrector change max 2e-2, growing to a long-time
blow-up). The discrete divergence operator sees a spurious **max |div| = 1.44** on
that div-free field, located entirely on the **fine cells at the high-side
interface faces** (single block: exactly 0.0). Root cause: the slaved fine
sub-face normal velocity is filled by **piecewise-constant injection** of the
covering coarse face, which drops the coarse→fine tangential variation of the
normal velocity; the fine cell's `1/dx_fine` amplifies that O(h) mismatch to O(1).
Linear (tangential) prolong of that slaved face cuts the divergence 26× on a frozen
field (1.44 → 0.056) — the right *direction* — but every attempt to use it in the
relaxation either blew up the contraction or corrupted the velocity. This note
explains the missing structure.

## 1. The discrete projection is an SPD system; the operators must be adjoint

Continuous projection: `u = u* − dt ∇φ`, with `∇·u = 0`, so `∇²φ = ∇·u*/dt`.

Discrete MAC analogue, with face-velocity DOFs `u`, cell-pressure DOFs `p`, a
**divergence** operator `D` (faces → cells) and a **gradient** operator `G`
(cells → faces):

```
u = u* − dt·G p,     D u = 0
  ⇒  (D G) p = D u*/dt           (the pressure Poisson system)
```

For this to be a *true* (orthogonal) projection the discrete Helmholtz
decomposition must hold, which requires **G = −Dᵀ** in the grid inner products
(face areas / cell volumes). Then `L = D G = −D Dᵀ` is **symmetric negative
(semi-)definite** — SPD up to the overall sign — and two things follow:

1. The corrected `u = u* − dt G p` is **exactly discretely divergence-free**
   (`D u = 0`) at convergence.
2. The correction `dt·G p` is L²-orthogonal to the divergence-free subspace — no
   spurious energy is injected.

If `G ≠ −Dᵀ` (the operators are **not adjoint**), `L` is non-symmetric, the
"projection" is **oblique**, and it can *create* divergence/energy rather than
remove it. That is precisely the instability: an oblique interface projection
pumps the corner pressure mode. (Almgren–Bell–Colella–Howell–Welcome 1998;
Martin–Colella 2000; Guittet et al. 2015 — the composite/synchronization
projection is built around keeping `D` and `G` exact adjoints across levels.)

**This is the lens for everything below: a fix is correct iff it keeps `G = −Dᵀ`.**

## 2. The 2:1 interface as a slaved DOF (the prolong operator P)

At a 2:1 face the fine block carries 4 sub-faces (2×2 tangential) sharing one
coarse face. In mobydiff's convention the **coarse (high-side) block owns** the
shared face as its interior low face; the fine block holds the 4 sub-faces as a
**slaved** high-side halo. So the fine sub-face normal velocities are not
independent DOFs — they are a linear map of the coarse face DOF(s):

```
u_F^face = P · u_C^face            (P = the tangential prolong, 1 coarse → 4 fine)
```

- **Injection (current):** `P` copies the one covering coarse value to all 4
  sub-faces. Columns of `P` are indicator-of-coarse-face; each fine sub-face reads
  exactly one coarse face.
- **Linear prolong (wanted):** `P` is bilinear tangential interpolation from the
  covering coarse face and its tangential neighbours; each fine sub-face is a
  weighted blend. For conservation the weights must be **mean-preserving**:
  `mean_4(P u_C) = u_C` over the covering face, i.e. the 4 fine fluxes sum to the
  coarse flux. (Standard bilinear face prolong satisfies this; verify the actual
  E3 weights do — §6.)

Now the operators:

- **Divergence `D`.** A fine interface cell's divergence uses its slaved sub-face
  `u_F^face = P u_C`. So `D` restricted to the interface **contains `P`**:
  `(D u)_F = … + s · (P u_C) · d1x_F`  (s = ±1 for the face side).
- **Gradient `G = −Dᵀ`.** Taking the transpose, the coarse face DOF receives
  `Pᵀ` of the adjacent fine-cell pressures: the coarse face correction is
  `(G p)_C^face = … − Pᵀ (p_F of the 4 fine cells) · (volume weights)`. The
  same `P` that spreads the coarse velocity to the fine cells must, transposed,
  **gather the fine pressures back to the coarse face**.

Injection's `Pᵀ` is the *sum* over the 4 fine cells (restriction by addition);
mean-preserving linear `Pᵀ` is the conservative restriction plus the tangential
weights. **The transpose is not optional** — it is what makes the projection
orthogonal.

## 3. Why the three bolt-ons failed (all break `G = −Dᵀ`)

- **(1) Linear-prolong the slaved velocity every sweep.** This changes `D` (the
  divergence now reads the linear-prolong face) but leaves `G` as the injection-
  based reconstruction, *and* re-derives the slaved face from the still-relaxing
  coarse interior each sweep — a non-stationary, non-adjoint coupling. Oblique +
  non-stationary ⇒ the SOR amplifies (corrector change 2e-2 → 4.48).

- **(2) Pre-freeze + per-sweep linear prolong.** Same adjoint break, same result.

- **(3) Frozen-δ reflux** (`divCorr = div(linProlong) − div(injection)` added to
  the RHS). This corrects **only the divergence `D`** (and only its constant part,
  frozen). `G` is untouched, so `G ≠ −Dᵀ`. The SOR relaxes a system whose RHS and
  operator disagree: it drives the *injection* divergence to `−δ` (a large value),
  not the true divergence to 0, and the velocity it reconstructs from the
  unchanged `G` is inconsistent — corrector change → 1.8, `divCorr` compounding
  1.5 → 4 → 45 across substeps. **A RHS-only correction can never restore
  adjointness.**

The unifying lesson: **the linear prolong `P` must enter `D` and `G` together.**

## 4. The mobydiff recast (ifGrad + the missing tangential coupling)

mobydiff's composite sweep (`pressure_solver.f90 redblack_sweep`) already encodes
the **normal-direction** interface coupling correctly:

- Reconstruction (lines 312–323):
  `q_face = qs_face − dt·μ·ifGrad·(Δp_above − Δp_below)`,
  with `ifGrad = 1/gap`, `gap = (dx_fine+dx_coarse)/2 = 1.5 dx_fine`, and
  `Δp = p − pStart`. This is the normal part of `G`: the slaved face responds to
  the coarse-minus-fine pressure across the gap.
- Divergence (line 302) reads `q_face` back; the diagonal `denom` (lines 289–300)
  carries the matching `μ·ifGrad·d1(P)` sensitivity — i.e. for the *normal*
  direction the code is already an adjoint pair (`D` reads `q_face`, `G` writes it,
  same `ifGrad`), which is why single-direction band-refined channels are stable.

What is missing is the **tangential** structure of `P`:

1. **In `D`:** the slaved sub-face must be the *linear-prolong* of the coarse
   face — i.e. the reconstruction's `qs_face` and its pressure response must carry
   the tangential blend `P`, not the covering-coarse constant. Concretely the
   slaved face for fine cell `(i=nb, j, k)` should read coarse pressures at the
   covering coarse cell **and its tangential neighbours**, with the `P` weights,
   so `Δp_above` is the *interpolated* coarse pressure at the fine sub-face
   location, not the single covering coarse cell.
2. **In `G = −Dᵀ`:** the coarse owner face must gather the 4 fine-cell pressure
   residuals with `Pᵀ` (the same weights, transposed) when it relaxes — so a fine
   cell's pressure update feeds back to the coarse face (and hence to the coarse
   cell's divergence) through the conservative restriction.

So the interface row of the SPD system, for a fine cell adjacent to a high-side
face, becomes (schematically):

```
denom_F · Δp_F  =  div_F
   div_F = (other faces) + [ qs_face − dt·μ·ifGrad·( Σ_t w_t p_C,t − p_F ) ] · d1x_F
   denom_F = (other faces) + dt·μ·ifGrad·d1x_F · (1 + Σ_t w_t²)     ← P enters the diagonal
```

where `Σ_t w_t p_C,t` is the **tangentially interpolated coarse pressure** at the
fine sub-face (the `P`-weighted stencil, `Σ_t w_t = 1` by mean-preservation), and
the coarse cell's own row symmetrically receives `w_t·(…)` from each covering fine
cell (the `Pᵀ`). The diagonal picks up `Σ w_t²` from the interpolation, keeping the
matrix SPD and diagonally dominant (so SOR still contracts — *this* is how we get
the linear-prolong accuracy without the attempt-(1) blow-up: the coupling is
**stationary and symmetric**, baked into the operator, not re-interpolated from a
moving field).

## 5. The SOR sweep, consistent form

Red-black Gauss-Seidel / SOR on `L Δp = b`, `b = D u*/dt`:

```
for each colour, each cell C:
    r_C   = b_C − (L Δp)_C            # residual = current divergence/dt
    Δp_C += ω · r_C / denom_C         # SOR update
    # write back the faces C owns/holds, from the UPDATED pressures:
    #   interior faces : q ± ω r_C · d1 · μ                 (unchanged)
    #   normal interface face : q = qs − dt μ ifGrad (Δp_above − Δp_below)
    #   NEW tangential part   : Δp_above is the P-interpolated coarse pressure;
    #                           the coarse owner face simultaneously gathers
    #                           this cell's Δp_F with weight w_t (the Pᵀ term)
```

Key points for the implementation:

- The **residual `r_C` must be computed with the linear-prolong slaved face**
  (the `P`-interpolated `q_face`), because `D` contains `P`. This is what the
  frozen-δ tried to approximate — but here it is exact and *paired with the
  matching `G`*, so it converges to `D u = 0` instead of to a wrong fixed point.
- The **diagonal `denom_C` must include the `Σ w_t²` term** so the update is
  consistent with the off-diagonal `Pᵀ` coupling. Omitting it (frozen-δ) is the
  adjoint break.
- The coupling is **stationary**: `P` (the geometric weights `w_t`) is fixed by
  the grid, computed once; only the *pressures* it multiplies change across sweeps.
  No re-exchange of a moving velocity inside the loop (that was attempt-1's sin).

## 6. Implementation plan

1. **Get the `P` weights** for the normal-velocity tangential prolong from the E3
   gather (`comm.f90 entry_gather_map`, `lLin` dims). Confirm they are
   mean-preserving (`Σ_t w_t = 1`); if not, symmetrize/renormalize for
   conservation (the SPD derivation needs it).
2. **Reconstruction (`redblack_sweep` 312–323):** replace the single covering
   coarse `Δp` in `Δp_above` with the `P`-weighted tangential blend of coarse
   cell pressures at the fine sub-face location.
3. **Diagonal (`denom` 289–300):** add the `Σ w_t²` factor on the interface faces.
4. **Coarse owner row (`Pᵀ`):** when the coarse cell relaxes, gather each covering
   fine cell's pressure with `w_t` into the shared-face term — the transpose of
   step 2. This is the genuinely new coupling and the part most worth a careful,
   small first step (start with the flat face, defer edges/corners).
5. **Drop** the per-sweep velocity linear-prolong and the frozen-δ; the
   correctness now lives in the operator, not in a halo refill.
6. **Corners/edges:** the divergence map showed the corner is worst. The corner
   fine cell has two slaved interface faces; the `P` stencils there are 2D blends
   of coarse pressures — implement after the flat face is verified.

## 6b. Machinery audit (2026-06-22) — what exists vs. what the fix needs

A full read of the transfer machinery (`comm.f90 gather_taps`/`entry_gather_map`,
`pressure_solver.f90`) pins down exactly why this is a new-operator job, not a
flag flip:

- **Existing linear prolong is the plain (3/4, 1/4) interpolation** (`gather_taps`
  ~line 1556): a fine cell-centred halo = 3/4·(covering coarse) + 1/4·(adjacent
  coarse). Face-staggered own-dim: coincident → 1 tap, midpoint → (1/2,1/2).
- **It is velocity-only.** `lin = merge(1, 0, c%linProlong .and. var /= VAR_P)`
  (lines 1208, 1335) — **VAR_P is excluded**, pressure always injects. So the
  reconstruction's `Δp_above` (the halo pressure) is piecewise-constant today.
- **It is NOT mean-preserving.** The four fine sub-faces over one coarse face do
  not average back to the coarse value (each is its own 3/4-1/4 blend). For a
  *projection* that is fatal: non-mean-preserving prolong ⇒ the coarse flux ≠
  Σ fine fluxes ⇒ global mass is not conserved.
- **Restrict is the uniform average** (`copy_local_scalar_entries` ~line 1268:
  `val/(c1*c2*c3)`; `entry_gather_map` gc=2). In the area-weighted inner product
  the uniform average **is exactly the adjoint of injection** (`average = Pᵀ_inj`).
  So **injection-prolong ↔ average-restrict is already a consistent adjoint pair**
  — which is precisely why the current scheme is *stable* (just O(1)-inaccurate at
  the interface). It is a correct projection of the *wrong* (injected) velocity.

Therefore the adjoint-consistent fix is **not** "turn on linear prolong"; that
would pair a linear `P` with the average `Pᵀ_inj` (non-adjoint) and a
non-conservative `P` (mass loss) — exactly the failure modes seen in attempts 1–3.

**What must be built (all four, together):**
1. A **mean-preserving** tangential prolong `P_mp` for the projection. The clean
   construction: `f_m = c0 + δ_m·(c_+ − c_−)/2` with the covering coarse `c0` and
   its two tangential neighbours `c_±`, offsets `δ_m = ±1/4` (coarse units). This
   is O(h²) AND `mean_m(f_m) = c0` (conservative). (NOT the stock 3/4-1/4.)
2. Its **exact transpose** restrict `P_mpᵀ` (not the uniform average): `c0`
   gathers the covering fine cells, and the neighbour coarse cells `c_±` gather
   the same fine cells with `±δ_m/2`. This is a genuinely new restrict stencil.
3. **Pressure** must use `P_mp`/`P_mpᵀ` at the interface (remove the `var/=VAR_P`
   exclusion *for the projection path only* — the momentum path keeps stock
   prolong), so `Δp_above` in the reconstruction is the tangential blend.
4. The **predictor slaved velocity `qs`** must be filled by `P_mp` (a
   `composite_qs_setup` fed from a mean-preserving prolong), and the diagonal
   `denom` must add the `Σ δ_m²` self-coupling so the SPD operator stays
   symmetric and diagonally dominant (SOR contraction preserved).

This is a self-contained but non-trivial unit: new `P_mp`/`P_mpᵀ` gather stencils
(device kernels + entry weights in `comm.f90`), a projection-only pressure-prolong
path, and the `denom` term. Build flat-face first (1 tangential neighbour pair),
verify each regression metric, then generalize to the 2D corner blend. Until all
four land together the system is non-adjoint and reproduces the blow-up — partial
application is worse than injection (empirically: attempts 1–3).

## 6c. Exact operator formulas (derived 2026-06-22; Piece 1 validated)

Conventions: high-side COARSE owns the face (its interior low face); the FINE
block (physHigh==FACE_COARSE) slaves its high halo. 1D tangential shown; the flat
face is the 2D product, the corner the full 2D blend.

**Piece 1 — mean-preserving prolong `P_mp` (VALIDATED, D side).** Slaved fine
sub-face m over coarse face c0, tangential neighbours c_±:
```
u_F,m = c0 + delta_m*(c_+ - c_-)/2,   delta_m = -1/4 (lower fine), +1/4 (upper)
```
mean over the covered fine = c0 (conservative). Edge reads clamp to the injected
halo (holds the right coarse per cell). Divergence of the exact Beltrami:
injection 1.44 -> 0.021 (68x), beats the non-conservative 3/4-1/4 linear (0.056).
Prototype: `apply_mp_interface_velocity` (main.f90, host, MOBY_RKLINPROLONG).

**Piece 2 — transpose restrict `P_mp^T` (G side, the crux).** The coarse low-halo
pressure the coarse owner reconstructs against:
```
dp_C(0,J) = mean(covered fine) + (1/16)*(p_up - p_low)_{C-1} - (1/16)*(p_up - p_low)_{C+1}
```
i.e. the average restrict PLUS the fine pressure *gradient* inside the two
tangential-neighbour coarse cells (the (3/4-1/4)-prolong analogue would have its
own transpose; mean-preserving gives the symmetric 1/16 form above). This
gradient is destroyed by averaging, so `P_mp^T` cannot be a post-exchange patch on
the coarse halo -- it must be a NEW restrict in the comm.f90 gather: each coarse
halo cell gathers its covering fine (weight 1, area-scaled to the mean) plus the
fine of its tangential-neighbour coarse cells with weights +-1/16 (the delta/2
transpose, area-scaled). 2D flat face: the product of two 1D transposes; corner:
the 2D blend. This is the substantial remaining unit.

**Piece 3 — none.** The fine couples to the coarse pressure THROUGH the prolonged
velocity (Piece 1), not via a coarse->fine pressure prolong. Confirmed by the
D = P_mp, G = -D^T structure: only the velocity prolong and the pressure restrict
carry the interface coupling.

**Piece 4 — denom (tractable).** Fine interface-cell diagonal gains the P_mp P_mp^T
self term:
```
denom += (sum_dims delta^2 / 2) * dt_gamma * mu * ifGrad * d1(VAR_P)
```
delta=1/4 -> +1/32 per tangential dim (~3% per dim). Keeps the operator SPD and
diagonally dominant so the SOR contraction is preserved (the property attempts
1-3 lacked).

**Integration order:** (i) Piece 1 -> device kernel, applied after each projection
velocity exchange + in composite_qs_setup; (ii) Piece 2 -> new restrict in the
gather (the hard part); (iii) Piece 4 -> one sweep line; (iv) full regression. All
four must land together: with D=P_mp but G,denom unmatched the operator is
non-adjoint and the relaxation diverges (empirically attempts 1-3). Status: Piece 1
validated in isolation (divergence dump); Pieces 2+4 derived; Piece 2's gather
restrict is the remaining substantial implementation.

## 6d. Implementation status (2026-06-22, this session)

All operator machinery is implemented behind `MOBY_MP=1` (off by default; the
production path is the unchanged injection/average composite projection, verified
**bit-exact** to the prior baseline on refined Beltrami and single-level
poiseuille, CPU + GPU). `tools/beltrami_regression.py` drives the three tests.

What landed (in `comm.f90 gather_taps`, a TENSOR product of 1D operators so the
prolong P and restrict P^T are exact transposes; `gather_point` widened to 6
taps; per-entry `lMpt/sMpt` mark interface-FACE-restrict tangential dims; runtime
flags `mpInterface`, `mpFineOwns`):

- **Piece 1 (velocity MP prolong, gather mode 1)** — VALIDATED. `MOBY_RKLINPROLONG`
  test A (divergence of the exact field) drops **1.44 -> 0.0214** (the bilinear
  tensor form, ~ the additive prototype's 0.021). This is the clean D-side win.
- **Piece 2 (owner-pressure transpose restrict P_p^T, mode 2)** + **Piece 3
  (coarse-pressure MP prolong into the fine slave halo, mode 3)** + the in-situ
  reconstruction reading them, so D = P_p, G = -D^T. The fine block OWNS its near
  face (`mpFineOwns`: the velocity prolong fills only the deep nb+2 halo, the
  in-situ reconstruction owns nb+1) so the relaxed slaved velocity is never
  re-injected between sweeps. Two-pass per-colour exchange (Pass A injection makes
  the fine-fine pressure halos current; Pass B the MP operators) for determinism.
- **Piece 4 (denom self term, `MP_SELF`, runtime `MOBY_MPSELF`)** — present.

**Why it is NOT yet production-ready (the remaining work, both confirmed
empirically):**

1. **Stability.** The MP prolong's wider tangential stencil (covering +/- 1/8
   neighbours) makes the cross-level coupling, frozen per sweep in mobydiff's
   block-iteration, lose diagonal dominance: the corrector change GROWS with niter
   at sor=1.5 (1.8e-2 -> 8e-2 -> 0.47 over niter 10/50). Lowering sor only slows
   it; Piece 4 on the fine-cell diagonal (even MP_SELF up to 4) barely moves it,
   so the destabilising mode is NOT the fine-cell diagonal -- it is the corner
   (next point) and/or a coupled coarse-fine mode the per-level damping misses.
   Diagnostic ladder used: enabling only the prolong (average restrict) still
   diverges (non-adjoint, expected); the in-situ-velocity vs gather-velocity
   double-handling gave a WRONG fixed point (fixed by removing the direct
   velocity prolong from the relaxation, mode 1 now only feeds the diagnostic and
   the deep halo).
2. **Corners.** The flat-face P_p/P_p^T are tensor products of 1D operators that
   CLAMP at the patch tangential edges; the residual corrector change is now
   localised there (max ~5e-2 at the corners vs the flat-face L2 ~1e-3, comparable
   to baseline). The corner needs the genuine 2D blend (doc 6 step 6), still TODO.

So flat faces are adjoint and ~accurate, but the corner blend and the relaxation
stability for the wider stencil are the substantial remaining unit -- exactly the
scope the original plan flagged ("build flat-face first... corners after"). The
next session should consider §6e (the recommended next direction) before more MP
tuning, since it removes the instability at its source rather than damping it.

## 6e. RECOMMENDED next direction — interface pressure agglomeration

The whole-session learning is an **asymmetry**: of the two 2:1 orientations, only
ONE manufactures the O(1) divergence, and understanding why points at a cleaner
fix than the MP-prolong route.

**The asymmetry (why only one interface is the problem).** Ownership is by
STORAGE, not by level: the high-coordinate block owns the shared face (holds it
as its interior low face `q(1)`); the low-coordinate block holds the slave halo
`q(nb+1)`. So:
- **Orientation A — fine below, coarse above.** Transfer is coarse->fine =
  PROLONG; the FINE cell's divergence reads a prolonged face. Prolongation is an
  EXPANSION: it must INVENT the sub-grid tangential variation. Injection invents
  "flat" (all 4 sub-faces = the covering `c0`), wrong by O(h_coarse); the fine
  cell amplifies that by 1/h_fine -> **O(1) spurious divergence** (the 1.44).
- **Orientation B — coarse below, fine above.** Transfer is fine->coarse =
  RESTRICT; the COARSE cell reads the average and the FINE cells own their faces.
  Restriction is a COMPRESSION: averaging is the exact, conservative coarse
  representation (no invention), and it lands the residual on the coarse side
  (gentle 1/h_coarse). **No cell carries a spurious local divergence.**

Crucially injection is mean-preserving: the four bad fine divergences SUM TO ZERO
(global mass conserved), so the error is purely LOCAL -- which is why the smooth-
flow / global-mass gates are blind to it and only a local-divergence metric on a
tangentially-varying field (Beltrami) exposes it. And the reason B is stable while
A's accurate fix is not: averaging is a CONTRACTION (||restrict|| <= 1), but the
mean-preserving prolong is an EXPANSION (||P|| ~ 1.25, the +/-1/8 slope), and the
per-sweep-frozen block iteration amplifies expansions. The prolong CANNOT be made
contracting without discarding the accuracy it adds -- so damping (Piece 4 family)
is the only MP route, and it must target the coupled coarse<->fine P*P^T mode, not
the fine diagonal (which is why MP_SELF up to 4 barely moved it).

**The fix: make A behave like B by agglomerating the interface PRESSURE.** The
projection currently drives EACH fine interface cell's divergence to zero; since
those residuals are conservative in the mean, enforce `div = 0` on the coarse MEAN
of the covering 2x2(x2) fine cells instead -- which injection already satisfies --
so the corrector has nothing spurious to chase: no pressure dipole, no corner
mode, no blow-up. Concretely:
- The 2x2(x2) fine cells touching the interface KEEP their fine velocity faces but
  SHARE ONE (coarse) pressure unknown, relaxed against their SUMMED divergence
  (Galerkin coarsening `P_agg^T L P_agg`: SPD-preserving, a contraction -> stable
  by construction, conservative by construction). Pressure is coarse at the
  interface, velocity stays fine -- the hybrid real composite/AMR projections use.
- **Do NOT agglomerate the full cells**: merging control volumes just RELOCATES
  the 2:1 interface one fine layer inward (orientation A again). Pressure-only.

**Tradeoff and the both-worlds combo.** Pure agglomeration is a STABILITY fix: the
slaved velocity stays injection (O(h) at the interface), so the bounded error sits
near injection level (no blow-up, but not 2nd-order). To get BOTH: feed the
predictor the MP-prolonged slave (Piece 1, already validated, has the tangential
variation) AND agglomerate the interface pressure. The velocity correction is then
the shared coarse gradient (constant across the 4), which does NOT remove the MP
variation, and the projection only constrains the mean, so it never fights it ->
plausibly stable AND O(h^2). This is the most promising thing to try.

**Cheap pre-test before building it.** Temporarily make the projection NOT correct
the fine interface faces individually (enforce only the coarse-mean divergence
there) and run the t=8 / t=40 Beltrami: if the blow-up disappears, the mechanism
is confirmed and the full pressure agglomeration is worth implementing. Verify SPD
+ conservation with `MOBY_HALO_AUDIT` and the global-mass-residual gate; watch that
the grouped relaxation keeps exact nb-/rank-independence (the grouping must be
derived identically on every rank, like the existing Z-order entry enumeration).

## 6f. Implementation status — flat-face pressure agglomeration (2026-06-23)

The agglomeration of §6e is implemented behind **`MOBY_AGGLOM=1`** (off by default;
the production composite path is bit-exact to before — verified on refined Beltrami
A/B, CPU+GPU). It targets ONLY orientation A (the dirty side, `physHigh==FACE_COARSE`,
a fine block below a coarse one), which §6e identifies as the sole spurious source.

**What it does.** `pressure_solver.f90`:
- `agg_skip` (declare-target): the volume sweep's skip predicate. True for a PURE
  flat-face orientation-A cell — its sole interface is one `FACE_COARSE` high face
  and its covering 2x2 tangential group does not run into a tangential interface.
  It is true for a block's redundant HALO copy of a neighbour's interface cell too
  (the copy sits at `idx==hi(d)` in the interface direction; only the tangential
  index is the halo), so the copy is skipped and cannot contaminate the shared face
  the owner keeps quiescent — this halo-awareness is essential (a narrow interior-
  only skip leaves the halo copies chasing the local spurious divergence and the
  corrector got WORSE, B 1.0e-3 -> 1.7e-3).
- `agglomerate_sweep(dInt, gcolor)`: relaxes the pure flat groups on the dInt
  interface. The 2x2 covering fine cells share one pressure DOF, updated by a single
  `phi = -sor * (summed divergence) / (summed diagonal)` and applied to all four
  (internal tangential faces get +phi and -phi and cancel; the slaved face is
  reconstructed from the shared pressure). `sum_m denom_m` over-estimates the true
  Galerkin diagonal `P^T L P` (it double-counts the internal faces), which only
  under-relaxes — the iteration still converges to summed-div = 0. Groups are
  coloured by covering-coarse parity so same-colour groups never share a face;
  called for dInt = 1,2,3 each colour. Inert single-level (no interfaces).

**Results (Beltrami, niter=200 unless noted).**
- A (div-of-exact): **1.44 unchanged** — correct. Agglomeration does not change the
  per-cell divergence operator; it stops the corrector from CHASING the (real, but
  mean-zero) local divergence.
- B (corrector change on the exact field): single-direction slabs **7.7e-4 -> 1.5e-4
  (5x)**; 3D patch **1.0e-3 -> 4.2e-4 (2.5x)**. Direct measure of the projection's
  spurious interface activity — the user-visible "oscillation" — substantially
  reduced. (The 3D patch keeps a 4.2e-4 floor because its edges/corners are still
  per-cell; a slab, all flat face, reaches 1.5e-4.)
- C (whole flow t=8): **3.60e-2 vs baseline 3.74e-2** — at the injection floor, as
  §6e predicts (pure agglomeration is a STABILITY fix; the slaved velocity stays
  injection / O(h), so the bounded error sits near injection level, not 2nd order).
- Trajectory (max|u| excess over the exact decay, niter=10): baseline ~2.2% rising
  to ~3.8% by t=18; agglomeration ~0.9% and flat to t=13 — the interface artifact
  is roughly HALVED and its growth delayed. Past t~15 the residual growth is driven
  by the per-cell edges/corners and converges back toward baseline.

**The remaining hard part — edges/corners.** A cell on >1 interface (patch edge =
2, corner = 3) is left per-cell. Agglomerating them by the naive MULTIPLICATIVE
per-direction scheme (process the x-, then y-, then z-groups, overwriting the
corner pressure each time) is UNSTABLE — it grows the trajectory from t=0
(1.011 -> 1.058 by t=18, worse than baseline) even though it improves the single-
substep B (to 2.4e-4). The corner needs a genuinely coupled / symmetric (SPD) 2D-3D
blend, exactly the corner work flagged in §6 step 6 and §6d; it is NOT a multiplicative
sweep. Until then the flat-face agglomeration is the stable, validated increment.

## 6g. Edge/corner agglomeration is a PRESSURE fix that REGRESSES velocity (2026-06-23)

Following the §6f handout's 3-step plan, the edge/corner agglomeration was carried
through damped multiplicative -> symmetric multiplicative -> **additive Schwarz**
(`MOBY_AGGLOM_OMEGA=<w>`; unset = the validated flat-only path). The decisive
finding is that **all three schemes give the SAME niter=10 velocity trajectory**,
and that trajectory is WORSE than flat-only -- so the regression is intrinsic to
agglomerating the edge/corner pressure at all, not to how the relaxation is
ordered.

Refined-patch, niter=10, ratio = max|u|/(2*exp(-0.01t)) to t=20 (GPU, nofma):

| scheme | trajectory peak | Test B (refined) |
|---|---|---|
| off (no agglom) | 1.038 | 1.034e-3 |
| **flat-only (validated)** | **1.043** | 4.16e-4 |
| edges, damped multiplicative (w=0.25-1.0) | 1.054 | 3.60e-4 |
| edges, symmetric multiplicative | 1.054 | -- |
| edges, **additive Schwarz** (w=1.0) | **1.054** | 3.79e-4 |

Key controls: the single-block (no-interface) trajectory is FLAT at 0.999 to t=20
(growth is purely interface, not bulk niter=10 under-convergence); niter=30 off =
niter=10 off (not a convergence-rate artefact); the trajectory peak is INSENSITIVE
to w and to mult-vs-additive while Test B and interface_decay pn DO respond
(additive lowers both, as an SPD scheme should). So additive Schwarz works as
designed -- it removes the stale per-direction face reconstruction of the
multiplicative sweep (Test B 4.16e-4 -> 3.79e-4, fully stable, all flat/off/
single-level gates bit-exact) -- yet it does not move the velocity trajectory.

**Why (root cause).** A FLAT interface cell has exactly ONE slaved face, and its
agglomerated pressure IS the coarse-face pressure that face must match -- so flat
agglomeration is both stable AND accurate. An EDGE/CORNER cell has 2-3 slaved
faces feeding 2-3 DIFFERENT coarse neighbours; a single agglomerated pressure DOF
(the multiplicative overwrite, the additive sum, or even a Galerkin averaging-P
value) can match at most one of them, so the others reconstruct from the wrong
coarse pressure and carry an O(h) velocity error. That error is NOT mean-zero
(unlike the spurious divergence the agglomeration removes), so it shows directly in
max|u|. Leaving the edge per-cell keeps each face tied to its own pressure gradient
-- worse Test B but better velocity. The two goals genuinely conflict at the
multi-face cells.

**Consequence.** Pure pressure agglomeration cannot fix edges/corners on the
velocity metric -- Galerkin PtLP (step iv) would make the operator SPD but keeps
the single-DOF-per-cell limitation, so it would not beat flat on velocity either.
The velocity at edges needs the slaved faces to carry their own coarse->fine
tangential variation, i.e. the §6e BOTH-WORLDS combo: agglomerate the pressure for
stability AND feed the predictor the MP-prolonged slave velocity (Piece 1) so each
of the 2-3 faces reconstructs from a face-specific value, not the shared coarse
mean. That is the recommended next direction; until it lands, **flat-only
agglomeration is the production operating point** and the edge/corner additive
path stays opt-in (a correct, SPD, gate-passing scaffold for the both-worlds work,
the right pressure treatment if only the pressure metrics matter). Implementation:
`agglomerate_accumulate`/`agglomerate_apply` + the `aggDP` accumulator in
`pressure_solver.f90`; the projection loop runs the broad volume-sweep skip, three
snapshot accumulations and a red-black apply once per iteration.

## 7. Verification (mandatory regression, see memory `beltrami-interface-regression`)

- **Pressure-correction-alone (exact-overwrite, `MOBY_RKEXACT`):** corrector
  change to the exact div-free input must drop from 2e-2 toward ~0 (single block
  stays 0.0).
- **Divergence-of-exact (`MOBY_RKDIV` + `tools/plot_div_slice.py`):** max|div|
  must drop from 1.44 toward the linear-prolong floor (~0.05) and below.
- **Whole Beltrami t=8:** refined L2 must fall from 3.7e-2 toward the single-block
  1.4e-4; then the long t=40 run must stay bounded (no blow-up).
- **Single level + nb-independence:** channel nb=4 vs default and chanp 1-vs-8
  ranks must stay **bit-exact vs Phase 2** under `-Mnofma` (the interface row
  changes must be inert where there is no interface, and rank-independent).
- **Interface-decay gate** and the uniform-flow-through-patch exactness.

## Sources

- Almgren, Bell, Colella, Howell, Welcome, *A Conservative Adaptive Projection
  Method for the Variable Density Incompressible Navier–Stokes Equations*, JCP 1998.
- Martin, Colella, *A Cell-Centered Adaptive Projection Method for the
  Incompressible Euler Equations*, JCP 2000.
- Guittet, Theillard, Gibou, *A stable projection method for the incompressible
  Navier–Stokes equations on non-graded adaptive grids* (octree MAC), JCP 2015 —
  https://arxiv.org/pdf/2306.09957 (related nodal-projection follow-up).
- PhysBAM adaptive incompressible flow (multitude of Cartesian grids):
  https://physbam.stanford.edu/papers/stanford2013-02.pdf
