# Numerics review (2026-09-26, head `6f2c6c1`)

Scope: the momentum predictor, the passive-scalar transport (including the
conjugate branch), the pressure projection, the physical boundary conditions,
and the prepare/solve split. The RANS/LES/IDDES kernels were read only where
the momentum or scalar path calls them. Files read in full: `moby_solve.f90`,
`moby_prepare.f90`, `step.f90`, `pressure_solver.f90`, `boundary.f90`, the
transport half of `scalar.f90`, the metric builders in `init.f90` /
`blocks.f90`, the exchange API of `comm.f90`, `update_ibm_mu` /
`set_ibm_coeff` in `ibm.f90`, `docs/numerical-methods.md`,
`docs/prepare_solve_strategy.md`, `docs/conjugate/conjugate_ibm.tex` and the
STATUS header of `docs/next_session_conjugate.md`. Nothing was run except two
small Python checks (appendix). Every "should" below is a recommendation, not
a change; nothing in the tree was modified.

The nine questions are answered in order. Section 9 collects the findings
ranked by how much they matter, with the file and line of each.

---

## 1. Should `moby_prepare` do ALL preprocessing and `moby_solve` only solve?

**Yes, and the code is closer to that than it looks.** Today the solver still
has three init paths (`moby_solve.f90:94-120`): analytic `refine_body` (classify
inline), analytic `remove_solid` (classify inline), and plain. The file path
(`[ibm] coeff_file`) *also* rebuilds the leaf table from the masks and
cross-checks the file row by row (`read_ibm_coeff_file`, `ibm.f90:420`;
`check_block_table` in `field_hdf5.c`). And `moby_prepare` refuses two things
the solver accepts: a body-free case and an unset `[blocks] nb`
(`moby_prepare.f90:90-93`). So there are two grid/leaf-table builders in the
solve path (ini-driven and file-driven) plus a cross-check whose only purpose is
to guard the redundancy.

What "prepare does everything" would look like:

- **The case file becomes the single source of truth for the grid and the leaf
  table**: node lines per level, the leaf table (`blocks`), face kinds, the
  masks (kept only for diagnostics), and, when there is a body, the
  coefficient / `dwall` tiles. The solver reads the table and never rebuilds
  it; the `[grid]`/`[blocks]` ini keys become prepare inputs. The restart
  cross-check then compares two files that were both written by the code, not a
  file against an ini-derived rebuild.
- **Prepare loses its two refusals**: body-free cases produce a case file with
  no coefficient datasets; `nb` becomes mandatory (the one-block-per-rank
  layout `DIST_RANKBOX` disappears with it, which also removes the "nothing to
  narrow" caveat of `rdenomBlocks` and one whole branch of `build_block_metrics`).
- **What must stay in the solver** because it depends on the rank count:
  the Z-order split, the exchange entries (`init_block_exchange`), the boundary
  point lists, the projection's static metric tables, the device maps. All of
  these are cheap and rank-local.
- **Load balancing fits this split cleanly**: prepare knows the geometry, so it
  can write one *weight per leaf* (body block: `rdenom` recomputed every
  substage; cut cells: the 3× conjugate rate; trip blocks; refined blocks cost
  the interface kernels) without knowing the rank count. The solver replaces
  the closed-form `zorder_start/count` by a weighted prefix sum over the same
  Morton order. Results stay rank-count independent because the split only
  moves work, never arithmetic. Autotuning (`nb`, GPU mapping) is a separate
  axis and already lives in `tools/moby_tune.sh`.

Costs and risks:

- Every run needs a prepare step. For a channel that is seconds; a
  `moby_solve --prepare-if-missing` convenience (case file absent or its
  stored ini-hash stale) keeps the tutorials one-command.
- The analytic IBM coefficients are computed on the device in the solver and on
  the host in prepare (libm ulps). CLAUDE.md already declares the CPU prepare
  canonical; making it the only path retires the GPU `set_ibm_coeff` kernel and
  one bit-exactness reference moves once.
- The 7-case suite's inline analytic cases (wavy walls, Beltrami) become
  prepare+solve pairs; the P0 gates already prove they are bit-exact either way.

Net: roughly −300 lines in `moby_solve.f90`/`ibm.f90`/`blocks.f90` (the inline
classify dispatch, the legacy coefficient-file reader, the rank-box layout, the
leaf-table cross-check), one fewer concept ("which of three init paths am I
on"), and a natural home for load balancing. I would do it.

---

## 2. Correctness of the momentum, scalar and pressure equations

### 2.1 Momentum (`step.f90:130-399`)

**Time integration.** `rk_alpha/beta/gamma` (`step.f90:28-30`) are the Wray /
Spalart-Moser-Rogers low-storage RK3: α = (8/15, 5/12, 3/4), β = (0, −17/60,
−5/12), γ = α+β = (8/15, 2/15, 1/3). Correct. The predictor uses the *old*
pressure gradient scaled by `dt_gamma` and the projection uses `dt_gamma`
(`pressure_solver.f90:656`), i.e. the incremental (Le & Moin) form. Correct.

**Convection.** Divergence form on pair sums (arithmetic means), with the
skew correction `+0.25 u_i (Σ pair-sum differences)` (`step.f90:241-247`).
Checked: the pair-sum difference is 2× the stencil divergence at the u point,
so `0.25·u·2·div = ½ u div u`, and `div − ½ u div u = ½(div + adv)`. The
interpolations are exact at the cell centre and first-order at a stretched
face (the face is not the midpoint of two centres); that is the standard
choice that preserves skew symmetry on non-uniform grids (Morinishi's
non-uniform weights are not used, and the residual is documented). Correct as
a second-order energy-neutral scheme.

**Pressure gradient / divergence.** `d1x(i,VAR_U) = 1/(x_c(i) − x_c(i−1))`
and `d1x(i,VAR_P) = 1/(x_f(i+1) − x_f(i))` (`init.f90:585-590`). The
projection operator `D·diag(mu)·G` is then self-adjoint in the cell-volume
inner product. Correct and SPD.

**Diffusion — finding F1 (`init.f90:592-598`).** The three-point stencil
`lapM = 2/(hm(hm+hp))`, `lapP = 2/(hp(hm+hp))` is used for every variable in
every direction. For the *face-staggered* direction of a component (u in x) it
is exactly the conservative flux form, because the centre-to-centre distance
is (hm+hp)/2. For the *cell-centred* directions (u in y and z) it is **not**:
the conservative form is `lapP = 1/(hp·Δy_j)`, and `Δy_j ≠ (hm+hp)/2` on a
stretched line (geometric ratio r: they differ by (r−1)²/4r). Consequences:

- the molecular viscous term is not discretely conservative on stretched
  grids: the wall-normal momentum balance does not telescope to the two wall
  stresses; the SGS correction (`step.f90:508-511`, flux form on
  `turb%inv_dy`) and the scalar diffusion (flux form on `sc%invDy`) *are*
  conservative, so the three diffusion operators are mutually inconsistent;
- the operator is not symmetric, which matters the day diffusion goes implicit
  (section 7): a Jacobi/Chebyshev solve needs the SPD flux form.

Truncation order is unaffected on smoothly stretched lines (the difference is
O(h·dh)). The fix is two lines in `slice_grid_direction` and changes results
at the truncation level on every stretched case (the channels), so it needs
its own re-validation, not a bit-exactness gate.

**IBM penalization.** `qs *= mu` after the whole substage increment including
the pressure term (`step.f90:252-253`), and the projection correction is
`× mu` on every face (`pressure_solver.f90:704-711`). Inside the body
(mu → 1e-30) the velocity stays zero; in a graded cut cell the *whole* substage
update is penalized, which is consistent. `oldrhs` stores the unpenalized RHS,
so the scheme is Luchini's centre-weight correction integrated with implicit
Euler on the λ term instead of his exact factor B(λΔt). Stable for any Δt, but
first order in time at cut cells — this is the `TODO` at the top of
`moby_solve.f90`. Upgrading to the exact factor is two coefficients instead of
one: `qs = (B·q + Δ)/(λΔt + B)` with `B = λΔt/(e^{λΔt} − 1)`; solid cells
(B → 0) still give exactly zero. Cheap and worth doing (finding F7).

**Time-step limiter (`step.f90:695-709`).** `cflmax` is compared against the
*maximum* single-component CFL, not the sum over directions. RK3 with
second-order central convection is stable for `Σ_d CFL_d ≤ √3`; the tutorials
run `cflmax = 0.8`, so the safety margin is real only while one direction
dominates (channels, boundary layers). Near a stagnation region or a nose,
where two components are comparable, the sum can reach 1.6. Not a bug, but
`cflmax` should be documented as a per-direction number, or the limiter
should use the sum (finding F8).

### 2.2 Passive scalar (`scalar.f90:2772-3368`)

Called before `momentum` with the start-of-substage velocity and this
substage's `nut` — the right choice, and the reasoning in the caller
(`moby_solve.f90:364-378`) is correct. RK3 memory (`oldrhs`) is the momentum's
scheme verbatim. Convection: divergence form on the six face velocities with
arithmetic-mean face values; the optional `skew`/`advective` subtraction uses
the same pair sums. Diffusion: flux form on `invDx`, the face diffusivity
`1/(Re Pr) + ½(ν_t,L + ν_t,R)/Pr_t`. IBM dirichlet mode: the same implicit
`mu_s` as momentum with `coef_p/Pr`, which is the scalar's own Luchini λ
because the stored coefficient carries 1/Re. All correct.

One design asymmetry worth stating: the momentum was moved to skew form
because divergence form is energy-neutral only for a *discretely*
divergence-free advecting field, which a 12-iteration projection never
delivers. The scalar default (`convMode = SC_CONV_DIV`, `scalar.f90:148`) is
the form the momentum team found unstable at 2:1 interfaces. It has been
gated in the refined scalar cases (`ibmwavyr`) and it is conservative, so this
is a documented trade, not a defect — except in the conjugate mode, where it
interacts with the cut-face masking (section 5, finding F2).

### 2.3 Pressure projection (`pressure_solver.f90:213-344`)

Each iteration: residual `phi = −ω div(q)/denom`, optional Chebyshev combine,
phi halo exchange (interface-row restrict), `p += phi/dt_γ`, `q_face +=
(phi_L − phi_R)·d1f·mu`, then a halo refresh. The velocity is corrected
every iteration, so the divergence of the *current* velocity is the residual:
errors do not accumulate across steps (the projection is self-correcting),
which is what makes a fixed small `niter` viable. The operator is
`D diag(mu) G` with `G = −Dᵀ` in the volume inner product, including the
outlet pair (2·d1f in the diagonal, d1f against the mirrored ghost) and the
2:1 composite (2/3, 4/3) metrics. SPD as claimed.

**Chebyshev recurrence.** `alpha/beta/gamma` (`pressure_solver.f90:300-311`)
are Saad's Algorithm 12.1 with `α = 2ρ/δ`: verified term by term. Correct.

**Chebyshev bounds — finding F3 (`pressure_solver.f90:185-189`).** `lmax` is
set to the Gershgorin bound 2.0 *exactly*, and the checkerboard mode of the
Jacobi-preconditioned Poisson operator attains λ = 2 exactly on a uniform
periodic grid. The Chebyshev residual polynomial satisfies |T_k(−1)| = 1, so
`|P_k(2)| = 1/T_k(d/c) ≈ 1/(1 + k²·lmin)`. Computed for N = 288:

| bounds | k | \|P(2)\| | \|P(1.9)\| | \|P(10·lmin)\| |
|---|---|---|---|---|
| lmin auto, lmax = 2.0 | 12 | 0.989 | 0.637 | 0.889 |
| lmin auto, lmax = 2.2 | 12 | 0.477 | 0.093 | 0.899 |
| plain Jacobi ω = 0.8 | 12 | 0.002 | 0.007 | ≈1.000 |

So the Chebyshev projection as configured **does not damp the 2-Δx divergence
mode at all**; it relies on physical viscosity to remove it. That is stable
in periodic channels (the mode is not excited there) and it is exactly the
class of failure recorded for the boundary layer: "chebyshev + niter 6 +
Dirichlet-p outlet + dt ~ 0.5 unstable, 2-Δx pressure mode, plain Jacobi
niter 6 stable" (CLAUDE.md, B0). The standard remedy is `lmax = 1.1 × bound`
(what every multigrid Chebyshev smoother does), or one damped-Jacobi sweep
appended to the polynomial. This is a one-line config experiment
(`cheb_lmax = 2.2`) on the B0 case before any code changes. I rate the causal
link as likely but untested.

The table also says something about what 12 Chebyshev iterations buy: the
smooth modes are reduced by ~1 % per iteration (|P| ≈ 0.89–0.99 for λ ≤
10·lmin) — the projection is a smoother, not a solver, at any `niter` a
production run can afford. The residual divergence that the skew form was
introduced to tolerate is therefore permanent by design. The block structure
makes a two-level geometric multigrid (coarse level = each leaf at nb/2; the
2:1 restrict/prolong entries already exist in `comm.f90`) the natural next
step if the projection is ever to converge rather than smooth. That is a
project, not a fix; noted here because it bounds what every other pressure
tweak can achieve.

---

## 3. Boundary conditions: per RK stage, per projection iteration?

**Per RK stage: required.** The tangential ghost (`ghost = 2v − interior`)
and the pressure ghost are functions of the interior, and the interior changes
every substage. With time-independent boundary data (every current case) the
value written is the same; with time-dependent inlets the predictor should
see the data at the substage time. Applying the end-of-step BC to the
intermediate velocity `u*` is the standard incremental-projection choice and
costs O(Δt) in the wall tangential slip; second order overall.

**Per projection iteration (`pressure_solver.f90:326`, `:1021`): not needed
for any consistent face kind.** Inside the loop the only consumer of halo or
ghost data is the divergence, which reads the cell's own six faces (indices
1..nb+1) and the `+axis` normal-component halo that `sync_divergence_halos`
refreshes. Going through `apply_bc` (`boundary.f90:670-720`) face kind by face
kind:

| face kind inside the loop | what `apply_bc` writes | read by the loop? |
|---|---|---|
| Dirichlet normal velocity (wall, inlet) | the same pinned face value | no change (idempotent) |
| outlet normal velocity (BC_OUTFLOW) | skipped (`ofc = 0`) | — |
| tangential Dirichlet/Neumann ghosts | refreshed from the corrected interior | no (momentum reads them next substage) |
| pressure ghosts | copy / mirror of interior | no (the predictor never reads a physical p ghost: `momentum_face_start` skips the wall face, the outlet high face is not predicted) |
| **Neumann NORMAL velocity** (`boundary.f90:686-696`) | `q(nb+1) = q(nb) + dn·value` | **yes: it changes the last cell's divergence** |

So the per-iteration call is dead work for every face kind except the last
row, and the last row is a **consistency defect (finding F4)**: with the
normal face copied from the interior face every iteration, the last cell's
x-faces move together, its divergence does not respond to its own `phi`
(true diagonal 0 in that direction while `compute_rdenom` counts `d1f·mu`),
and the coupling to the neighbouring cell is one-sided. The projection
operator is then not symmetric at that row, which is what neither the
Chebyshev bounds nor the SPD argument cover. The only ini using it is
`tutorials/sailplane/input.ini:66` (`x_max_u_type = neumann`, with Neumann p
everywhere — the legacy "leaky outflow"); it runs red-black at `niter = 20`,
which tolerates it. The A0 outlet patch (`x_max_patch = outlet`) is the
consistent replacement and should be adopted there; after that, a Neumann
type on a normal velocity component should be a config error, or be honoured
only in the `outflow_copy` call.

Once F4 is closed, `apply_bc` can leave the projection loop entirely: one
call after the last `jacobi_apply`, before the final full exchange (the
exchange's tangential extension copies a neighbour's ghosts into edge halos,
so the BC write must precede it — that ordering is already what the loop
does). That removes `(niter − 1)` launches per substage — 33 per step at
`niter = 12`, i.e. roughly the 4 % `projection_bc` bucket. Same for the
red-black path (per colour).

---

## 4. A unified ghost approach with BC "providers"

> **STATUS (2026-09-27): done in REDUCED form** -- section 10 step 6. The
> three kernels became one affine write `dst = w*src + C` from tables
> resolved at init (velocity, pressure and scalar rows as columns of the
> same tables), `apply_scalar_bc_q` is gone; the rows did NOT move into the
> exchange entries, for the reasons measured/read there (the same-level copy
> entries read physical ghosts and pack runs first, so a BC launch is needed
> before the exchange in any design).

The exchange already *is* a provider model: every halo point is a weighted
gather from source rows given by per-dimension affine maps
(`entry_gather_map`, `comm.f90:1052`), and the same kernel serves same-level
copies, restrictions, prolongations and the blended pressure ghost. A physical
boundary ghost is a degenerate case of exactly that map: one source (the
interior cell, `ga = −1` reflected index or `ga = 1`), one weight (−1 for
Dirichlet-mirror, +1 for Neumann/copy) and a per-point constant (`2·value`
or `dn·value`). The gather kernel has weights (`lWp`) but no constant term, so
the change is: one constant array indexed by point, filled from
`pointBcValue`. Then:

- `apply_bc`, `apply_scalar_bc` and `apply_scalar_bc_q` collapse into entries
  (three kernels and two point-list formats become zero kernels and one
  format); `bc%pointFace/slot/i/j/k` become ordinary entry boxes;
- corner and edge ghosts at physical faces are defined by construction (today
  they come from the neighbour block's `apply_bc`'d ghost through the
  tangential extension, which is correct but only because of call ordering);
- the divergence round, the copy-only round and the full round each pick the
  entries they need, so the BC "cost" follows the same prefix/compaction
  logic as everything else, and the ordering hazard of section 3 disappears.

What cannot be a ghost: the **normal velocity on a physical face** is a face
DOF on the boundary, not a halo cell. Two facts make it cheap anyway: (a) the
predictor and the projection never write a pinned face, so a Dirichlet normal
value needs to be written once — at init, restart, or when its value changes
— and not every substage; (b) the outlet face is corrected by the projection
and copied once per substage by `outflow_copy`. So the face writes reduce to
one tiny kernel per substage for outflow faces only.

Would it be faster? Marginally: it removes ~40 launches per step at `niter =
12` (the 3–4 % BC buckets) and puts the BC points into kernels that run
anyway. The real gain is one abstraction instead of three and a corner
treatment that is right by construction. I would do it after the F4 fix and
after the projection-loop `apply_bc` is gone, in that order, so each step is
gated bit-exact.

---

## 5. The conjugate heat transfer treatment

**The scheme is sound where it is one-dimensional, and the design note is
honest about the rest.** The single face coefficient `dm/(w/κ_L +
(1−w)/κ_R + R_c·dm/h)` on `w = φ_L/(φ_L − φ_R)` is the Patankar series
resistance; it enforces temperature and flux continuity by construction,
reduces to Luchini's centre-weight correction as κ_s → ∞ and to zero flux as
κ_s → 0, and the obliquity lemma (both perpendicular distances shortened by
the same cosine) is exact for a plane and O(h·curvature) otherwise. The C3
fluid-fraction capacity `f + (1−f)C_s` with a plane reconstruction from ∇φ is
the correct volume average and is what makes transients second order (gated
1.99/2.00 against a 1.71/1.09 control). Conservation to round-off follows
from writing everything as face fluxes. The dropped tangential term is
measured, not assumed, and the C2 verdict (worse than dropping it at high
contrast, because the discrete `s_t` does not vanish when the true one does)
is the right reading of the numbers. The time-step analysis is correct: the
interface is not stiff (R ≥ h/max κ), the binding limit is the solid
diffusivity and the Gershgorin factor of a cut cell.

**Finding F2 — the cut-face convective mask is inconsistent with continuity
in the default divergence mode (`scalar.f90:2969-2985`, `:3002-3009`).**
Convection is masked on every solid-node face *and every cut face*. A cut
face whose staggered node lies on the fluid side (w > ½) carries a genuine
fluid velocity `u_f ≈ δ·∂u/∂n` (the graded penalization is designed to give
exactly that). Dropping its flux leaves the fluid cut cell with the other five
fluxes, whose sum is `−u_f`, so in divergence form a uniform scalar changes at
the rate `s·u_f/h ≈ s·(δ/h)·∂u/∂n` in every such cell. In wall units with the
first cell at h⁺ the ratio of this spurious rate to the cell's own diffusive
coupling is roughly `0.02·Pr·h⁺²` (δ/h averaged over cut positions ≈ 1/8,
∂u⁺/∂y⁺ = 1, six-face coupling `6/(Re Pr h²)`), i.e. one to a few per cent at
h⁺ ≈ 1–2. The
`skew` mode halves it; only `advective` (`f = 1`, the divergence built from
the same masked velocities, `:3002-3009`) removes it — and then the same
quantity reappears as a conservation defect, which is the honest place for
it, since the mass that crosses that face turns inside the fluid sliver of the
solid-centred cell and the grid cannot represent that. It is first order and
local either way, but the default should be the bounded, uniform-preserving
form. The CHT tutorials (`tutorials/cht/pipe/*.ini`, `channel/*.ini`) set no
`[scalar] convection` key, so they run divergence mode; the flat channel is
immune (the wall lies on a face, all cut-face nodes are solid), the pipe is
not. No gate sees it: conservation budgets are satisfied by construction and
the uniform-scalar gates run without a curved conjugate body. Test: pipe
flow, `θ ≡ 1` initial condition, conjugate on, divergence vs advective;
the fluid cut cells drift in the first and not in the second. Then decide the
default for conjugate scalars (I would force `advective`) and re-check the
pipe Nusselt numbers.

Smaller points:

- `ν_t` is excluded from cut faces. Correct with a resolved wall (WALE →
  0); the config error for wall functions guards the other case. Fine.
- The band material (F2 shell) and the `CONJ_MIN_COSINE` grazing guard
  (`w = ½` for near-tangent arms) are reasonable; the guard's threshold is a
  tuning constant whose sensitivity is not recorded.
- The volumetric source is weighted by `vfrac` on both sides. Correct.

---

## 6. Performance-neutral simplifications of the numerics logic

Ordered by (simplification gained)/(risk):

1. **Drop `apply_bc` from the projection loop** (after F4): −33 launches per
   step, −1 ordering hazard, and the red-black path loses its per-colour BC
   call too. Bit-exact once no Neumann-normal face exists.
2. **Stop rewriting pinned normal faces every substage**: they are written at
   init/restart and never modified by anything; the per-substage write is a
   no-op in disguise. Bit-exact.
3. **One BC mechanism instead of three** (`apply_bc`, `apply_scalar_bc`,
   `apply_scalar_bc_q`): section 4. Bit-exact by construction if the
   arithmetic is `2v − q` and `q + dn·v` as today.
4. **One init path instead of three** (section 1): the inline classify
   dispatch, the legacy global-layout coefficient reader
   (`ibm.f90:454-473`), the rank-box layout, the leaf-table cross-check.
5. **The three `interface_correct` kernels → one** (already noted in
   CLAUDE.md as ~1 ms/step at 16 ranks): one launch over a face list instead
   of three plane loops.
6. **Split the scalar kernel by mode.** `scalar_transport_kernel` carries the
   S1/S2/S3/S5a/C1/C2/F2 branches in one body with ~60 privates; the
   `jacobi_apply` register lesson says the dormant branches cost occupancy
   even when not taken. A conjugate kernel and a plain kernel, each half the
   size, is simpler to read and probably faster; bit-exact if the expressions
   are moved verbatim.
7. **`qs → q` copy** (`step.f90:382-397`): a full three-component copy per
   substage (6 doubles/cell of traffic). Avoidable by double-buffering `q`
   with a substage index (`q(:,:,:,:,:,cur)`), but the device-mapped flat-array
   convention makes a descriptor swap awkward and the exchange kernels index
   `q` directly. Worth measuring, not obviously worth doing.
8. **Chebyshev `cheb_combine`** could be fused into `jacobi_compute_phi`
   (one launch, one fewer read of `phi`); bit-exact if the operation order is
   kept.

Not recommended: fusing the cross-level copy into the same-level kernel
(register pressure on every round, already measured) and any
`maxregcount` games (already measured to backfire).

---

## 7. Optional implicit diffusion with the existing iterative machinery

**Feasible, but the iterative route only pays for moderate stiffness; the
stiff cases want a per-block line solve. And the strongest practical case is
the conjugate solid, not the momentum.**

Where the explicit diffusion limit actually binds today (appendix A):

- `tutorials/channel_kmm180`: natural stretching with `natural_dyw_plus =
  0.05` gives a first cell of 0.053 wall units. The Péclet limiter then
  caps `dt` at 7.7e-6, forty times below the ini's own 3.125e-4 and 230×
  below the CFL limit. Either the wall spacing is ten times finer than a
  wall-resolved channel needs (dy⁺ ≈ 0.5–1 is usual, which relaxes the limit
  100–400×) or this case is exactly the one an implicit wall-normal treatment
  is for. Check what the archived runs actually stepped at before choosing.
- The CHT pipe cases with `solid_rhocp = 0.0625` (α_s/α_f = 16): the solid's
  diffusive limit is 16× the fluid's and sets `dt` for the whole run; the
  fluid cut cells are not affected (the harmonic mean keeps κ_face ≤ 2 at
  κ_s = 1).
- IBM/refined airfoil cases: convection-limited by a wide margin; implicit
  diffusion buys nothing.

Options, with what each costs:

**(a) Chebyshev-Jacobi Helmholtz solve, a twin of `pressure_projection`.**
The operator `I − dt_γ·ν·L` is SPD (once F1 puts the momentum Laplacian in
flux form) and its Jacobi-preconditioned spectrum lies in `[1/(1+2σd), 2]`
with `σ = ν dt/h²_min`, `d` = dimensions. It fits the exchange and kernel
shapes exactly (~250 lines reusing `exchange_scalar_halos` and the metric
tables), and the pressure update needs no rotational correction at second
order in velocity (pressure carries an O(ν dt) wall error either way; Kim &
Moin's form). But the iteration count scales with √(1+6σ): at σ ≈ 3–5 (a
5–10× relaxation of the explicit limit) ~20–25 iterations reach 1e-4 on the
increment, comparable to today's whole projection; at σ ≈ 100 (the KMM wall
cell at its ini `dt`) it is hopeless. So this route buys about a factor 5 in
`dt` on diffusion-limited cases at roughly the projection's cost. Not
attractive for momentum.

**(b) Per-block line-implicit in the stretched direction only (wall-normal
y), block-Jacobi across block faces.** Exact tridiagonal solves per (i,k)
line inside each block (one thread per line, nb = 32–64 unknowns), halo
values lagged, one halo exchange and optionally a second pass. The coupling
across a y block face is weak precisely where it matters (the stiff cells are
in the wall block, the block face sits nb cells away where dy is large). This
is the Kim–Moin/Le–Moin approximate factorization with the factorization
error O(Δt²) inside RK3, restricted to one direction, and it is the only
method that handles σ ≈ 100. Cost: a tridiagonal kernel (~150 lines), one
extra exchange per substage, and the F1 flux-form Laplacian so the solve is
consistent with the explicit x/z parts. Moderate complexity, contained in
`step.f90` plus one kernel.

**(c) Implicit solid conduction for conjugate scalars, fluid explicit.** The
solid operator is static (coefficients never change), has no convection,
lives on a minority of cells, and its Chebyshev bounds can be computed once.
A Jacobi/Chebyshev Helmholtz solve restricted to solid cells with the fluid
values as lagged Dirichlet data at cut faces removes the `C_s` limit entirely
and keeps the fluid at its own explicit limit (which is where the harmonic
mean already caps the cut-face coefficient). This is the note's own option 3
(`conjugate_ibm.tex` §"What does constrain the time step") and it is the
cheapest of the three to build (the operator is a subset of the existing
transport kernel's diffusion line) with the largest measured payoff on the
CHT campaign. I would do this one first.

None of the three is a large code change; (b) and (c) do not complicate the
projection at all. Prerequisite for all: F1.

---

## 8. Higher-order convection (4th-order central, QUICK) for momentum/scalar

**Consistency with the pressure equation.** The projection enforces `D u = 0`
with the second-order `D`; the predictor can use any convective operator
without breaking that. What a 4th-order convective term cannot do is raise
the *formal* order: the pressure gradient `G = −Dᵀ` stays second order, so the
scheme is second order with reduced dispersion/aliasing error — a real
benefit, but not "fourth order". Full fourth order (Morinishi 1998) needs the
wide-stencil `D` and `G` together, a wide-stencil Poisson operator, and a
consistent interface transfer; that is not compatible with the 2:1 machinery
as built and I would not attempt it. Energy neutrality is preserved only if
the 4th-order skew form is built with the matching 4th-order interpolation and
difference pair; QUICK is upwind-biased and not energy-neutral by design, so
it re-opens the interface instability that skew closed. QUICK/TVD is
appropriate for RANS scalars (today first-order upwind) and for bounded
passive-scalar fronts, not for DNS/LES momentum.

**The blocker is halo depth, not the pressure.** Every candidate stencil
needs a second halo layer, and `comm.f90` is single-layer everywhere: entry
boxes, gather maps, the 2:1 restrict/prolong rows, the divergence subset,
`phiIfaceRow`. That is the same increment S5b/TVD and the RANS van Leer have
been waiting for; the CLAUDE.md log already carries its measured motivation
(the SD7003 γ front). So the honest plan is: 2-layer halo first (comm-only,
gated bit-exact because nothing reads the new layer), then per-cell order
switching.

**Order reduction where the stencil cannot reach.** Fourth order across a 2:1
interface would need a two-layer prolongation of at least quadratic accuracy
to keep its order; near a physical wall and inside the IBM band the wide
stencil reads pinned/penalized values. Standard practice is a mask (a band of
one cell around interfaces, physical faces and the near-body band that
`init_ibm_band` already builds) that selects the second-order flux, applied as
`F = F2 + m·(F4 − F2)` so the masked cells are bit-exact with today.

**Cost.** The momentum kernel already moves ~17 doubles/cell and is
occupancy-limited; a 4th-order version roughly doubles the loads. Expect a
sizeable predictor slowdown unless the per-component split that CLAUDE.md
lists as the next momentum idea is done first.

Recommendation: yes for scalars (TVD/QUICK-class, after the halo increment,
starting with the RANS scalars where the gain is measured); no for momentum
beyond an experimental 4th-order central toggle, because the pressure caps the
order and the interface machinery caps the stencil.

---

## 9. Problems found, ranked

| # | severity | what | where |
|---|---|---|---|
| F2 | correctness (local, first order) — **FIXED 2026-09-27** (section 10 step 2: conjugate scalars default to `advective`, explicit `divergence` refused; measured drift 4.1e-5 → 0.0; conservation gate re-baselined to the O(h) leak) | conjugate cut-face convective mask inconsistent with continuity in divergence mode; uniform scalar not preserved at fluid-node cut faces of curved bodies; CHT pipe tutorials ran this mode | `scalar.f90:2969-2985`, `:3002-3009` |
| F3 | ~~stability (likely the B0 outlet instability)~~ **RETRACTED 2026-09-27 as a cause** (section 10 step 1: the B0 instability reproduces on the B0-era binary with e-fold 33–37 t.u. and `lmax = 2.2` grows FASTER, 27–32 t.u.; `main` is stable at niter 6 with either bound). What stands is the polynomial fact: the auto `lmax` leaves λ = 2 undamped (\|P₁₂(2)\| = 0.99 vs 0.48 at 2.2), a property of the smoother with no measured consequence. No code change. | `pressure_solver.f90:185` |
| F4 | consistency (SPD lost at one row) — **FIXED 2026-09-27** (section 10 step 3: sailplane on the outlet patch; a Neumann normal velocity is a config error) | Neumann normal-velocity faces rewritten inside the projection loop; was used by the sailplane outflow | `boundary.f90 resolve_face_bcs`, `tutorials/sailplane/input.ini` |
| F1 | consistency/conservation on stretched grids — **FIXED 2026-09-30** (section 10 step 8: flux form on the variable's own control volume; `stencil_test` 48/48; turb180, Blasius and KMM180 re-validated) | molecular viscous term non-conservative in cell-centred directions; inconsistent with SGS and scalar diffusion; non-symmetric | `init.f90:592-598` |
| F5 | dead work — **FIXED 2026-09-27** (section 10 step 4: one `apply_bc` per projection, both solvers; max_abs 0 at production flags; bc bucket −91.6 %) | `apply_bc` in the projection loop and pinned-face rewrites each substage | `pressure_solver.f90` projection loops |
| F6 | configuration — **FIXED 2026-09-27** (section 10 step 5: archived run confirmed at dt 7.696e-6; ini to `natural_dyw_plus = 0.5`, dt 3.125e-4 measured) | `channel_kmm180` first cell 0.053 wall units → Péclet-capped `dt` 7.7e-6 vs the ini's 3.125e-4 | `tutorials/channel_kmm180/input.ini`, README |
| F7 | accuracy (first order in time at cut cells) | implicit-Euler penalization factor instead of Luchini's exact B(λΔt) | `ibm.f90:1307-1350`, `step.f90:252` |
| F8 | documentation/safety — **DONE 2026-09-27** (section 10 step 5: documented in configuration.md; init-time print of the worst-case directional sum) | `cflmax` compared against the max component, RK3 limit is on the sum (√3) | `step.f90 get_timestep_rates` |
| F9 | design | the scalar default (divergence) is the form the momentum abandoned; fine while gated, but it is why F2 bites | `scalar.f90:148` |

Suggested order of work: F3 (one config experiment — DONE, retracted), F2 (one test run, then
a default — DONE), F4 (sailplane to the outlet patch, config error — DONE), F5 (DONE) and the
section-4 unification (DONE in reduced form, step 6), F1 (its own validation round on the stretched cases),
F7, then the section-7 (c) implicit solid if the CHT campaign is to be re-run
at `niter = 12`.

### Plan item: F1, the flux-form viscous stencil (added 2026-09-27; DONE 2026-09-30, section 10 step 8 MEASURED)

**Do the two forms coincide in this code?** Only where `Δ = (hm+hp)/2`, and
that is guaranteed in exactly one place: a component's face-staggered
direction, because `cell_center_at` (`init.f90:722-730`) is by construction
the midpoint of two nodes. In a component's cell-centred directions the
control volume is the cell itself, `Δ_j = node(j+1) − node(j)`, while
`(hm+hp)/2 = (Δ_{j−1} + 2Δ_j + Δ_{j+1})/4`; these agree only on a `uniform`
line. Every stretched distribution (`natural`, `geometric`, `tanh`, `cosine`,
`blayer`, `nodes_file`) breaks it. Measured on the KMM180 natural line: the
code's `lapP` is 0.913× the flux weight at the first wall cell (local ratio
1.66), the volume-weighted operator is 3.2 % asymmetric there, and the
Poiseuille wall-normal budget misses the telescoped wall fluxes by 2.8e-5
relative. On a geometric line with r = 1.10 the coefficients are 0.23 % off
everywhere and the budget error is 2.3e-3.

**The change** (`slice_grid_direction`, `init.f90:592-598`): one formula for
every variable,

```
lapM(i,var) = d1(i,var)/hm
lapP(i,var) = d1(i,var)/hp
lap0(i,var) = -(lapM(i,var) + lapP(i,var))
```

where `d1(i,var)` is the inverse width of the variable's OWN control volume,
already computed on the lines above (`1/(x_c(i) − x_c(i−1))` for the
staggered direction, `1/(x_f(i+1) − x_f(i))` otherwise). In the staggered
direction this is mathematically the old `2/(hm(hm+hp))`; keep the old
expression textually there if a bit-exact uniform-grid gate is wanted,
otherwise accept a few ulps.

**What it buys**: the molecular momentum diffusion, the SGS correction and
the scalar diffusion then share one face gradient at every face; the
wall-normal operator becomes symmetric in the cell-volume inner product
(the precondition for any implicit diffusion, section 7); the integrated
momentum balance telescopes exactly to the wall stresses. Accuracy stays
second order (the flux form trades pointwise exactness on quadratics for
conservation).

**How to gate it**: it is a numerics change, not a refactor. Uniform-grid
cases (Beltrami, min_channel, the IBM/airfoil cases) must be bit-exact if the
staggered expression is kept textually, else ≤ 1e-14. Stretched cases move
at truncation level: re-run `validation/rans_sst` (natural y), the
turbulentBoundaryLayer Blasius gate (geometric y, `compare_blasius.py`:
expect θ/H/du errors to move by less than their recorded bands), and a short
KMM180 statistics comparison against its archived reference. A conservation
gate is new and cheap: on a laminar Poiseuille channel, the column sum of the
discrete viscous term must equal the two wall fluxes to round-off.

**When**: with the section-7 implicit-diffusion work, or in the one-by-one
re-validation pass of the stretched-grid cases, whichever comes first. Not
in isolation, because the validation it needs is the same run either way.

---

## Appendix A — the two numbers computed

KMM180 natural line (`natural_wall_coordinate`, ny = 136, jb = 16, dyw = 0.05):
y⁺_max = 179.4, first spacing 0.0526 wall units, centre spacing 3.79.
Péclet rate (1/Re)/Δy²_min = 6.50e4 → `dt ≤ 0.5/rate = 7.7e-6`; the ini's
3.125e-4 corresponds to σ_wall = 20.3. CFL at 3.125e-4 with u⁺ = 20 is 0.14;
the CFL-limited step would be 1.75e-3.

Chebyshev residual polynomial `P_k(λ) = T_k((d−λ)/c)/T_k(d/c)` with
`lmin = (2/3) sin²(π/288) = 7.9e-5`: the table in section 2.3. Plain damped
Jacobi factor `(1 − ωλ)^k` at ω = 0.8, λ = 2, k = 12: 2.2e-3.

---

## 10. Execution sequence, one item at a time (added 2026-09-27)

Each step is self-contained, ends with a gate, and names what it depends on.
"Bit-exact" means `max_abs 0` on the 7-case suite (+ the scalar legs) at
production flags when the arithmetic source is untouched, at `nofma` when an
expression moves (`compile.sh cpu_nofma/gpu_nofma`, `submit_nofma_gate.sh`).

**Step 0 — reference set.** Cut CPU+GPU `nofma` binaries of `moby_solve` and
`moby_prepare` from the current COMMIT (not a working tree), record the hash
in a PROVENANCE file, as `~/s5c_ref_binaries/` was. Half an hour. Everything
below compares against it.

> **MEASURED (2026-09-27).** DONE. `~/numrev_ref_binaries/` cut from commit
> `6db60ef` (clean tree, `git status` empty), nvhpc 25.9, on ISTM-io2:
> `moby_solve`/`moby_prepare` × {cpu, gpu} × {production, nofma}, eight
> binaries with md5s in the PROVENANCE file. It supersedes `~/s5c_ref_binaries`
> for work on `main`, whose binaries predate the skew lockdown and the
> `[scalar] convection` key and cannot run today's inis.

**Step 1 — F3, Chebyshev `lmax` (experiment, no code).** On the boundary-layer
case that produced the B0 instability (`tutorials/turbulentBoundaryLayer`,
`accel = chebyshev`, `niter = 6` and `12`, `dtmax = 0.5`), run the default
bounds against `cheb_lmax = 2.2` and log `max|div|` and `max|pn|` per step
over ~100 time units. Prediction: the default e-folds at ~36 t.u., the 2.2 run
does not. Control: `min_channel` with both settings, divergence residual
history only (the field moves at truncation level, which is expected). If
confirmed: set the auto bound to `1.1 × 2.0` (`pressure_solver.f90:185`),
note in `numerical-methods.md` that `cheb_lmax` is a safety margin, and
re-gate the Chebyshev cases at the physics level (every Chebyshev ini moves,
by construction). One session of runs, one line of code.

> **MEASURED (2026-09-27). PREDICTION CONTRADICTED — no code change; the
> causal claim in section 9 is retracted (see F3 there).** The B0 case was
> recovered from git (`4f819fc^:tutorials/turbulentBoundaryLayer/topbc_outlet/
> {template,blasius2d}.ini` + `make_blasius_ic.py`; 384×160×4, Blasius inlet,
> outlet x_max AND y_max, `dtmax = 0.5`, 4000 steps = 2000 t.u., snapshots
> every 100 steps; `max|div|` = interior max of the staggered divergence,
> `rms2dx` = the 2-Δx content of the z-mean u row at y ≈ 23 as in the
> original `wiggle_rms.py`), and run on TWO binaries on the local RTX 3060:
> `main` at `e126d1a` and the B0-era binary `5d44eb4` (the commit whose README
> records the instability; it still had the `cheb_lmax` key).
>
> | binary | run | max\|div\| t=300 / 700 / 2000 | rms2dx t=300 / 700 / 2000 | verdict |
> |---|---|---|---|---|
> | `5d44eb4` (B0 era) | cheb, niter 6, lmax auto (2.0) | 2.7e-5 / **0.98** / 1.03 | 3.5e-5 / **0.092** / 0.146 | e-folds, saturates at O(1) div |
> | `5d44eb4` | cheb, niter 6, **lmax 2.2** | 3.9e-5 / **0.99** / 1.01 | 3.8e-5 / **0.043** / 0.050 | e-folds, same schedule |
> | `5d44eb4` | cheb, niter 12, lmax auto | 2.1e-6 / 8.9e-7 / 1.4e-8 | 3.0e-5 / 3.1e-5 / 3.1e-5 | stable |
> | `5d44eb4` | plain Jacobi, niter 6 | 8.7e-5 / 1.1e-4 / 6.8e-5 | 4.0e-5 / 2.2e-5 / 3.3e-5 | stable |
> | `e126d1a` (main) | cheb, niter 6, lmax auto | 2.8e-5 / 1.3e-5 / 5.9e-6 | 3.4e-5 / 3.1e-5 / 3.1e-5 | stable |
> | `e126d1a` | cheb, niter 6, lmax 2.2 | 2.3e-5 / 1.2e-5 / 9.3e-6 | 3.2e-5 / 2.9e-5 / 3.2e-5 | stable |
> | `e126d1a` | cheb, niter 12, lmax auto | 1.1e-5 / 9.5e-7 / 2.7e-8 | 2.6e-5 / 3.1e-5 / 3.1e-5 | stable |
> | `e126d1a` | cheb, niter 12, lmax 2.2 | 2.7e-6 / 1.6e-6 / 2.2e-8 | 3.0e-5 / 3.1e-5 / 3.1e-5 | stable |
> | `e126d1a` | plain Jacobi, niter 6 | 1.1e-4 / 1.3e-4 / 4.2e-5 | 3.8e-5 / 2.5e-5 / 3.9e-5 | stable |
>
> Growth in the window t = 350–600 on the B0 binary: e-folding time **32.6
> t.u. (max|div|) / 36.9 (rms2dx)** at the default bound — the recorded "~36
> t.u." reproduced — and **26.8 / 31.9 t.u. at `lmax = 2.2`**: the enlarged
> bound does not damp the mode, it grows slightly FASTER. So the mode is not
> the undamped λ = 2 checkerboard of the Chebyshev polynomial (which 2.2 would
> have cut from |P| 0.99 to 0.48 per projection); what removes it is `niter`
> (12 is stable on both binaries, as the B0 README's niter 60 was) and, on
> `main`, whatever changed between `5d44eb4` and now — `main` is stable at
> niter 6 with either bound. The one identified physics change in that span is
> the momentum convection form: `5d44eb4` has no skew form at all (`git grep
> skew 5d44eb4 -- src` is empty; the port is `df89591`, 2026-07-23) and ran
> DIVERGENCE convection, which `main` can no longer run (error stop). That is
> a correlation, not a proof; the retraction rests on the 2.2 run alone.
> Control (`tutorials/min_channel/input_gpu.ini`, 200 steps, per-block
> divergence residual, GPU): lmax 2.0 vs 2.2 read max|div| 1.60e-3 vs
> 1.62e-3 at step 200 (rms 3.1e-4 vs 3.4e-4), interleaved over the whole
> history — indistinguishable, as the polynomial table predicts for the
> smooth modes the residual is made of. The recovered case, inis and the
> history/fit scripts are not committed (the B0 case was retired on
> 2026-09-26); the recipe above rebuilds them in minutes.

**Step 2 — F2, conjugate convective mask (test, then a default).** Cheapest
curved conjugate case in `validation/conjugate/` (the cylinder-type oblique
case, or the pipe at reduced resolution), `θ ≡ 1` initial condition, flow on,
conjugate on: run `[scalar] convection = divergence` and `= advective`, report
`max|θ − 1|` restricted to fluid cut cells after ~1000 steps. Prediction:
O(1e-2) drift in divergence, round-off in advective. Then: in
`validate_conjugate_config`, resolve an UNSET `convection` to `advective` when
any scalar is conjugate, and reject an explicit `divergence` with a message
that names this finding. Gates: the 9-case scalar suite bit-exact (nothing
non-conjugate changes); C1–C3 re-run, expecting ONE gate to move — the
`Σ C θ dV` drift with the flow on becomes an O(h) cut-cell leak instead of
1e-16, and its recorded band must be re-baselined, with a note that the leak
is the same quantity the divergence form put into the field. Finally re-run
the pipe Nusselt comparison against the Neuhauser reduction. One to two
sessions, the pipe run on a remote GPU.

> **MEASURED (2026-09-27). PREDICTION CONFIRMED in kind; default changed.**
> (a) Test on `validation/conjugate/wavy.ini` (the oblique analytic wavy wall,
> 32×32×8, Re 100, `forcing_x = 1`, conjugate κ_s = 5, C_s = 2) with
> `initial = solid_init = 1.0`, 1000 steps of `dt = 2e-4` from rest, CPU.
> `max|θ − 1|` over the 256 FLUID CUT cells (centre in the fluid, ≥ 1 of 6
> neighbour centres in the solid, analytic marker): **divergence 4.09e-5
> (rms 1.35e-5), skew 2.05e-5 (exactly half), advective 0.0** — and 0.0 over
> every fluid and solid cell. The divergence-form drift grows with the
> spin-up (2.0e-5 / 2.9e-5 / 3.5e-5 / 4.1e-5 at steps 250 … 1000 while
> max|u| = 0.053 / 0.106 / 0.157 / 0.204), i.e. ∝ the cut-face velocity as
> the mechanism says. The review's "O(1e-2)" was sized for a developed DNS
> thermal layer at h⁺ ≈ 1–2; the qualitative prediction (drift in divergence,
> round-off in advective, half in skew) is what was tested and it holds.
> (b) `scalar.f90`: `convSet` records whether `[scalar] convection` was
> written; `validate_conjugate_config` resolves an UNSET key to advective
> when any scalar is conjugate (printed at init), refuses an explicit
> `divergence` with a message naming F2, and warns on explicit `skew`.
> Non-conjugate scalars keep the divergence default (early return).
> Verified: the key-less wavy run prints the resolution and is
> `max_abs 0` (un/vn/wn/pn/theta) against the explicit-advective run;
> explicit divergence exits 1 with the F2 message. No ini in the tree sets
> the key on a conjugate case (every CHT tutorial and conjugate gate now
> runs advective).
> (c) Gates vs step 0 (`~/numrev_ref_binaries`, nofma): 9-case scalar suite
> **9/9 max_abs 0 CPU and 9/9 GPU**; 7-case suite **7/7 max_abs 0 CPU**
> (`turbles`/`turbslab`/`detles`/`les_ibm_refine` needed their generated
> inputs linked from `mobydiff.scalar`/`mobydiff.bl`; not committed).
> C1: exactly ONE gate moved, as predicted — the insulated-box `Σ C θ dV`
> drift **1.22e-16 → 2.235e-11 relative** (1.286e-12 absolute, 100 steps),
> the advective cut-cell leak; re-baselined as a leak-magnitude band
> (`--tolerance 1e-9`, README (3) rewritten, `run_gates_c1.sh` comment).
> All other C1 groups PASS incl. determinism 1 == 4 ranks / CPU == GPU at
> tolerance 0. C2 **all PASS** (75 PASS lines, 0 FAIL, 42 min), C3 **all
> PASS** (15/0, 13 min) — their gates are diffusion-only fixed points or
> flow-on budgets that already carried the leak, so nothing else moved.
> (d) Pipe Nusselt comparison: leg D of the 2026-09-20 campaign re-launched
> on istmcorax (RTX 5090) at 17:10 from the SAME settled state
> (`p_pr_settle2_30000.h5`, found with the whole campaign in
> `mobydiff.scalar/validation/conjugate/`) and the SAME ini (168 000 steps,
> `niter = 6` Chebyshev — kept so the binary is the ONLY change; the niter-12
> re-measurement belongs to the one-by-one pass), new binary
> `build_gpu_corax`, run dir `~/pipe_rerun_f2/` (`p_f2_stat.*`); log confirms
> "resolved to advective". **DONE 2026-09-28** (ended 17:28, 168 000 steps at
> 0.520 s/step, no NaN). Post-processed with the tutorial's own recipe
> (`pipe_stats.py` on the 15 snapshots and the plane statistics,
> `make_caches.py`, `plot_pipe.py`, `pipe_numbers.py`, `variance_budget.py`),
> old → new on the production grid: interface heat c0 9.760 → 9.761 (every
> scalar within 0.04 %), q_w 0.24854 → 0.24856, ⟨θ⟩⁺ axis +0.1 → +0.7 %, θ′⁺
> axis −6.2 → −6.0 %, θ′⁺ (y⁺ < 50) +0.7 → −0.1 %, interface signature at y⁺
> 7.25 moved ≤ 0.004 toward the reference (c4 0.903 → 0.907, mbc 0.862 →
> 0.866, isof 1.209 → 1.207); the ONE systematic move is θ′ through the
> solid: 0.798/0.638/0.528/0.420/0.385/0.374 → 0.802/0.633/0.521/0.404/
> 0.365/0.353 vs Neuhauser 0.757/0.591/0.490/0.388/0.353/0.343 — the deep
> excess +9 → +3 %, and renormalised at 0.1 d the decay reads −2.0 … +0.4 %
> where it read +0.1 … +7.7 %. Direction as F2 predicts (the drift was a
> spurious source at the fluid cut faces, where the solid is forced), but
> ONE realisation: the velocity statistics F2 cannot touch moved by
> comparable amounts (u_z′ core +3.1 → −3.4 %, ⟨u_r′u_z′⟩ −1.2 → +1.0 %,
> u_b 1.0075 → 1.0092), which bounds the claim. `tutorials/cht/pipe/asset/`
> accumulators + figures 1–7 refreshed, README production tables updated
> with a provenance note, `pipe_report.html` keeps the campaign and gains a
> §11 addendum with the old/new table.

**Step 3 — F4, Neumann normal velocity.** (a) `tutorials/sailplane/input.ini`:
replace the `x_max_{u,v,w}_type = neumann` rows by `x_max_patch = outlet`
(Dirichlet p at the outlet replaces the all-Neumann pressure; README updated).
(b) `resolve_face_bcs`/`validate`: a `BC_NEUMANN` on the normal component of
a non-periodic face is a config error pointing at the outlet patch. Gates:
the sailplane 1-step legacy bit-exact leg of `run_gates_big.sh` is
re-baselined (the BC changed on purpose); sailplane stable 200 steps with the
DEFAULT Jacobi/Chebyshev solver (it needed red-black before, which is itself
a data point); `validation/freestream` unchanged. Half a session.

> **MEASURED (2026-09-27). DONE.** (a) `tutorials/sailplane/input.ini`: the
> five `x_max_{u,v,w}_type = neumann` / `x_max_p_*` rows replaced by
> `x_max_patch = outlet` (README section added). (b) `boundary.f90
> resolve_face_bcs`: after the patch resolution, a `BC_NEUMANN` on the
> NORMAL component of any non-periodic face is an `error stop` whose message
> names the face and points at the outlet patch; verified on a pois_io
> variant without the outlet declaration (exit 1, the new message). On a
> declared outlet an explicit `x_max_u_type = neumann` still trips the
> earlier "contradicts the declared patch type" error first, which is what
> the freestream config gate greps for. The only in-tree user of a Neumann
> normal velocity was the sailplane. Gates: sailplane prepared case file
> (nb 10, 2 ranks, 72 s) vs the committed legacy coefficient file, 1 step,
> production CPU, SAME outlet ini both sides: **un/vn/wn/pn max_abs 0**
> (re-baselined in the sense that both sides now carry the outlet).
> Sailplane 200 steps on the DEFAULT projection, GPU (RTX 3060, 18 M cells):
> plain damped Jacobi `sor 0.8, niter 20` **1.72 s/step, bounded**
> (max|u| 1.33 / 2.00 / 2.71 / 1.98 at steps 50/100/150/200, max|p| 400–600,
> no NaN) and Chebyshev `niter 12` **1.50 s/step, bounded** (max|u| 2.00 /
> 2.00 / 1.00 / 1.71, max|p| 550–1120) — so the case no longer NEEDS
> red-black; the earlier blow-up on the default solver was the `sor = 1.5`
> Jacobi confound (CLAUDE.md, 2026-09-26), not the BC. `solver = redblack`
> is kept in the ini as the faster choice. `validation/freestream/run_gates.sh
> all` at production CPU: oblique exact, pois PASS, 1 == 4 ranks EXACT (both),
> declared == inferred wall EXACT, both contradiction rows error-stop —
> status 0, unchanged.

**Step 4 — F5, `apply_bc` out of the projection loop.** Depends on step 3.
`pressure_projection`: one `apply_bc` after the last `jacobi_apply` /
`interface_correct`, before the final full exchange; same for the red-black
loop (once, after the last colour). Gate: bit-exact at production flags on
the suite, `validation/freestream` (outlets) and `validation/redblack_interface`
1 and 4 ranks — every write the removed calls made was idempotent, so
`max_abs 0` is the expected result, and anything else means a face kind was
missed. Measure `proj_timing: bc` before/after (expect ≈ −(niter−1)/niter of
the bucket). Half a session.

> **MEASURED (2026-09-27). DONE.** `pressure_projection`: the per-iteration
> `apply_bc` is gone; ONE call at `iIter == nIter` after `jacobi_apply` /
> `interface_correct` and before the final full exchange (comment states why
> nothing between iterations reads the ghosts, and that F4 removed the one
> non-idempotent write). `redblack_projection`: once, after the last colour of
> the last iteration, same position. Gates at PRODUCTION flags against the
> step-0 production binaries (`~/numrev_ref_binaries/moby_solve_{cpu,gpu}`),
> all **max_abs 0**: 7-case suite CPU 1 rank 7/7, CPU 4 ranks 7/7, GPU 7/7
> (incl. every RANS scalar); `validation/freestream` pois_io (outlets, 200
> steps) / oblique / lamboseen pre-vs-new 4/4 datasets each, and
> `run_gates.sh all` status 0 unchanged; `validation/redblack_interface`
> refined_channel pre-vs-new 1 rank and 4 ranks 4/4 datasets, and its own
> driver ALL PASS. Profile (`[output] profile = true`, GPU RTX 3060, 200
> steps, niter 12), `proj_timing: apply_bc` before → after:
>
> | case | calls | bc s/step | Δ | projection s/step | step s/step |
> |---|---|---|---|---|---|
> | min_channel (2:1, Chebyshev) | 7200 → 600 | 9.59e-4 → 8.07e-5 | **−91.6 %** | 36.7 → 35.7 ms (−2.8 %) | 46.0 → 41.6 ms |
> | freestream pois_io (outlets) | 7200 → 600 | 6.38e-4 → 5.38e-5 | **−91.6 %** | 4.30 → 3.71 ms (−13.7 %) | 4.55 → 3.97 ms |
>
> −(niter−1)/niter = −91.7 % predicted; the bucket landed on it. The step
> deltas beyond the bucket are single-run scatter on a shared workstation
> GPU, not a claim.

**Step 5 — F8 and F6, documentation-level.** `cflmax` documented as a
per-direction number in `configuration.md`, plus an init-time print of the
worst-case sum over directions (cheap, from the same reduction). KMM180: read
the archived run's `cfl` prints in the sibling checkout to see what it
actually stepped at; if the limiter cut it to ~8e-6, set `natural_dyw_plus`
to a wall-resolved value (0.5) and the caps accordingly, and record the
decision in the tutorial. An hour each. Step 5b decides whether the momentum
half of step 11 has a customer.

> **MEASURED (2026-09-27). DONE.** (a) F8: `docs/configuration.md` `cflmax`
> row now says it bounds `dt` by the largest SINGLE-COMPONENT Courant rate
> and that the RK3 limit is on the SUM (√3); `get_timestep_rates` gained an
> optional `sum_rate` (max over cells of Σ_d |u_d|/Δx_d, a second reduction
> variable in the same kernel — `rates` untouched) and
> `update_timestep_limits` prints it ONCE, on the initial field:
> `cfl: max single-component Courant rate X; worst-case SUM Y = r × the
> component max; cflmax × ratio = …`. Inert: 7-case suite at production
> flags vs the step-0 binaries, CPU 1 rank 7/7 AND GPU 7/7 max_abs 0. On min_channel and KMM180 the initial field is streamwise
> (ratio 1.00), so the print is informational there; it exists for the 3D
> cases. (b) F6: the archived KMM180 restart's metadata
> (`mobydiff.scalar/tutorials/channel_kmm180/channel_kmm180_restart.h5`,
> 2026-08-03) reads `dt = 7.696034892725926e-06`, `cfl = [0.0037, 0.5]`,
> `step 650000` at `t = 5.0` — **the Péclet limiter had cut the step to
> 7.7e-6, as Appendix A computed (7.7e-6), for 650 k steps per 5 h/u_τ.**
> `tutorials/channel_kmm180/input.ini` now sets `natural_dyw_plus = 0.5`
> (comments + a new README record the decision). 100-step GPU runs of the
> ini, both spacings: 0.5 → `dt = 3.125e-4` (dtmax binds; cfl 0.1306,
> Péclet 0.2268, dy_wall⁺ 0.498); 0.05 → `dt = 7.696035e-06` (Péclet 0.5
> binding, cfl 0.0032) = the archived value to every digit. 40.6× fewer
> steps. **Step 5b's verdict for step 11: the production KMM case is NOT
> viscous-limited any more at the wall-resolved spacing (Péclet 0.23 vs
> cfl 0.13, both under dtmax); the momentum half of step 11 has no customer
> here** — a customer would be a case that must keep dy_wall⁺ ≪ 0.5.

**Step 6 — Q4, one BC mechanism (refactor).** Depends on step 4. (a) Add a
per-point constant to the gather (`lC`/`rC`, zero for every existing entry).
(b) `init_block_exchange` emits physical-face ghost entries from
`faceBcType`/`pointBcValue`: source = the adjacent interior cell (or its
same-rank halo copy — ordering: BC entries run after the same-level copies,
so edge/corner ghosts read a refreshed halo), weight −1 (Dirichlet mirror) or
+1 (Neumann/copy), constant `2v` or `dn·v`; scalars in `q` ride the same
entries; the projection's outlet phi mirror becomes the same entry with
weight −1 and constant 0. (c) Delete `apply_bc`'s ghost branch,
`apply_scalar_bc_q` and `apply_scalar_bc`; keep one tiny kernel that writes
pinned normal faces at init/restart and the outflow copy per substage.
Gate: bit-exact at production flags — `2v − q` and `−q + 2v` are the same
IEEE operation, and the corner ghosts are the same expression on the same
values; `min_channel` 1 == 2 == 3 == 4 ranks EXACT is the load-bearing leg.
Two to three sessions.

> **MEASURED (2026-09-27). DONE, with a DELIBERATE REDUCTION of scope: the
> boundary rows are resolved into ONE affine write and one kernel form, but
> they stay a boundary.f90 kernel and did NOT become exchange entries.** Why
> the entry version was not taken, found while reading comm.f90 before
> writing any code: (1) the same-level COPY entries READ physical ghosts --
> `entry_boxes` extends a face entry into the ghost row (index 0 / nb+1)
> wherever the combined edge neighbour is absent, i.e. at every physical
> wall, and `pack_entries` runs FIRST (to overlap the messages), so the BC
> writes must complete in a launch of their own before pack whatever list
> they sit in; "zero kernels" is unreachable and the launch count is
> unchanged either way. (2) Adding a per-point constant to the existing
> gather changes the arithmetic path of every halo point (`x + 0.0` flips
> the sign of a negative zero), so the existing entries would need a branch
> anyway. (3) The 2:1 cross-level entries do NOT extend into physical ghost
> rows (`interface_boxes`), so the ordering hazard of section 3 is confined
> to same-level copies, is stated at `apply_bc`, and is what the 1 == 4
> ranks legs gate. (UPDATE 2026-10-01: they extend now, as same-level
> entries do -- not extending left the normal velocity of a physical HIGH
> face unwritten next to a level jump; `comm.f90 candidate_boxes`,
> `docs/next_session_outlet.md` O4. The ordering rule of `apply_bc` then
> covers every transfer. `FACE_OUTFLOW` below is gone the same day: the
> outlet face is predicted.) What WAS done (boundary.f90, scalar.f90): the BC type of
> every q variable is resolved at init into `bcKind(var,face)` (GHOST /
> FACE / FACE_OUTFLOW / NONE), `bcW(var,face)` and `bcC(var,point)`
> (`resolve_affine_rows`; the scalar columns come from `init_scalar` via
> `set_scalar_bc_rows`, the tables are sized `dns%nVar`), and every write
> is `dst = w*src + C`: Dirichlet ghost w=-1 C=2v, Neumann ghost w=+1 C=dn v,
> pinned face w=0 C=v, outflow face w=+1 C=dn v (=0, predictor call only).
> `apply_bc(blk, bc, vars, outflow_copy)` is generic over a variable list
> (`scalar_sync` passes the scalar columns), `apply_scalar_bc_q` is DELETED,
> `apply_scalar_bc` (standalone RANS/phi arrays) is the same line with its
> per-face mode resolved per call; neither kernel maps `blk%x/y/z` or the
> type table any more and neither has a type branch (the three kernels were
> 110 + 62 + 69 = 241 lines, the two are 76 + 48 = 124; boundary.f90 as a
> whole 891 -> 917 with the resolver and its comments). Bit-exactness
> argument: `2v - q` and `(-1)*q + 2v` are the same IEEE operation, `0*q + v`
> is `v` for every finite q, and the Neumann `dn*v` moved out of the kernel
> (it used to be contractible into `fma(dn, v, q)`), so the formal gate is
> nofma and the only rows that can move at PRODUCTION flags are nonzero
> Neumann data -- present only in `tutorials/cht/channel` (heat-flux walls),
> in no suite. Gates vs `~/numrev_ref_binaries`, all **max_abs 0**: 7-case
> suite nofma CPU 1 rank 7/7, CPU 4 ranks 7/7, GPU 7/7; 9-case scalar suite
> CPU 9/9, GPU 9/9 (uniform3 = scalar Dirichlet inlets + Neumann outlet);
> extra legs -- freestream pois_io (parabola inlet, outlet, walls, 200
> steps), oblique, lamboseen, `redblack_interface/refined_channel`, and
> `rans_inlet/inlet_channel` (the SCALAR_BC_VALUE mode, k/omega inlet
> ghosts) -- CPU 1 rank 5/5, CPU 4 ranks 5/5, GPU 5/5; and the 7-case suite
> at PRODUCTION flags CPU 7/7, as predicted. **61 comparisons, every one
> max_abs 0** (incl. every RANS scalar and nut). MEASUREMENT LANDMINE: the
> suite drivers key their output prefix on `MODE`, so two suites with the
> same MODE running concurrently in the same directories delete each
> other's snapshots (one leg read NO OUTPUT and was rerun under its own
> label); the drivers now say so.

**Step 7 — Q1, prepare does everything (refactor).** Independent of step 6.
(a) `moby_prepare` accepts body-free cases and writes node lines + leaf
table + face kinds always. (b) The solver reads the leaf table and stops
rebuilding it; the row-by-row cross-check stays behind a flag for one
release. (c) The solver requires a case file; `moby_solve --prepare` (or
auto-prepare when the file is absent or its stored ini hash is stale) keeps
the tutorials one command. (d) Delete the inline classify dispatch
(`moby_solve.f90:94-120`), the legacy global-layout coefficient reader, the
rank-box layout (`nb` becomes mandatory). (e) Later: a per-leaf weight column
and the weighted Morton split. Gates: P0/P1/P1b prepare gates; the 7-case
suite run FROM prepared files bit-exact vs the inline reference of step 0
(the wavy/Beltrami P0 gates already prove this for the analytic path); the
GPU `set_ibm_coeff` kernel is retired with the CPU prepare canonical. Two to
three sessions.

> **HANDOUT (2026-09-29): `docs/next_session_prepare_everything.md`** —
> increments 7-0 … 7-4 with their gates, the `nb` rule (settled 2026-09-29:
> the key stays optional — unset means one block per rank, chosen by the
> builder from its rank count and STORED in the case file, i.e. today's
> rank box as a case property; only the rank-box LAYOUT code goes), the 23
> live inis without `nb` (unchanged under the rule) and the 9
> inline-analytic-IBM cases that become prepare+solve pairs, and the
> inherited landmines. Not started.

> **MEASURED (2026-09-29, increments 7-0 … 7-2; production flags, ref =
> `~/step7_ref_binaries` at `d2249b1`).** 7-0: eight binaries + PROVENANCE.
> 7-1: `moby_prepare` accepts body-free and nb-less cases; the nb rule is
> `derive_block_nb` (blocks.f90), the derived size rides the file as
> `block_nb_auto`/`block_nb_ranks`. 7-2: the solver takes block size, leaf
> table and node lines from the case file (`read_case_layout`,
> `init_block_set_from_table`) and rebuilds nothing; the coefficient/dwall
> readers lost the row-by-row self-comparison; a restart-snapshot row
> cross-check was ADDED (`fdm_h5_check_block_table` — the one the handout
> said "stays" only ever checked dataset extents). Gates
> (`validation/prepare/run_gates_step7.sh` + the existing drivers), all
> `max_abs 0`: min_channel and beltrami/slab_y (explicit nb) prepared on
> 1 == 4 ranks IDENTICAL files, prepared(p)+solved(r) == inline ref(r) for
> p, r ∈ {1, 4} (8 legs); wf180_y30 (8×6×8, nb-less) prepares to nb 8 6 8 on
> 1 rank and 4 3 8 on 4 (dims 2 2 1 — an ODD nb, accepted: the red-black
> `colorOffset` is `sum(origin) mod 2`, continuous across any block face,
> so parity only matters for refinement, which still needs an explicit nb),
> DIFFERENT tables and SAME fields on 1 and 4 ranks, and the 1-rank file on
> 4 ranks stops with "rank owns no blocks" as rule 5 says; conduction
> (4×16×4, scalar `s1`) the same at nb 4 16 4 / 2 8 4; wavy (analytic body,
> pinned dims) 1 == 4 identical + solve == inline; GPU twin of all five
> from the CPU-prepared files vs the GPU reference: 16/16. P0 22/22, P1 STL
> 16/16, block_nb 7/7 (CPU). 7-case suite 1 and 4 ranks, 9-case suite 1 rank
> (its `turbles`/`turbslab` pin `[mpi] dims`, so a 4-rank run of that suite
> is not a valid invocation), GPU 7- and 9-case: every leg `max_abs 0` except
> `les_ibm`, whose committed `ibm_coeff.h5` is a LEGACY global-layout file —
> the handout named the sailplane as the last such user; this is a second —
> gated after the layout reader learned to report a legacy file as
> layout-less (below). FOUND: `[mpi] dims` in a prepare input describes
> the SOLVE (the P0 wavy inis pin `1 1 1`), so `comm_cart_dims` takes a
> `for_solve` flag and prepare uses `product(dims)` as the rank count;
> `block_levels` means "number of levels" in mobygeom files and "finest
> level" in prepare files — the reader derives the count from the table.
>
> **MEASURED (2026-09-29, increment 7-3).** One builder, two entry points:
> `src/modules/prepare.f90` `prepare_case(input, case_file, c)`;
> `moby_prepare` is a thin main and `moby_solve` calls it IN-PROCESS when
> the case file is absent (`--prepare` forces it). `[case] file` is the one
> user-facing name (default `<field_prefix>.case.h5`); `[ibm] coeff_file`
> is an alias for one release (a note is printed; both set = config
> error). STALENESS is an INPUT ECHO, not a hash: `case_input_echo`
> (config.f90) prints every value the file is a function of, one
> `key = value` per line (grid, lengths, re, periodicity, nb or `auto`,
> refine keys and boxes, remove_solid/keep_buried, ibm_enabled, STL list +
> transform or the analytic wall parameters, the scalar column); the
> builder stores it (`case_inputs` attr, readable with h5dump), the solver
> rebuilds it from its ini and the first differing line stops the run with
> both values printed. Files without the attr (mobygeom, pre-7-3) fall back
> to the explicit checks. The committed legacy `les_ibm/ibm_coeff.h5` was
> converted ONCE, number for number, with the temporary `moby_prepare
> --convert-legacy` into `ibm_coeff_case.h5` (committed, 15 MB): a solve
> from it equals the 7-0 reference's solve from the legacy file at
> `max_abs 0` (un/vn/wn/pn/nut); `channel_ibm.ini` reads it, the legacy
> file stays for `measure_nut.py`. Manual gates on min_channel: a run with
> no case file prepares `<prefix>.case.h5` and reproduces the inline
> reference at `max_abs 0`; a second run reads it; `--prepare` rebuilds it;
> `re` 180 → 190 in the ini stops with `file has: re = 1.8…E+002 / ini
> has: re = 1.9…E+002`; `[case] file` gives the same fields; the alias
> prints its note; both keys error-stop. Suites, every run now
> auto-preparing, all `max_abs 0` vs 7-0: step-7 driver CPU 23/23 + P0
> 22/22 + STL 16/16 + block_nb 7/7 (65), 7-case 1 + 4 ranks (16), 9-case
> (10), GPU step-7 driver + 7-case + 9-case (33).
>
> **MEASURED (2026-09-29, increment 7-4 — the deletions; −1097 / +279
> lines over 17 files).** GONE: the inline classify dispatch in
> `moby_solve.f90` and the file branches of `classify_active_mask` /
> `classify_refinement_masks` (+ the mask readers `read_block_active` /
> `read_block_masks` / `read_mask_window` and their C functions); the legacy
> global-layout coefficient reader (`fdm_h5_read_ibm_coeff`,
> `read_ibm_coeff_legacy`, the 7-3 converter); `DIST_RANKBOX` / `distMode`
> and every branch keyed on it (`init_block_set` requires nb, `face_kind`,
> `build_block_metadata/metrics`, `resolve_neighbors`, `owner_of_origin`,
> the rank-box peer walk, `rank_box` / `local_range`, `dns%localSize`,
> `set_serial_local_size`); the DEVICE `set_ibm_coeff` and its device
> helpers (`bisection`, `add_neighbor_coeff`, the `declare target` lines) —
> the host twin is now THE kernel, named `set_ibm_coeff`, over the
> indicator; the solver-side analytic dwall branches in rans.f90 /
> scalar.f90 (the file always carries dwall). Every run reads a case file;
> a file without a blocks table (legacy layout) is refused naming the fix.
> The sailplane tutorial declares its STL + `[blocks] nb = 10` + `[case]
> file` and prepares on first run; its legacy `sailplane_ibm_coeff.h5` is
> retired (git history, last at the 7-3 commit); `run_gates_big.sh`'s solve
> leg is now prepared-file 1 == 2 ranks. GATES (production flags vs 7-0,
> UNCHANGED inis): step-7 driver CPU 10 cases (the five of 7-1/7-2 + the
> live analytic-IBM cases `validation/conjugate/wavy`, scalar `ibmwavy` /
> `ibmwavyr` (refine_body), rans_geometry `wavy` / `wavy_refine`) + P0
> 22/22 + STL 16/16 + block_nb 7/7 = 86/86; 7-case suite 1 + 4 ranks 16/16;
> 9-case 10/10; freestream (oblique / Poiseuille in-out / Lamb–Oseen /
> 1 == 4 ranks / CPU vs GPU / config twins), red-black 2:1 (1 == 4 ranks
> EXACT), RANS inlet (1 == 4 EXACT), min_channel 1 == 2 == 3 == 4 ranks
> EXACT — all `max_abs 0`. GPU: body-free cases vs the GPU reference and
> the analytic cases vs the CPU reference 47 legs `max_abs 0`, plus
> 7-case + 9-case + red-black GPU; the four "failures" — `ibmwavy` /
> `ibmwavyr` GPU vs CPU-ref — are `theta` alone at **1.1e-16** and are
> NOT this step: the new GPU binary solving from the prepared file equals
> the OLD GPU binary's inline solve at `max_abs 0` on every field
> (velocities, pn, theta, phi — the analytic coefficients even came out
> identical between the device and host kernels on these cases), and the
> old GPU binary shows the same 1.1e-16 against the CPU reference: it is
> the scalar transport's CPU-vs-GPU production-flag arithmetic, which
> predates step 7. So the "one expected move" of the handout (device libm
> ulps) did not materialise on any gated case, and there is no measured
> GPU deviation attributable to 7-4.
> FOUND BY THE FREESTREAM DRIVER: the 7-3 default name `<prefix>_case.h5`
> matched every driver's snapshot glob `<prefix>_*.h5` (the Lamb–Oseen
> checker read the case file as a snapshot); the default is now
> `<field_prefix>.case.h5`, re-gated: freestream 9/9 (reflected fraction
> 2.230e-2, the niter-12 value), 7-case CPU 7/7 + GPU 7/7, red-black GPU,
> min_channel step-7 legs 5/5.
>
> **REVERSED (2026-09-30, at the user's request): no in-process prepare.**
> "Keep prepare / solve separate without repetitions." The builder was one
> routine with two callers, so no code was duplicated, but the solver
> re-read the ini and re-built the grid inside `prepare_case`, applied the
> restart metadata twice, carried the whole preprocessing in its binary and
> answered "which rank layout" in two places. Now: `moby_solve` never
> builds a case file — a missing one stops with the exact `mpirun -n N
> moby_prepare input.ini` printed, a stale one with the key named;
> `moby_prepare input.ini` defaults its output to the solver's name through
> the ONE resolver `case_file_name` (config.f90); `--prepare` is gone. The
> gate drivers (25 scripts, ~60 launch sites) call
> `tools/prepare_if_missing.sh <ranks> <solver> <ini>` before a solve: it
> resolves the file the solver will read and prepares it ONLY when absent —
> a blind prepare would have recomputed the committed `ibm_coeff_case.h5`
> and the zero-force twins the multilevel_body/refine2d inis point at — with
> the same rank count as the solve and the moby_prepare next to the solver
> (none beside a pre-step-7 reference: it builds inline). The tutorials
> gained the explicit prepare line. FOUND: the P1 convention (prepare ini
> names the STL, the solve ini drops it) made the input echo report the STL
> line as stale on every conjugate/CHT/NACA case; `compare_case_inputs` is
> now key-based with the GEOMETRY RULE — an ini naming no geometry source
> accepts the file's stl_*/wall_* lines, one naming STL files must match
> them — and `h5same.py --ignore-attr case_inputs` lets the S3 "coef_p_blocks
> is the only change" gate ignore the echo, which legitimately records
> `scalar_coef`. Re-gated with the drivers: step-7 driver + P0 + STL +
> block_nb 92/92, 7-case CPU 1 + 4 ranks 16/16, 9-case 10/10, freestream /
> red-black / RANS inlet / min_channel 1 == 2 == 3 == 4 (12/12), conjugate C1
> + scalar S3 88 legs (S3's cylinder leg SKIPPED: its campaign restart is
> not on this machine, as before), GPU 47 + the four known 1.1e-16 `theta`
> legs. S4 is not runnable here either (needs `turbles_67600.h5` from the S2
> LES campaign) and was not re-gated.

**Step 8 — F1, flux-form viscous stencil (numerics change).** The plan entry
in section 9. Gates as stated there; run it in the same session as the
stretched-case re-validation so the RANS channel and Blasius numbers are
re-measured once, not twice. One session plus remote GPU reruns.
**MEASURED (DONE 2026-09-30, `docs/next_session_after_step7.md` STEP 8).**
`slice_grid_direction` builds the flux form on the variable's own control
volume: `lapM = 1/(hm·width)`, `lapP = 1/(hp·width)` in a cell-centred
direction (`width` = the cell, computed once beside `d1`), the Taylor
spelling `2/(hm(hm+hp))` kept textually in the face-staggered direction
where the two coincide. New unit test `src/test_stencil.f90`
(`build_cpu/stencil_test`): for uniform / geometric r = 1.10 / cosine /
tanh lines, all three directions, all four variables, the telescoping
identity holds to ≤ 2.7e-15 and the volume-weighted symmetry to ≤ 6.5e-14
(48/48); the informational column puts the OLD weights 2.4 % / 33 % / 4.2 %
off the flux weights on those lines and 0 on uniform and staggered ones.
Suites vs `~/step7b_ref_binaries` at nofma, CPU AND GPU digit for digit:
uniform-line cases move at ROUND-OFF only (les_ibm ± refine u ≤ 8.2e-14,
pn ≤ 1.5e-12; Beltrami 1.0e-14 / 4.9e-14; det/conserve 1.1e-14 / 2.9e-13;
lam30t 1e-37; uniform3, conduction, prsweep, wf180_y30 EXACTLY 0) — the
cell width and the centre-to-centre spacings round differently, so the
"bit-exact on uniform lines" half is ulps, not zero, as the plan allowed;
stretched cases move at truncation level over 20 steps (min_channel u
7.7e-4, turb180 1.2e-3, turbsst 9.3e-4, the turbulent LES channels
turbles/turbslab/detles 2.2-2.4e-2 on a chaotic field). Converged RANS
turb180 (natural line, GPU, 132 565 steps, ONE binary lineage, baseline
run the same morning at niter 12): T2 gate PASS both before and after —
u_tau 1.0008 → **1.0000**, log-law max dev 0.049 → 0.050 (tol 0.06), U+
centreline 18.16 both (DNS 18.20); the converged field moved 7.8e-3 in u
(0.05 % of the velocity scale), 1.7e-1 in omega (of 8.5e4). Blasius
precursor (laminar, one-sided natural y, recovered from `4f819fc^`; GPU,
t = 2000 from the analytic IC): gate numbers unchanged to the printed digit
between the pre-F1 and F1 binaries (theta 1.33 %, H 0.45-0.46 % at
niter 12), fields 3.4e-5 apart in u. KMM180 developed statistics (natural
y, dyw+ 0.5; 20 t.u. windows after a 30 t.u. transient from the archived
restart, pre-F1 and F1 in parallel on istmcetus): every profile within the
sampling scatter -- U+ centreline 18.515 / 18.510, u'+ peak 2.652 / 2.623,
-<u'v'>+ 0.7192 / 0.7193, u_tau 1.0000 both, both ~1 % from KMM DNS
(`tutorials/channel_kmm180/asset/f1_2026-09-30/`). Nothing owed.

**Step 9 — F7, exact penalization factor (numerics change at cut cells).**
`update_ibm_mu` produces two factors, `muA = B/(λΔt + B)` on `q` and
`muB = 1/(λΔt + B)` on the substage increment, `B = λΔt/expm1(λΔt)` with
`B = 0` for `λΔt > 700`; the predictor applies them; the projection keeps
`muB`. Body-free: `λ = 0` gives `B = 1`, both factors exactly 1.0, so
bit-exact by construction. Body cases: cylinder Re 40 `C_D` and the `les_ibm`
law of the wall re-measured; a `dt` halving study on the cylinder drag should
show the cut-cell time error drop faster than before. One session.

**MEASURED (DONE 2026-10-01, `docs/next_session_after_step8.md` item 3).**
Implemented in the equivalent closed form `x = λΔt`: `ibm%mu` =
`(1 − e^{−x})/x` (= `muB`, series below x = 0.1) and the state factor
`e^{−x}` (= `muA`), `ibm.f90 penal_incr_factor` / `penal_state_factor`.
DEVIATION from the text above, for bit-exactness and cost: the fused
predictor is NOT touched (`qs = (q + Δ)·mu`); the missing `(e^{−x} − mu)·q`
is a separate kernel over the BODY BLOCKS only
(`step.f90 add_penalization_state_correction`), and `update_ibm_mu` visits
the same list — so a body-free rank runs neither, by construction, and a
body case pays an `exp` only where there is a body.
Gates: unit test `penalization_test` (factors vs `expm1` references ≤ 6e-16,
exact at x = 0 and x = 1e28); nofma suites vs `~/step8_ref_binaries`, CPU
and GPU — the 5 body-free flow cases and all 9 scalar cases `max_abs 0`
(28/28), `les_ibm` ± refine move (1e-3 in u after 20 steps), as intended;
with a body: 1 == 4 ranks `max_abs 0`, CPU == GPU `max_abs 0`, solid-cell
velocity 3e-28. THE ORDER, measured where it can be
(`validation/penalization/`, `du/dt = f − λu` in a uniform periodic box):
implicit Euler 1.23e-2 → 1.45e-3 over three halvings (order 1.04, 1.03,
1.01), exact factor ≤ 5.6e-17 at every dt. Body cases: cylinder Re 40
`C_D` 1.69234524 → 1.69234525 (same steady fixed point); `les_ibm`
developed statistics t = 5..25 paired old/new — U+(29) 13.08 / 13.10, bulk U
15.065 / 15.066, differences of opposite sign at the two walls (sampling
scatter). **The predicted "cut-cell time error drops faster" on the cylinder
is NOT observable: it was never visible.** On the Re 100 shedding flow the
two factors differ by ≤ 1.1e-5 in velocity at every dt from 5e-3 to
6.25e-4, three orders below the dt-to-dt differences near the body — a body
at rest keeps its cut cells quasi-steady. And the dt study is not a clean
convergence test anyway: the stored pressure pollutes within one time unit
at the Dirichlet-p outlet (Chebyshev niter 60; `validation/cylinder/README.md`),
identically for both binaries. NOT done, same first-order factor: the
Dirichlet SCALAR penalization (`scalar.f90` `mus = 1/(1 + dt_γ coef_p/Pr)`).

**REVISED THE SAME DAY (2026-10-01, `docs/next_session_after_step9.md`
item 1): the exponential is replaced by the AMPHIBIOUS rational form, for
momentum AND the Dirichlet scalar.** With `P3 = 1 + x + x²/2 + x³/6`:
`state = 1/P3`, `incr = (1 + x/2 + x²/6)/P3` (`ibm%mu`), and
`state − incr = −x(1/2 + x/6)/P3` from its own function
(`penal_state_minus_incr`), so the correction kernel no longer reads
`ibm%mu`. The structure above is unchanged (untouched predictor + the
body-blocks pass). The scalar: `ss = (s0 + Δ)·incr + (state − incr)·s0 +
(1 − state)·s_body` at `x = dt_γ coef_p/Pr`, in cells with `coef_p ≠ 0`
only. A one-time `min(coef) ≥ 0` check at the case-file read replaces
AMPHIBIOUS's per-cell clamp (P3 has a root at −1.596).
Why: with a time-varying right-hand side the RK3 coupling sets the error and
the two factors are indistinguishable (handout table), the exponential cost
an `exp` per DOF and substage, and a rational function calls no libm on the
device. MEASURED:
- `penalization_test`: the three functions vs exact rational references
  ≤ 2e-16, `state + x·incr = 1` ≤ 1.1e-16, exactly (1, 1, 0) at x = 0,
  `1 − state` exactly 1 in a solid cell, positive and monotone over
  1e-8 … 1e30, `(state − e^{−x})/(x⁴/24) → 1`.
- `validation/penalization/` is now a THIRD-ORDER gate with a scalar twin:
  velocity 1.01e-4 → 2.61e-7 over three halvings (orders 2.77, 2.88, 2.94 —
  the pre-registered row), scalar 7.46e-5 → 1.92e-7 (the same orders; implicit
  Euler 9.06e-3 → 1.07e-3, order 1); `λ = 200` reaches the steady fixed point
  to ≤ 2.9e-15; CPU and GPU give the same 16 digits.
- nofma suites vs `~/step9_ref_binaries`, CPU and GPU: the 5 body-free flow
  cases and the 9 scalar cases `max_abs 0` (28/28). Body cases move at
  O(x⁴): `les_ibm` 2.5e-9 in u after 20 steps (the exponential had moved it
  1e-3 from implicit Euler), `refine_body` 5.0e-7; CPU == GPU `max_abs 0` on
  both, and on the wavy-wall body case with scalars (1 == 4 ranks too).
- Cylinder Re 40 control-volume `C_D`, exponential / rational:
  1.6923452515 / 1.6923452519 (the sampled series differ by ≤ 1.0e-8).
- Scalar body gates (`validation/scalar/README.md`, last section): solid cell
  == body value exactly, heat diagnostic 2.7e-15, the `ibmwf180` /
  `ibmwf1000` closed-form budgets 3.7e-15 / 4.4e-15; the two energy-budget
  residuals fall 3× (3.9e-4 → 1.2e-4, 4.5e-4 → 1.5e-4), as a time error
  should.
- Cost, `les_ibm` (256/640 body blocks), 400 steps, ms/step, implicit Euler /
  exponential / rational — RTX 3060: `ibm_mu` 0.54 / 1.06 / 0.40, step
  24.01 / 24.74 / 24.01; A6000: `ibm_mu` 0.257 / 0.439 / 0.212, step
  13.35 / 13.72 / 13.47. PRE-REGISTERED `ibm_mu ≤ 0.6 ms` and the step within
  1 % of the implicit-Euler binary: met (+0.0 % and +0.9 %). **It was NOT met
  by the first form of the functions**, which wrote the cubic coefficient as
  `x/6.0d0`: a second fp64 divide per DOF that the compiler may not turn into
  a multiply (0.72 ms and +2.5 % on the 3060). The 1/6 is now a constant
  multiply; one divide per evaluation. What remains over implicit Euler is the
  separate state-correction pass (`momentum` +0.17 ms on the A6000), which
  the exponential version carried too.
Reference set: `~/step9b_ref_binaries`.

**Step 10 — Q7(c), implicit solid conduction for conjugate scalars.** Does
NOT depend on step 8 (the scalar operator is already flux form). Design: a
Chebyshev-Jacobi Helmholtz solve over solid + cut cells only, coefficients
static (the conjugate face diffusivities), Gershgorin bounds computed once at
init, fluid values lagged as Dirichlet data at cut faces, Crank–Nicolson per
substage so the C3 transient order gates stay 2.00 (implicit Euler would make
the solid first order). `scalar_conjugate_peclet_rate` then drops the solid
cells' `C_s` term and keeps the fluid cut-cell share. Gates: every C1–C3
gate reproduced (steady states identical, transient orders 2.00), the pipe
`rhocp = 0.0625` case at the FLUID's `dt`. Two to three sessions.

**Step 11 — Q7(b), line-implicit wall-normal momentum diffusion.** Only if
step 5b shows a viscous-limited production case. Depends on step 8. Per-block
tridiagonal per (i,k) line in y, block-Jacobi across y block faces with one
extra exchange, approximate-factorization error O(Δt²) inside RK3. Gate:
laminar Poiseuille exact at σ = 20; developed KMM statistics at the larger
`dt` against the archived reference. Two sessions.

**Step 12 — Q8, second halo layer, then bounded scalar convection.** Lowest
priority, largest change. (a) `comm.f90` two-layer entries, gather maps, 2:1
rows and the divergence subset; nothing reads layer 2 yet, so the gate is
bit-exact. (b) TVD/van Leer for the RANS scalars (the SD7003 front is the
measured motivation), then optionally for passive scalars. (c) A 4th-order
central toggle for momentum with the interface/body/wall mask, as an
experiment only. Several sessions.

Dependency graph: 0 → {1, 2, 5, 7, 9, 10} in any order; 3 → 4 → 6;
5b → 11; 8 → 11; 12 last. Steps 1–5 together are about three sessions and
close every correctness finding; 6–7 are the simplifications; 8–12 are
numerics changes and features, each with its own re-validation.
