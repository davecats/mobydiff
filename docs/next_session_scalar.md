# Passive scalars (branch `scalar`) — implementation plan

STATUS: **S0 + S1 + S2 + S3 + S4 LANDED AND GATED (S0–S2 2026-08-03, S3 and
S4 2026-08-04, branch `scalar`).** S5 (thermal wall function, TVD scalar
convection, Boussinesq) is NOT started.

What landed (the plan below is unchanged except where noted):

- **S4** — statistics and tooling. New `src/modules/scalar_stats.f90` owning
  `scalar_stats_type`: per-row accumulators of SEVEN columns per scalar —
  `<s>`, `<s²>`, `<u_c s>` and, on the cell's LOW and HIGH y face, the
  convective flux `<v s>` and the TOTAL flux `J = <v s − D ∂s/∂y>` built with
  the TRANSPORT KERNEL's own face diffusivity (`1/(Re Pr) + ν_t/Pr_t(face)`,
  the same `eddy_diffusivity`, the same FACE_CLOSED / adiabatic-body masks).
  rms is `sqrt(<s²> − <s>²)`, the resolved turbulent flux is the convective
  column, and a WALL row's `J` is the exact discrete wall flux — which is what
  makes `theta_tau` and the Nusselt number exact rather than reconstructed.
  Plus the immersed body's HEAT RELEASE in the cancellation-free form, and
  `tools/compare_fields.py` dataset discovery.
  - **DEVIATION from §9 (deliberate): the statistics are a SOLVER-LEVEL
    facility driven by `[scalar]` keys and called from `moby_solve.f90`, not a
    case component "alongside channel_stats.f90".** The case `after_step`
    interface carries neither `sc` (Pr, Pr_t, the wall modes) nor `turb%nut`,
    and the same statistics have to serve the channel, the boundary layer AND
    body cases (the heated cylinder runs the `airfoil` case). Extending the
    case interface would have touched every case for no gain. What IS taken
    from `channel_stats`/`bl_stats` is their row machinery, verbatim: the
    per-level `lvlOff` tables and the `(x,y)` y-fastest flattening, and their
    HDF5 writers/readers — so S4 adds NO new C code and the existing readers
    work (`nstat` = 7·nScalar; stat `s` of scalar `is` is column
    `7(is−1)+s`).
  - Config, all in `[scalar]`, all off by default:
    `stats_sample_interval`, `stats_write_interval`, `stats_file`,
    `stats_layout = profile|plane`, `heat_interval`, `heat_file`.
    `profile` = wall-normal rows x-z averaged, one file per refinement level
    (`name.h5`, `name_l1.h5`, …) — the channel form; `plane` = rows of the
    global `(x,y)` plane z averaged — the boundary-layer form. Accumulators
    continue from the file on restart (the `channel_stats` recipe).
  - **The Nusselt/heat diagnostic uses the cancellation-free form the S3
    FINDING forced**: `heat_interval` writes, per scalar, the flux across
    every staircase face separating a solid cell from a fluid one PLUS the
    penalization delivered into the graded fluid cells — never
    `∫coef_p (s_body − s) dV` alone, which sees only ~63 % of the heat. Each
    interior face is visited once, as the LOW face of the cell that owns it,
    so blocks and ranks never double count; FACE_CLOSED faces are skipped and
    lose nothing (both their sides are solid by the removal criterion). An
    `adiabatic` scalar reports EXACTLY zero — it exchanges nothing with the
    body by construction — which the gate checks as a positive control on a
    manifestly non-zero field.
  - `tools/scalar_stats.py` (NEW) is the production reader: `profile`
    (mean/rms/turbulent flux, both wall fluxes, `theta_tau`, the total-flux
    constancy, Nusselt), `plane` (per-station wall flux and `Nu_x`), `heat`
    (the runtime file → Nusselt, with a literature band).
    `tools/compare_fields.py` with no dataset arguments now compares the
    datasets present in BOTH files (rank-of-`un` selects field datasets, so
    the `blocks` table and the node lines drop out; canonical `un vn wn pn`
    first, then scalars/`nut`/RANS variables alphabetically).
  - **Bit-exactness is again by construction**: with the intervals off (the
    default) nothing is allocated and no kernel is called; the only edits to
    existing files are the `[scalar]` key handlers, two `use` lines and three
    call sites in `moby_solve.f90`. Gated both ways — `run_bitexact.sh` /
    `run_bitexact_s3.sh` at max_abs 0, AND a statistics-ON vs statistics-OFF
    twin whose fields are bit-identical.
  - The runtime heat file is written at FULL double precision (`ES24.16`):
    `ES16.8` put a 1e-9 floor under the comparison with the independent
    Python transcription, which agrees to 1.9e-15 once the digits are there.

- **S3** — the immersed body. `ibm%coef` gains its cell-centred (VAR_P)
  column ONLY when `[scalar]` is configured (`init_ibm(ibm, blk,
  cell_centred)`), so every scalar-free case keeps its memory, its kernels
  and its case-file layout exactly. `set_ibm_coeff` / `set_ibm_coeff_host`
  were already generic in `var`: the guard now reads
  `var > ubound(ibm%coef,4)`, which admits `VAR_P` exactly when it exists.
  The analytic path adds one `set_ibm_coeff(dns, blk, ibm, VAR_P)` in
  `moby_solve.f90` and in `moby_prepare.f90` (the host twin there for STL
  geometry); the file path reads the NEW, OPTIONAL `coef_p_blocks` case-file
  dataset. The transport kernel gained the two wall modes and `dt_gamma`.
  - **`dirichlet`** (default): `ss = ss*mu_s + (1 - mu_s)*ibm_value` with
    `mu_s = 1/(1 + dt_gamma coef_p/Pr)` formed inline — no `mu` array per
    scalar — and `oldrhs` keeping the UNpenalized rhs, the momentum
    predictor's structure verbatim. The stored coefficient carries the
    `1/Re` scaling, so `coef_p/Pr` is the scalar's own `1/(Re Pr)` and ONE
    cell-centred array serves every scalar at its own `Pr`. Diffusive fluxes
    are not masked.
  - **`adiabatic`**: no penalization; the convective AND diffusive flux is
    masked on each of the six p-cell faces whose staggered velocity
    coefficient exceeds `SOLID_FACE_THRESHOLD` (rans.f90's `solw/sole/...`
    test for `k`). The mask is symmetric across a face, so the flux form
    still telescopes and the fluid conserves `∫s dV` exactly.
  - **DEVIATION from §4 (deliberate):** only `ibm%coef` gains the VAR_P
    extent — **`ibm%mu` keeps `VAR_U:VAR_W`**. The scalar penalization
    factor is Prandtl-dependent and formed inline, so a VAR_P `mu` plane
    would be a third more of a full field array, never read. `update_ibm_mu`
    is unchanged (it loops `VAR_U, VAR_W`).
  - **CORRECTION to landmine (9) / §4's `keep_buried` note:**
    `[blocks] keep_buried` is honoured **only in the `refine_body` branch**
    (`classify_refinement_masks`); a single-level case with the default
    `remove_solid` needs `[blocks] remove_solid = false` instead. The
    heated-cylinder gate sets both. For the THERMAL body integral this turns
    out to be belt-and-braces rather than load-bearing (unlike the A2
    penalization FORCE): a removed buried core consists of deep-solid cells,
    which hold the body value to the last bit and therefore contribute
    exactly zero to `∫coef_p (s_body − s) dV`.
  - **`coef_p_blocks`** is a 4-D per-leaf tile dataset shaped exactly like
    `dwall_blocks`; both now go through one pair of C helpers
    (`read_leaf_tiles` / `case_append_leaf_tiles`) parameterised by the
    destination stride and component-plane offset. The three coefficient
    entry points gained an `n_comp` argument = the Fortran array's component
    extent; the DATASETS are unchanged, so a case file for a scalar-free run
    is byte-identical to what P0–P3 produced (gated). With scalars on and no
    `coef_p_blocks` the solver hard-errors naming the fix.
  - **FOUND WHILE GATING:** that hard error is not hypothetical — it caught
    the inherited S1 `uniform3` gate, whose zero-force twin
    (`validation/multilevel_body/ibm_coeff_ml3_zero.h5`) predates S3.
    `make_uniform_twin.py` now writes a ZEROED `coef_p_blocks` alongside the
    zeroed `coef_blocks` (the scalar analogue of "the body exerts no force"
    is "the body penalises no scalar"). The twin is a GENERATED artifact, not
    a tracked file (`*.h5` is gitignored; `validation/multilevel_body/
    setup.sh` builds it), so the generator is the fix — a stale local copy is
    repaired by re-running `make_uniform_twin.py`, NOT by the "re-run
    moby_prepare with [scalar]" the solver's error message suggests, which is
    the right advice only for a real case file. Any pre-S3 coefficient file needs the same treatment or a
    re-run of `moby_prepare` with `[scalar]`.
  - **FINDING that changed the Nusselt gate (§9's S3 bullet is superseded):
    the A2 penalization integral does NOT transpose to a Dirichlet scalar.**
    `F = ∫coef·u dV` works for the FORCE only because `u_body = 0` and the
    stored velocity keeps a ~1e-26 residual, so `coef·u` is O(1). A Dirichlet
    scalar's stored value is the body value TO THE LAST BIT (that is gate
    (n)), so `coef_p (s_body − s)` evaluates to `1e28 × 0 = 0` in every solid
    cell — while the cell is really re-heated every substage by exactly the
    flux it loses to its fluid neighbours. Measured on the wavy case: solid
    cells contribute `0.000000e+00`, the graded fluid band `3.74e-02`, truth
    `6.95e-02`, i.e. the integral sees 63 % of the heat. The body heat
    release is therefore measured as **staircase-interface flux + graded-cell
    penalization** (`check_scalar_ibm.py surface`), which is
    cancellation-free, and validated by a full discrete energy budget: on a
    case with no boundary flux at all it equals `d/dt ∫s dV` to **3.9e-4**.
    Note this also means the `keep_buried` landmine is inert for a THERMAL
    body integral in both directions — the deep-solid cells the rule
    protects are exactly the ones the integral cannot see anyway.
  - **Bit-exactness has TWO tiers here.** Without a body (`dns%ibm_enabled`
    false) `useIbm` is false, `adiab` is `.false.` for every scalar, the six
    masks collapse to the FACE_CLOSED flags and the penalization statement
    is not entered — the S2 arithmetic survives BY CONSTRUCTION
    (`run_bitexact_s3.sh`). With a body whose coefficients are all zero (the
    `uniform3` twin) the statement IS entered and exactness rests on the
    IEEE identity `x*1.0 + 0.0*v = x` — the same class of argument as the
    IDDES `fd_force = 0` blend identity, and measured, not assumed.
  - **The heated-cylinder gate needed the `pn`-drift recipe** (a second
    instance of the A2 caveat, with a new trigger): the committed steady
    Re = 40 restart carries `|pn| ~ 1.2e3` of velocity-neutral mode whose
    spurious `∇p` is only self-consistent with the coefficient field that
    grew it. Restarting it as-is on the freshly PREPARED coefficients — which
    differ from the retired mobygeom ones only in the near-grazing envelope —
    kicked the forces to `|C_L| ~ 1.7e2`, oscillating with no decay. Zeroing
    `pn` fixes it completely AT THE PRODUCTION `niter = 6` (`C_D` back to 1.69
    within 40 steps, `|C_L| < 6e-4`); the `niter = 60` rebuild phase the A2
    momentum cross-check needed is unnecessary here because the thermal
    balance never reads the pressure.

- **S2** — the turbulent closure. The transport kernel's face diffusivity is
  now `D_face = 1/(Re Pr) + ½(ν_t,L + ν_t,R)/Pr_t(face)`, reading `turb%nut`
  — the ONE blended field LES (WALE/Smagorinsky), RANS (SST) and IDDES all
  write, so one code path covers all three. `Pr_t` is the per-scalar constant
  (`[scalar.N] prt`) or the **Kays–Crawford** correlation (`prt_model =
  kays`), implemented as the pure `!$omp declare target` function `prt_kays`
  and unit-tested host-side by `src/test_scalar.f90` (CMake target
  `scalar_test`) against an mpmath transcription. `[rans] wall_treatment =
  wall_function` together with scalars is a hard config error
  (`validate_turbulence_values`) — a thermal wall function is S5. §8 is
  finished: `get_timestep_rates` scales the eddy Peclet rate by
  `ν_eff = ire/Pr_min + ν_t/Pr_t,min` (`scalar_min_prt`), gated on
  `nScalar > 0` so both factors are exactly 1 otherwise.
  - **`prt_kays` numerics** (comment in scalar.f90): written directly, the
    correlation subtracts `(C Pe_t)²[1 − exp(−x)]` from a term that cancels
    it to leading order, so it loses all precision at large `Pe_t` and
    overflows `a²` past `Pe_t ~ 1e150`. The small-`x` branch evaluates the
    analytically equivalent series `1/Pr_t = (1/Prt_inf)[1 − x/3! + x²/4! −
    …]`, which also makes the `Pe_t → ∞` limit exact by inspection; the
    branches join at `x = 1/2` to their own slope.
  - **DEVIATION (deliberate, one line):** `scalar_transport` moved from
    "immediately before `momentum`" to "after the turbulence block, still
    before `momentum`" in `moby_solve.f90`. It still reads the
    start-of-substage `q` (nothing in between touches it — that was the S1
    reason for the placement), but it now sees THIS substage's `nut` with
    current halos, the same one the momentum predictor uses, instead of the
    previous substage's.
  - **ADDED beyond the plan (an S0 defect found while building the gate
    cases):** a scalar's explicit `<face>_type` key on a **declared** `wall`
    patch was a hard config error, because S0 applied the A0 contradiction
    rule to walls. §5 of this plan says "PATCH_WALL → Neumann 0 (adiabatic)
    **unless the ini gives a type/value**", and an isothermal wall is the
    canonical scalar BC, so `resolve_scalar_row` gained a `strict` flag: a
    wall honours an explicit type, inlet and outlet keep the strict rule (an
    inlet must impose a value, an outlet must be zero-gradient).
  - **Bit-exactness is by construction on BOTH counts.** `count = 0` allocates
    no scalar and calls no kernel (S0/S1's argument, re-gated). A scalar run
    with turbulence OFF keeps the S1 arithmetic byte-for-byte: every face
    diffusivity is initialised to the molecular value and the eddy term sits
    behind `if (useNut)`, so nothing is added — not even `+ 0.0`. When no
    model is active `turb%nut` does not exist at all, and the kernel is handed
    a 1-cell `sc%nutNone` dummy to map (the `turb_type` "1-cell dummies,
    uniform device maps, all accesses model-guarded" idiom). **No turbulent
    diffusivity exists at all in a DNS run** — every face keeps `1/(Re Pr)`,
    `eddy_diffusivity` is never called — and the dead branch was MEASURED, not
    assumed: interleaved S1-vs-S2 nofma runs on a DNS channel + scalar give
    S2/S1 = 0.9998 (CPU 4 ranks) and 1.0038 (GPU), both inside the noise floor
    of the `count = 0` control (0.9962 / 1.0068), whose code path is
    identical. Numbers in `validation/scalar/README.md`.

S4 gates (all in `validation/scalar/`, README records the commands and
numbers; `run_gates_s4.sh [stats|accum|plane|levels|restart|det|noeffect|heat|adia|cyl|tools]`):

| gate | result |
|---|---|
| the solver's rows vs the SNAPSHOT's rows, all seven columns, one sample (LES channel, 2 scalars) | **2.9e-14** / 2.8e-14 (constant `Pr_t` / Kays) |
| the same, four accumulated samples | **2.6e-14** / 2.7e-14 |
| the `plane` (boundary-layer) layout | **4.8e-16** / 4.5e-16 |
| the 2:1 wall-band case: per-LEVEL row tables (16/48 core rows, 64/96 band rows) | **2.6e-14** (level 0), **2.1e-13** (level 1) |
| restart continuation: 2 + 2 samples == one 40-step run | **0.000e+00** (and the fields are bit-identical) |
| statistics 1 rank == 4 ranks / CPU == GPU (atomics + allreduce reorder the sums) | 2.3e-14 / 9.4e-16 |
| statistics ON vs OFF: the fields | max_abs **0** |
| body heat release vs `check_scalar_ibm.py surface` (20 samples, both terms) | **2.8e-15** |
| the same, energy budget with the SOLVER's own `Q` (10-step window) | **4.5e-04** |
| an `adiabatic` scalar's heat columns, on a manifestly non-zero field | **0.000e+00** exactly |
| heated cylinder Re 40: the RUNTIME Nusselt number | **3.3653** (S3 post-processed 3.3655, Churchill–Bernstein 3.35) |
| `compare_fields.py` with no dataset arguments | discovers `un vn wn pn nut theta theta_kc` |
| `[scalar] count = 0` vs the S3 binaries, 7-case suite, nofma | max_abs **0**, CPU **and** GPU |
| the S1 + S2 + S3 scalar cases vs the S3 binaries, nofma | max_abs **0** |
| every S1 / S3 gate re-run with the S4 binary (19/19 body gates; `scalar_test` ALL PASS) | identical numbers to every digit |

LANDMINE for the re-runs: `run_gates.sh`'s `det` group compares CPU vs GPU at
TOLERANCE 0, so it must be given the nofma binaries
(`BIN=…/build_cpu_nofma/moby_solve GBIN=…/build_gpu_nofma/moby_solve`).
With the default FMA-contracted builds it reports ~1e-14 on the velocities and
6.7e-16 on the scalar and the suite ends in "FAILURES" — the documented
contraction difference (README S0/S1 (f)), not a regression.

S3 gates (all in `validation/scalar/`, README records the commands and
numbers; `run_gates_s3.sh [solid|conserve|balance|prep|missing|refine|det|cyl]`):

| gate | result |
|---|---|
| `dirichlet`: solid cell == `ibm_value` (analytic, restarted, `refine_body`, file path) | **0.000e+00** in all four, on 704 / 704 / 5984 / 704 solid cells |
| the analytic solid set vs the prepared `coef_p` solid set | **identical** (704 cells) |
| `adiabatic`: `∫φ dV` with a body, 200 steps | drift **0.000e+00** (of 1.25e-01) |
| `adiabatic`: cells sealed on all six staggered faces | 624 cells, max\|Δφ\| **0.000e+00** |
| body heat release vs `d/dt ∫θ dV` (no boundary flux), 10-step window | **3.9e-04** (200-step window: 3.1 %, the trapezoid's `O(dt²)`) |
| the A2 penalization integral on the same case | sees **63 %** of the heat — see the FINDING above |
| prepared case file vs the `[scalar]`-stripped twin | dataset-identical apart from `coef_p_blocks` |
| `coef_p_blocks` vs an independent transcription of the graded formula | **1.4e-16** (single level), **2.2e-16** (level 1 of `refine_body`) |
| prepare 1 rank == 4 ranks | identical file |
| solve from the case file == the inline analytic solve (single level / 2 levels) | max_abs **0** incl. both scalars |
| `[scalar]` + a case file without `coef_p_blocks` | rejected, message names the fix |
| heated cylinder Re 40, Pr 0.71: Nu (body-local, t = 105→120) | 3.520 → 3.415 → 3.380 → **3.366**, last drift 0.42 % — Churchill–Bernstein **3.35**, numerical band 3.2–3.5 |
| the same, vs the INDEPENDENT Gauss/CV border flux at 3 box sizes | **0.90 % / 1.98 % / 3.20 %** (the A2 box-size signature) |
| 1 rank == 4 ranks, CPU == GPU (scalar + IBM), analytic path, file path, +LES | **EXACT** (max_abs 0, incl. `nut`) |
| `[scalar] count = 0` vs the S2 binaries, 7-case suite, nofma | max_abs **0**, CPU **and** GPU |
| scalar WITHOUT a body vs the S2 binaries, 9 scalar cases (S1 + S2), nofma | max_abs **0**, CPU **and** GPU |
| every S1 gate re-run with the S3 binary | identical numbers to every digit |

S2 gates (all in `validation/scalar/`, README records the commands and
numbers; `run_gates_s2.sh [kays|wferr|sst|les|band|det]`):

| gate | result |
|---|---|
| Kays–Crawford unit test, every branch + both limits | 29 tabulated values to **1e-13** rel; `Pe_t = 1e300` → 0.85 exactly; branches join at `x = 1/2`; monotone over 12 decades |
| steady SST, resolved walls: `theta` vs the `nut`-integral prediction | **0.096 %** of the wall difference; flux constancy **5.2e-08**; `θ⁺/(Pr y⁺) = 1.0000` |
| the same with `prt_model = kays` | **0.094 %**; flux constancy **2.6e-08**; `Pr_t` spans **0.8832 … 1.7000** in-kernel (its `Pe_t → ∞ / → 0` limits), wall flux 7.4 % below the constant-`Pr_t` twin |
| LES channel Re_τ 180, Pr 0.71: `θ⁺/U⁺` | first cell **0.7233** vs `Pr` 0.71 (1.9 %); log layer **0.8558** vs `Pr_t` 0.85 (0.7 %) |
| the same, `θ⁺` vs Kader over the wall layer y⁺ ∈ [1,35] | mean **6.8 %**, max 20.7 % — at y⁺ 6.3, on Kader's own blend kink |
| the same, resolved `<v'θ'>/J` | 0.07 / 0.52 / 0.80 / 0.92 / 0.94 at y⁺ 5 / 15 / 30 / 60 / 120, peak **0.95 θ_τ u_τ** |
| the same, constant-flux residual (convergence + conservation) | **3.3 %** of `θ_τ`; `θ_τ` drifts **1.0 %** between the halves of the window |
| 2:1 wall-band channel: NO spurious scalar band | `θ'_rms` interface/core ratio excess **+0.12 %**; no row within six of the interface off the control by >1.1 % |
| `wall_treatment = wall_function` + scalars | rejected, explicit message |
| 1 rank == 4 ranks, CPU == GPU (LES + scalar) | **EXACT** (max_abs 0 incl. `nut` and `theta`) |
| `[scalar] count = 0` vs the S1 binaries, 7-case suite, nofma | max_abs **0**, CPU **and** GPU |
| scalar with turbulence OFF vs the S1 binaries, nofma | max_abs **0** on 5 scalar cases |
| every S1 gate re-run with the S2 binary | identical numbers to every digit |

LANDMINE for the LES gates (cost an hour): the `../channel_interface/les`
campaign's **final** `*_50001.h5` files carry an O(1e6) velocity-neutral
pressure mode (the niter = 6 pn-drift family, CLAUDE.md A2) — restarting on
one blows the run up in a single step. Every mid-run snapshot has `|pn| ~ 10`;
the gate cases restart from `*_49600.h5`.

METHOD note for the LES gates: the thermal field relaxes on
`tau = L²/(pi² D_eff) ~ 6 t.u.`, so reaching a stationary mean from any IC is
~3 tau of wall clock in which nothing is measured. `make_theta_ic.py` seeds
the Reynolds analogy of the run's own mean velocity, and `retarget_theta.py`
then replaces the mean by the profile the run's own MEASURED total
diffusivity implies while keeping every fluctuation. Neither is circular with
the gate: both derive from the run itself, not from Kader, and stationarity
is measured afterwards (constant-flux residual + first-half/second-half
`θ_τ`), not assumed.

- **S0** — layout/config/io/halos/BCs. Scalars are extra variables of `blk%q`
  (`dns%nVar = NVAR + dns%nScalar`, `VAR_S0+is` in `q`, `SCR_S0+is` in
  `qs`/`oldrhs`, both in init.f90); new `src/modules/scalar.f90` owns
  `scalar_type` (allocatable per-scalar arrays, no `MAX_SCALARS` anywhere) +
  `enter_/exit_scalar_data`; `[scalar]` / `[scalar.N]` parsed straight into it
  (`apply_scalar_config`, grown on demand); `comm%activeVars` is allocatable
  sized `dns%nVar` and the send/recv buffers are sized with it;
  `fdm_h5_write_field`/`fdm_h5_read_field` take a variable count + a packed
  32-char name table (`field_var_names`, `FDM_VAR_NAME_LEN`) plus a
  per-variable `found` flag (`read_force_file` still passes 4);
  `apply_scalar_bc_q` (boundary.f90) is the var-indexed twin of `apply_bc`'s
  cell-centred branch; patch-derived BC defaults resolve in `init_scalar`
  (inlet ⇒ Dirichlet `inlet`, wall/outlet/generic ⇒ Neumann 0, explicit keys
  win, contradiction = hard error).
- **S1** — the transport kernel: one fused `scalar_transport` (collapse(4),
  inner scalar loop, one GPU launch for any N), 2nd-order central
  divergence-form convection on the p cell's own face velocities, `[flow]
  convection = skew` honoured (`conv_div − s·(div u)|stencil`), face-flux
  diffusion on the module's own `invD*` tables, molecular diffusivity only,
  momentum RK3 structure verbatim, no upwind option; `scalar_finish` (qs→q,
  ghosts, ONE batched `exchange_halos` over the scalar vars) after
  `pressure_projection`. Peclet limit scaled by `1/Pr_min` (gated on
  `nScalar > 0`; the `ν_t/Pr_t` half of §8 belongs to S2, where ν_t first
  enters the scalar diffusivity).
- **DEVIATION from §2 of this plan (deliberate, one line):**
  `scalar_transport` is called **immediately BEFORE `momentum(...)`**, not
  after it. The plan's own rationale — "the advecting velocity is the previous
  substage's projected (divergence-free to solver tolerance) field, the same
  field momentum advects with" — is only true before `momentum`, because
  `momentum` ends by copying `qs → q` AND the velocity halos are refreshed
  only by the exchange that follows it. Calling after would advect the scalar
  with the non-solenoidal predicted velocity and with stale halo values. Both
  kernels read the same `q` and write disjoint `qs`/`oldrhs` slots, so the
  order is otherwise free; everything stays outside the projection.
- **ADDED beyond the plan:** the diffusive flux is masked at `FACE_CLOSED`
  faces (blocks removed inside the immersed body hold a zeroed halo; the
  convective flux already vanishes there because the face velocity is pinned).
- **NOTE on `convection = skew` for a scalar** (implemented as specified, flagged
  here): `conv_div − s·(div u)` is the ADVECTIVE form — it preserves a uniform
  scalar exactly for ANY advecting field but gives up the divergence form's
  exact global conservation. The skew-SYMMETRIC form (what the momentum kernel
  uses, neutral in `∑s²`) would subtract HALF of it, and would leave
  `½ s div u` on a uniform field. Flip the factor if `∑s²` neutrality is
  wanted instead; every gate above runs the default divergence form.

Gates (all in `validation/scalar/`, README records the commands and numbers):

| gate | result |
|---|---|
| `[scalar] count = 0` bit-exact, 7-case suite, nofma, CPU **and** GPU | max_abs **0** on every dataset |
| uniform scalar through a 3-level refined layout (faces/edges/corners) | max\|s − c\| = **0.0** (both scalars), velocities still 0.0 |
| global conservation, periodic box, Beltrami flow | drift **−6.8e-19** of max\|s\|·V over 200 steps |
| pure conduction, linear profile on a stretched grid | exact fixed point (CPU **0.0**, GPU 1.1e-16) |
| conduction with a source, stretched grid, ny 16/32/64 | order **1.996 / 1.999** |
| advection–diffusion MMS (uniform oblique flow), n 16/32/64 | order **1.99 / 2.00** |
| Pr sweep 0.1 / 1 / 10, laminar channel vs the analytic series | max dev **1.0e-5 / 1.4e-4 / 1.6e-3**, wall-flux ratio 1.0011 / 1.0004 / 1.0041 |
| 1 rank == 4 ranks, CPU == GPU (refined transport case) | **EXACT** (max_abs 0 incl. the scalar) |
| S0 restart round-trip (dataset present / absent) | exact read / warn + reinit |

LANDMINE found while implementing (cost an hour): a local variable named
`nVar` in `blocks.f90` SHADOWS the use-associated `NVAR` parameter (Fortran is
case-insensitive), so `allocate(q(..., nVar, ...))` silently allocated a
zero-size variable dimension and every kernel wrote outside `q`. The local is
now `nQ`.

---

## Original plan (2026-08-03)

Goal: a user-selectable number of passive scalars `s_1..s_N` (N = 0 default),
each with its own Prandtl/Schmidt number, stored at the PRESSURE point (cell
centre), transported by the flow, coupled to the LES SGS model and to the RANS
SST model through a Reynolds analogy (constant `Pr_t` or the Kays–Crawford
correlation), with its own immersed-boundary coefficients evaluated at the cell
centre. Minimally invasive: ONE solution array, ONE output file, and
bit-exactness with the current solver when `N = 0` BY CONSTRUCTION.

---

## 0. The one design decision everything follows from

**Scalars are extra variables of `blk%q`.** `NVAR` (currently the parameter 4 =
u,v,w,p) becomes the runtime count `dns%nVar = 4 + nScalar`; scalar `is` lives
at `VAR_S0 + is` with `VAR_S0 = 4`.

Why this is the cheap route, not the expensive one:

1. **The metric arrays do not grow.** `blk%x/d1x/lapX*` etc. are dimensioned
   `(...,NVAR,nBlocks)` because the four variables sit at four different
   staggered positions. A cell-centred scalar sits EXACTLY at the pressure
   position, so every scalar kernel reads the `VAR_P` column
   (`blk%d1x(i,VAR_P,b)`, `blk%lapX0(i,VAR_P,b)`, `blk%x(i,VAR_P,b)`). The
   metric arrays keep their `VAR_U:VAR_P` extent forever.
2. **The 2:1-interface halo exchange comes for free and is already correct.**
   `copy_local_entries` / `pack_entries` (comm.f90) treat any variable index
   `> 3` as cell-centred: all three dims take the `lGC == 2` restrict average
   (8-cell fine→coarse average), the prolong injects the covering coarse cell,
   and the `(2·coarse + fine)/3` ghost blend (`lWp`/`lWpDst`) is applied — the
   pressure treatment, which is precisely what a cell-centred advected scalar
   wants. The interface-normal skip `var /= c%lNrm(e)` can never match a scalar
   index (lNrm ∈ {0,1,2,3}), so scalar halos are always written. Scalars ride
   the SAME MPI message as u,v,w,p — no extra rounds, no new pack/unpack code.
3. **Output and restart are one file, one collective write**: `fdm_h5_write_field`
   already receives the whole `blk%q`; it needs a variable COUNT and a name
   table instead of the hardcoded `4` / `{"un","vn","wn","pn"}`.
4. **`N = 0` is bit-exact by construction**: `dns%nVar == 4` reproduces every
   allocation shape exactly, the scalar kernels are never called, `apply_bc`'s
   `do var = VAR_U, VAR_P` loop is untouched, and no extra exchange is issued.
   This is the same argument the body-force increment used (separate correction
   kernel, fused predictor byte-for-byte untouched) — not a `+0.0` argument.

Rejected alternative: separate `scalar_type` arrays in the style of
`sst%k/omg` with `exchange_scalar_halos` per scalar. It needs no comm changes
but costs one MPI round per scalar per substage, duplicates the io path, and
contradicts the "one array / one file" requirement.

### Scratch and RK storage

- `blk%q(0:nb+1,...,1:nVar,nBlocks)` — u,v,w,p,s_1..s_N.
- `blk%qs(0:nb+1,...,1:NVEL+nScalar,nBlocks)` — the substage scratch.
- `blk%oldrhs(1:nb,...,1:NVEL+nScalar,nBlocks)` — RK3 low-storage.

LANDMINE to document in `block_set_type`: the scalar index differs between the
arrays — `q` slot `4 + is`, `qs`/`oldrhs` slot `3 + is` (p has no scratch). Two
inline mappings, defined once in init.f90:

```fortran
integer(C_INT), parameter :: VAR_S0 = NVAR      ! q:      VAR_S0 + is
integer(C_INT), parameter :: SCR_S0 = NVEL      ! qs/rhs: SCR_S0 + is
```

The alternative (waste one `qs` plane so the indices coincide) costs ~25 % of
`qs`+`oldrhs` on the GPU and is not worth it.

### Per-scalar state lives in its own derived type — NO fixed bound

There is **no `MAX_SCALARS`**. `src/modules/scalar.f90` defines `scalar_type`,
which owns flat contiguous ALLOCATABLE arrays sized `nScalar` and maps them to
the device once in `enter_scalar_data` / `exit_scalar_data` — the pattern
`turb_type`, `sst_type` and `bodyforce_type` already follow (CLAUDE.md:
"Derived types own flat contiguous allocatable arrays; map them to the device
once"). Nothing per-scalar goes into `dns_type` (which is mapped as a flat
value type and would force a fixed bound), and nothing goes into
`boundary_type`:

```fortran
type :: scalar_type
    integer(C_INT) :: n = 0                       ! number of scalars
    real(C_DOUBLE), allocatable :: pr(:), prt(:)  ! (n)
    integer(C_INT), allocatable :: prtModel(:)    ! (n) constant | kays
    real(C_DOUBLE), allocatable :: source(:), init(:), ibmValue(:)
    integer(C_INT), allocatable :: ibmMode(:)     ! (n) dirichlet | adiabatic
    integer(C_INT), allocatable :: bcType(:,:)    ! (n, NFACES)
    real(C_DOUBLE), allocatable :: bcValue(:,:)
    logical,        allocatable :: bcTypeSet(:,:), bcValueSet(:,:)   ! host only
    real(C_DOUBLE), allocatable :: invDx(:,:), invDy(:,:), invDz(:,:)! (0:nb+1, nBlocks)
    character(len=32), allocatable :: name(:)     ! host only (io labels)
end type scalar_type
```

`dns_type` gains only two plain integers: `dns%nScalar` (so `blocks.f90` can
size `q`/`qs`/`oldrhs`) and the derived `dns%nVar = NVAR + dns%nScalar`.

Config parsing writes straight into `scalar_type`: `read_runtime_config`
already takes `turb`, `les`, `ps` and `bc`, so `sc` joins that list, and the
`[scalar.N]` handler grows the arrays on demand (host-only init code — a
grow-by-copy on each new section index, or a two-pass count-then-fill; either
is fine at ini-parse cost). `name`, `bcTypeSet`, `bcValueSet` are never read in
a kernel and stay host-side.

---

## 2. Governing equation and discretisation

For each scalar (non-dimensional, `Re` the solver's Reynolds number):

```
∂s/∂t + ∇·(u s) = ∇·( (1/(Re·Pr) + ν_t/Pr_t) ∇s ) + q_src
                  − coef_p·(s − s_body)                    [IBM, Dirichlet mode]
```

Staggered layout makes this simpler than the momentum equation: the p-cell's
six faces carry exactly the six staggered velocity components, so the
divergence form needs no velocity interpolation and is discretely conservative.

**Convection** (2nd-order central, divergence form), per cell `(i,j,k)` of block `b`:

```
conv = ( q(i+1,j,k,VAR_U,b)·½(s_i + s_{i+1}) − q(i,j,k,VAR_U,b)·½(s_{i−1} + s_i) )·d1x(i,VAR_P,b)
     + (y-term with VAR_V, d1y(j,VAR_P,b)) + (z-term with VAR_W, d1z(k,VAR_P,b))
```

`[flow] convection = skew` must apply to the scalar too: `conv_skew = conv_div
− s·(div u)|stencil`, the divergence built from the SAME face velocities. The
momentum lesson (docs/next_session_skew_convection.md, the C11 interface
instability) is that the divergence form is only neutral for exactly
divergence-free advection, which the incremental projection and the 2:1
interface halos never grant.

**No upwind option.** First-order upwind exists in exactly one place today —
rans.f90's k/ω/γ/R̃e_θt transport (lines 1641-1665, 1768-1783) — as a documented
deviation awaiting the TVD/van-Leer increment, and the project stance is that
global upwinding is off-limits in an energy-conserving code
(docs/next_session_airfoil.md:104). Passive scalars therefore get central /
skew only; when the deferred TVD increment lands (it needs a second upwind halo
cell, the blocker shared with the RANS scalars) it serves both from the same
machinery. Until then, sharp-front over/undershoot is a KNOWN limitation, not a
configuration choice.

**Diffusion** (face-flux form, variable diffusivity):

```
D_face = 1/(Re·Pr) + ½(ν_t,L + ν_t,R)/Pr_t(face)
flux   = D_face·(s_R − s_L)·inv_d(face)        [0 if the face is masked]
diff   = Σ_dim (flux_+ − flux_−)·d1{x,y,z}(·,VAR_P,b)
```

`inv_d` = inverse centre-to-centre distance at the p-position. These tables live
in `turb_type` (`turb%inv_dx(:,VAR_P,:)`) but `init_turbulence` returns early
when the model is `none`, so a DNS run has none. The scalar module builds its
OWN three small tables `(0:nb+1, nBlocks)` at init (mirror of the turbulence
ones, ~30 lines) — allocated only when scalars are on, no module coupling, no
branch in the kernel.

**Time integration**: the momentum RK3, verbatim in structure:

```
ss = s + dt_alpha·rhs + dt_beta·oldrhs_s
ss = ss·mu_s + (1 − mu_s)·s_body            [IBM Dirichlet mode only]
oldrhs_s = rhs
```

with the implicit penalization factor `mu_s = 1/(1 + dt_gamma·coef_p/Pr)`
computed inline (one FMA + one divide per cell per scalar — no `mu` array per
scalar). Inside the body `coef_p = SOLID/Re = 1e30/Re`, so `mu_s → 0` and
`ss → s_body` exactly, with no `dt` restriction (ibm.f90's implicit form).

**Placement in the substage** (moby_solve.f90 main loop), all of it OUTSIDE the
projection:

```
update_ibm_mu / bodyforce / turbulence(nut) / exchange nut   [unchanged]
momentum(...)                                                [unchanged]
scalar_transport(...)          <-- NEW: reads q (velocities + scalars) + nut, writes qs scalar slots
apply_bc / exchange velocities / pressure_projection         [unchanged]
scalar_finish(...)             <-- NEW: qs→q copyback, scalar ghosts, exchange_halos(scalar vars)
```

`scalar_transport` and `momentum` both read the start-of-substage `q`, so their
order is irrelevant and the advecting velocity is the previous substage's
projected (divergence-free to solver tolerance) field — the same field momentum
advects with. One extra exchange per substage when scalars are on (batched over
all scalars); merging it into the projection's final exchange is a later
optimisation and deliberately NOT done now (the projection is bit-exactness-
critical).

All scalars are handled by ONE fused kernel (`collapse(4)` over b,k,j,i with an
inner `do is = 1, nScalar`), so the GPU sees one launch regardless of N.

---

## 3. Turbulence coupling

The scalar kernel reads `turb%nut` — the single blended field that LES (WALE /
Smagorinsky), RANS (SST) and IDDES all write. One code path therefore covers
all three models:

- **LES**: `D_sgs = ν_t/Pr_t` with constant `Pr_t` (default 0.85 for the SGS
  turbulent Prandtl number; typical LES values 0.4–0.9 — the default is
  configurable per scalar). `ibm_aware` already zeroes `ν_t` in solid cells and
  WALE gives `ν_t → 0` at the wall, so no extra masking is needed.
- **RANS (k-ω SST) / IDDES**: identical expression, Reynolds analogy. `Pr_t`
  either constant or from the **Kays–Crawford** correlation, evaluated
  pointwise from the local eddy-viscosity ratio:

  ```
  Pe_t  = (ν_t/ν)·Pr
  1/Pr_t = 1/(2·Pr_t∞) + 0.3·Pe_t/sqrt(Pr_t∞) − (0.3·Pe_t)²·(1 − exp(−1/(0.3·Pe_t·sqrt(Pr_t∞))))
  Pr_t∞ = 0.85
  ```

  Implemented as a pure `!$omp declare target` function next to the `lm_*`
  correlations in the scalar module and unit-tested against tabulated values
  the way `src/test_transition.f90` tests the LM correlations (add the cases to
  that driver or a new `src/test_scalar.f90`). Guard `Pe_t → 0` (returns
  `2·Pr_t∞`… clamp to `Pr_t∞` at the molecular limit) and `Pe_t → ∞`.

- **RANS wall functions**: `[rans] wall_treatment = wall_function` together
  with scalars is a **hard config error** in this phase — a thermal wall
  function (Kader/Jayatilleke P-function) is a separate, separately-validated
  increment (S5). Resolved-wall RANS works: `ν_t → 0` at the wall leaves the
  molecular diffusivity, which is the correct near-wall limit.

- **Buoyancy is explicitly out of scope.** The hook is obvious and cheap later:
  a Boussinesq source is `bodyforce.f90`'s `custom` hook reading `q(...,VAR_S0+is,...)`.

---

## 4. Immersed boundary: cell-centred coefficients

`set_ibm_coeff` / `set_ibm_coeff_host` (ibm.f90) are ALREADY generic in `var` —
they read `blk%x(ix,var,b)` and the graded second-order stencil from the same
column. The only thing blocking `VAR_P` is the guard
`if (var < VAR_U .or. var > VAR_W) error stop`. So:

- `ibm%coef` / `ibm%mu` allocate `VAR_U:VAR_P` (4 components) when scalars are
  on, `VAR_U:VAR_W` otherwise (unchanged memory for every existing case).
- Analytic IBM (inline): one extra `call set_ibm_coeff(dns, blk, ibm, VAR_P)`
  in moby_solve.f90 / the prepare path.
- The stored coefficient carries the `1/Re` scaling; the scalar uses
  `coef_p/Pr`, i.e. `1/(Re·Pr)` = the scalar diffusivity — the consistent
  scaling for the graded sharp-interface formula `Σ((d0−d)/d)/d0²`. One
  cell-centred coefficient array therefore serves ALL scalars, each with its
  own Pr.

**Two body wall modes per scalar** (`[scalar.N] ibm_wall`):

- `dirichlet` (default): penalization toward `ibm_value`. Diffusive fluxes into
  the body are NOT masked (the solid cell holds the wall value and delivers the
  flux). First-order staircase accuracy, exactly as for the velocity.
- `adiabatic`: no penalization; instead the six p-cell faces are masked when the
  corresponding staggered velocity coefficient is solid — the `solw/sole/...`
  test rans.f90 already uses for k. Convective AND diffusive flux masked ⇒
  fluid cells conserve the scalar exactly (zero flux through the body).

**File-based IBM (STL / `moby_prepare`)**: the cell-centred tiles are a NEW,
OPTIONAL case-file dataset `coef_p_blocks`, written by `moby_prepare` when the
prepare ini declares scalars (one extra `set_ibm_coeff_host(..., VAR_P, inside)`
+ an append function mirroring `fdm_h5_case_append_coef`). Existing case files
stay byte-identical and readable; the solver hard-errors with an explicit
"re-run moby_prepare with [scalar]" message if scalars are on and the dataset
is absent. `dwall_blocks`, masks, `block_active` are untouched.

LANDMINE (inherited from A2/A3): if a body-integrated flux diagnostic
(Nusselt from `∫coef_p·(s − s_body) dV`) is added later, it needs
`[blocks] keep_buried` for the same reason the penalization forces do — a
removed buried core absorbs loading outside the coefficient bookkeeping.

---

## 5. Boundary conditions

`apply_bc`'s non-normal branch (boundary.f90:697-719) is EXACTLY the
cell-centred treatment a scalar needs: Dirichlet via the ghost mirror
`ghost = 2·value − interior`, Neumann via `ghost = interior + dn·value`.

Do NOT extend `apply_bc`'s `do var = VAR_U, VAR_P` loop: `apply_bc` is called
`nIter` times per substage from INSIDE the projection loop, and it is a
bit-exactness-critical kernel. Instead:

- the scalar BC tables live in `scalar_type` (`sc%bcType(n,NFACES)`,
  `sc%bcValue`), NOT in `boundary_type`: `bc` keeps its current fixed
  `VAR_U:VAR_P` shape, its device maps stay trivial, and the restart
  `bc_type`/`bc_value` attributes keep their 24-entry layout, so restart
  metadata is byte-unchanged. `bc` is still the owner of the boundary POINT
  lists, which the scalar BC kernel reuses.
- a NEW `apply_scalar_bc_q(blk, bc, sc)` kernel — the mechanics of
  `apply_bc`'s cell-centred branch over the same `bc` point lists, indexing
  `blk%q(...,VAR_S0+is,...)`. (The existing generic `apply_scalar_bc` takes a
  standalone `s(0:,0:,0:,1:)` array and cannot be fed a strided `q` slice
  portably under OpenMP target mapping — hence a `var`-indexed twin rather than
  a refactor.) Called once per substage from `scalar_finish`.
- **patch-derived defaults**, set-if-unset in the `resolve_face_bcs` slot, so
  the A0 face concept extends naturally:
  - `PATCH_WALL` → Neumann 0 (adiabatic) unless the ini gives a type/value
  - `PATCH_INLET` → Dirichlet with the scalar's `inlet` value
  - `PATCH_OUTLET` → Neumann 0 (zero gradient; `BC_OUTFLOW` is a
    face-staggered concept and does not apply to a cell-centred scalar)
  - `PATCH_GENERIC` / unset → Neumann 0
  Explicit `[scalar.N] <face>_type/_value` keys win, and a contradiction with
  the declared patch is a hard config error (same rule as A0).

---

## 6. Configuration surface

Numbered sections, parsed like the existing `[grid.x]` axis sections
(`apply_grid_axis_value` / `grid_axis_index` is the pattern to copy):

```ini
[scalar]
count = 2                 ; optional; otherwise inferred as the highest [scalar.N]

[scalar.1]
name        = theta       ; HDF5 dataset name (default s1); validated against the
                          ; reserved set {un,vn,wn,pn,nut,k,omega,gamma,rethetat,fd,x,y,z,blocks}
pr          = 0.71        ; molecular Prandtl/Schmidt number
prt         = 0.85        ; constant turbulent Prandtl number
prt_model   = constant    ; constant | kays
initial     = 0.0         ; uniform initial value
init_profile = uniform    ; uniform | linear_y   (channel hot/cold-wall convenience)
source      = 0.0         ; constant volumetric source
ibm_wall    = dirichlet   ; dirichlet | adiabatic
ibm_value   = 1.0
y_min_type  = dirichlet   ; per-face overrides, same key style as [boundary]
y_min_value = 1.0
y_max_type  = dirichlet
y_max_value = 0.0
inlet       = 0.0         ; value seeded at PATCH_INLET faces
```

Storage: `scalar_type`'s allocatable arrays (§0), grown as `[scalar.N]`
sections are parsed. Sections must be contiguous from 1; a gap, or a `count`
that disagrees with the highest section index, is a hard config error.

---

## 7. IO — one file

- `fdm_h5_write_field` gains `n_var` and a packed name table; `block_stride =
  var_stride*n_var`; the loop writes `n_var` datasets. With `n_var == 4` and
  the existing names the output is byte-identical.
- `fdm_h5_read_field` likewise, plus a per-variable `found` flag: a restart
  file without a scalar dataset (older run, or a newly added scalar) →
  warn + initialise from `[scalar.N] initial` — the precedent set by the RANS
  named-scalar restart.
- `read_force_file`'s `qtmp` must keep passing `4` explicitly (it reads a
  velocity-layout file, not the scalar-extended one).
- No new append call, no second file: scalars are part of the same
  collective `write_field`, in the same `blocks`-table layout. `nut`, `k`,
  `omega`, … keep their existing append path unchanged.
- `tools/compare_fields.py`: `FIELDS` becomes "the datasets present in both
  files" (or accepts the scalar names on the command line — it already takes a
  `datasets` positional).

---

## 8. Time-step limits

- CFL: unchanged (same velocities).
- Peclet: `dns%peclet_rate` is `max over grid of ν/h²` — multiply by
  `max(1, 1/min_s Pr_s)` in `precompute_peclet_rate`, and in
  `get_timestep_rates` use `ν_eff = ire/Pr_min + ν_t/Pr_t,min`. Both are gated
  on `nScalar > 0`, so a scalar-free run keeps the exact current limit.
- Watch the A3 lesson: with LES/RANS, the explicit eddy-diffusion limit is what
  actually binds in fine cells; a low-Pr scalar (liquid metals, Pr ≪ 1) tightens
  it further. Report the binding limit in the terminal line as today.

---

## 9. Increments and gates

Each increment ends with the standard bit-exactness gate: `count = 0` must be
max_abs 0 (nofma, CPU AND GPU) versus the pre-increment binaries on the
standard 7-case suite (min_channel, les_ibm ± refine_body, Beltrami y-slab,
turb180, wf180_y30, lam30t).

**S0 — layout + config + IO + halos + BCs (no transport).**
`dns%nVar`, q/qs/oldrhs extents, comm `activeVars` sizing, the C io var count,
config parsing, initial condition, `apply_scalar_bc_q`, the scalar exchange.
The scalar is initialised and carried but not advanced.
Gates: bit-exact off; a constant scalar survives write→restart→write
byte-identically; 1 rank == 4 ranks EXACT; CPU == GPU EXACT; a scalar declared
on a 2:1-refined case has correct halos (manufactured-linear-field check on
every exchange-written halo cell, the `MOBY_HALO_AUDIT` method).

**S1 — transport kernel (molecular, no IBM, no turbulence).**
Gates:
- *uniform-scalar preservation*: `s ≡ const` stays constant to 0.0 in an
  arbitrary flow, including through a 3-level refined patch and across
  edges/corners (the scalar analogue of the uniform-flow interface gate);
- *global conservation*: `d/dt ∫s dV` = boundary flux to round-off in a
  periodic box (adiabatic walls);
- *pure conduction*: no flow, Dirichlet walls → exact steady linear profile;
  order 2 on a stretched grid;
- *advection–diffusion MMS / Taylor–Green scalar*: order 2 in space;
- *Pr sweep* (0.1 / 1 / 10) in a laminar channel: thermal-BL thickness scaling
  and Nusselt vs the analytic value;
- 1 == 4 ranks EXACT, CPU == GPU EXACT.

**S2 — turbulent closure (`ν_t/Pr_t`, constant + Kays–Crawford).**
Gates: correlation unit test (tabulated values, every branch); turbulent
channel Re_τ 180, Pr 0.71 — mean `θ⁺` against Kader / Kawamura DNS, turbulent
heat flux profile; the same case under SST (resolved walls) — log-law slope;
2:1 wall-band channel: NO spurious scalar band across the interface (jump
ratios, the u'/v' band lesson from `validation/channel_interface/`);
`wall_treatment = wall_function` + scalars is rejected with a clear message.

**S3 — IBM cell-centred coefficients.** (DONE — see the STATUS header for
what actually landed and what this planned text got wrong.)
Analytic path + `moby_prepare` `coef_p_blocks` + both body wall modes.
Gates: solid-cell scalar == `ibm_value` exactly (Dirichlet mode); adiabatic
mode conserves `∫s dV` to round-off with a body in the domain; heated cylinder
(Re 40) Nu vs literature and vs a Gauss/CV border-flux cross-check (the A2
method); prepared case file identical to the previous one apart from the added
dataset (`h5same` minus `coef_p_blocks`); `refine_body` per-level coefficients
consistent.
SUPERSEDED IN ONE PLACE: "Nu from `∫coef_p (s − s_body) dV`, the A2 method"
does not work for a Dirichlet scalar — the solid cell's stored value IS the
body value bit-for-bit, so the product vanishes and ~37 % of the heat is
invisible. The gate measures the staircase solid/fluid face flux plus the
graded-cell penalization instead, validated by a full energy-budget closure
(3.9e-4) on a case with no boundary flux. Details in the STATUS header.

**S4 — statistics and tooling (optional but wanted for production).**
Channel/boundary-layer scalar statistics (mean profile, rms, turbulent flux,
Nusselt) alongside `channel_stats.f90`; `compare_fields.py` dataset discovery.

**S5 — deferred, each its own session:** thermal wall function
(Kader/Jayatilleke) for `wall_treatment = wall_function`; TVD/van-Leer scalar
convection (shared with the RANS scalars — blocked by the single halo layer);
Boussinesq buoyancy via the body-force custom hook.

---

## 10. Files touched (S0–S3)

| File | Change |
|---|---|
| `src/modules/init.f90` | `VAR_S0`, `SCR_S0`, `dns%nScalar`, `dns%nVar` (two integers — nothing per-scalar) |
| `src/modules/blocks.f90` | q/qs/oldrhs extents from `dns%nVar` |
| `src/modules/comm.f90` | `activeVars` allocatable (`dns%nVar`), `maxCount` sizing, the two bound checks (lines ~453, ~1629) |
| `src/modules/boundary.f90` | `apply_scalar_bc_q`, patch-derived defaults in `resolve_face_bcs` (tables live in `scalar_type`) |
| `src/modules/config.f90` | `[scalar]` / `[scalar.N]` parsing into `scalar_type` + validation |
| `src/modules/scalar.f90` (NEW) | `scalar_type` (per-scalar config, BC tables, metric tables), init/enter/exit, `scalar_transport`, `scalar_finish`, `prt_kays` |
| `src/modules/step.f90` | Peclet scaling only |
| `src/modules/ibm.f90` | `VAR_P` allowed in `set_ibm_coeff{,_host}`, coef/mu extent, `coef_p_blocks` read |
| `src/modules/io.f90`, `field_hdf5.c` | `n_var` + names in write/read field, `coef_p_blocks` case append |
| `src/moby_solve.f90` | two call sites + init/finalise |
| `src/moby_prepare.f90` | cell-centred coefficient pass + append |
| `src/test_scalar.f90` (NEW) | Kays–Crawford unit test |
| `tools/compare_fields.py` | dataset discovery |
| `validation/scalar/` (NEW) | the S1–S3 gate cases and drivers |

---

## 11. Landmines (read before writing code)

1. `comm.f90:117` `activeVars(NVAR)` becomes an ALLOCATABLE sized `dns%nVar`
   (`comm_init` already receives `dns`), and `comm.f90:453-454`'s `maxCount`
   send/recv buffer sizing must use `dns%nVar` too — otherwise scalar halos
   silently truncate. Same for the two bound checks at `comm.f90:1629/1634`.
2. The `q` vs `qs` scalar index offset (`VAR_S0+is` vs `SCR_S0+is`).
3. `field_hdf5.c` hardcodes `4` in THREE places (write loop, `block_stride`,
   the read path) plus `read_force_file`'s `qtmp`.
4. Do NOT add scalars to `apply_bc` (called `nIter`× inside the projection) or
   to the projection's exchange var lists.
5. Central convection is unbounded: expect over/undershoot at sharp fronts
   until the shared TVD increment lands. There is deliberately no upwind
   fallback (§2).
6. Scalars take the PRESSURE 2:1 transfer (blended prolong ghost), while the
   RANS scalars take plain injection — an intentional difference; the S2
   interface-band gate is what justifies it.
7. Penalization Dirichlet at the immersed body is first-order/staircase, like
   the velocity; the graded (second-order) coefficient mitigates it. Do not
   expect better than the velocity's near-body accuracy.
8. The advecting velocity is divergence-free only to the projection tolerance
   (niter = 6 production runs) — use `[flow] convection = skew` for the scalar
   in the same runs where it is used for momentum.
9. `keep_buried` if a body-integral flux diagnostic is added (A2 rule).

---

## 12. Next-session prompt (S5 — thermal wall function, TVD, Boussinesq)

Paste this to start the next implementation session. S5 is THREE independent
increments; do ONE per session.

> Implement increment **S5a/S5b/S5c** of `docs/next_session_scalar.md` (pick
> one) on branch `scalar` — read the doc in full first: its STATUS header
> records exactly what S0–S4 landed, every deviation, every finding and every
> gate number; also read the "Active work" section of CLAUDE.md for the
> project conventions, and `validation/scalar/README.md` for the gate
> machinery you inherit.
>
> S0–S4 are DONE and gated: scalars are extra variables of `blk%q`, the fused
> `scalar_transport` kernel carries central/skew convection, the turbulent
> closure `1/(Re Pr) + nut/Pr_t` (constant or Kays–Crawford) and BOTH
> immersed-body wall modes, `moby_prepare` writes `coef_p_blocks`, and
> `scalar_stats.f90` accumulates the channel/boundary-layer statistics and the
> body heat release while the solver runs.
>
> - **S5a — thermal wall function.** `[rans] wall_treatment = wall_function`
>   together with `[scalar]` is currently a HARD CONFIG ERROR
>   (`validate_turbulence_values`); lift it by adding the Kader/Jayatilleke
>   P-function to the scalar's wall treatment, mirroring T3's momentum wall
>   functions (`rans.f90`: the wall-cell `nut` and the `omega` blend are the
>   pattern, and the T3 gates in `validation/rans_sst/` are the shape of the
>   validation). Gate against the resolved-wall scalar channel the way T3
>   gated `wf180_y30` against `turb180`, and keep the resolved path bit-exact.
> - **S5b — TVD/van-Leer scalar convection.** Shared with the RANS scalars
>   (`rans.f90`'s documented first-order-upwind deviation) and blocked by the
>   SINGLE halo layer — a second upwind cell is needed and a block-edge
>   fallback would break the nb/rank-independence that Phase 1 established.
>   The measured motivation is in CLAUDE.md (the SD7003 γ front, 104 level-4
>   cells) — this increment is about the halo, not about the scalars.
> - **S5c — Boussinesq buoyancy.** The hook is `bodyforce.f90`'s `custom`
>   path reading `q(...,VAR_S0+is,...)`; the design note is §3 of this plan.
>   Deliberately out of scope until someone needs it.
>
> Conventions that are not negotiable: build both paths with the module
> loaded; always `mpirun`; save the S4 nofma binaries outside the tree BEFORE
> touching any source; `run_bitexact.sh` / `run_bitexact_s3.sh` must stay at
> max_abs 0; re-run `run_gates.sh`, `run_gates_s2.sh`, `run_gates_s3.sh` and
> `run_gates_s4.sh`; new gates in `validation/scalar/` with a README block;
> never declare an increment done with a failing build or an ungated result.

---

### Historical: the S4 prompt (consumed 2026-08-04)

> Implement increment **S4** of `docs/next_session_scalar.md` — passive-scalar
> STATISTICS AND TOOLING — on branch `scalar` (read the doc in full first: its
> STATUS header records exactly what S0/S1/S2/S3 landed, every deviation, every
> finding and every gate number; also read the "Active work" section of
> CLAUDE.md for the project conventions, and `validation/scalar/README.md` for
> the gate machinery you inherit).
>
> S0-S3 are DONE and gated: scalars are extra variables of `blk%q`, the fused
> `scalar_transport` kernel carries central/skew convection, the turbulent
> closure `1/(Re Pr) + nut/Pr_t` (constant or Kays-Crawford) and BOTH immersed-
> body wall modes (`dirichlet` penalization toward `ibm_value`, `adiabatic`
> face masking), `moby_prepare` writes the optional `coef_p_blocks` case-file
> dataset, and `[scalar] count = 0` / a scalar run without turbulence / a scalar
> run without a body are each bit-exact CPU+GPU against the previous
> increment's binaries.
>
> **S4 — statistics and tooling.** Section 9 of the plan: channel and
> boundary-layer scalar statistics (mean profile, rms, turbulent flux, Nusselt)
> alongside `channel_stats.f90`, and `tools/compare_fields.py` dataset
> discovery (`FIELDS` becomes "the datasets present in both files"). Read
> `validation/scalar/check_scalar_turb.py` first: its `channel` analysis is the
> post-processing form of exactly the statistics the solver should accumulate,
> and the S2 gate numbers it produced are what the in-solver version must
> reproduce.
>
> BEFORE YOU START, read the S3 FINDING in the STATUS header about the
> penalization integral: a Dirichlet body's heat release CANNOT be measured as
> `int coef_p (s_body - s) dV` (a solid cell holds the body value to the last
> bit, so the product is 1e28 x 0 = 0 and ~37 % of the heat is invisible). If
> S4 adds a runtime Nusselt statistic, it must use the cancellation-free form
> `check_scalar_ibm.py surface` uses — the staircase solid/fluid face flux plus
> the graded-cell penalization — or the control-volume balance.
>
> Conventions that are not negotiable:
> - Build both paths (`./compile.sh cpu && ./compile.sh gpu`) with the module
>   **loaded** (`module load toolkits/nvhpc/25.9`) — see the environment
>   landmine in Section 13; always launch through `mpirun`, even on one rank.
> - **Save the S3 nofma binaries outside the tree BEFORE touching any source.**
> - `[scalar] count = 0` must stay bit-exact (`run_bitexact.sh`, 7-case suite,
>   max_abs 0, CPU AND GPU), and `run_bitexact_s3.sh` (scalar without a body,
>   9 cases) must stay max_abs 0 — that pair is also the cheap, sharp form of
>   "every S1/S2/S3 gate still reads the same number". Re-run `run_gates.sh`,
>   `run_gates_s2.sh` and `run_gates_s3.sh`.
> - New gate cases and drivers go in `validation/scalar/`, and the README grows
>   an S4 block with what was run and measured.
> - Never declare an increment done with a failing build or an ungated result.
>
> Stop after S4's gates and report. S5 (thermal wall function, TVD/van-Leer
> scalar convection, Boussinesq buoyancy) is a separate session; do not start
> it. Update this document's STATUS header and `validation/scalar/README.md`.

---

### Historical: the S3 prompt (consumed 2026-08-04)

> Implement increment **S3** of `docs/next_session_scalar.md` — the
> cell-centred immersed-boundary coefficients for the passive scalars — on
> branch `scalar` (read the doc in full first: its STATUS header records
> exactly what S0/S1/S2 landed, every deviation and every gate number; also
> read the "Active work" section of CLAUDE.md for the project conventions,
> and `validation/scalar/README.md` for the gate machinery you inherit).
>
> S0 + S1 + S2 are DONE and gated: scalars are extra variables of `blk%q`
> (`VAR_S0+is`; `SCR_S0+is` in `qs`/`oldrhs`), `src/modules/scalar.f90` owns
> `scalar_type` and the fused `scalar_transport` kernel (2nd-order central
> divergence-form convection, face-flux diffusion with `D_face = 1/(Re Pr) +
> ½(ν_t,L+ν_t,R)/Pr_t(face)` reading `turb%nut`, constant `Pr_t` or
> Kays–Crawford, momentum RK3), `scalar_finish` runs after
> `pressure_projection`, and BOTH `[scalar] count = 0` AND a scalar run with
> turbulence off are bit-exact CPU+GPU on their suites.
>
> **S3 — the immersed body.** Section 4 of the plan is the design and it was
> written from a complete read of `ibm.f90`; the short version:
>
> - `set_ibm_coeff` / `set_ibm_coeff_host` are ALREADY generic in `var` (they
>   read `blk%x(ix,var,b)` and the graded second-order stencil from the same
>   column) — the only thing blocking the pressure position is the guard
>   `if (var < VAR_U .or. var > VAR_W) error stop`. Allow `VAR_P`, and give
>   `ibm%coef` / `ibm%mu` the extent `VAR_U:VAR_P` **only when scalars are on**
>   so every existing case keeps its current memory and its current field
>   layout exactly.
> - Analytic IBM (inline): one extra `set_ibm_coeff(dns, blk, ibm, VAR_P)` in
>   `moby_solve.f90` and in the prepare path.
> - The stored coefficient carries the `1/Re` scaling, so the scalar uses
>   `coef_p/Pr` = its own `1/(Re Pr)` — one cell-centred coefficient array
>   serves ALL scalars, each with its own `Pr`.
> - Two body wall modes per scalar (`[scalar.N] ibm_wall`, already parsed into
>   `sc%ibmMode` by S0 but unused): `dirichlet` (default) penalises toward
>   `ibm_value` with the implicit factor `mu_s = 1/(1 + dt_gamma coef_p/Pr)`
>   computed inline (no `mu` array per scalar), diffusive fluxes NOT masked;
>   `adiabatic` masks BOTH the convective and the diffusive flux on each of
>   the six p-cell faces whose staggered velocity coefficient is solid — the
>   `solw/sole/...` test `rans.f90` already uses for `k`.
> - File-based IBM: a NEW, OPTIONAL case-file dataset `coef_p_blocks`, written
>   by `moby_prepare` when the prepare ini declares scalars (one extra
>   `set_ibm_coeff_host(..., VAR_P, inside)` + an append mirroring
>   `fdm_h5_case_append_coef`). Existing case files must stay byte-identical
>   and readable; the solver hard-errors with an explicit "re-run moby_prepare
>   with [scalar]" message when scalars are on and the dataset is absent.
>   `dwall_blocks`, the masks and `block_active` are untouched.
>
> Conventions that are not negotiable:
> - Build both paths (`./compile.sh cpu && ./compile.sh gpu`) with the module
>   **loaded** (`module load toolkits/nvhpc/25.9`) — see the environment
>   landmine in Section 13; always launch through `mpirun`, even on one rank.
> - **Save the S2 nofma binaries somewhere outside the tree BEFORE touching
>   any source**: they are the reference for every bit-exactness gate below.
> - `[scalar] count = 0` must stay bit-exact (`run_bitexact.sh`, 7-case suite,
>   max_abs 0, CPU AND GPU) and a scalar run WITHOUT an immersed body must be
>   bit-exact vs the S2 binaries (`run_bitexact_s1.sh`, same idea). Re-run
>   `run_gates.sh` and `run_gates_s2.sh` — every S1 and S2 gate must still pass
>   with the same numbers (none of them has a body).
> - New gate cases and drivers go in `validation/scalar/`, and the README's
>   results section grows an S3 block with what was run and measured.
> - Never declare an increment done with a failing build or an ungated result.
>
> S3 gates (all must pass and be recorded): solid-cell scalar == `ibm_value`
> EXACTLY in `dirichlet` mode; `adiabatic` mode conserves `∫s dV` to round-off
> with a body in the domain; a heated cylinder (Re 40) Nusselt number against
> literature AND against a Gauss/CV border-flux cross-check (the A2 method in
> `validation/cylinder/`); the prepared case file dataset-identical (`h5same`)
> to the previous one apart from the added `coef_p_blocks`; `refine_body`
> per-level coefficients consistent; 1 rank == 4 ranks EXACT and CPU == GPU
> EXACT on a scalar + IBM case.
>
> LANDMINES for this increment specifically: (1) `set_ibm_coeff_host` is the
> host twin of the device kernel and carries a KEEP IN LOCKSTEP comment —
> whatever you change in one, change in the other, and the prepare-vs-inline
> gates are what catch a divergence; (2) a body-integrated flux diagnostic
> (Nusselt from `∫coef_p (s − s_body) dV`) needs `[blocks] keep_buried` for
> exactly the reason the penalization FORCES do (a removed buried core absorbs
> loading outside the coefficient bookkeeping — CLAUDE.md A2); (3)
> penalization Dirichlet at the body is first-order/staircase like the
> velocity, and the graded coefficient only mitigates it — do not expect
> better than the velocity's near-body accuracy.
>
> Stop after S3's gates and report. S4 (scalar statistics + `compare_fields.py`
> dataset discovery) and S5 (thermal wall function, TVD/van-Leer scalar
> convection, Boussinesq buoyancy) are separate sessions; do not start them.
> Update this document's STATUS header and `validation/scalar/README.md` with
> what landed and what each gate measured.

---

## 13. Environment and workflow notes (from the S2 session)

- **`./compile.sh` without the nvhpc module loaded destroys the build cache.**
  It re-runs `cmake` and re-discovers MPI; with no module the wrappers resolve
  to `/usr/bin/mpif90` (gfortran module files), which nvfortran then rejects
  with `Corrupt or Old Module file .../gfortran-mod-15/openmpi/mpi_f08.mod`.
  Recovery is simply to `module load toolkits/nvhpc/25.9` and re-run
  `./compile.sh`, but the failure looks like a source error and is not one.
  `validation/scalar/compile_nofma.sh` reads the matching `build_cpu` /
  `build_gpu` cache, so a clobbered cache breaks the nofma builds too.
- **A new CMake target must be added to `SOLVER_TARGETS`**, not just given an
  `add_executable` — that list is what carries `-O2`, `-Mpreprocess`, the MPI
  and HDF5 include paths and the per-target module directory. A target left
  out of it fails with "Label field of continuation line is not blank" (the
  preprocessor never ran) and "cannot open source file hdf5.h".
- **Do not restart a channel case from a run's FINAL snapshot.** The niter = 6
  projection accumulates a large velocity-neutral mode in the stored `pn` (the
  pn-drift family, CLAUDE.md A2), and the final write catches it at its worst:
  `../channel_interface/les`'s `*_50001.h5` files carry `|pn| ~ 1e6` where every
  mid-run snapshot has `|pn| ~ 10`, and restarting on one blows the run up in a
  single step. Restart from a mid-run snapshot, or zero `pn` first.
- **`pkill -f <pattern>` matches the shell running it** when the pattern also
  appears in the command line, so it kills the caller before (or instead of)
  the target. Kill by PID from `ps -C moby_solve`, and beware that two runs
  writing the same `field_prefix` with different `field_interval`s produce
  files with identical `(step, t)` labels — one silently overwrites the other
  and the statistics silently mix.
- **The local GPU is shared.** Budget the wall clock: the 64x48x64 LES channel
  ran at 0.2-0.4 s/step with another job resident, i.e. ~7 min per 0.4 t.u.
  snapshot, and the thermal field of a scalar channel relaxes on
  `tau = L²/(pi² D_eff) ~ 6 t.u.` — see the METHOD note in the STATUS header
  for how `make_theta_ic.py` + `retarget_theta.py` cut that down.

---

### Historical: the S2 prompt (consumed 2026-08-03)

> Implement increment **S2** of `docs/next_session_scalar.md` — the turbulent
> closure of the passive scalars — on branch `scalar` (read the doc in full
> first: its STATUS header records exactly what S0/S1 landed, the three
> flagged deviations, and every gate number; also read the "Active work"
> section of CLAUDE.md for the project conventions, and
> `validation/scalar/README.md` for the gate machinery you inherit).
>
> S0 + S1 are DONE and gated: scalars are extra variables of `blk%q`
> (`VAR_S0+is`; `SCR_S0+is` in `qs`/`oldrhs`), `src/modules/scalar.f90` owns
> `scalar_type` and the fused `scalar_transport` kernel (2nd-order central
> divergence-form convection, molecular face-flux diffusion `D = 1/(Re·Pr)`,
> momentum RK3), `scalar_finish` runs after `pressure_projection`, and
> `[scalar] count = 0` is bit-exact CPU+GPU on the 7-case suite.
>
> **S2 — turbulent scalar diffusivity.** Add the eddy part to the face
> diffusivity in `scalar_transport`: `D_face = 1/(Re·Pr) + ½(ν_t,L + ν_t,R)/
> Pr_t(face)`, reading `turb%nut` — the ONE blended field LES (WALE /
> Smagorinsky), RANS (SST) and IDDES all write, so one code path covers all
> three models. `Pr_t` is either the per-scalar constant (`[scalar.N] prt`,
> default 0.85) or the **Kays–Crawford** correlation (`prt_model = kays`,
> already parsed into `sc%prtModel` by S0 but unused): implement it as a pure
> `!$omp declare target` function in scalar.f90 next to the `lm_*`
> correlations' style, with the `Pe_t → 0` and `Pe_t → ∞` limits guarded, and
> unit-test it against tabulated values the way `src/test_transition.f90`
> tests the LM correlations (add `src/test_scalar.f90` or extend that driver).
> `[rans] wall_treatment = wall_function` together with scalars must be a HARD
> CONFIG ERROR in this phase (a thermal wall function is S5). Also finish §8:
> `get_timestep_rates` must use `ν_eff = ire/Pr_min + ν_t/Pr_t,min` in the
> Peclet rate (S1 did the molecular half only, in `precompute_peclet_rate`).
>
> Conventions that are not negotiable:
> - Build both paths (`./compile.sh cpu && ./compile.sh gpu`, module
>   `toolkits/nvhpc/25.9`); always launch through `mpirun`, even on one rank.
> - **`[scalar] count = 0` must stay bit-exact** and a scalar run with
>   turbulence OFF must be bit-exact vs the S1 binaries: prove both with
>   `validation/scalar/compile_nofma.sh` builds and
>   `validation/scalar/run_bitexact.sh` (7-case suite, max_abs 0, CPU AND GPU).
>   Re-run `validation/scalar/run_gates.sh` — every S1 gate must still pass
>   with the same numbers (they all run turbulence-free).
> - New gate cases and drivers go in `validation/scalar/`, and the README's
>   results section grows a S2 block with what was run and measured.
> - Never declare an increment done with a failing build or an ungated result.
>
> S2 gates (all must pass and be recorded): the Kays–Crawford correlation
> unit test (tabulated values, every branch, both limits); a turbulent channel
> at Re_τ 180 with Pr 0.71 under LES — mean `θ⁺` against Kader / Kawamura DNS
> and the turbulent heat-flux profile; the same case under SST with resolved
> walls — log-law slope; a 2:1 wall-band-refined channel showing **NO spurious
> scalar band across the interface** (jump ratios, the u'/v' band lesson in
> `validation/channel_interface/README.md`); `wall_treatment = wall_function`
> + scalars rejected with a clear message; 1 rank == 4 ranks EXACT and
> CPU == GPU EXACT on a scalar+turbulence case.
>
> Stop after S2's gates and report. S3 (cell-centred IBM coefficients +
> `moby_prepare coef_p_blocks`) is a separate session; do not start it. Update
> this document's STATUS header and `validation/scalar/README.md` with what
> landed and what each gate measured.

### Historical: the S0 + S1 prompt (consumed 2026-08-03)

> Implement passive scalars in mobydiff, increments **S0** and **S1** of
> `docs/next_session_scalar.md` (read it first, in full — it is the design and
> it was written from a complete read of the affected code; also read the
> "Active work" section of CLAUDE.md for the project conventions). Branch
> `scalar`, already checked out and identical to `boundaryLayer` in `src/`.
>
> **S0 — layout, config, io, halos, BCs, no transport.** Scalars become extra
> variables of `blk%q` (`dns%nVar = NVAR + dns%nScalar`; scalar `is` at
> `VAR_S0+is` in `q`, `SCR_S0+is` in `qs`/`oldrhs`). New
> `src/modules/scalar.f90` with `scalar_type` owning ALLOCATABLE per-scalar
> arrays (no fixed bound anywhere — no `MAX_SCALARS`) plus
> `enter_/exit_scalar_data`; `[scalar]` / `[scalar.N]` config parsed straight
> into it; `comm%activeVars` allocatable sized `dns%nVar` (fix the `maxCount`
> buffer sizing at comm.f90:453-454 and the bound checks at 1629/1634);
> `fdm_h5_write_field` / `fdm_h5_read_field` take a variable count + name table
> (field_hdf5.c hardcodes `4` in three places: the write loop, `block_stride`,
> and the read path — and `read_force_file` must keep passing 4);
> `apply_scalar_bc_q` in boundary.f90 (a `var`-indexed twin of apply_bc's
> cell-centred branch, boundary.f90:697-719 — do NOT extend apply_bc's own
> `do var = VAR_U, VAR_P` loop, it runs inside the projection); patch-derived
> BC defaults in the `resolve_face_bcs` slot; initial condition + restart
> (absent dataset ⇒ warn + reinit, the RANS named-scalar precedent).
> The scalar is carried and written but not advanced.
>
> **S1 — the transport kernel**, molecular diffusivity only (no `nut`, no IBM
> coefficients yet): one fused kernel `scalar_transport` called from
> moby_solve.f90 right after `momentum(...)`, plus `scalar_finish` (qs→q
> copyback, ghosts, one batched `exchange_halos` over the scalar vars) after
> `pressure_projection` returns. Everything OUTSIDE the projection. 2nd-order
> central divergence-form convection using the p-cell's own face velocities,
> honouring `[flow] convection = skew` (`conv_div − s·(div u)|stencil`);
> face-flux diffusion on the scalar module's own `invD*` tables; the momentum
> RK3 structure verbatim. NO upwind option.
>
> Conventions that are not negotiable:
> - Build both paths (`./compile.sh cpu && ./compile.sh gpu`, module
>   `toolkits/nvhpc/25.9`); always launch through `mpirun`, even on one rank.
> - **`[scalar] count = 0` must be bit-exact by construction** — identical
>   allocation shapes, no kernel called, no extra exchange, apply_bc untouched.
>   Prove it with `-Mnofma` / `-gpu=nofma` builds and `tools/compare_fields.py`
>   (max_abs 0) against the pre-change binaries on the standard 7-case suite
>   (min_channel, les_ibm ± refine_body, Beltrami y-slab, turb180, wf180_y30,
>   lam30t), CPU AND GPU. Do this at the END of S0 and again at the end of S1.
> - New gate cases and drivers go in `validation/scalar/` with a README
>   recording what was run and measured.
> - Never declare an increment done with a failing build or an ungated result.
>
> S1 gates (all must pass and be recorded): uniform scalar `s ≡ const`
> preserved to 0.0 through a 3-level refined patch including edges/corners
> (the scalar analogue of the uniform-flow interface gate); `d/dt ∫s dV` =
> boundary flux to round-off in a periodic box; pure conduction between
> Dirichlet walls reproduces the exact linear profile, order 2 on a stretched
> grid; an advection–diffusion MMS (or Taylor–Green scalar) at order 2; a
> Pr sweep (0.1 / 1 / 10) in a laminar channel; 1 rank == 4 ranks EXACT;
> CPU == GPU EXACT.
>
> Stop after S1's gates and report. S2 (turbulent closure, `nut/Pr_t` +
> Kays–Crawford) and S3 (cell-centred IBM coefficients + `moby_prepare`
> `coef_p_blocks`) are separate sessions; do not start them. Update this
> document's STATUS header with what landed and what each gate measured.
