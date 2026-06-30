# Next session — code cleanup (remove the MOBY_* diagnostic hooks + dead gated paths)

Branch `claude/jacobi-interface`. The block-refinement + 2:1-interface + LES + IBM
work is validated and merged into this branch (see CLAUDE.md "Active work" and the
`validation/channel_interface/{les,les_ibm,core_patch,developed}` campaigns). What
remains is **debt removal**: the development/diagnostic scaffolding accumulated
while debugging the interface and projection is still in the hot files. This
session strips it so the production code is clean and readable (a stated CLAUDE.md
coding convention: "the code must stay very easy for humans to read").

This is a **pure refactor**: every standard case must stay **bit-exact** vs the
pre-cleanup binary (see Verification). The env hooks are no-ops when their variable
is unset, so removing them must not change any default run.

## The hooks (enumerate fresh: `grep -rn "MOBY_" src/`)

19 distinct `MOBY_*` env-var hooks, 63 references across `src/main.f90` (43),
`pressure_solver.f90` (11), `step.f90` (4), `comm.f90` (3), `init.f90` (1),
`flow/generic_flow.f90` (1). Two categories — handle them differently:

### A. Pure diagnostics — DELETE outright (hook + the buffers/branches they drive)
These only observe or dump; unset = inert. Remove the env read, the gated code, any
dedicated buffers (e.g. `divBuf`, RHS-capture slots), and the comment blocks.
- `MOBY_STEPDIV` — per-step divergence monitor (`print_step_divergence`, `divBuf`).
- `MOBY_DIVDUMP` — writes D·u before/after projection as companion field files.
- `MOBY_RHSDUMP` — dumps the discrete momentum RHS L_h(u) into velocity slots.
- `MOBY_TERMDUMP=<var>` — dumps each momentum term of a component separately.
- `MOBY_MANUF=<amp>` — adds amp·grad(phi) manufactured perturbation.
- `MOBY_KEBAL`, `MOBY_KESKEW` — kinetic-energy-balance / skew diagnostics.
- `MOBY_PHASETIME` — per-phase wall-time accumulation (the `les_timing:` style lines).
- `MOBY_HALO_AUDIT` — manufactured-linear-field halo-exchange audit (pre-loop).
- `MOBY_RESLOG` — projection residual logging.
- `MOBY_PROJONLY`, `MOBY_PREDONLY` — run only projection / only predictor (debug
  isolation). Removing these restores the unconditional full step.

### B. Algorithmic toggles — COLLAPSE to the production branch, do NOT just delete
Each gates a real code path with a *validated production winner* (now a config key
or the documented default). Replace the hook with its production branch hardwired,
and delete the losing branch. Confirm the production value against CLAUDE.md before
collapsing — get this wrong and you change physics, not just remove debt.
- `MOBY_CHEB`, `MOBY_CHEB_LMAX`, `MOBY_CHEB_LMIN` — Chebyshev-Jacobi acceleration.
  PRODUCTION = `[pressure] accel = chebyshev` (config `ps%cheb`); the env vars are
  dev overrides of lmax/lmin. Keep the config path; remove the env overrides (or
  fold sane defaults in). `pressure_solver.f90:75-81`.
- `MOBY_PHIINTERP` — phi-interpolation variant in the projection. Decide the
  production branch (the validated interface transfer) and hardwire.
- `MOBY_VELINJECT` — velocity prolong-injection toggle in `comm.f90:179`. PRODUCTION
  = const-1/2 inject (the `interface_constant_half` default). Hardwire.
- `MOBY_IFFILT=<alpha>` — coarse-interface-band tangential filter coefficient.
- `MOBY_NORECON` — disables the deep-halo cubic reconstruction. PRODUCTION = the
  const-1/2 default already skips the cubic (see CLAUDE.md "constant-1/2 interface
  is the DEFAULT"), so reconcile this hook with that and remove the dead branch.
- `MOBY_PHIINTERP`/`MOBY_IFFILT`/`MOBY_NORECON`/`MOBY_VELINJECT` are interface-
  transfer levers from the 2:1 debugging; cross-check each against the RESOLVED
  decisions in CLAUDE.md (`interface_constant_half=true`, `momentum_reflux=false`)
  and the `interface-validation-suite` / `interface-normal-treatment` memories so
  the surviving branch is exactly the validated one.

Order: do category A first (mechanical, low-risk, big readability win), verify
bit-exact, commit. Then category B one hook at a time, each followed by a bit-exact
check, so a physics regression is bisectable to a single hook.

## Watch for
- `main.f90` parks diagnostic buffers in the pressure/velocity SLOTS for output
  (DIVDUMP/RHSDUMP companions) — when you remove the hook, also remove the slot
  parking and any `write_field` companion calls, or you leave dangling writes.
- Removing `MOBY_STEPDIV` deletes `print_step_divergence` + `capture_divergence` +
  `divBuf` allocation/mapping — check the GPU enter/exit data maps too.
- Don't break the `validation/.../run_*.py` drivers that set `MOBY_STEPDIV=1` etc.
  for monitoring (grep `validation/ tools/` for `MOBY_`); update or drop those uses.
- Keep `compile.sh cpu && compile.sh gpu` green after every commit.

## Verification (this is a pure refactor — bit-exact is the gate)
Build BOTH sides `-Mnofma` (CPU) / `-Mnofma -gpu=nofma` (GPU); compare with
`tools/compare_fields.py` (CLAUDE.md Verification):
- `tutorials/channel_kmm180` (or a short channel), `tutorials/wavychannel`
  (analytic IBM), `tutorials/sailplane` (file IBM) — bit-exact vs pre-cleanup.
- channel `nb=4` (many blocks) bit-exact vs pre-cleanup.
- 2:1 interface: a refined channel + the Beltrami interface regression
  ([[beltrami-interface-regression]] memory) — bit-exact.
- LES+IBM: `validation/channel_interface/les_ibm` short run still stable + gates
  hold; `validation/channel_interface/les` short run unchanged.
- Never declare done with a failing build or a non-bit-exact case.

## NEXT-SESSION PROMPT

> Read `docs/next_session_cleanup.md` and CLAUDE.md. Branch `claude/jacobi-interface`.
> The block-refinement / 2:1-interface / LES / IBM work is validated and done; now
> remove the development debt. Strip the 19 `MOBY_*` diagnostic/dev hooks (63 refs
> across main.f90, pressure_solver.f90, step.f90, comm.f90, init.f90,
> generic_flow.f90) and the dead gated paths. Category A (pure diagnostics:
> STEPDIV/DIVDUMP/RHSDUMP/TERMDUMP/MANUF/KEBAL/KESKEW/PHASETIME/HALO_AUDIT/RESLOG/
> PROJONLY/PREDONLY) — delete outright incl. their buffers and slot-parking. Category
> B (algorithmic toggles: CHEB*/PHIINTERP/VELINJECT/IFFILT/NORECON) — collapse each
> to its VALIDATED production branch (config `accel=chebyshev`,
> `interface_constant_half=true`, `momentum_reflux=false`; cross-check the
> interface-validation-suite / interface-normal-treatment memories) and delete the
> losing branch; do NOT just remove the hook. This is a PURE REFACTOR: every standard
> case (channel, wavychannel, sailplane, channel nb=4, a refined-interface channel,
> the Beltrami interface regression, the LES+IBM short run) must stay BIT-EXACT vs
> the pre-cleanup binary built `-Mnofma`/`-gpu=nofma` (`tools/compare_fields.py`). Do
> category A first (verify+commit), then category B one hook at a time (verify after
> each so a regression bisects to one hook). Keep both CPU and GPU builds green. Make
> a plan first; execute after.
