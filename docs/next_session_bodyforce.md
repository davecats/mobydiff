# Next session — optional spatially-varying volumetric (body) force

Branch `claude/jacobi-interface`. Goal: add a **config-gated, spatially-varying
volumetric force** `f(x)` to the momentum equation, on top of the existing
constant per-direction forcing. Off by default → **bit-exact** with the current
solver. This is a feature addition (not a refactor); the gate is: disabled ==
byte-identical, enabled == physically correct on a manufactured check.

## What already exists (do NOT remove or duplicate)

`[flow] forcing_x/y/z` (`dns%forcing(1:3)`, default 0) is a **constant** body
force per direction, added in the momentum predictor RHS:

```
src/modules/step.f90  (subroutine momentum), for each component e.g. u:
    rhsu = ( -0.25*(advection) + forcing(VAR_U) + ire*(diffusion) )
    blk%qs(...) = blk%q(...) + dt_alpha*rhsu + dt_beta*oldrhs - dt_gamma*dpx
    blk%qs(...) = blk%qs(...) * ibm%mu(...)          ! IBM mask (solid -> 0)
```

The channel flow case wires `forcing_x = 1` as the mean pressure gradient. Keep
`forcing_*` exactly as is. The new term is **additive and independent**: the RHS
becomes `... + forcing(VAR_U) + bforce_u + ...`.

## Design (decided with the user)

- **Own module + own datatype.** New file `src/modules/bodyforce.f90` with a
  `bodyforce_type` deriving-type that OWNS its flat contiguous array (per the
  CLAUDE.md convention: derived types own flat arrays mapped to the device once
  in `enter_*_data`/`exit_*_data`). Do NOT bolt the array onto `block_set_type`.
  ```
  type :: bodyforce_type
      logical              :: enabled = .false.
      integer              :: source  = SRC_NONE   ! profile | file | custom
      real(C_DOUBLE), allocatable :: f(:,:,:,:,:)  ! (1:nb,1:nb,1:nb,NVEL,nBlocks)
      ! profile params / file path as needed
  end type
  ```
  Shape mirrors `blk%oldrhs` (interior 1:nb, per-component `NVEL`, trailing block
  index) — the force is read only at interior predicted faces (`i>=uStartX`
  etc.), so no halos are needed. Each component lives at ITS staggered location
  (`f(:,:,:,VAR_U,:)` at the u-face, etc.), matching how `forcing`/`q` are used.
- **Three sources** (`[force] type =`):
  - `profile` — a built-in analytic form filled once at init at each component's
    staggered coordinate (`blk%x(i,VAR_U,b), blk%y(j,VAR_U,b), blk%z(k,VAR_U,b)`
    for u, and the V/W-staggered lines for v/w). Fortran has no expression
    parser, so expose a small set of NAMED profiles with parameters (start with
    what you need; e.g. `sine`: `amp`, wavenumber `k`, direction), NOT a general
    string evaluator. Document that anything beyond the named set uses `custom`.
  - `file` — read `fx/fy/fz` from an HDF5 field into `f` at init, reusing the
    block-table read path in `io.f90` (`read_field` machinery; the force file is
    laid out like a velocity field, un/vn/wn = fx/fy/fz at the staggered points).
  - `custom` — the user fills `f` **inside the time loop** each step/substage
    (time-dependent forcing, controllers, actuators). See the API below.
- **Config `[force]` section** in `config.f90` (`apply_config_value`, add a
  `case ("force")` block) + fields on `dns` (or better, parse straight into the
  `bodyforce_type` via a small `apply_force_value`, mirroring how LES/MPI values
  are dispatched). Keys: `enabled` (bool), `type` (profile|file|custom),
  per-component profile params, `file` (path for type=file). Default disabled.

## Where it plugs in

1. `main.f90`: declare `type(bodyforce_type) :: bf`; after the block set + grid
   exist, `call init_bodyforce(bf, dns, blk, g)` (allocates `f`, fills it for
   `profile`/`file`); `call enter_bodyforce_data(bf)` (device map); in the RK
   substage loop, for `custom`, `call update_bodyforce(bf, blk, dns, g,
   dns%t_current)` BEFORE `momentum(...)`; `call exit_bodyforce_data(bf)` in the
   teardown block with the other `exit_*_data` calls.
2. `momentum(...)` (step.f90): add an argument for the force (pass `bf` or, to
   keep the kernel lean, `bf%f` + `bf%enabled`). Add `+ fbu` to `rhsu` (and
   `fbv`/`fbw` to `rhsv`/`rhsw`), read at `(i,j,k,VAR_U,b)` etc. Add `bf%f` to
   the kernel `map(to: ...)`. The `*ibm%mu` mask already zeroes the force in
   solid cells — no extra IBM masking needed (note this; it is the desired
   behaviour: no body force inside the immersed body).

### The bit-exactness knob (get this right)

When `enabled=.false.` the momentum RHS must be **byte-identical** to today. Two
viable routes:
- **Optional argument.** `momentum(..., bforce)` with `bforce` present only when
  enabled; guard the `+ fb*` term with `if (present(bforce))`. Cleanest for
  bit-exactness but the term sits inside one big fused `!$omp target` loop — a
  compile-time `present()` branch there is fine (the compiler specializes), but
  verify.
- **Always-add a zeroed array.** Simpler, but `rhsu = (... + forcing + 0.0 +
  ...)` must be shown bit-exact vs `(... + forcing + ...)`. With `-Mnofma` and no
  reassociation `a + 0.0 == a` exactly, so this usually holds — but it is a
  DECISION to verify empirically on the full suite, not assume.

Recommend the optional-argument route (provably no-op when off); fall back to
always-add only if the two-kernel-variant cost is unacceptable, and then prove
bit-exact.

## The `custom` API to expose (from bodyforce.f90)

Let the user write the force each step without editing the solver core:
- `bf%f` is **public** (the user writes into it), plus the block geometry they
  need is already on `blk` (`blk%x/y/z`, `blk%nb`, `blk%nBlocks`, `blk%origin`,
  `blk%level`).
- `bodyforce_zero(bf)` — device-side zero of `f`.
- `bodyforce_update_to_device(bf)` / `_from_device(bf)` — `!$omp target update`
  wrappers so a host-side fill reaches the GPU (and vice-versa).
- A single **user hook** `update_bodyforce(bf, blk, dns, g, t)` in bodyforce.f90
  with a clearly-marked "USER CODE HERE" body (the research-code idiom): for
  `custom`, the user edits this one routine to fill `bf%f` (host or a `!$omp
  target` kernel) as a function of `blk%x/y/z`, `blk%q` (the current field, for
  controllers), and `t`. For `profile`/`file` it is a no-op (already filled).
  Document that a controller (e.g. hold bulk velocity) reduces `blk%q` over the
  domain (use `comm_allreduce_sum`) then writes a uniform `f`.

## Verification

- **Disabled (the hard gate):** bit-exact vs pre-feature `-Mnofma`/`-gpu=nofma`
  on the full suite (`tools/compare_fields.py`, max_abs 0): min_channel, the
  les_ibm channel + refine_body, the Beltrami y-slab. Reuse the harness from the
  cleanup session (scratchpad `run_suite.sh` / `compare_suite.sh` pattern, or
  rebuild it: 4 short cases, CPU + GPU).
- **Enabled — manufactured balance:** pick `f` that sustains a known steady
  solution and check the field stays on it. Simplest: on a periodic box with
  `initial = uniform` and re=const, a uniform `f` reproduces the existing
  `forcing_*` behaviour — assert `type=profile` constant == `forcing_*` to the
  same value gives the SAME trajectory (cross-check against the const path).
  Then a spatially-varying `f = -nu*lap(u_target) + grad(...)` that holds a
  chosen `u_target` (e.g. a cosine profile) steady to truncation order.
- **`custom`:** a per-step controller that holds bulk velocity constant on a
  channel — assert bulk `U` converges to target and stays.
- **GPU parity:** enabled case CPU==GPU to round-off (the mask branch and the
  `f` map).
- Build BOTH `compile.sh cpu` and `compile.sh gpu` green after each step.

## Watch for

- `field_type` is GONE — the solver state is `block_set_type`; put the array in
  `bodyforce_type`, not a resurrected field type.
- Staggering: fill/read each component at ITS face coordinate, not cell centres.
- IBM: the `*ibm%mu` mask zeroes the force in solid — intended; don't re-mask.
- `custom` device coherence: if the user fills `f` on the host, they MUST
  `bodyforce_update_to_device` before `momentum`; provide and document it.
- Do not touch the removed interface config (const-1/2 is unconditional now) or
  reintroduce `momentum_reflux` (see CLAUDE.md "Production-config lockdown").

## NEXT-SESSION PROMPT

> Read `docs/next_session_bodyforce.md` and CLAUDE.md. Branch
> `claude/jacobi-interface`. Add a **config-gated, spatially-varying volumetric
> body force** `f(x)` to the momentum predictor, KEEPING the existing constant
> `[flow] forcing_x/y/z`. New module `src/modules/bodyforce.f90` owning a
> `bodyforce_type` (its own flat `f(1:nb,...,NVEL,nBlocks)` array, device-mapped
> via `enter_/exit_bodyforce_data`, per the CLAUDE.md GPU-data convention). New
> `[force]` config: `enabled` (default false), `type = profile|file|custom`.
> `profile` = built-in named analytic forms filled at init at each component's
> staggered coords; `file` = read fx/fy/fz from an HDF5 field via the io
> block-table path; `custom` = expose `bf%f` (public) + `bodyforce_zero`,
> `bodyforce_update_to_device/_from_device`, and a single user hook
> `update_bodyforce(bf, blk, dns, g, t)` the user edits to fill the force inside
> the RK loop each step. Apply as `rhsu += bf%f(i,j,k,VAR_U,b)` (and v/w) in
> `momentum`, before the `*ibm%mu` mask (which correctly zeroes force in solid).
> DISABLED must be **bit-exact** vs the pre-feature `-Mnofma`/`-gpu=nofma` binary
> on the full suite (min_channel, les_ibm channel + refine_body, Beltrami y-slab;
> `tools/compare_fields.py`, max_abs 0) — use the optional-argument route so the
> term is a provable no-op when off, or prove `+0.0` bit-exactness if you go
> always-add. ENABLED: verify a manufactured `f` holds a known solution steady,
> a constant `f` matches the `forcing_*` path, and CPU==GPU. Keep both builds
> green. Make a plan first; execute after.
