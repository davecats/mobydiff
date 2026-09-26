# Validation & verification

`mobydiff` is verified at two levels: **correctness** against exact solutions and reference
data, and **regression** (bit-exact reproduction across refactors). The `validation/`
directory holds the reference cases and drivers; `tools/` holds the checkers.

**[`validation/README.md`](../validation/README.md) is the index of every case** — what each
exercises and which driver runs it; each case directory has its own README with commands
and recorded results.

## Reference flows (examples)

| Case | Directory | Exact solution / reference | Checker |
|------|-----------|----------------------------|---------|
| Beltrami / ABC flow | `validation/beltrami/` | Analytic 3D Beltrami (2π-periodic cube), uniform and through 2:1 slabs/patches | `tools/check_beltrami.py` |
| Inflow/outflow Poiseuille, free stream | `validation/freestream/` | Periodic reference, exact oblique free stream | `run_gates.sh` |
| Turbulent channel + 2:1 interface | `validation/channel_interface/` | Uniform-resolution DNS reference | channel stats tools |
| k-ω SST / transition / wall functions | `validation/rans_sst/` | DNS channel profiles, laminar parabola | `run_gates.sh`, `check_gates.sh` |
| Passive scalars, conjugate heat transfer | `validation/scalar/`, `validation/conjugate/` | Closed forms, MMS, energy budgets | per-directory drivers |

The Beltrami case exercises the core discretization and time integration against a
closed-form solution and is the fastest sanity check (`[flow] initial = beltrami`). The
`channel_interface` suite is the turbulence validation of the block refinement and the 2:1
interface (flat interfaces, edges/corners, with LES and with IBM walls), each compared
against a uniform-resolution reference run.

## Refinement-specific checks

Two properties are checked to round-off for the block-refinement and interface machinery:

- **Uniform-flow preservation** — a spatially uniform velocity field advected through a
  refined patch must be preserved exactly (max deviation 0.0), since every consistent
  interface transfer is exact for a constant field. This isolates transfer bugs from
  physics.
- **Global mass conservation** — the total divergence residual must stay at round-off
  (≈ 1e-20 relative to the velocity scale) with a refinement patch present.
- **Interface stability** — white noise on a refined patch must decay everywhere
  (`validation/interface_decay/`).

## Bit-exact regression for refactors

A change advertised as a **pure refactor** must reproduce the pre-refactor output bit-for-bit.
Because default FMA contraction introduces 1–2 ulp differences for arithmetically identical
source, both the reference and the candidate are built with FMA disabled:

- CPU: `-Mnofma` (`./compile.sh cpu_nofma`)
- GPU: `-Mnofma -gpu=nofma` (`./compile.sh gpu_nofma`)

and compared with `tools/compare_fields.py` or `tools/h5maxdiff` (max_abs 0) on the
velocity components, pressure and every model field (`nut`, RANS scalars, passive scalars
named explicitly). The comparison is run on the CPU **and** the GPU path, so the two
backends are also checked against each other.

The standard regression suite (listed in `validation/README.md`) covers:

- a minimal channel with blocks + a 2:1 interface + Chebyshev acceleration (1 and 4 ranks),
- the Beltrami y-slab interface case,
- a channel with file-based IBM and WALE LES (with and without `refine_body`),
- k-ω SST with resolved walls, wall functions and transition, and
- passive-scalar transport, with and without an immersed body.

## Rules of thumb

- Never declare work done with failing builds or unverified results.
- Build the CPU path as well as the GPU path — the CPU build is the debugging reference.
  Host/device staleness bugs are invisible on the CPU: a GPU-only suite run is required.
- For a physics change, check the relevant reference flow; for a refactor, check bit-exactness.
