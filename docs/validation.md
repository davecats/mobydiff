# Validation & verification

`mobydiff` is verified at two levels: **correctness** against exact solutions and reference
data, and **regression** (bit-exact reproduction across refactors). The `validation/`
directory holds the reference cases and drivers; `tools/` holds the checkers.

## Reference flows

| Case | Directory | Exact solution / reference | Checker |
|------|-----------|----------------------------|---------|
| Taylor–Green vortex | `validation/taylor_green/` | Analytic 2D TGV (2π-periodic) | `tools/check_tgv.py` |
| Beltrami / ABC flow | `validation/beltrami/` | Analytic 3D Beltrami (2π-periodic cube) | `tools/check_beltrami.py` |
| Poiseuille channel | `validation/poiseuille/` | Laminar parabolic profile | `tools/check_parabolic_channel.py` |
| Turbulent channel + 2:1 interface | `validation/channel_interface/` | Uniform-resolution DNS reference | channel stats tools |

The Taylor–Green and Beltrami cases exercise the core discretization and time integration
against closed-form solutions and are the fastest sanity checks. The Poiseuille case checks
the wall-bounded forced channel against the discrete steady parabola. The
`channel_interface` suite is the turbulence validation of the block refinement and the 2:1
interface (flat interfaces, edges/corners, and with LES), each compared against a
uniform-resolution reference run.

## Refinement-specific checks

Two properties are checked to round-off for the block-refinement and interface machinery:

- **Uniform-flow preservation** — a spatially uniform velocity field advected through a
  refined patch must be preserved exactly (max deviation 0.0), since every consistent
  interface transfer is exact for a constant field. This isolates transfer bugs from
  physics.
- **Global mass conservation** — the total divergence residual must stay at round-off
  (≈ 1e-20 relative to the velocity scale) with a refinement patch present.

## Bit-exact regression for refactors

A change advertised as a **pure refactor** must reproduce the pre-refactor output bit-for-bit.
Because default FMA contraction introduces 1–2 ulp differences for arithmetically identical
source, both the reference and the candidate are built with FMA disabled:

- CPU: `-Mnofma`
- GPU: `-Mnofma -gpu=nofma`

and compared with `tools/compare_fields.py` on the velocity components and pressure
(`un vn wn pn`). The comparison is run on the CPU **and** the GPU path, so the two backends
are also checked against each other.

A representative regression set covers:

- a minimal channel with blocks + a 2:1 interface + Chebyshev acceleration,
- a channel with file-based IBM and WALE LES (with and without a 2:1 interface),
- geometry-driven `refine_body`, and
- the Beltrami y-slab interface regression.

## Rules of thumb

- Never declare work done with failing builds or unverified results.
- Build the CPU path as well as the GPU path — the CPU build is the debugging reference.
- For a physics change, check the relevant reference flow; for a refactor, check bit-exactness.
