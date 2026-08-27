# `jacobi_apply` — the mu planes, A6000, 2026-08-27

Same-day A/B of "do not read `ibm%mu` when there is no immersed body" against
the HEAD binary (`2e74cc4`), both default flags, both `PROFILE=1`, 400 cold
steps, 138.4 M cells, istmcetus GPU 1.

| layout | reference | candidate | gain |
|---|---|---|---|
| `nb` unset (control) | 1.17811 | 1.02871 | **12.7 %** |
| **`nb = 64 44 48`** (production) | 1.21922 | **1.09031** | **10.6 %** |

Both runs are drift-clean: over the last 100 steps the cumulative average rises
by 0.00014 (reference) and 0.00018 (candidate) s/step — the same slope, so the
machine state is comparable and the delta is the change.

The runtime diagnostic columns are identical between the two binaries
(`L2_div 1.84481025E-05`, `Linf 1.00341584E+00` at step 400) on both layouts.

## Where the time went (`rect_jacobi`, s/step)

| step phase | ref | cand | delta |
|---|---|---|---|
| momentum | 0.3503 | 0.3500 | −0.1 % |
| **ibm_mu** | 0.0401 | 0.0000 | **−100 %** |
| bodyforce | 0.0595 | 0.0594 | −0.2 % |
| apply_bc | 0.0014 | 0.0014 | — |
| vel_exchange | 0.0104 | 0.0104 | — |
| **projection** | 0.7461 | 0.6578 | **−11.8 %** |
| io_stats | 0.0114 | 0.0113 | — |
| **total** | **1.2192** | **1.0903** | **−10.6 %** |

| projection phase | ref | cand | delta |
|---|---|---|---|
| sweep (`jacobi_compute_phi`) | 0.2813 | 0.2663 | **−5.3 %** |
| **apply (`jacobi_apply`)** | 0.4074 | 0.3339 | **−18.1 %** |
| phi_exchange | 0.0228 | 0.0228 | — |
| vel_exchange | 0.0256 | 0.0257 | — |
| apply_bc | 0.0074 | 0.0075 | — |
| setup | 0.0015 | 0.0016 | — |

Nothing regressed in absolute terms, on either layout.

## The measurement that chose the fix

`ncu` on the reduced case (`~/.moby_prof/rect_small.ini`: production ini with
nx/8, nz/4, same 64×44×48 block shape; ncu locks the SM clock to 1.40 GHz, so
absolute durations are ~30 % above production and compute-bound kernels are
penalised relative to memory-bound ones — read the *ratios*).

| kernel | duration | DRAM | B/cell | regs | occ |
|---|---|---|---|---|---|
| `apply` k1 (`p += phi/dt`) | 160.2 µs | **91.1 %** | 24.6 | 59 | 57.9 % |
| `apply` k2 (face corrections) | 585.3 µs | **84.0 %** | 82.9 | 86 | 38.7 % |
| `compute_phi` | 762.6 µs | 45.7 % | 58.8 | 128 | 31.5 % |

Both `apply` kernels sit at 84–91 % of peak DRAM throughput, so **hypothesis 2
of the handout (register/occupancy pressure) is dead for this kernel** — 38.7 %
occupancy still reaches 84 % of peak. The only lever is bytes. The measured
bytes per cell match the hand model to 3 %:

| item | B/cell | share of `apply` |
|---|---|---|
| `q` u/v/w read+write | 48 | 45 % |
| **`ibm%mu` (3 planes)** | **24** | **22 %** |
| `q` p read+write | 16 | 15 % |
| `phi` read (once per launch) | 16 | 15 % |

`ibm%mu` is the largest removable item, and in every non-IBM run it is
identically 1: `dns%ibm_enabled = .false.` makes both `set_ibm_coeff` branches
zero `coef` and `cycle`, `read_ibm_coeff_file` is guarded by the same flag, and
`mu = 1/(1+dt·0) = 1`. Multiplying by exactly 1.0 is the identity in IEEE, so
skipping the loads is bit-exact by construction, not by tolerance.

## After the change (same `ncu` case)

| kernel | duration | DRAM | B/cell | regs |
|---|---|---|---|---|
| `apply` k1 | 160.0 µs (—) | 91.5 % | 24.6 | 59 |
| `apply` k2 | **464.2 µs (−20.7 %)** | 74.2 % | **58.0** (−24.9) | 86 |
| `compute_phi` | 746.0 µs (−2.3 %) | 26.7 % | 34.8 (−24.0) | 122 |

Exactly the predicted 24 B/cell came off each kernel. `apply` k2 converted it
into time almost one-for-one because it was DRAM-bound; `compute_phi` did not,
because at 45.7 % DRAM and 82 % SM throughput it is **compute/latency-bound** —
its 41 % traffic cut bought 2–13 %. That asymmetry is the whole reason the
handout insisted on profiling before choosing.

## Consequence for the next increment

`apply` k2 has fallen off the roofline (84 % → 74 % DRAM). It now has headroom
again, so the remaining `apply` levers are the ones that were pointless before:

- **Fuse the two launches.** `phi` is read once per launch over the whole
  interior; one fused kernel drops 8 of the remaining ~83 B/cell. The two
  launches write disjoint outputs (`VAR_P` vs `VAR_U/V/W`) and no cell reads
  another cell's `q` output, so the fusion is a pure scheduling change.
- **Split the rare high-face work into its own kernel.** k2 carries 86
  registers and 38.7 % occupancy for six `face_grad_corr` calls, three of which
  fire only at `i == nx` / `j == ny` / `k == nz` on interface or outlet faces.

`compute_phi` (0.266 s/step, 24 % of the step) is now the better target of the
two, and it is a *different* problem: 122 registers, 31.5 % occupancy, 82 % SM
throughput. Instruction count and register pressure, not traffic.

## Caveat

The gain applies to runs without an immersed body — every channel, boundary
layer, Beltrami and turbulence-model case, including the production boundary
layer. IBM runs (naca, sd7003, sailplane, les_ibm) take the same path as before
and are bit-identical; they are unchanged, not slower.
