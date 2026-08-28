# `jacobi_apply` — what limits it, A6000, 2026-08-27

The `ncu` characterisation of the projection kernels, plus one optimisation that
was implemented, gated, measured and then **REVERTED by decision** — see the
last section. The characterisation is the part that survives; it is what should
drive the next increment.

## What limits each kernel

`ncu` on the reduced case (`~/.moby_prof/rect_small.ini`: the production ini with
nx/8, nz/4, same 64×44×48 block shape; ncu locks the SM clock to 1.40 GHz, so
absolute durations run ~30 % above production and compute-bound kernels are
penalised relative to memory-bound ones — read the *ratios*).

| kernel | duration | DRAM | B/cell | regs | occupancy |
|---|---|---|---|---|---|
| `apply` k1 (`p += phi/dt`) | 160.2 µs | **91.1 %** | 24.6 | 59 | 57.9 % |
| `apply` k2 (face corrections) | 585.3 µs | **84.0 %** | 82.9 | 86 | 38.7 % |
| `compute_phi` | 762.6 µs | 45.7 % | 58.8 | 128 | 31.5 % |

**`jacobi_apply` is DRAM-bound.** Both launches sit at 84–91 % of peak, and k2
gets there at 38.7 % occupancy — so the handout's hypothesis 2
(register/occupancy pressure) is dead for this kernel, and the only lever is
bytes.

**`jacobi_compute_phi` is not.** 45.7 % DRAM against 82 % SM throughput, 128
registers, 31.5 % occupancy: it is compute/latency-bound. Do not carry a
conclusion from one kernel to the other — the reverted experiment below cut the
same 24 B/cell from both and they responded completely differently.

## The byte budget of `apply`

Measured bytes per cell match a hand model to 3 % (24.6 vs 24, 82.9 vs 80),
which makes the budget actionable rather than indicative:

| item | B/cell | share of `apply` |
|---|---|---|
| `q` u/v/w read+write | 48 | 45 % |
| `ibm%mu` (3 planes) | 24 | 22 % |
| `q` p read+write | 16 | 15 % |
| `phi` read (once per launch) | 16 | 15 % |

## The reverted experiment: skipping the mu planes

`ibm%mu` is identically 1 whenever there is no immersed body — `ibm_enabled =
.false.` makes both `set_ibm_coeff` branches zero `coef` and `cycle`,
`read_ibm_coeff_file` is guarded by the same flag, so `mu = 1/(1+dt·0) = 1`. A
warp-uniform `useIbm` scalar in `jacobi_compute_phi` / `jacobi_apply` skips the
loads, and `update_ibm_mu` can be skipped outright. Multiplying by exactly 1.0
is the IEEE identity, so it was bit-exact by construction.

Same-day A/B against the HEAD binary (`2e74cc4`), both `PROFILE=1`, 400 cold
steps, 138.4 M cells, istmcetus GPU 1. Both runs drift-clean (cumulative average
rising 0.00014 / 0.00018 s/step over the last 100 steps — the same slope) and
the runtime diagnostic columns identical between binaries.

| layout | reference | candidate | gain |
|---|---|---|---|
| `nb` unset (control) | 1.17811 | 1.02871 | 12.7 % |
| `nb = 64 44 48` (production) | 1.21922 | 1.09031 | **10.6 %** |

| projection phase | ref | cand | delta |
|---|---|---|---|
| sweep (`jacobi_compute_phi`) | 0.2813 | 0.2663 | **−5.3 %** |
| apply (`jacobi_apply`) | 0.4074 | 0.3339 | **−18.1 %** |
| `ibm_mu` (step phase) | 0.0401 | 0.0000 | −100 % |

Gates all passed: 7-case suite `max_abs 0` on CPU **and** GPU (the four channel
cases and Beltrami exercised the new branch, `les_ibm` ± refine the unchanged
IBM path), `validation/block_nb/run_gates.sh` CPU + GPU, 1 rank == 4 ranks.

**The `ncu` prediction held on both counts**, and that is the transferable
result: `apply` k2 was DRAM-bound and converted the 24 B/cell into time almost
one-for-one (−20.7 % in the ncu case, B/cell 82.9 → 58.0); `compute_phi` took
the *same* cut (58.8 → 34.8 B/cell) and gained 2–13 %, because traffic was never
its limiter.

**Reverted (2026-08-28), by decision, not because of a gate.** The 2:1 interface
is mostly going to be used together with the IBM — near-body refinement is the
whole point of it — so the body-free path is not the path the code will spend
its time in. A second branch through the two hottest projection kernels is not
worth carrying for it. Recover the change from history if a body-free campaign
ever justifies it.

## Consequence for the next increment

Because the mu cut is gone, `apply` k2 is back at 84 % of peak DRAM and its
byte budget is again the one tabulated above. The levers that remain, in the
order the profile supports:

- **`jacobi_compute_phi` is the better target** (0.281 s/step, 23 % of the
  step), and it is a *different* problem from `apply`: 128 registers, 31.5 %
  occupancy, 82 % SM throughput at 46 % DRAM. Instruction count and register
  pressure, not traffic — the six `face_grad_denom` calls and the division are
  where to look, not the array accesses.
- **Fuse the two `apply` launches.** `phi` is read once per launch over the
  whole interior; fusing drops 8 of 107 B/cell (~7 % of `apply`, ~2.5 % of the
  step). The launches write disjoint outputs (`VAR_P` vs `VAR_U/V/W`) and no
  cell reads another cell's `q` output, so it is a pure scheduling change.
- **Split the rare high-face work out of `apply` k2.** It carries 86 registers
  and 38.7 % occupancy for six `face_grad_corr` calls, three of which fire only
  at `i == nx` / `j == ny` / `k == nz` on interface or outlet faces. Worth
  little while k2 is pinned at 84 % of peak, but it is the natural companion to
  the fusion, which would otherwise raise the fused kernel's register count.
