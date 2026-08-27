# Per-direction `[blocks] nb` — Phase 1 gates

`[blocks] nb = 64 44 48` sets the block edge per direction; `nb = 16` still
broadcasts. The halo tax is `(1+2/nb_x)(1+2/nb_y)(1+2/nb_z)`, so a direction
that buys no refinement — y under `refine_dims = xz` — should carry the largest
`nb` the interface placement tolerates.

```bash
./run_gates.sh cpu
./run_gates.sh gpu
```

## What is gated, and why these gates

**Single level: nb-independence.** The block decomposition is bookkeeping, so
tiling the same grid differently must give bit-identical fields. That makes
non-cubic `nb` *self-gating*: no reference implementation is needed, the cubic
run **is** the reference. Gated at `--tolerance 0`:

| gate | result |
|---|---|
| `nb = 32 16 8` == `nb = 8` | PASS (max_abs 0) |
| `nb = 16 8 4` == `nb = 8` | PASS (max_abs 0) |
| `nb = 8` == `nb` unset (one block per rank box) | PASS (max_abs 0) |
| `nb = 32 16 8` on 1 rank == 4 ranks | PASS (max_abs 0) |

**2:1 interfaces: uniform-flow preservation**, not a cross-layout comparison —
see the finding below. A constant field is preserved exactly by any
*consistent* set of transfer operators, so a per-direction indexing mistake
anywhere in the entry generation, the gather maps or the ghost blend shows up
immediately as a nonzero deviation. `uniform_rect.ini` (quadtree `xz`) and
`uniform_rect_xyz.ini` (octree) run uniform oblique flow `(0.9397, 0.3420, 0.2)`
through a 3-level patch at `nb = 8 4 8`, 50 steps:

| case | levels | max deviation | pn spread |
|---|---|---|---|
| `uniform_rect` (xz) | 448 / 248 / 32 | **0.0** | **0.0** |
| `uniform_rect_xyz` (xyz) | 448 / 504 / 64 | **0.0** | **0.0** |

The tiling is deliberately fine enough that levels 0, 1 **and** 2 all survive,
so both l0–l1 and l1–l2 interfaces are exercised. A coarser `nb` lets the 2:1
smoothing swallow level 0 entirely (at `nb = 16 8 16`: 280 leaves, all refined)
and the gate silently loses half its coverage — `run_gates.sh` therefore checks
the level histogram, not just the deviation.

## FINDING: the 2:1 interface is NOT nb-independent

Measured 2026-08-27. Two layouts refining the **identical** cell range
(verified with `leaftable_test`: 32768 base cells, y 0..31 refined in both)
produce fields differing by ~1e-4:

| comparison | max_abs (un / pn) after 20 steps |
|---|---|
| cubic `nb = 8` vs cubic `nb = 4`, same refined region | 8.1e-05 / 9.1e-05 |
| non-cubic `nb = 16 16 8` vs `nb = 32 32 4` | 6.7e-05 / 5.6e-05 |

**This predates per-direction nb** — the first row uses the plain scalar form
and reproduces to 9 significant digits on a pre-Phase-1 binary
(8.1000457953e-05 vs 8.1000457967e-05, the residue being nofma vs fma). It is
not a Phase-1 regression, and it is consistent with the historical record:
every earlier phase gated "channel `nb = 4` **without refinement** bit-exact",
never with.

So the CLAUDE.md statement that results are "EXACTLY independent of nb and rank
count" holds for the **single-level lattice** (where the redundant halo-layer
sweep makes it exact) and for **rank count always**, but not across a 2:1
interface, where a block boundary lying on the interface plane is not
equivalent to one lying away from it. Refined layouts differ at the level of
the interface truncation error, which is the expected order — but they are not
bit-identical, and no gate should assume they are.
