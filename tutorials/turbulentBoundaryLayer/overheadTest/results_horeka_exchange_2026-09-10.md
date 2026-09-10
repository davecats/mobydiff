# The exchange at scale is a per-ROUND cost, not a per-cell or per-point one

<!-- DRAFT: sections 3 onward are filled from job 5139461 (dev_accelerated) and
     the 16-rank follow-up. Sections 1-2 are arithmetic on committed logs and are
     final. -->

Phases 1 and 2 of `docs/next_session_2to1_performance.md`. The handout set the
target as "the refined case does 82 % of the single-level twin's exchange work
while carrying 44 % of its cells — that factor ~1.9 is the 2:1 interface's cost
at scale, it is device-local, and it is the target", and pre-registered three
branches: cross-level volume, copy volume, or a per-point cost difference
between op kinds.

**None of the three is what the data says.** The exchange at 16 ranks is bought
almost entirely by the *number of exchange rounds* multiplied by a fixed cost per
kernel launch, and by nothing about the mesh at all.

## 1 — The 1.9x is not extra points, and this was already in the repository

`exchange sizes` is printed by every profiled run and has been since `be48d44`.
Summing its two numbers (local copy points, which never become messages, and
send points, which do) over the whole job:

| config | cells | exchange points | points per cell |
|---|---|---|---|
| `rect_jacobi` (single level) | 138.41 M | 15.390 M | 1.112e-1 |
| `refined_yp82_rect_jacobi` | 60.56 M | 6.425 M | 1.061e-1 |
| `refined_big_rect_jacobi` | 242.22 M | 25.786 M | 1.065e-1 |

**0.954x and 0.957x.** The refined cases exchange slightly FEWER points per cell
than the single-level twin, and the total is independent of rank count (the same
halo points, some of them promoted from a device copy to a message). So the
handout's first two branches cannot fire: there is no extra volume of either
kind to find, and "1.9x exchange per cell" is **1.9x per exchange POINT**.

## 2 — Where it does come from: a fixed cost per round

Fitting microseconds per exchange ROUND against points per rank, over the rank
counts of one configuration (`results_scaling/new/`, job 5139351; the 1-rank runs
are excluded — they have no peers, so `pack`/`unpack` are identically zero, and
their copy kernel is an order of magnitude larger than at any other rank count):

| config | pack fixed | unpack fixed | local_copy fixed | local_copy ns/pt | worst residual |
|---|---|---|---|---|---|
| `base_jacobi` (34.5k–574k send pts/rank) | **95.9 us** | **86.8 us** | — (no local copies) | — | 2.6 us |
| `rect_jacobi` | 101.5 | 93.6 | **56.1** | 0.098 | 1.3 us |
| `refined_yp82_rect_jacobi` | 98.4 | 85.5 | **114.5** | 0.099 | 0.9 us |
| `refined_big_rect_jacobi` | 103.6 | 97.5 | **109.8** | 0.101 | 1.6 us |
| `refined_yp82_rect_redblack` | 97.3 | 85.1 | **117.9** | 0.099 | 1.7 us |

Residuals of 1–3 us on numbers of 60–120 us, across a 17x span in points: the
affine split is not being forced onto the data. Reading it:

- **The marginal point costs the same everywhere.** 0.098–0.101 ns in every
  configuration — single level or refined, `nb = 16` or `nb = 64 44 48`. The 2:1
  interface carries **no per-point penalty**, which retires the handout's third
  branch as it was written.
- **The fixed part is 250–300 us per round**, and there are 39 rounds in every
  step of every configuration ever measured here. That is 9.8 ms/step for
  `rect_jacobi` and 11.6 ms/step for `refined_yp82` — **14 % and 28 %** of their
  16-rank steps.
- At 16 ranks the two configurations spend **the same absolute device-local
  exchange time**, 13.24 against 13.29 ms/step, while carrying 138.41 M and
  60.56 M cells. *That identity is the 1.9x.* The refined case is not doing more
  work; it is paying an identical, mesh-independent bill on 44 % of the cells.
- **The refined copy intercept is double the single-level one** (114.5 / 109.8 /
  117.9 against 56.1) at an identical slope. The candidate is structural: a full
  exchange launches the cross-level kernel as a SECOND kernel, in the 21 of 39
  rounds per step that carry one.

Everything above is arithmetic on logs that were already committed. What the new
runs add is (a) the same fit inside ONE allocation, (b) the cross-level kernel
timed in its own bucket instead of inferred from a difference, (c) the op split
that says how many points are cross-level at all, and (d) A0 at 8 and 16 ranks.
The predictions written before those runs are in
`horeka/exchange/PREREGISTERED.md`.

## 3 — Phase 2: A0 repeated at 8 ranks. It fails again, harder.

A0 (2026-09-02) inserted a compute-bound target kernel between the `Isend`/
`Irecv` posts and the `MPI_Waitall` at **2** GPU ranks and found that nothing was
hidden. Its own text asked for a repeat on a genuinely many-rank machine before
the conclusion was carried, because the message sizes, peer counts and transports
all change. Repeated here on 2 nodes, 8 ranks, 100 steps, job 5139461 — a
5.7 ms kernel against waits of 49 and 75 us, i.e. **76–116x** the quantity it
would have to hide:

| case | probe us/round | `mpi_wait` us/round, no probe | with probe | | s/step no probe | with probe |
|---|---|---|---|---|---|---|
| `refined_yp82_rect_jacobi` 8x2 | 5701.5 | 49.1 | **70.7** | **1.44x** | 0.060138 | 0.283078 |
| `rect_jacobi` 8x2 | 5700.4 | 75.0 | **79.4** | **1.06x** | 0.109144 | 0.331868 |

**The wait goes UP, never down.** And the step time rises by 0.2229 s = 39 x
5.7 ms exactly: the probe kernel is serialised into the step in full, with no
part of it overlapping anything. There is no asynchrony to exploit here at all —
not a partial one, not a small one.

So the 2-rank conclusion carries: **the transfer happens inside `MPI_Waitall`
and nowhere else, and P2 (overlap) is closed** on this stack. The scope caveat
the 2-rank probe attached to itself is now discharged rather than inherited.
The mild degradation (6 % and 44 %) is the same contention the 2-rank probe saw
(9 %).

