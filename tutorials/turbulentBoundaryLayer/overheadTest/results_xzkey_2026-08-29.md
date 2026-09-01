# The xz leaf-key change, measured. istmcetus, 2026-08-29 / 09-01

A/B of the y-tile-to-low-bits leaf key (`4f04f6b`) against `2e06bef`, both
built with production flags via `compile.sh`, refined production case
(`refined_yp82_rect_jacobi.ini`, 448 leaves, 60.56 M cells).

## GPU, 2 ranks — the target configuration (2026-09-01, clean machine)

400 steps a side, back to back, **zero foreign compute apps before AND after
each run**, drift flat (last 75 steps: ref 0.32358 → 0.32394, cand 0.31202 →
0.31210 s/step).

| | ref | cand | |
|---|---|---|---|
| **step** | **0.323940** | **0.312171** | **−3.63 %** |

The runtime diagnostic lines are **identical** between the two binaries at
every output step (`L2_div 2.93355685E-05`, `Linf 1.00242065E+00` at step 400)
— the renumbering is physics-invariant on the real production case at 2 ranks,
independently of the field comparison already done.

| exch bucket (s/step) | ref | cand | delta | % of step |
|---|---|---|---|---|
| pack | 0.005345 | 0.004938 | −0.000407 | −0.13 % |
| mpi_post | 0.000097 | 0.000097 | — | — |
| **mpi_wait** | **0.021820** | **0.003003** | **−0.018817** | **−5.81 %** |
| unpack | 0.005143 | 0.004376 | −0.000768 | −0.24 % |
| local_copy | 0.021436 | 0.024587 | +0.003151 | +0.97 % |
| **total exchange** | **0.053842** | **0.037001** | **−0.016841** | **−5.20 %** |

`mpi_wait` falls **7.3×**, from 6.74 % of the step to 0.96 %; the whole exchange
from 16.6 % to 11.9 %. Inside the projection, `phi_exchange` −2.16 % of the step
and `vel_exchange` −1.42 %.

**Why 3.63 % and not the ~6 % predicted.** The traffic prediction was exact
(20.3× predicted, 20.5× measured), but two things eat the time conversion:

- `local_copy` **rises 0.97 % of the step** — the interface traffic did not
  vanish, it moved onto the device, and a local copy is cheap but not free;
- the **volume kernels get 0.37 % slower** (`sweep` +0.32 %, `apply` +0.05 %).
  The new leaf order changes which blocks are adjacent in memory, and the block
  loop pays a little for it. That gives back 7 % of the exchange saving.

Net −3.63 %, of which the exchange contributes −5.20 % and the rest is given
back. The lesson for the next byte-reduction estimate: **a byte model predicts
traffic, not time** — budget for where the bytes land and for second-order
effects on the kernels that read them.

## CPU, 4 ranks (2026-08-29, contended machine)

Run while the GPUs were held by a foreign job; kept because the bucket result
is unambiguous and it exposed a mechanism the GPU run also shows.

### Traffic — exact, and independent of anything running on the machine

| 4 CPU ranks | old key | new key | |
|---|---|---|---|
| peer send points | 1 390 900 | **138 000** | 10.1× |
| per-rank send min / max | 101 300 / 963 600 | **23 000 / 46 000** | asymmetry 9.5× → 2.0× |
| peers per rank | 3 | **2** | |
| scalar (nv=1) round | 11.127 MB | **1.104 MB** | |
| full (nv=4) round | 44.509 MB | **4.416 MB** | |
| local copy points | 5 033 780 | 6 286 680 | +25 % |

The interface traffic moved off the wire and into host-local copies, as
designed, and the per-rank send volume went from a 9.5× spread to 2×.

## Timing — read the buckets, not the step

| exch bucket (s/step) | old key | new key | factor |
|---|---|---|---|
| pack | 1.2170 | 1.1372 | 1.1× |
| mpi_post | 0.0021 | 0.0020 | 1.0× |
| **mpi_wait** | **7.7716** | **0.0853** | **91×** |
| unpack | 1.8182 | 0.3481 | 5.2× |
| local_copy | 5.3243 | 8.5593 | 0.6× (rises, as designed) |
| **total exchange** | **16.1333** | **10.1319** | **1.59×** |

`mpi_wait` falls from **7.9 % of the step to 0.09 %**. The exchange as a whole
falls from 16.5 % to 11.1 %.

**Do not quote the step time.** Both sides ran at ~98 / 91 s/step against the
**18.77 s/step this exact configuration measured on a quiet machine**
(`results_multirank_2026-08-28.md`) — a ~5.2× inflation from the foreign
GPU job's host-memory traffic. Both sides saw the same environment, so the
−6.7 % step delta is indicative, but it is not the clean number and the clean
number is not in hand.

The bucket result does not depend on that: a **91× collapse in one bucket**,
with the byte count independently known to have fallen 10×, is far outside
anything contention can manufacture.

## What this corrects about the CPU exchange

`results_exchange_diag_2026-08-28.md` read the CPU exchange as *latency /
progress*-bound: 0.46 GB/s effective, nothing progressing outside the
`Waitall`. That was only half right. A 10× byte reduction produced a **91×**
wait reduction — strongly superlinear, so per-round latency was never the
binding term.

The likely mechanism is the one the size report makes visible: under the old
key one rank sent 963 600 points while another sent 101 300, a **9.5×
straggler**, and `MPI_Waitall` pays for the slowest peer. The new key
equalises that to 2×. So a large part of what looked like transfer cost was
**exchange load imbalance**, invisible in the compute-side load balance (which
was 1.000 throughout, since every rank owns the same number of leaves).

That is worth carrying: **leaf-count balance is not exchange balance.** The
size report's `send pts/rank min … max` is the number that shows it, and it
should be read on any new decomposition.

## Both platforms, side by side

| | CPU 4 ranks | GPU 2 ranks |
|---|---|---|
| peer send points | 10.1× fewer | 20.5× fewer |
| `mpi_wait` | 91× lower | 7.3× lower |
| `mpi_wait` share of step | 7.9 % → 0.09 % | 6.74 % → 0.96 % |
| total exchange share | 16.5 % → 11.1 % | 16.6 % → 11.9 % |
| step | not usable (contended) | **−3.63 %** |

The CPU's 91× against the GPU's 7.3× is the straggler effect: the old key gave
the CPU decomposition a 9.5× send-volume spread across 3 peers, and `Waitall`
pays for the slowest. The GPU case has one peer per rank, so there is no
straggler to remove — only volume — and its reduction is correspondingly the
one the byte model predicts.
