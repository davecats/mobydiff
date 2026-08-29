# The xz leaf-key change, measured. CPU 4 ranks, istmcetus, 2026-08-29

A/B of the y-tile-to-low-bits leaf key (`4f04f6b`) against `2e06bef`, both
built with production flags via `compile.sh`, refined production case
(`refined_yp82_rect_jacobi.ini`, 448 leaves, 60.56 M cells), 4 CPU ranks,
20 steps a side, back to back.

**The GPU A/B is still not done**: istmcetus' two GPUs were held continuously by
another user's `nekrs` job (85–93 %, ~14.6 GB each) through a 6-hour wait
window, and it is the only 2-GPU host available. A 12-hour waiter is armed.

## Traffic — exact, and independent of anything running on the machine

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

## Still open

- The GPU 2-rank A/B (waiter armed, 12 h). Expected to be larger than the CPU
  result in relative terms: at 2 GPU ranks the exchange was 17 % of the step
  with `mpi_wait` 7.2 %, against a much faster compute baseline.
- Whether the reduced `unpack` (5.2×) and raised `local_copy` (+61 % in time,
  +25 % in points) net out the same way on GPU, where local copies run at
  device bandwidth rather than host.
