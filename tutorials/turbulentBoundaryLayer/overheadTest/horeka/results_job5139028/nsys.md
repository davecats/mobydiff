# Timeline probe — is the wait device work, and is it every round?

No MPI events: OpenMPI's Fortran mpi_f08 bindings reach the C layer as
PMPI_*, so nsys's MPI_* interception never fires. The exchange's kernel
signature carries the same information: `pack -> copy_local -> [gap =
post + Waitall] -> unpack`.

**READ SHAPE, NOT MAGNITUDE.** Tracing inflates the wait unevenly
(measured: base 120 -> ~690 us, rect 738 -> ~1255 us), so no absolute
value here may be quoted against an untraced run. Ratios WITHIN one
traced run, and the distribution shape, are what this file is for.


## `base_jacobi_r8_N2` — 2 node(s): hkn0507 hkn0526

| rank | span ms | **GPU busy** | idle ms | gaps>50us | median | p90 | share of idle | rounds | mean wait | p50 | top 1 % | H2D n | P2P ms | copy ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 3973 | **69.6 %** | 1207 | 2065 | 73 us | 148 | 75 % | 588 | 599.4 us | 441.1 | 18 % | 60218 | 154.5 | 0.0 |
| 1 | 3973 | **69.5 %** | 1208 | 2069 | 75 us | 147 | 75 % | 588 | 616.8 us | 455.1 | 18 % | 60218 | 157.8 | 0.0 |
| 2 | 3976 | **70.3 %** | 1179 | 1864 | 68 us | 118 | 72 % | 588 | 631.1 us | 474.8 | 17 % | 62500 | 155.3 | 0.0 |
| 3 | 3979 | **70.2 %** | 1185 | 1880 | 68 us | 118 | 72 % | 589 | 645.0 us | 482.0 | 17 % | 62559 | 161.0 | 0.0 |
| 4 | 3978 | **69.2 %** | 1227 | 1994 | 83 us | 179 | 75 % | 588 | 559.1 us | 454.9 | 14 % | 60260 | 157.8 | 0.0 |
| 5 | 3964 | **69.1 %** | 1222 | 1998 | 74 us | 165 | 74 % | 587 | 601.8 us | 453.5 | 18 % | 60107 | 159.5 | 0.0 |
| 6 | 3981 | **70.3 %** | 1183 | 1777 | 67 us | 131 | 72 % | 589 | 605.9 us | 470.1 | 17 % | 62574 | 165.9 | 0.0 |
| 7 | 3981 | **70.1 %** | 1189 | 1781 | 68 us | 127 | 72 % | 589 | 621.4 us | 473.8 | 17 % | 62568 | 160.5 | 0.0 |

- GPU busy 69.1–70.3 % across ranks.
  **The GPU is idle through the wait — the stall is NOT unfinished
  device work.** A rank blocked on a CUDA stream would show the GPU busy.
- The slowest 1 % of rounds hold 14–18 % of the wait.

## `rect_jacobi_r4_N1` — 1 node(s): hkn0507

| rank | span ms | **GPU busy** | idle ms | gaps>50us | median | p90 | share of idle | rounds | mean wait | p50 | top 1 % | H2D n | P2P ms | copy ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 7006 | **83.8 %** | 1132 | 753 | 82 us | 209 | 63 % | 1138 | 229.3 us | 208.4 | 2 % | 79829 | 19.7 | 412.1 |
| 1 | 7008 | **83.7 %** | 1144 | 572 | 95 us | 403 | 63 % | 1138 | 244.3 us | 208.9 | 4 % | 79839 | 43.5 | 410.7 |
| 2 | 7005 | **83.8 %** | 1135 | 547 | 92 us | 321 | 62 % | 1138 | 231.6 us | 208.4 | 3 % | 79829 | 40.8 | 411.7 |
| 3 | 7006 | **83.4 %** | 1159 | 640 | 112 us | 347 | 64 % | 1138 | 248.4 us | 207.1 | 3 % | 79829 | 20.9 | 413.3 |

- GPU busy 83.4–83.8 % across ranks.
  **The GPU is idle through the wait — the stall is NOT unfinished
  device work.** A rank blocked on a CUDA stream would show the GPU busy.
- The slowest 1 % of rounds hold only 2–4 % of the
  wait: it is EVERY round, so a per-round mean is a fair description.

## `rect_jacobi_r8_N2` — 2 node(s): hkn0507 hkn0526

| rank | span ms | **GPU busy** | idle ms | gaps>50us | median | p90 | share of idle | rounds | mean wait | p50 | top 1 % | H2D n | P2P ms | copy ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 5016 | **56.0 %** | 2205 | 1923 | 368 us | 1552 | 83 % | 1054 | 1235.5 us | 1197.2 | 3 % | 73951 | 18.3 | 184.8 |
| 1 | 5016 | **55.8 %** | 2215 | 2177 | 359 us | 1370 | 83 % | 1054 | 1270.0 us | 1238.6 | 2 % | 73972 | 28.2 | 182.9 |
| 2 | 5016 | **56.3 %** | 2192 | 1396 | 555 us | 1848 | 82 % | 1055 | 1304.1 us | 1607.3 | 2 % | 74015 | 32.2 | 183.8 |
| 3 | 5016 | **60.1 %** | 2003 | 362 | 56 us | 117 | 30 % | 1055 | 1312.5 us | 1674.7 | 2 % | 75731 | 13.6 | 183.9 |
| 4 | 5019 | **57.9 %** | 2113 | 949 | 71 us | 178 | 32 % | 1056 | 1318.8 us | 1658.4 | 2 % | 155939 | 14.3 | 184.2 |
| 5 | 5010 | **56.0 %** | 2205 | 2274 | 325 us | 1342 | 83 % | 1053 | 1263.5 us | 1534.2 | 2 % | 73878 | 29.1 | 183.2 |
| 6 | 5006 | **56.2 %** | 2187 | 2177 | 350 us | 1175 | 82 % | 1053 | 1236.5 us | 1542.5 | 2 % | 73863 | 30.2 | 183.3 |
| 7 | 5020 | **55.8 %** | 2220 | 2022 | 345 us | 1398 | 83 % | 1055 | 1223.8 us | 1035.1 | 2 % | 74042 | 18.3 | 184.4 |

- GPU busy 55.8–60.1 % across ranks.
  **The GPU is idle through the wait — the stall is NOT unfinished
  device work.** A rank blocked on a CUDA stream would show the GPU busy.
- The slowest 1 % of rounds hold only 2–3 % of the
  wait: it is EVERY round, so a per-round mean is a fair description.

## `rect_jacobi_r8_N4` — 4 node(s): hkn0507 hkn0526 hkn0534 hkn0604

| rank | span ms | **GPU busy** | idle ms | gaps>50us | median | p90 | share of idle | rounds | mean wait | p50 | top 1 % | H2D n | P2P ms | copy ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4193 | **66.2 %** | 1415 | 886 | 194 us | 242 | 73 % | 1044 | 304.2 us | 240.4 | 7 % | 73233 | 13.7 | 182.9 |
| 1 | 4192 | **65.6 %** | 1441 | 897 | 199 us | 404 | 74 % | 1043 | 338.8 us | 295.8 | 6 % | 73185 | 14.8 | 181.1 |
| 2 | 4193 | **65.8 %** | 1436 | 1230 | 138 us | 259 | 74 % | 1044 | 315.7 us | 300.4 | 3 % | 73240 | 16.0 | 182.1 |
| 3 | 4193 | **65.8 %** | 1432 | 943 | 184 us | 328 | 73 % | 1044 | 313.4 us | 303.6 | 3 % | 73240 | 13.7 | 181.8 |
| 4 | 4177 | **65.7 %** | 1433 | 942 | 209 us | 353 | 74 % | 1040 | 324.5 us | 333.1 | 3 % | 72950 | 14.2 | 181.0 |
| 5 | 4177 | **65.7 %** | 1430 | 1118 | 163 us | 273 | 74 % | 1040 | 320.2 us | 312.3 | 3 % | 72950 | 16.3 | 181.2 |
| 6 | 4195 | **65.6 %** | 1444 | 940 | 209 us | 393 | 74 % | 1045 | 331.8 us | 307.8 | 3 % | 73311 | 13.2 | 181.8 |
| 7 | 4195 | **65.8 %** | 1434 | 1068 | 203 us | 264 | 74 % | 1045 | 319.0 us | 285.3 | 3 % | 73311 | 15.9 | 183.5 |

- GPU busy 65.6–66.2 % across ranks.
  **The GPU is idle through the wait — the stall is NOT unfinished
  device work.** A rank blocked on a CUDA stream would show the GPU busy.
- The slowest 1 % of rounds hold only 3–7 % of the
  wait: it is EVERY round, so a per-round mean is a fair description.
