# 2:1 interface — validation cases

Curated benchmarks for the energy-conserving (constant-1/2) 2:1 interface work
(`docs/interface_band_handout.md`). Ordered fast → slow. Build first:

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
./compile.sh cpu && ./compile.sh gpu
```

| # | case | what it checks | runtime |
|---|------|----------------|---------|
| 1 | **interface benchmark** | (A) stability: constant-1/2 BOUNDED vs baseline BLOW-UP; (B) interface banding: u'/v' rms excess (tracked, to reduce) | ~5 min (GPU) |
| 2 | **KEBAL energy gate** (slab) | convective KE production at the interface band | ~5 s (CPU) |
| 3 | **bit-exact no-interface** | single-level run unchanged by the interface code | ~1 min |
| 4 | **bit-exact multi-rank** | 1-rank == 2-rank (x-split) | ~1 min (GPU) |
| 5 | **developed stats** | time-averaged mean + Reynolds stresses, interface signatures | ~7-14 h |
| 6 | **full turbulence validation** | reference vs refined_y110/y55, spectra + profiles | hours-days |

## 1. Stability benchmark (the headline regression — THIS session's case)

The ~250-step refined channel, checking BOTH (A) stability and (B) the interface
banding. Self-checking on stability; banding reported as tracked metrics.

```bash
cd validation/channel_interface/interface_benchmark
python3 run_benchmark.py --arch gpu --ranks 1        # energy (bounded) + baseline (blow-up)
python3 run_benchmark.py --arch gpu --ranks 1 --energy-only   # just the stable case
```
(A) Stability PASS = energy `div_max < 0.5` (measured ~0.09, settles) AND baseline
blows up (`div_max -> 1e11` ~step 200). Discriminates the fix from the bug.
(B) Banding (energy final field, current default values, LOWER = better):
`u' excess x1.45-1.56` at the interfaces (streak band) and `v' kink 0.24-0.31`
(the wall-normal step). These are the numbers a better interface treatment must
reduce -- the constant-1/2 default stabilises the band but does not remove it.

## 2. KEBAL energy gate (Beltrami slab, V&V metric)

```bash
MOBY_KEBAL=1 mpirun -n 1 ./build_cpu/main validation/beltrami/slab_y_diag.ini | grep SKEW
```
Reports convective KE production split band vs interior, DIV and clean SKEW form.
Uniform (no interface) interior SKEW = round-off (2e-14); slab band SKEW = +0.13
(the defect). Toggle `MOBY_NORECON` / `MOBY_VELINJECT` / `MOBY_KESKEW` to see each
V&V lever cut the band energy (0.127 -> 0.054 -> 0.045 -> 0.036).

## 3. Bit-exact no-interface (single-level regression)

A single-level (no 2:1) run must be bitwise identical with the interface code
present — the constant-1/2 default and the skew/KEBAL paths are inert without an
interface. Build BOTH sides `-Mnofma` and compare with `tools/compare_fields.py`
(any single-level case, e.g. `validation/beltrami/uniform.ini`). Expect max_abs = 0.

## 4. Bit-exact multi-rank (x-split)

1-rank vs 2-rank must be identical (the y wall-bands stay whole under an x-split):

```bash
# run benchmark.ini (or any refined ini) at -n 1 and -n 2 in separate dirs, then:
python3 tools/compare_fields.py run1/channel_field_*.h5 run2/channel_field_*.h5
```
Expect max_abs = 0 for u,v,w,p (verified for constant-1/2 default AND MOBY_KESKEW).
Split x or z (`[mpi] dims = N 1 1`), NOT y. The solver auto-pins rank→GPU.

## 5. Developed-flow time-averaged statistics

The converged statistics (the single snapshot at t≈0.08 is NOT converged — see the
handout's v' interface-normal step observation). Two-leg driver + post-processor:

```bash
cd validation/channel_interface/developed
python3 run_developed.py --arch gpu --ranks 2            # constant-1/2
python3 run_developed.py --arch gpu --ranks 2 --skew     # + interface_skew
python3 ../../../tools/plot_channel_stats.py stats.png \
    runs/default/stats/channel_stats.h5:constant-1/2 \
    runs/skew/stats/channel_stats.h5:+KESKEW
```
See `developed/README.md`. Compare against a uniform-fine reference to quantify the
v' interface-normal step and whether the residual transfer defect remains.

## 6. Full turbulence validation (existing)

`./run_validation.sh [gpu|cpu] [nranks]` — reference (uniform 256x128x256) vs
refined_y110 / refined_y55, transient + long stats; post with
`tools/channel_interface_validation.py`.
