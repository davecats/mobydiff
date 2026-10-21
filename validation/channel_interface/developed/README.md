# Developed 2:1 wall-band channel — time-averaged statistics

Run a developed turbulent channel (Re_tau = 180) with 2:1-refined wall bands and
collect **time-averaged** statistics, to quantify the interface signatures cleanly
(the single-snapshot probe at t≈0.08 is not converged — see
`docs/history/next_session_edges_les.md`). The energy-conserving
**constant-1/2 interface** keeps the refined channel stable; it is hardwired
(the `[blocks] interface_constant_half` toggle was removed 2026-07-01).

> **Note (2026-09-26).** This directory records the reflux ON-vs-OFF study and
> an `interface_skew` variant. Both toggles (`[blocks] momentum_reflux`,
> `interface_skew`) and the `MOBY_KESKEW` hook were REMOVED once the study
> settled it: the solver now always runs the reflux-OFF, const-1/2 interface,
> i.e. what the runs below call `reflux_off`. The driver options that set those
> keys (`run_developed.py --skew/--no-reflux`, `run_reflux_study.sh`) are gone;
> the recorded conclusions below are unchanged. (`plot_interface_validation.py`,
> which plotted `runs/reflux_off` against `runs/reflux_on`, was removed on
> 2026-09-26; git history.)

Two phases (the stats accumulator starts fresh in phase 2, so the discarded
transient does not pollute it):

| phase | ini | time | stats |
|-------|-----|------|-------|
| 1 transient (discard) | `transient.ini` | t = 0 .. 5 | OFF |
| 2 statistics | `developed.ini` | t = 5 .. 25 | ON → `channel_stats.h5` (+`_l1`) |

At `dtmax = 3.125e-4` that is ~16 000 + ~64 000 steps (≈0.6 s/step on one RTX-class
GPU → ~14 h total; ~half on 2 GPUs).

## Re-validated 2026-10-21: the first statistics at `niter = 12`

RE-VALIDATED 2026-10-21, commit 532ef54 (banner), production flags, istmcetus
A6000, 1 rank, transient 16 000 steps (35 min) + statistics 64 000 steps
(0.120 s/step, 2.1 h): `run_developed.py --ranks 1 --name cheb12` from
`IC.h5` (the interpolated KMM180 restart). No number had been recorded for
this case before (its reflux study below is settled history). Settings
changed for the rerun: `cflmax = 1.0` (dtmax 3.125e-4 binds, so the step is
unchanged), `sor = 0.8` kept EXPLICIT (ignored by Chebyshev, but the IC
carries the KMM180 run's red-black `sor = 1.5` and the restart check of
2026-10-21 refuses a defaulted key that disagrees with the file).
`analyse_developed.py` (this directory) prints the table; both levels'
statistics files combined, walls folded:

| quantity | cheb12 (Chebyshev 12) | KMM / Moser Re_tau 180 |
|---|---|---|
| U+ centreline | 18.39 | 18.3 (+0.5 %) |
| bulk U+ | 15.82 | -- |
| u'+ peak | 2.630 | 2.66 (-1.1 %) |
| v'+ peak | 0.839 | 0.84 |
| w'+ peak | 1.082 | 1.08 |
| -<u'v'>+ peak | 0.736 | 0.73 (+0.9 %) |

Interface band at the 2:1 jump (y = 0.5986 and its mirror; the jump of each
profile between the last fine and the first coarse row, against the
variation over the two rows on either side -- a ratio near 1 is a smooth
profile): U+ 0.87, u'+ 0.52, v'+ 0.12, w'+ 0.75, -<u'v'>+ 0.41. No band.
PASS: within ~1 % of KMM on every quantity, like the single-level KMM180
tutorial run of the same day (`tutorials/channel_kmm180/README.md`).

### The red-black twin: the R1 admission test (2026-10-21)

`run_developed.py --ranks 1 --redblack` (`[pressure] solver = redblack`,
`niter = 6`, `sor = 1.5`, the same two legs from the same IC; istmcetus
A6000, the other GPU): the test `validation/redblack_interface` left open --
the developed Re_tau 180 band channel with red-black SOR across the 2:1
wall bands against the Chebyshev solution of the same case.

| quantity | cheb12 | redblack 6 | KMM |
|---|---|---|---|
| U+ centreline | 18.39 | 18.37 | 18.3 |
| bulk U+ | 15.82 | 15.82 | -- |
| u'+ peak | 2.630 | 2.618 | 2.66 |
| v'+ peak | 0.839 | 0.839 | 0.84 |
| w'+ peak | 1.082 | 1.084 | 1.08 |
| -<u'v'>+ peak | 0.736 | 0.736 | 0.73 |
| interface-band ratios (U+, u'+, v'+, w'+, -<u'v'>+) | 0.87 0.52 0.12 0.75 0.41 | 0.87 0.52 0.06 0.76 0.42 | -- |
| statistics leg, s/step (A6000, 1 rank) | 0.120 | 0.143 | -- |

The two projections give the same developed statistics to the sampling
scatter of a 20 t.u. window (differences <= 0.5 %) and neither has an
interface band: red-black SOR across a 2:1 interface WORKS (decision 1 of
the handout: it is admissible as a production projection on refined cases).
It is not the faster one here: 19 % more per step than Chebyshev 12 on this
band-refined grid (the per-colour cross-level `phi` exchange R1 adds; the
single-level pipe of `tutorials/cht/pipe` reads the opposite, red-black 6
at 0.630 s/step against Chebyshev 12 at 0.715 on the RTX 5090). The omega
stability limit with the interface present is between 1.6 and 1.7
(`validation/redblack_interface/run_r1_items.sh cpu`, 200 steps on this
grid): 1.5 has margin.

## Quick path: `run_developed.py`

One command runs both legs (generates `IC.h5` if missing, runs transient → stats,
prints the plot command):

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
./compile.sh gpu
cd validation/channel_interface/developed
python3 run_developed.py --arch gpu --ranks 2              # -> runs/default
python3 run_developed.py --arch gpu --ranks 2 --refine-dims xz   # -> runs/xz (quadtree bands)
# -> runs/<name>/stats/channel_stats.h5 (+ _l1)
python3 ../../../tools/plot_channel_stats.py stats.png \
    runs/default/stats/channel_stats.h5:constant-1/2
```

`--ranks N` sets the x-decomposition (`dims = N 1 1`) and `mpirun -n N`. The manual
steps below are equivalent.

## 0. Build (on the run machine)

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
./compile.sh gpu        # build_gpu/moby_solve
```

## 1. Initial condition

Interpolate the bundled KMM180 restart onto the refined grid (once):

```bash
python3 tools/make_channel_restart.py --mode refined --band-cells 24 \
    --source tutorials/channel_kmm180/channel_kmm180_restart.h5 \
    --out validation/channel_interface/developed/IC.h5
```

## 2. Phase 1 — transient (writes the t=5 restart)

```bash
cd validation/channel_interface/developed
mkdir -p run_transient && cd run_transient
mpirun -n 2 ../../../../build_gpu/moby_solve ../transient.ini
# final field is channel_field_<laststep>.h5 -- copy it up as the developed restart:
cp channel_field_*.h5 ../transient_t5.h5
cd ..
```

## 3. Phase 2 — developed run with statistics

`developed.ini` restarts from `transient_t5.h5`. **Run it in a FRESH directory**
(no pre-existing `channel_stats.h5`) so the accumulator starts at zero:

```bash
mkdir -p run_developed && cd run_developed
mpirun -n 2 ../../../../build_gpu/moby_solve ../developed.ini
# -> channel_stats.h5 and channel_stats_l1.h5 (flushed every stats_write_interval
#    and at the end). runtimedata.txt has the running bulk quantities.
cd ..
```

## 4. Post-process

```bash
python3 tools/plot_channel_stats.py channel_stats.png \
    run_developed/channel_stats.h5:constant-1/2
# overlay several runs: append more <stats.h5>:<label> arguments
```

Produces mean U+ vs y+ (law of the wall), u'/v'/w' rms, and −⟨u'v'⟩, fine near the
walls + coarse in the core, with the 2:1 interface marked.

## Reflux ON-vs-OFF study (SETTLED: reflux removed; historical record)

The u'/v' interface BANDS are a momentum-reflux artifact (the reflux injects the
fine-side resolved Reynolds-stress flux into the under-resolved coarse interface
cell; see `docs/history/next_session_edges_les.md`). The study (driver
`run_reflux_study.sh`, removed 2026-09-26 with the key; in git history) ran the
developed two-leg stats for reflux ON and OFF to test the trade: reflux OFF
removes the band but the reflux existed to conserve the MEAN interface flux
(-<u'v'>). Outcome (CLAUDE.md): reflux OFF removes the band at no cost to the
mean profile / Reynolds shear vs the uniform-fine reference, so OFF is now the
only behaviour.

`run_reference.sh gpu N` runs the **uniform-fine reference** (256x128x256, single
level, two-leg, same t-window/params) — the ground truth the refined reflux_on /
reflux_off runs are compared against (the refined fine level shares this grid
bitwise). It generates `runs/reference/REF_IC.h5` on first use and prints the
3-way overlay command (reflux-on / reflux-off / uniform). Tested on 1 GPU
(~0.29 s/step, fits in 6 GB at 256^3); much faster on a multi-GPU box.

## Multi-rank / multi-GPU notes

- The solver pins each rank to a GPU automatically (node-local rank → device,
  `comm.f90`), so `mpirun -n 2` uses both GPUs with no extra flags.
- **Decompose x (or z), not y** — `[mpi] dims = 2 1 1` in the inis. The 2:1 wall
  bands are in y; an x/z split keeps every interface whole on its rank (no
  cross-level exchange across a rank boundary). y-splits are untested for the
  refined path.
- Verified: 1-rank vs 2-rank (dims 2 1 1) is **bit-identical** (u,v,w,p), both for
  the constant-1/2 interface and (historically) with the removed `interface_skew`.
- `dims` must divide the block lattice: nb=8 → 16 blocks in x, so 2/4/8/16 ranks
  in x are all valid.
