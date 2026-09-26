# Developed 2:1 wall-band channel — time-averaged statistics

Run a developed turbulent channel (Re_tau = 180) with 2:1-refined wall bands and
collect **time-averaged** statistics, to quantify the interface signatures cleanly
(the single-snapshot probe at t≈0.08 is not converged — see
`docs/next_session_edges_les.md`). The energy-conserving
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
cell; see `docs/next_session_edges_les.md`). The study (driver
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
