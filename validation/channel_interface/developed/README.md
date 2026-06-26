# Developed 2:1 wall-band channel — time-averaged statistics

Run a developed turbulent channel (Re_tau = 180) with 2:1-refined wall bands and
collect **time-averaged** statistics, to quantify the interface signatures cleanly
(the single-snapshot probe at t≈0.08 is not converged — see
`docs/interface_band_handout.md` "Session 2026-06-26b"). The energy-conserving
**constant-1/2 interface is the default** (`[blocks] interface_constant_half`),
which keeps the refined channel stable.

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
python3 run_developed.py --arch gpu --ranks 2              # constant-1/2 default
python3 run_developed.py --arch gpu --ranks 2 --skew       # + interface_skew
# -> runs/<name>/stats/channel_stats.h5 (+ _l1)
python3 ../../../tools/plot_channel_stats.py stats.png \
    runs/default/stats/channel_stats.h5:constant-1/2 \
    runs/skew/stats/channel_stats.h5:+KESKEW
```

`--ranks N` sets the x-decomposition (`dims = N 1 1`) and `mpirun -n N`. The manual
steps below are equivalent.

## 0. Build (on the run machine)

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
./compile.sh gpu        # build_gpu/main
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
mpirun -n 2 ../../../../build_gpu/main ../transient.ini
# final field is channel_field_<laststep>.h5 -- copy it up as the developed restart:
cp channel_field_*.h5 ../transient_t5.h5
cd ..
```

## 3. Phase 2 — developed run with statistics

`developed.ini` restarts from `transient_t5.h5`. **Run it in a FRESH directory**
(no pre-existing `channel_stats.h5`) so the accumulator starts at zero:

```bash
mkdir -p run_developed && cd run_developed
mpirun -n 2 ../../../../build_gpu/main ../developed.ini
# -> channel_stats.h5 and channel_stats_l1.h5 (flushed every stats_write_interval
#    and at the end). runtimedata.txt has the running bulk quantities.
cd ..
```

To also test the optional skew band correction, set `interface_skew = true` in
`developed.ini` (or `MOBY_KESKEW=1`), run in a separate directory, and overlay.

## 4. Post-process

```bash
python3 tools/plot_channel_stats.py channel_stats.png \
    run_developed/channel_stats.h5:constant-1/2
# overlay two runs:
python3 tools/plot_channel_stats.py channel_stats.png \
    run_dev_default/channel_stats.h5:constant-1/2 \
    run_dev_skew/channel_stats.h5:+KESKEW
```

Produces mean U+ vs y+ (law of the wall), u'/v'/w' rms, and −⟨u'v'⟩, fine near the
walls + coarse in the core, with the 2:1 interface marked.

## Reflux ON-vs-OFF study (the current open question)

The u'/v' interface BANDS are a momentum-reflux artifact (the reflux injects the
fine-side resolved Reynolds-stress flux into the under-resolved coarse interface
cell; see `docs/next_session_orientation_asymmetry.md`). `run_reflux_study.sh`
runs the developed two-leg stats for reflux ON (default) and reflux OFF and prints
the overlay command, to test the trade: reflux OFF removes the band but the reflux
exists to conserve the MEAN interface flux (-<u'v'>), so check whether reflux OFF
degrades the mean profile / Reynolds shear vs the uniform-fine reference.

```bash
cd validation/channel_interface/developed
./run_reflux_study.sh gpu 2          # arch nranks; runs reflux_on + reflux_off
# (single case: run_developed.py --no-reflux)
```

## Multi-rank / multi-GPU notes

- The solver pins each rank to a GPU automatically (node-local rank → device,
  `comm.f90`), so `mpirun -n 2` uses both GPUs with no extra flags.
- **Decompose x (or z), not y** — `[mpi] dims = 2 1 1` in the inis. The 2:1 wall
  bands are in y; an x/z split keeps every interface whole on its rank (no
  cross-level exchange across a rank boundary). y-splits are untested for the
  refined path.
- Verified: 1-rank vs 2-rank (dims 2 1 1) is **bit-identical** (u,v,w,p), both for
  the default constant-1/2 and with `interface_skew`/`MOBY_KESKEW`.
- `dims` must divide the block lattice: nb=8 → 16 blocks in x, so 2/4/8/16 ranks
  in x are all valid.
