# LES &harr; IBM coupling validation (off-grid plane-wall channel)

**STATUS: VALIDATED (2026-06-30) with converged developed statistics (t=5..25).**
All gates pass on a single grid AND across the 2:1 interface; see `ibm_les_profiles.png`.

LES (WALE) is validated for **grid-aligned** channel walls across block refinement
and the 2:1 interface (`../les/`). The one untested piece was the **LES&harr;IBM
coupling** &mdash; the `ibm_aware` solid-cell `nut` masking in `src/modules/les.f90`
(a cell is forced `nut=0` iff any of its 6 staggered faces has `|ibm%coef|>1e20`).
This directory exercises it on an immersed wall that does **not** lie on grid nodes.

No solver code was changed: the existing `ibm_aware` mask + WALE `sd2` operator
already give the physically correct `nut&rarr;0` into the IBM wall with no spurious
band, on a single grid **and** across the 2:1 interface.

## The test case

A plane-wall channel whose two flat walls are the **file-based IBM** (built from
two wall-slab STLs via `tools/mobygeom.py`), deliberately placed **off grid node**:

| | |
|---|---|
| grid | `64 x 80 x 64`, **uniform** in all directions, `nb=8` |
| domain | `lx=4&pi;`, `ly=2.5`, `lz=2&pi;` (`dy=0.03125`) |
| walls | `y=0.259375` and `y=2.259375` &mdash; mid-cell (8.3&middot;dy, 72.3&middot;dy), off node AND off cell-centre |
| fluid gap | exactly **2.0** (half-height 1) &rarr; Re&tau;&asymp;180 with `re=180`, `forcing_x=1` |
| solid | ~8 cells thick each side &rarr; fully-solid cells (`coef=1e30/Re`) + one band cell per wall |
| LES | WALE, `momentum_reflux=false`, `interface_constant_half=true` |
| IC | KMM180 developed field y-shifted into the gap, solid zeroed |

Cases (`run_ibm_les.py`):
- **a_wale** &mdash; single level, IBM wall, WALE LES (the coupling under test).
- **b_none** &mdash; single level, IBM wall, LES off (control for the mean flow).
- **c_refine** &mdash; `refine_body` at the IBM wall (2:1-interface &times; IBM &times; LES);
  body-driven block-table refines the wall bands+buffer to level 1, coarse core
  (y&isin;[0.75,1.75]), buried-solid blocks removed.

## How to run

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
./compile.sh gpu                       # build_gpu/main must carry the nut output (this branch)
cd validation/channel_interface/les_ibm
# all prerequisite .h5 files are committed (or rsync this dir); to rebuild them: ./setup.sh
MP=/opt/.../hpcx-2.25.1/ompi/bin/mpirun   # or just "mpirun" on another host
python3 run_ibm_les.py --arch gpu --case all --mpirun "$MP"   # a_wale, b_none, c_refine
python3 measure_nut.py   --run runs/a_wale/stats              # gates 1-2
python3 ibm_les_stats.py                                       # gates 3-4 -> ibm_les_profiles.png
```
Defaults: `--t-transient 5 --t-average 20 --snap-interval 800` (the developed-stats
campaign). `measure_nut.py`/`ibm_les_stats.py` need only numpy + h5py.

## Gates & results

Converged over **51 snapshots, t=5..25** (a_wale / b_none / c_refine ran to ~25 600 steps):

| # | gate | result |
|---|------|--------|
| 1 | solid-cell `nut==0` (hard) | **PASS** &mdash; `max\|nut\|` over all solid cells = `0.0`, all 51 snapshots |
| 2 | no spurious wall-`nut` spike | **PASS** &mdash; `nut/&nu;` band cell 0.011 &rarr; core 0.222, **band/core = 0.05** (physical `nut&rarr;0`, not a spike) |
| 3 | mean-U law of the wall | **PASS** &mdash; a_wale, b_none AND the grid-aligned `../les/` LES channel collapse on `2.44 ln y&plus; +5` (U&plus;=13.1 vs 13.2 at y&plus;=29; sublayer 1.20 vs 1.13); WALE bulk U=15.09 &gt; no-LES 14.69 (SGS raises the log-layer U toward the reference); resolved stresses symmetric across both IBM walls |
| 4 | `nut` step (not band) across the 2:1 interface | **PASS** &mdash; the fine wall bands carry lower `nut` (halved filter width); at y=0.75/1.75 `nut` steps up to the coarse core by **2.03&times;** (toward the `delta^2` ratio; the resolved strain is also lower on the coarse side), a sharp **smooth step with NO overshoot/band** |
| 5 | stability 1000+ steps | **PASS** &mdash; all three cases ran ~25 600 steps to t=25, `div` bounded, mass ~1e-15, no NaN |
| 6 | bit-exact no-LES / no-IBM | **N/A** &mdash; no solver code changed; paths byte-identical by construction |
| 7 | CPU == GPU (masking branch) | **PASS** &mdash; case (a) 10 steps agree to 4.6e-14 (FMA round-off; the offloaded `ibm_aware` branch is correct on GPU) |

**Key finding:** the prime suspect &mdash; a spurious `nut` spike at the IBM band
cells where the SGS model reads the IBM-forced velocity drop as resolved strain
&mdash; does **not** materialise. WALE's `sd2` operator kills `nut` at the IBM band
just as it does at a grid-aligned wall, and the `nut` step across the 2:1 interface
that coincides with the IBM wall is the physical filter-width step, not a band. No
band-aware `nut` damping is needed; **no solver code change**. See `ibm_les_profiles.png`.

## Files

- `channel_ibm.ini` (a/b), `channel_ibm_refine.ini` (c).
- `grid.h5`, `wall_{lo,hi}.stl`, `ibm_coeff.h5` (single level), `IC.h5` &mdash; committed.
- `ibm_coeff_blocks.h5`, `IC_refine.h5` (case c, larger) &mdash; gitignored; `./setup.sh`
  regenerates them (needs the geometry venv) or rsync them with the directory.
- `make_walls_stl.py`, `make_ibm_ic.py`, `setup.sh` &mdash; data generators.
- `run_ibm_les.py` &mdash; campaign driver. `measure_nut.py` &mdash; gates 1-2.
  `ibm_les_stats.py` &mdash; gates 3-4. `RESUME_STATUS.md` &mdash; session handoff notes.

## Notes / gotchas
- The IBM is **implicit** (`mu = 1/(1+dt&gamma;&middot;coef)`, `ibm.f90:506`; `qs *= mu`)
  &mdash; no dt restriction; standard CFL governs.
- `mobygeom` needs `--grid-file` (from `mobygrid`); `read_restart_metadata`
  overwrites grid/BC/re/forcing/`ibm_enabled` from the IC attrs (minted via cold-start).
- les.f90 solid rule: cell solid iff **any** of its 6 staggered faces has `|coef|>1e20`
  (so the cell straddling the wall is masked; the first full-`nut` cell is ~1.5 cells in).
