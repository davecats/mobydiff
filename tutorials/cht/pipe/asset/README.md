# `asset/` — the reference data and the comparison

Everything here reproduces the comparison in [`../README.md`](../README.md)
**without the run outputs and without the 11.8 GB published archive**.

```bash
./pipe_numbers.py                 # every table in the README
./plot_pipe.py                    # figures 1-7  -> figures/
./plot_pipe.py --blocks           # figure 8
./compare_neuhauser.py velocity pipe_prod_snaps.npz
./compare_neuhauser.py thermal  pipe_prod_statsD.npz --scalar c0 --wall-flux 0.24854
```

## The reference

`neuhauser_profiles.npz` (250 kB) holds the radial profiles of all seven cases
plus the ⟨θ′²⟩ budget terms of c0, reduced from

> J. Neuhauser, *Conjugate Heat Transfer in Turbulent Pipe Flows with
> Non-Uniform Heating Effects at Two Prandtl Numbers*, NekRS v26.0,
> DOI 10.35097/26za20q32xsz43yk

The published bag is 11.8 GB and holds a 6.8 GB zarr of (time, r, phi) moments
for every case, Prandtl number and boundary condition. The comparison uses a
handful of 226-point radial profiles, so the tutorial ships the reduction.

**It is checkable, not merely convenient.** `extract_neuhauser.py` rebuilds
the cache from the archive:

```bash
tar xf '10.35097-26za20q32xsz43yk(1).tar'
tar xf .../data/dataset/cht_short/joined_datasets.interp.zarr.tar -C /somewhere
NEUHAUSER_ZARR=/somewhere/joined_datasets.interp.zarr ./extract_neuhauser.py
```

and figures 1–5 come out **byte-identical** whether `plot_pipe.py` reads the
cache or the archive. `compare_neuhauser.ref_case` prefers the cache and falls
back to the zarr automatically.

Conventions that must hold on both sides of every comparison, and do:

* wall units from **their** definitions — u_τ from d⟨w⟩/dr as r → R from the
  fluid side, y⁺ = (R − r)/δ_v, δ_v = ν/u_τ;
* θ_τ = q_w/u_τ with q_w from the **interface-heat balance** |Q|/(2πRL),
  never from a near-wall slope (it is differenced across the cut cells);
* the mean temperature referenced to the **interface value** on both sides —
  theirs is wall-referenced, ours absolute, so the raw levels differ by an
  additive constant;
* their radial grid (226 points, dr down to 0.04 wall units) is interpolated
  onto **our** bin centres, never the reverse.

## What is in here

| file | |
|---|---|
| `neuhauser_profiles.npz` | the reference: 7 cases × radial profiles + the c0 budget |
| `extract_neuhauser.py` | rebuilds it from the published archive |
| `pipe_prod_statsD.npz` | our plane statistics, binned into r (all 7 scalars, fluid and solid) |
| `pipe_prod_snaps.npz` | our snapshot statistics, binned into r (velocity + scalars) |
| `pipe_prod_heat.txt` | the solver's interface-heat diagnostic over the statistics leg |
| `grid_residual.npz` | the z-averaged cross-sections behind Fig. 6b,c |
| `pipe_prod_plane.npz` | one z-plane of one production snapshot (Fig. 7) |
| `pipe_prod_blocks.npz` | the leaf table and node lines (Fig. 8) |
| `make_caches.py` | rebuilds the last three from the run outputs |
| `compare_neuhauser.py` | the profile-by-profile comparison, and the shared reference loader |
| `pipe_numbers.py` | every table in the README |
| `plot_pipe.py` | the figures |
| `variance_budget.py` | the ⟨θ′²⟩ budget, both sides — **needs snapshots**, see below |
| `figures/` | fig1…fig8 |
| `pipe_report.html` | the full report: derivations, figure discussion, findings |

`variance_budget.py` is the one script that cannot run from the shipped
assets: it differentiates the instantaneous field, so it needs the snapshots
and the case file.

```bash
./variance_budget.py /path/to/p_pr_stat_1[89]*.h5 \
    --case /path/to/pipe_prod.h5 --qw 0.24854
```
