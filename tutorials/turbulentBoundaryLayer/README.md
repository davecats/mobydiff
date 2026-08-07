# Turbulent zero-pressure-gradient boundary layer (ZPG TBL)

A spatially-developing incompressible **ZPG turbulent boundary layer** DNS,
non-dimensionalised by the inlet displacement thickness (Re_δ*,0 = 450, U∞ = 1),
validated against the **SIMSON pseudo-spectral reference** (`passivewall.hdf5`,
Schmitt/KIT). This is the final, best-resolved configuration that emerged from an
extended trip/resolution/solver study — the full history is in
**[`tests_record.md`](tests_record.md)**.

## The case shipped here

| ingredient | value | why |
|---|---|---|
| grid | 4096 × 176 × 192 (138 M) | Δx⁺≈4 ≈ Δz⁺≈3.7, Δy⁺_wall≈0.2, Δy⁺_max≈4 (in δ₉₉) |
| x-grid | geometric, stretch 1.5 | Δx⁺ ~ uniform as u_τ falls downstream |
| y-grid | blayer (wall-clustered + freestream coarsening) | resolve the BL, keep the tall domain affordable |
| domain | 750 × 100 × 32 δ*₀ | ly = 100 δ*₀ (very tall) so the top pins a true ZPG |
| trip | Schlatter–Örlü, `trip_amp = 0.03` | gentle trip (matches the reference's development) |
| convection | skew-symmetric | energy-neutral under the incremental projection |
| **pressure solver** | **red-black SOR, niter = 6, sor = 1.5** | stable at low niter on this outlet case + ~1.8× faster than Chebyshev-Jacobi niter=12 |

## Result (converged, ~4000 t.u. of averaging)

At Re_θ = 677, vs the SIMSON spectral reference:

| quantity | this DNS | SIMSON | Δ |
|---|---|---|---|
| c_f | 0.00464 | 0.00471 | **−1.6 %** |
| H | 1.501 | 1.507 | **−0.4 %** |
| u′/v′/w′_rms peak | 2.69 / 1.02 / 1.31 | 2.63 / 1.02 / 1.30 | ~1–2 % |
| −u′v′ peak | 0.874 | 0.882 | −0.9 % |

Mean U⁺(y⁺) collapses onto the reference through the sublayer and log region; the
Reynolds stresses overlay the spectral profiles; the shape factor sits on the
reference. The residual ~1.6 % in c_f is the intrinsic 2nd-order-FD-vs-spectral
floor (see `tests_record.md` for the full decomposition of the original ~5 % gap
into tripping + streamwise resolution + this floor).

![c_f, H, U+ and Reynolds stresses vs SIMSON](assets/figures/passivewall_compare.png)

## Layout

```
README.md                 this file
tests_record.md           the full study: every case and what it showed
production.ini            phase 1 (re-equilibration), red-black solver
production_stats.ini      phase 2 (statistics accumulation)
production_stats.h5       the converged span+time statistics  [local, ~115 MB]
restart_field.h5          a developed field (IC + grid + snapshot)  [local, 4.4 GB]
reproduce.py              regenerate every figure from the above
assets/figures/           the PNGs (committed)
assets/postpro/           the post-processing scripts + reference data
    passivewall.hdf5      SIMSON spectral reference  [local, ~103 MB]
    tbl_uncontrolled.mat  earlier spectral reference (mean only)
    ref_schlatter_orlu_Re670.prof   KTH Re_θ=677 profile
```

Files marked *[local]* are too large for git (`restart_field.h5`,
`passivewall.hdf5` and `production_stats.h5` all exceed GitHub's 100 MB limit);
keep them in place to reproduce the figures.

## Reproduce

**Figures** (fast, from the shipped statistics):
```bash
python3 reproduce.py            # -> assets/figures/*.png
```

**The DNS** (multi-day GPU run):
```bash
module load toolkits/nvhpc/25.9
mpirun -n 1 ./build_gpu/moby_solve production.ini        # phase 1 from restart_field.h5
mpirun -n 1 ./build_gpu/moby_solve production_stats.ini  # phase 2, restart from the phase-1 end field
```
