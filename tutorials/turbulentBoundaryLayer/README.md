# Zero-pressure-gradient turbulent boundary layer

Spatially developing ZPG turbulent boundary layers, DNS on the CaNS reference
grid (650 × 100 × 26 δ*₀, 3200 × 384 × 136, Re_δ*,0 = 450, tripped Blasius
inflow, conservative convection, red-black SOR projection).

| directory | case |
|---|---|
| [`smooth/`](smooth/README.md) | smooth wall: the code comparison against SIMSON, CaNS and AMPHIBIOUS at Re_θ ≈ 677; the outlet zone; the trip / resolution study behind the case (`smooth/tests_record.md`); the MATLAB export (`smooth/assets/matlab`) |
| [`rough_macdonald/`](rough_macdonald/README.md) | the same layer over an immersed egg-carton roughness after MacDonald et al. (JFM 2016), k⁺ ≈ 10, λ⁺ = 113, running through the outlet; Hama roughness function against the smooth case. Starts from the smooth case's developed field |
| `overheadTest/` | the performance / exchange campaign run on this case (HoreKa job scripts, profiles) -- not a tutorial |

Each tutorial directory carries its inis (from scratch: stage 0 / phase 1 /
phase 2), a `reproduce.py` that regenerates the committed reference data
(`assets/mobydiff/xyz_3200_384_136/data.nc`) and figures from the statistics,
and a provenance note. The statistics files, the last restart field and the
snapshots stay local (`.gitignore`; each over GitHub's 100 MB limit).
