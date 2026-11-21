# Zero-pressure-gradient turbulent boundary layer

Spatially developing ZPG turbulent boundary layers, DNS on the CaNS reference
grid (650 × 100 × 26 δ*₀, 3200 × 384 × 136, Re_δ*,0 = 450, tripped Blasius
inflow, conservative convection, red-black SOR projection).

| directory | case |
|---|---|
| [`smooth/`](smooth/README.md) | smooth wall: the code comparison against SIMSON, CaNS and AMPHIBIOUS at Re_θ ≈ 677; the outlet zone; the trip / resolution study behind the case (`smooth/tests_record.md`); the MATLAB export (`smooth/assets/matlab`) |
| [`rough_macdonald/`](rough_macdonald/README.md) | the same layer over an immersed egg-carton roughness after MacDonald et al. (JFM 2016), k⁺ ≈ 10, λ⁺ = 113, running through the outlet; Hama roughness function against the smooth case. Starts from the smooth case's developed field |
| [`smooth.RANS/`](smooth.RANS/README.md) | the smooth layer as a steady 2D wall-resolved RANS (k-ω SST, Spalart–Rumsey free stream): laminar to the DNS trip (k pinned to 0), k sources where the DNS trips and transitions, turbulent downstream; c_f within 3–6 % of the DNS from Re_θ 450, Re_θ(x) within 2 %; the freestream-ω, transition-model and trip calibration |
| [`apg.RANS/`](apg.RANS/README.md) | the smooth.RANS layer under an adverse pressure gradient at constant Clauser parameter β = +0.5 from x = 200 (Re_θ 500 → 1075, no separation): the top outlet's pressure is a table found by an optimisation loop around the solver (`top_optimise.py`); the procedure, the Griffin edge on the local free-stream P0, the RANS → DNS β ratio (r ≈ 0.92) |
| [`fpg.RANS/`](fpg.RANS/README.md) | the favourable twin, β = −0.3 from x = 200 (Re_θ 400 → 700, K ≤ 1.3·10⁻⁶) |
| [`apg.DNS/`](apg.DNS/README.md) | the DNS of apg.RANS (Lz = 52, nz 272, the ZPG trip spectrum on the doubled span): the RANS top-pressure table × 1.087 gives β = +0.497 on x 200–550 (0.465 → 0.527, Re_θ 500 → 1020) over 27 000 t.u. on 8 A100; the RANS → DNS transfer measured (dP_e/dx exact, the SST error in τ_w: r ≈ 0.84 under the APG against 0.92 at ZPG), the RANS c_f 10–11 % low against the DNS |
| `overheadTest/` | the performance / exchange campaign run on this case (HoreKa job scripts, profiles) -- not a tutorial |

Each tutorial directory carries its inis (from scratch: stage 0 / phase 1 /
phase 2), a `reproduce.py` that regenerates the committed reference data
(`assets/mobydiff/xyz_3200_384_136/data.nc`) and figures from the statistics,
and a provenance note. The statistics files, the last restart field and the
snapshots stay local (`.gitignore`; each over GitHub's 100 MB limit).
