# Conjugate heat transfer — validation against reference DNS

Two campaigns, deliberately complementary, each validating mobydiff's
conjugate heat-transfer scheme (`[scalar.N] ibm_wall = conjugate`) against a
published DNS.

| | [`channel/`](channel/) | [`pipe/`](pipe/) |
|---|---|---|
| reference | Flageul et al. 2015, body-fitted DNS | Neuhauser, NekRS spectral-element DNS |
| interface | **flat, grid-aligned** — level-set weight w = ½ exactly | **curved** — arbitrary cut-cell weights, staircase wall |
| what it isolates | the conjugate **physics**: the face coefficient is the exact harmonic mean here | the **immersed-boundary approximation** on top of that physics |
| Re_τ, Pr | 149, 0.71 | 181, 0.71 |
| sweep | 6 scalars: κ_s and C_s each over four decades | 7 scalars: the effusivity bracket from isothermal to isoflux |
| headline | near-wall variance peak **+0.7 %**, full profiles 0.4–3.8 % | conjugate signature **≤ 2 %** across the whole bracket, conductivity pair **0.1 %** |

Read `channel/` first: it establishes that the scheme gets conjugate
turbulence right where the discretisation is exact, which is what makes the
pipe's residuals attributable to the curved interface rather than to the
physics.

Each directory holds what you need to **run** the case; its `asset/` holds the
reference data, the comparison scripts, the figures and the full report.

The feature-level gates behind both — the manufactured two-material slabs, the
C1–C3 convergence and conservation tests, the pipe geometry gates — live in
`validation/conjugate/`, not here.
