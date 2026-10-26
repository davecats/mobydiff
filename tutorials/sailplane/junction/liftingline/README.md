# Step 1: the lifting line (Route A) -- DONE 2026-10-26

The loading of the WHOLE FRUE V0 wing at the flight condition, from which the
truncated DNS box takes its one number, the fuselage-axis angle of attack
`alpha_box = alpha_flight - alpha_i(root region)` (handout
`docs/next_session_sailplane_junction.md`, sections 3 and 4.1). Decisions
applied: D2 C_L = 1.0, D3 Ncrit 9 and 12 as the band, D4 XFOIL built from
source (`../assets/xfoil/build.sh`). Tools: `../stl_sections.py` (sections),
`../xfoil_polars.py` (polars), `../liftingline.py` (AeroSandbox 4.2.10
nonlinear lifting line, NeuralFoil 0.3.3 section data, polar inversion on the
XFOIL polars). Re = 3e5 on the root chord 0.831 m (U = 5.27 m/s at sea
level), Re(y) = 3e5 c(y)/0.831 per section. Fuselage NOT modelled.

## Section polars (`polars/`, XFOIL 6.99 source build, 240 panels)

| section y [m] | chord | Re | Ncrit | C_l(0) | C_l(3) | C_l,max (alpha) | Top_Xtr(3) |
|---|---|---|---|---|---|---|---|
| 0.5 (root, t/c 0.144) | 0.831 | 3.00e5 | 9 / 12 | 0.690 / 0.695 | 0.991 / 0.998 | 1.32 / 1.33 (10) | 0.43 / 0.47 |
| 1.0 | 0.768 | 2.77e5 | 9 / 12 | 0.641 / 0.590 | 0.951 / 0.928 | 1.19 (5.5) / 1.23 (6) | 0.64 / 0.65 |
| 2.0 | 0.759 | 2.74e5 | 9 / 12 | 0.644 / 0.601 | 0.953 / 0.937 | 1.20 (10) / 1.23 (6) | 0.64 / 0.66 |
| 3.0 | 0.749 | 2.70e5 | 9 / 12 | 0.639 / 0.592 | 0.951 / 0.920 | 1.21 (10) / 1.23 (6) | 0.65 / 0.65 |
| 5.0 | 0.690 | 2.49e5 | 9 | 0.632 | 0.946 | 1.21 (10) | 0.65 |
| 7.0 | 0.598 | 2.16e5 | 9 / 12 | 0.621 / 0.488 | 0.911 / 0.813 | 1.19 (10) / 1.23 (6.5) | 0.63 / 0.63 |
| 9.0 | 0.413 | 1.49e5 | 9 / 12 | 0.538 / 0.315 | 0.808 / 0.689 | 1.18 / 1.21 (7) | 0.62 / 0.65 |

alpha to the chord line (the wing sits at +2.3 deg to the fuselage axis).
The y = 5.0 m polar at Ncrit 12 did not converge above alpha = 0 (8 points;
not needed). The thickened root section (y < 0.75 m) transitions at 0.43 c
at 3 deg; the wing sections proper carry a laminar separation bubble with
transition at 0.64-0.66 c and the lower surface laminar to the trailing
edge. The flap (0.87 c, about 6 deg down) is what puts C_l(0) at 0.6-0.7.

## The lifting line at C_L = 1.0

`liftingline.py SECTIONS --cl 1.0 --ncrit N --polars polars`
(`strips_cl1.0_ncrit{9,12}.csv`: y, c, c_l, alpha_i downwash, alpha_i polar,
alpha_eff per strip, right half):

| quantity | value |
|---|---|
| S (wing area, root chord carried to y = 0), AR | 12.75 m2, 29.9 |
| alpha_fuselage at C_L 1.0 (NeuralFoil sections) | 2.24 deg (chord line 4.54) |
| C_L,max of the lifting line | 1.11 at alpha 4 (C_L 1.2 unreachable) |
| root-region c_l, y = 0.55-1.05 m | 1.043 (1.06 at the thickened root, 1.04 at y = 1) |
| c_l at the cut, y = 2.95 m | 1.022: **4.2 % below the root**, under 1 % at y < 0.9 m |
| alpha_i root region, downwash at the bound vortex | 0.45 deg |
| alpha_eff root region from the XFOIL polars, Ncrit 9 / 12 | 3.76 / 3.81 deg (chord line) |
| **alpha_box, polar route (XFOIL), Ncrit 9 / 12** | **1.46 / 1.51 deg** |
| alpha_box, downwash route (NeuralFoil alpha_flight) | 1.79 deg |

**Result: alpha_box = 1.5 deg** (fuselage axis, nose-up; `stl_sections.py
rotate STL 1.5 OUT.stl`), band 1.45-1.8 deg. The polar route is the one to
take: `alpha_eff(c_l,root) - 2.3` needs only the loading (which at AR 30 is
set by the planform, not by the section model) and the XFOIL polar, whereas
the downwash route inherits the lifting line's own alpha at C_L 1.0, which
NeuralFoil reads 0.3-0.4 deg high on this flapped section (C_l 4-5 % low).
The Ncrit band is worth 0.05 deg. The 0.3 deg between the routes is the
NeuralFoil / XFOIL inconsistency and is covered by the +-0.5 deg sensitivity
runs of step 2; if step 2 shows the junction is sensitive at that level, the
refinement is a lifting line on the XFOIL polars themselves (a 100-line
own solver, or AeroSandbox's `Airfoil.generate_polars`).

For C_L 0.8 (not taken): alpha_fuselage 0.08 deg, alpha_box about -0.35 deg,
loading drop to the cut 7.5 %.

## What this settles for the next steps

- The 3 c cut stands: 4.2 % loading drop at the cut, < 1 % over the DNS
  core; step 2 still runs the 4 c prepare once as the check.
- The fin-only STL loses nothing (no horizontal tail; the lifting line does
  not need one at fixed C_L).
- The first reference for the DNS: the free-wing section at c_l 1.04 has its
  bubble transition at about 0.63 c (XFOIL, alpha_eff 3.8 deg,
  interpolated), the thickened root section at about 0.40 c.
