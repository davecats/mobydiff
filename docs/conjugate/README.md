# Conjugate heat transfer at an immersed interface — derivation note

`conjugate_ibm.tex` is the long-form companion to
[`../next_session_conjugate.md`](../next_session_conjugate.md): the full
derivation, the limit checks, the accuracy argument and the sketches. The
session document is the actionable plan (increments, gates, files, landmines);
this one is why the scheme looks the way it does.

Build:

```bash
make            # -> conjugate_ibm.pdf   (needs pdflatex + latexmk + tikz)
make clean      # drop the LaTeX intermediates, keep the pdf
```

Contents:

**In a hurry?** Read §10 (can the problem be avoided?), then §6 (what to
implement), then §7. The rest is why those two terms are the right ones.

| section | content |
|---|---|
| 1 | CHT equations, the two material ratios `κ = k/k_f` and `C = ρc/(ρc)_f`, how to read the note, geometry and notation (Fig. 1) |
| 2 | The three source methods: the Luchini λ correction and its implicit time treatment (Fig. 2); EJIIM's jump-corrected differences; Cipelli's corner correction (COCO) — and why none alone suffices |
| 3 | The bridge — an explicit RK3 march collapses EJIIM's augmented system to a local evaluation |
| 4 | Derivation of the exact cut-face flux (Eq. 24), its limits (Table 1), why Dirichlet needs no tangential term and conjugate does (Fig. 4), accuracy, contact resistance |
| 5 | The COCO route: two complementary 180° wedges give the same flux without a Taylor expansion (Fig. 5a); the parameter count that makes a single precomputed λ impossible for CHT (Table 2); the two traps — vanishing/sign-changing ratio denominator, and point- vs face-anchoring (Fig. 6); conducting sharp corners (Fig. 5b) |
| **6** | **The DNS baseline — what to implement.** The obliquity lemma: the level-set fraction `w = φ_L/(φ_L−φ_R)` of the two signed distances is exactly the arm fraction the series resistance needs, at any orientation (Fig. 7). One face coefficient, no new stored data (Table 3); the closed-form error estimate and its run-time indicator; curvature and grazing arms; the escalation path back to Eq. 24 |
| 7 | The discrete scheme: face loop, convection masking, the tangential gradient (escalation only, Fig. 8), cut-cell capacity (Fig. 9) |
| 8 | Stability: the non-stiffness result (Fig. 10) and the solid-diffusivity time-step limit |
| 9 | What `dwall` is, and the two arrangements needed around it — no new datasets |
| 10 | Two cheaper routes to check first: grid-aligned interface; Robin thin-wall model |
| 11 | Alternatives considered and rejected |
| 12 | Validation ladder (the gates are stated in full in the session document) |

Source papers, in `literature/`: `luchini-ibm-2025.pdf` (JCP 539:114245),
`interface-method.pdf` (Wiegmann & Bube, SIAM J. Numer. Anal. 37(3):827–862)
and `cipelli-2025.pdf` (Cipelli, Chiarini, Quadrio, Gatti & Luchini, riblet
corner correction).
