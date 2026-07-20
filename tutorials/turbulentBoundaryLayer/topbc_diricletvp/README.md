# topbc_diricletvp — Dirichlet p + prescribed v (an over-specified top)

The top (`y_max`) sets **Dirichlet p = 0** (the outlet top's
freestream-pressure pin) together with the wall-normal v **prescribed**
from the Blasius entrainment and Neumann u — an attempt to pin ZPG AND
impose the entrainment at the same time. See the parent `../README.md` for
the shared setup and run recipe.

## Why this is over-specified — and what the solver does

At an incompressible boundary the normal velocity and the pressure are a
**conjugate pair**: you specify one or the other, not both. Prescribing v
(Dirichlet) already closes the problem there; adding Dirichlet p is
redundant. The solver resolves the redundancy by **ignoring the Dirichlet
p**:

- The projection activates the pressure pin (the outlet-style face
  correction, mirrored-phi ghost) ONLY for a declared `outlet` patch. On
  this generic face it treats a physical boundary as pinned-velocity /
  Neumann-phi (`face_grad` returns 0), regardless of `p_type`.
- The Dirichlet p only changes the top pressure GHOST via `apply_bc`
  (ghost = -interior instead of the Neumann copy), but that ghost is never
  read: the top v-face is pinned (not advanced by the momentum predictor),
  and the projection's divergence / velocity correction never touch it.

So `y_max_p_type = dirichlet` here is a **no-op**. This case is
**bit-identical to `../topbc_displacement`** (Neumann p) — verified with
`compare_fields.py` (max_abs 0), and the gate numbers match to every digit
(theta +0.3/+3.1/+6.9/+9.1%, U_e drifting to 1.011).

## Takeaway

You cannot pin the freestream pressure by setting `p_type = dirichlet` on a
non-outlet face — the pin lives in the `outlet` machinery, which in turn
frees v (the projection determines the outlet-face velocity). So the two
goals are mutually exclusive *within this face type*:

- want the ZPG pressure pin  -> declare `outlet` (v self-selects:
  `../topbc_outlet`, the best Blasius match);
- want the prescribed entrainment -> Dirichlet v + Neumann p
  (`../topbc_displacement`; freestream drifts).

A genuine "pin both" is not a boundary-condition choice but a solver
feature: it would need an outlet-style pressure pin that also accepts a
prescribed (rather than projection-determined) normal velocity — i.e. a
new patch type, not a config combination. Absent that, the **outlet top
remains the recommended ZPG boundary layer top**.

## Figures

- `blasius.png`, `blasius2d.png` — identical to `../topbc_displacement`.
