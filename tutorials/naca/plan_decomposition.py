#!/usr/bin/env python3
"""Decomposition planner for the Re 4e5 NACA campaign (2026-07-16 redesign):
computes per-level block/cell counts and draws the nested-box layout BEFORE
anything is generated or run.

Two configurations:
  A "requested": 8 refinement levels, finest box = profile + 0.5c all
    around at y+~2 resolution, coarser boxes at +1c (upstream/lateral)
    / +2c (downstream, wake skew) per level.
  B "proposal": same domain/nose/y+ target with 10 refinement levels
    (coarser base -> far field cheap), snug distance-band fine levels
    (refine_body-style, 1-block rings), and the requested generous
    wake-skewed boxes from L6 downward (0.5c at L6, +1c/level, x2
    downstream).

Cells include the span (ny = 8). Memory estimate: ~0.42 GB / M cells
(calibrated on the L5-3D 33M-cell run and the 121M-cell OOM at 49 GB).
"""
import numpy as np

DOM_X, DOM_Z = 126.0, 100.0
NOSE = (50.0, 50.0)
CHORD = 1.0
NY = 8
NB = 8
CF = 0.008
RE = 4.0e5
BODY_AREA = 0.0817          # NACA 0012
PERIM = 2.06                # section perimeter


def yplus(delta):
    utau = np.sqrt(CF/2.0)
    return 0.5*delta*utau*RE   # first cell CENTRE


def clamp(box):
    x0, x1, z0, z1 = box
    return (max(0.0, x0), min(DOM_X, x1), max(0.0, z0), min(DOM_Z, z1))


def area(box):
    x0, x1, z0, z1 = box
    return max(0.0, x1-x0)*max(0.0, z1-z0)


def report(name, levels, d0):
    """levels: list of (level_index, region_area, note). Coarser first."""
    print(f"\n=== {name} ===")
    print(f"base grid: {int(DOM_X/d0)} x {NY} x {int(DOM_Z/d0)} "
          f"(Delta0 = {d0:.4e} c = c/{1/d0:.0f})")
    print(f"{'lvl':>3s} {'Delta':>10s} {'y+(D/2)':>8s} {'plane Mcells':>12s} "
          f"{'blocks':>9s}  note")
    tot_plane = 0.0
    for k, A, note in levels:
        dk = d0/2**k
        cells = A/dk**2
        blocks = cells/NB**2
        tot_plane += cells
        print(f"{k:3d} {dk:10.3e} {yplus(dk):8.2f} {cells/1e6:12.3f} "
              f"{blocks:9.0f}  {note}")
    total = tot_plane*NY
    print(f"TOTAL: {tot_plane/1e6:.1f} M plane cells x ny={NY} = "
          f"{total/1e6:.1f} M cells  (~{total/1e6*0.42:.0f} GB device)")
    return total


# ---------------- Config A: as requested ----------------
LA = 8
d0A = CHORD/24.0            # Delta8 = c/6144 -> y+ ~ 2.06
x0, z0 = NOSE
boxA = {}
b = (x0-0.5, x0+CHORD+0.5, z0-0.56, z0+0.56)   # finest: profile + 0.5c
boxA[LA] = b
for k in range(LA-1, 0, -1):
    m = LA - k
    boxA[k] = clamp((b[0]-1.0*m, b[1]+2.0*m, b[2]-1.0*m, b[3]+1.0*m))
levelsA = []
levelsA.append((0, DOM_X*DOM_Z - area(boxA[1]), "far field (whole domain)"))
for k in range(1, LA):
    levelsA.append((k, area(boxA[k]) - area(boxA[k+1]),
                    f"box {tuple(round(v,2) for v in boxA[k])}"))
levelsA.append((LA, area(boxA[LA]) - BODY_AREA,
                f"box {tuple(round(v,2) for v in boxA[LA])} (0.5c margin)"))
totA = report("A: requested (8 levels, 0.5c finest box, +1c/+2c margins)",
              levelsA, d0A)

# ---------------- Config B: proposal ----------------
LB = 10
d0B = CHORD/6.0             # Delta10 = c/6144 -> same y+ ~ 2.06
dB = [d0B/2**k for k in range(LB+1)]
# fine levels 10..7: snug distance bands, ~3 blocks of their own level thick
# (touch + buffer + 2:1 ring), both sides of the surface
levelsB = []
band_blocks = 3.0
fine_areas = {}
for k in range(LB, 6, -1):
    t = band_blocks*NB*dB[k]                 # band thickness per side
    fine_areas[k] = PERIM*2*t + (BODY_AREA if k == LB else 0.0)*0.0
levelsB_f = [(k, fine_areas[k],
              f"snug band, ~{band_blocks:.0f} blocks/side "
              f"({band_blocks*NB*dB[k]:.4f} c)") for k in range(LB, 6, -1)]
# boxes from L6 down: 0.5c at L6, +1c per level upstream/lateral,
# +2c downstream (wake skew)
boxB = {6: (x0-0.5, x0+CHORD+0.5, z0-0.56, z0+0.56)}
for k in range(5, 0, -1):
    m = 6 - k
    b6 = boxB[6]
    boxB[k] = clamp((b6[0]-1.0*m, b6[1]+2.0*m, b6[2]-1.0*m, b6[3]+1.0*m))
levelsB.append((0, DOM_X*DOM_Z - area(boxB[1]), "far field (whole domain)"))
for k in range(1, 6):
    levelsB.append((k, area(boxB[k]) - area(boxB[k+1] if k < 5 else boxB[6]),
                    f"box {tuple(round(v,2) for v in boxB[k])}"))
levelsB.append((6, area(boxB[6]) - BODY_AREA - sum(fine_areas.values()),
                f"box {tuple(round(v,2) for v in boxB[6])} (0.5c margin)"))
levelsB += levelsB_f
totB = report("B: proposal (10 levels, snug y+2 bands, 0.5c box at L6, "
              "+1c/+2c margins below)", levelsB, d0B)

print(f"\nA / B cell ratio: {totA/totB:.0f}x")
print(f"largest available GPU: RTX 5090 32 GB ~ 75 M cells")

# ---------------- figure ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

fig, axs = plt.subplots(1, 3, figsize=(17, 6))
for ax, boxes, Lf, ttl in ((axs[0], boxA, LA, "A: requested (boxes L1..L8)"),
                           (axs[1], boxB, 6, "B: proposal (boxes L1..L6 + snug bands)")):
    ax.add_patch(Rectangle((0, 0), DOM_X, DOM_Z, fc="none", ec="k", lw=1.2))
    cmap = plt.cm.viridis
    for k in sorted(boxes):
        x0b, x1b, z0b, z1b = boxes[k]
        ax.add_patch(Rectangle((x0b, z0b), x1b-x0b, z1b-z0b,
                               fc="none", ec=cmap(k/Lf), lw=1.4))
    ax.plot([NOSE[0], NOSE[0]+1], [NOSE[1], NOSE[1]], "r-", lw=2)
    ax.set_xlim(-3, DOM_X+3); ax.set_ylim(-3, DOM_Z+3)
    ax.set_aspect("equal"); ax.set_title(ttl)
    ax.set_xlabel("x/c"); ax.set_ylabel("z/c")
ax = axs[2]
for k in sorted(boxB):
    x0b, x1b, z0b, z1b = boxB[k]
    ax.add_patch(Rectangle((x0b, z0b), x1b-x0b, z1b-z0b,
                           fc="none", ec=plt.cm.viridis(k/6), lw=1.6))
# indicate the snug bands
th = np.linspace(0, 2*np.pi, 200)
for k in range(7, LB+1):
    t = band_blocks*NB*dB[k]
    ax.plot(NOSE[0]+0.5+0.56*np.cos(th)*(1+t), NOSE[1]+0.075*np.sin(th)+0,
            lw=0.5, color=plt.cm.plasma((k-7)/4), alpha=0.0)  # schematic only
ax.plot([NOSE[0], NOSE[0]+1], [NOSE[1], NOSE[1]], "r-", lw=3)
ax.annotate("L7..L10 snug bands\n(~3 blocks/side each,\n0.031/0.016/0.008/0.004 c)",
            xy=(NOSE[0]+0.5, NOSE[1]+0.1), xytext=(NOSE[0]-4.5, NOSE[1]+4),
            arrowprops=dict(arrowstyle="->"), fontsize=9)
ax.set_xlim(NOSE[0]-8, NOSE[0]+16); ax.set_ylim(NOSE[1]-9, NOSE[1]+9)
ax.set_aspect("equal"); ax.set_title("B: zoom (L1..L6 boxes; bands at the body)")
ax.set_xlabel("x/c"); ax.set_ylabel("z/c")
fig.tight_layout()
fig.savefig("decomposition_plan.png", dpi=130)
print("wrote decomposition_plan.png")
