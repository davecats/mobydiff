#!/usr/bin/env python3
"""Write two watertight wall-slab STLs for the IBM plane-wall channel.

The fluid gap is [y_lo, y_hi]; everything below y_lo and above y_hi is SOLID.
mobygeom's default convention is inside-the-STL = solid, so each slab is a closed
box. The slabs are padded past the x/z domain edges and past the outer y boundary
so the solid fully covers the one-cell ghost layer (no winding ambiguity there).

Walls are placed mid-cell (off grid node AND off cell centre): y_lo = 8.3*dy,
y_hi = 72.3*dy with dy = ly/ny -- so neither the y-face nodes nor the cell centres
land on the wall. Run with the geometry venv:
    /home/davide/ibmc/bin/python make_walls_stl.py
"""
import trimesh

LX = 12.566370614359172   # 4*pi
LY = 2.5
LZ = 6.283185307179586    # 2*pi
NY = 80
DY = LY / NY              # 0.03125

Y_LO = 8.3 * DY          # 0.259375  (bottom fluid/solid interface)
Y_HI = 72.3 * DY         # 2.259375  (top fluid/solid interface)

PAD = 0.25               # extend past every domain face (> a few dy)


def slab(y0, y1):
    """Closed box spanning the full (padded) x,z and y in [y0, y1]."""
    ex = [LX + 2 * PAD, y1 - y0, LZ + 2 * PAD]
    box = trimesh.creation.box(extents=ex)
    box.apply_translation([LX / 2.0, 0.5 * (y0 + y1), LZ / 2.0])
    return box


def main():
    lo = slab(-PAD, Y_LO)         # bottom solid: y in [-PAD, 0.259375]
    hi = slab(Y_HI, LY + PAD)     # top solid:    y in [2.259375, 2.75]
    lo.export("wall_lo.stl")
    hi.export("wall_hi.stl")
    print(f"wall_lo.stl: y in [{-PAD}, {Y_LO}]  watertight={lo.is_watertight}")
    print(f"wall_hi.stl: y in [{Y_HI}, {LY + PAD}]  watertight={hi.is_watertight}")
    print(f"fluid gap   : y in [{Y_LO}, {Y_HI}]  width={Y_HI - Y_LO}")


if __name__ == "__main__":
    main()
