#!/usr/bin/env python3
"""Put the scalar fields of one snapshot into another's velocity field.

The fine grid's velocity is adapted first WITHOUT scalars (its time step is
then convection-limited, ~5x looser than the conjugate Peclet limit), which
leaves the adapted snapshot with no scalar datasets. This copies them back in
from the interpolated initial condition, giving one restart that carries the
adapted velocity and the interpolated temperature fields.

    ./merge_scalars.py <velocity snapshot> <scalar source> <out>
"""
import shutil, sys
import h5py

def main():
    vel, sca, out = sys.argv[1:4]
    shutil.copyfile(vel, out)
    with h5py.File(sca, "r") as s, h5py.File(out, "r+") as d:
        keep = {"un", "vn", "wn", "pn", "blocks", "x", "y", "z"}
        for k, v in s.items():
            if isinstance(v, h5py.Dataset) and v.ndim == 4 and k not in keep:
                if k in d:
                    del d[k]
                d.create_dataset(k, data=v[...])
                print(f"   {k}: copied")
    print(f"   wrote {out}")

if __name__ == "__main__":
    main()
