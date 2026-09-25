#!/usr/bin/env python3
"""Combine statistics files from consecutive runs into one window.

The solver's accumulator is cumulative from a run's FIRST step, and it stores
both `raw_sum` and `count`, so two consecutive legs combine exactly by summing
those and re-forming the mean -- not by averaging the two `profile` tables,
which would weight a short leg like a long one.

    ./combine_stats.py out.h5 leg1.h5 leg2.h5 [...]
"""
import shutil, sys
import h5py
import numpy as np

def main():
    out, first, *rest = sys.argv[1:]
    shutil.copyfile(first, out)
    with h5py.File(out, "r+") as d:
        rs, cnt = d["raw_sum"][...], d["count"][...]
        n0 = float(cnt.sum())
        for p in rest:
            with h5py.File(p, "r") as s:
                if s["raw_sum"].shape != rs.shape:
                    raise SystemExit(f"{p}: shape {s['raw_sum'].shape} != {rs.shape}")
                rs += s["raw_sum"][...]
                cnt += s["count"][...]
                d.attrs["step"] = int(s.attrs["step"])
                d.attrs["t_current"] = float(s.attrs["t_current"])
        d["raw_sum"][...] = rs
        d["count"][...] = cnt
        safe = np.where(cnt > 0, cnt, 1.0)
        d["profile"][...] = rs / safe[:, None]
        print(f"   {out}: {len(rest)+1} legs, sample weight {n0:.0f} ->"
              f" {float(cnt.sum()):.0f} (x{cnt.sum()/max(n0,1):.2f})")

if __name__ == "__main__":
    main()
