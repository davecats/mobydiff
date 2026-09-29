#!/usr/bin/env python3
"""Convert a 1-block block-table restart field to LEGACY GLOBAL-3D layout so a
MULTI-RANK run can read it (e.g. the 4-rank HoreKa job).

    make_global_restart.py <src_blocktable.h5> <dst_global.h5>

Our restart fields store un/vn/wn/pn as 4D block-table datasets
(nBlocksGlobal, nz, ny, nx) plus a `blocks` table. The field reader
(field_hdf5.c fdm_h5_read_field) dispatches PER VARIABLE on dataset rank: 4D ->
block-table path (requires the run's block layout to MATCH the file's, so a
multi-rank run on a 1-block file fails), 3D -> legacy global path
(read_global_dataset_blocks: one hyperslab per block by origin, works on ANY rank
count for a level-0 uniform grid). So drop the leading block axis (1 block = whole
domain => un[0] is exactly the global (nz,ny,nx) field) and omit the `blocks`
dataset (=> the reader takes the legacy path). Grid lines + attrs are copied
verbatim. Only the FIRST restart needs this; the multi-rank run's own output is
already block-table matching its layout.
"""
import sys
import h5py
import numpy as np

src_path, dst_path = sys.argv[1], sys.argv[2]
with h5py.File(src_path, "r") as src, h5py.File(dst_path, "w") as dst:
    for k, v in src.attrs.items():
        dst.attrs[k] = v
    for name in ("x", "y", "z"):
        dst.create_dataset(name, data=np.array(src[name]))
    for name in ("un", "vn", "wn", "pn"):
        d = src[name]
        assert d.ndim == 4 and d.shape[0] == 1, f"{name} not 1-block 4D: {d.shape}"
        dst.create_dataset(name, data=d[0])
        print(f"  {name}: {d.shape} -> {dst[name].shape}")
    # deliberately NO 'blocks' dataset -> reader takes the legacy global path
print(f"wrote {dst_path} (legacy global-3D, multi-rank readable)")
