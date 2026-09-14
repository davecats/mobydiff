#!/usr/bin/env python3
"""Read ncu's `--page raw` CSV and say what limits each kernel.

    collect_ncu.py <results_dir>

`--page raw` is WIDE: a preamble, then one header row and one row per profiled
launch with ~280 metric columns. (The long "Metric Name/Metric Value" shape is
`--page details`; assuming it cost job 5144931's first collector run.)

`jacobi_compute_phi` is the control -- the kernel already known to sit near its
limit -- so the others are read against it, not against a spec sheet. Minimum
traffic is counted from the source in doubles per cell, each array counted ONCE
per cell however many neighbour accesses it has (the +1 neighbour is the next
cell's own value): 5 for compute_phi (3 velocity + rdenom read, phi write),
3 for apply k1 (phi read, p read+write), 8 for apply k2 (phi + mu read,
3 velocity read+write), 4 for compute_rdenom (3 mu read, rdenom write).

The MIN keys are matched against the kernel name with its `F1L<line>` part
REMOVED: nvfortran names a kernel after the source line of its target construct,
so keying on the line number makes the table silently lose its `min` column the
moment anything above it in the file moves (it did, in job 5145100). The
trailing ordinal (`_14_`, `_16_`) is the construct's index in the file and is
stable as long as no target region is added before it.
"""
import csv
import re
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

MIN = {"jacobi_compute_phi": 5, "compute_rdenom": 4,
       "jacobi_apply__14": 3, "jacobi_apply__16": 8}


def unline(name):
    """Kernel name with nvfortran's source-line stamp removed."""
    return re.sub(r"__F1L\d+_", "__", name)

COLS = [("gpu__time_duration.sum", "us", 1e-3),
        ("dram__bytes_read.sum", None, None), ("dram__bytes_write.sum", None, None),
        ("l1tex__average_t_sectors_per_request_pipe_lsu_mem_global_op_ld.ratio", "ld sect/req", 1),
        ("l1tex__average_t_sectors_per_request_pipe_lsu_mem_global_op_st.ratio", "st sect/req", 1),
        ("lts__t_sector_hit_rate.pct", "L2 hit %", 1),
        ("gpu__dram_throughput.avg.pct_of_peak_sustained_elapsed", "DRAM %peak", 1),
        ("sm__throughput.avg.pct_of_peak_sustained_elapsed", "SM %peak", 1),
        ("sm__warps_active.avg.pct_of_peak_sustained_active", "occupancy %", 1),
        ("launch__registers_per_thread", "registers", 1)]


def cells_of(run):
    # The solver's stdout lands in ncu.csv, not ncu.log -- ncu writes its CSV to
    # stdout and its own chatter to stderr, so the "block refinement: N leaves"
    # line is interleaved with the metric table. Search BOTH; getting this wrong
    # silently falls back to the level-0 grid and understates every per-cell
    # figure by the refinement ratio (1.75x here).
    txt = "".join((run / f).read_text(errors="replace")
                  for f in ("ncu.log", "ncu.csv") if (run / f).is_file())
    m = re.search(r"block refinement:\s+(\d+)\s+leaves", txt)
    cfg = (run / "config.ini").read_text()
    nb = re.search(r"^\s*nb\s*=\s*(.+)$", cfg, re.M)
    if m and nb:
        d = [int(x) for x in nb.group(1).split()]
        d = d * 3 if len(d) == 1 else d
        return int(m.group(1)) * d[0] * d[1] * d[2]
    g = [re.search(rf"^\s*n{x}\s*=\s*(\d+)", cfg, re.M) for x in "xyz"]
    return None if not all(g) else int(g[0].group(1))*int(g[1].group(1))*int(g[2].group(1))


def main():
    res = Path(sys.argv[1])
    print("# What limits `jacobi_apply`\n")
    prov = res / "provenance.txt"
    if prov.is_file():
        print("```\n" + prov.read_text().strip() + "\n```\n")

    for run in sorted(res.glob("ncu_*")):
        f = run / "ncu.csv"
        if not f.is_file():
            continue
        rows = list(csv.reader(f.read_text(errors="replace").splitlines()))
        hi = next((i for i, r in enumerate(rows) if r and r[0] == "ID"), None)
        if hi is None:
            print(f"_{run.name}: no header row -- see ncu.log_\n")
            continue
        hdr = rows[hi]
        idx = {k: i for i, k in enumerate(hdr)}
        data = [r for r in rows[hi + 1:] if len(r) == len(hdr)]
        agg = defaultdict(list)
        for r in data:
            agg[r[idx["Kernel Name"]]].append(r)

        def m(rs, key):
            vals = []
            for r in rs:
                try:
                    vals.append(float(r[idx[key]].replace(",", "")))
                except (ValueError, KeyError):
                    pass
            return st.mean(vals) if vals else None

        cells = cells_of(run)
        print(f"## {run.name[4:]}" + (f"  ({cells/1e6:.2f} Mcell)" if cells else "") + "\n")
        head = ["kernel", "launches", "us"] + [c[1] for c in COLS[3:]]
        if cells:
            head[3:3] = ["doubles/cell", "min", "vs min"]
        print("| " + " | ".join(head) + " |")
        print("|" + "---|" * len(head))
        for k, rs in sorted(agg.items(), key=lambda x: -len(x[1])):
            short = re.sub(r"^nvkernel_pressure_solver_", "", k).rstrip("_")
            us = m(rs, "gpu__time_duration.sum")
            if us is None:
                continue
            cellsrow = [f"`{short}`", str(len(rs)), f"{us/1e3:,.0f}"]
            if cells:
                rd, wr = m(rs, "dram__bytes_read.sum"), m(rs, "dram__bytes_write.sum")
                d = (rd + wr) / cells / 8 if rd is not None and wr is not None else None
                key = next((kk for kk in MIN if kk in unline(short)), None)
                cellsrow += [f"{d:.2f}" if d else "-",
                             str(MIN.get(key, "-")),
                             f"**{d/MIN[key]:.2f}x**" if (d and key) else "-"]
            for c in COLS[3:]:
                v = m(rs, c[0])
                cellsrow.append(f"{v:,.2f}" if v is not None else "-")
            print("| " + " | ".join(cellsrow) + " |")
        print()


if __name__ == "__main__":
    main()
