#!/usr/bin/env python3
"""Read ncu's raw CSV and say what each kernel's memory traffic actually is.

    collect_ncu.py <results_dir>

The question is not "what is the peak bandwidth" but "how many bytes per CELL
does this kernel move, against how many it needs". `jacobi_compute_phi` is the
control: it is already known to run at its bandwidth limit, so its measured
bytes/cell is what a well-behaved kernel of this shape costs, and the others are
read against it.

Minimum traffic is counted from the source, in doubles per cell:
  compute_phi   3 velocity arrays read + rdenom read + phi write        = 5
  apply k1      phi read + p read + p write                             = 3
  apply k2      phi read + mu read + 3 velocity read + 3 velocity write = 8
Neighbour accesses are not counted separately: they land in the same cache lines
as the centre value for the i-direction and are served by L2 for j and k, so a
kernel that needs 8 and moves 8 is perfect, one that moves 16 is wasting half.
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

MIN_DOUBLES = {"jacobi_compute_phi": 5, "jacobi_apply__F1L470": 3,
               "jacobi_apply__F1L502": 8}
WANT = ("dram__bytes_read.sum", "dram__bytes_write.sum",
        "gpu__dram_throughput.avg.pct_of_peak_sustained_elapsed",
        "sm__throughput.avg.pct_of_peak_sustained_elapsed",
        "l1tex__average_t_sectors_per_request_pipe_lsu_mem_global_op_ld.ratio",
        "l1tex__average_t_sectors_per_request_pipe_lsu_mem_global_op_st.ratio",
        "lts__t_sector_hit_rate.pct",
        "sm__warps_active.avg.pct_of_peak_sustained_active")
SHORT = {"dram__bytes_read.sum": "DRAM read", "dram__bytes_write.sum": "DRAM write",
         "gpu__dram_throughput.avg.pct_of_peak_sustained_elapsed": "DRAM %peak",
         "sm__throughput.avg.pct_of_peak_sustained_elapsed": "SM %peak",
         "l1tex__average_t_sectors_per_request_pipe_lsu_mem_global_op_ld.ratio": "sect/req ld",
         "l1tex__average_t_sectors_per_request_pipe_lsu_mem_global_op_st.ratio": "sect/req st",
         "lts__t_sector_hit_rate.pct": "L2 hit %",
         "sm__warps_active.avg.pct_of_peak_sustained_active": "occupancy %"}


def cells_of(run):
    """Cells per rank: leaves x nb product, from the run's own log."""
    log = next((p for p in (run / "ncu.log", run / "run.log") if p.is_file()), None)
    txt = log.read_text(errors="replace") if log else ""
    m = re.search(r"block refinement:\s+(\d+)\s+leaves", txt)
    nb = re.search(r"^\s*nb\s*=\s*(.+)$", (run / "config.ini").read_text(), re.M)
    if not (m and nb):
        return None
    d = [int(x) for x in nb.group(1).split()]
    if len(d) == 1:
        d = d * 3
    return int(m.group(1)) * d[0] * d[1] * d[2]


def main():
    res = Path(sys.argv[1])
    print("# Why `jacobi_apply` costs 3x `jacobi_compute_phi` per cell\n")
    prov = res / "provenance.txt"
    if prov.is_file():
        print("```\n" + prov.read_text().strip() + "\n```\n")

    for run in sorted(res.glob("ncu_*")):
        csvf = run / "ncu.csv"
        if not csvf.is_file():
            continue
        rows = [r for r in csv.DictReader(
            l for l in csvf.read_text(errors="replace").splitlines()
            if l and not l.startswith("=="))]
        if not rows:
            print(f"_{run.name}: no metric rows_\n")
            continue
        kn = next(k for k in rows[0] if "Kernel Name" in k)
        mn = next(k for k in rows[0] if k.strip() in ("Metric Name", "Metric"))
        mv = next(k for k in rows[0] if k.strip() in ("Metric Value", "Value"))
        agg = defaultdict(lambda: defaultdict(list))
        for r in rows:
            try:
                agg[r[kn]][r[mn]].append(float(str(r[mv]).replace(",", "")))
            except ValueError:
                pass
        cells = cells_of(run)
        print(f"## {run.name[4:]}"
              + (f"  ({cells/1e6:.2f} Mcell/rank)" if cells else "") + "\n")
        print("| kernel | " + " | ".join(SHORT[w] for w in WANT) + " |")
        print("|" + "---|" * (len(WANT) + 1))
        for k, m in agg.items():
            short = re.sub(r"^nvkernel_pressure_solver_", "", k).rstrip("_")
            cellsline = []
            for w in WANT:
                v = m.get(w)
                cellsline.append(f"{sum(v)/len(v):,.2f}" if v else "-")
            print(f"| `{short}` | " + " | ".join(cellsline) + " |")
        if not cells:
            continue
        print("\n### Bytes per cell, against what the kernel needs\n")
        print("| kernel | DRAM B/cell | doubles/cell | minimum | **waste** |")
        print("|" + "---|" * 5)
        for k, m in agg.items():
            rd, wr = m.get("dram__bytes_read.sum"), m.get("dram__bytes_write.sum")
            if not (rd and wr):
                continue
            b = (sum(rd) / len(rd) + sum(wr) / len(wr)) / cells
            key = next((kk for kk in MIN_DOUBLES if kk in k), None)
            mind = MIN_DOUBLES.get(key)
            short = re.sub(r"^nvkernel_pressure_solver_", "", k).rstrip("_")
            print(f"| `{short}` | {b:,.1f} | {b/8:.1f} | "
                  + (f"{mind}" if mind else "-") + " | "
                  + (f"**{b/8/mind:.2f}x**" if mind else "-") + " |")
        print()


if __name__ == "__main__":
    main()
