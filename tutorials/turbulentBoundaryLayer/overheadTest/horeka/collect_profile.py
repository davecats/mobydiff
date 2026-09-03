#!/usr/bin/env python3
"""Collect the 2:1 timing/profiling matrix into one markdown report.

    collect_profile.py <results_dir> [> summary.md]

Reads every <results_dir>/<config>_n<ranks>/ produced by run_matrix.sh and
reports, per run: seconds per step, seconds per million cells (the only figure
comparable ACROSS configs, since the configs differ in size), the per-phase
breakdown, and the exchange traffic report.

Standard library only, on purpose: this runs on a compute node at the end of a
batch job and must not fail for a missing module.

The headline numbers it derives:

  * 2:1 overhead   -- refined_big_rect_jacobi vs rect_jacobi at equal rank
                      count, per cell. Same level-0 grid, same block shape,
                      same solver: the ratio is the refinement machinery alone.
  * block tax      -- rect_jacobi vs base_jacobi per cell (blocking, no
                      refinement).
  * strong scaling -- per config, s/step against its own 1-rank (or smallest)
                      run, with the phase that loses the time named.
"""

import re
import sys
from pathlib import Path

TIMING = re.compile(
    r"^(step_timing|proj_timing|exch_timing):\s+(\S+)\s+calls\s+\d+\s+nsteps\s+(\d+)"
    r"\s+seconds\s+(\S+)\s+seconds_per_step\s+(\S+)")
CHRON = re.compile(r"^timing:\s+nsteps\s+(\d+)\s+loop_seconds\s+(\S+)\s+seconds_per_step\s+(\S+)")
LEAVES = re.compile(r"block refinement:\s+(\d+)\s+leaves,\s+(\d+)\s+refined")
SIZES = re.compile(
    r"exchange sizes: peers/rank\(max\)\s+(\d+)\s+send pts/rank min\s+(\d+)\s+max\s+(\d+)"
    r"\s+total send pts\s+(\d+)\s+local copy pts\s+(\d+)")
MB = re.compile(r"exchange MB per round .*?scalar nv=1\s+(\S+)\s+copy-only nv=3\s+(\S+)\s+full nv=4\s+(\S+)")


def read_config(path):
    """nb / grid / solver / niter from the config that the run actually used."""
    cfg = {}
    if not path.exists():
        return cfg
    for line in path.read_text().splitlines():
        line = line.split(";")[0].strip()
        if "=" not in line:
            continue
        k, v = (s.strip() for s in line.split("=", 1))
        cfg[k] = v
    return cfg


def cells_of(cfg, leaves):
    """Total cells. With [blocks] nb the grid is a leaf lattice of equal-size
    blocks, so cells = leaves x nb_x nb_y nb_z -- NOT nx ny nz, which counts
    the level-0 grid only and undercounts every refined case."""
    try:
        nx, ny, nz = (int(cfg[k]) for k in ("nx", "ny", "nz"))
    except (KeyError, ValueError):
        return None
    nb = cfg.get("nb")
    if nb and leaves:
        p = [int(x) for x in nb.replace(",", " ").split()]
        if len(p) == 1:
            p = p * 3
        return leaves * p[0] * p[1] * p[2]
    return nx * ny * nz


def parse_run(d):
    log = d / "run.log"
    if not log.exists():
        return None
    r = {"name": d.name, "phases": {}, "leaves": None, "sizes": None, "mb": None,
         "sps": None, "nsteps": None}
    for line in log.read_text(errors="replace").splitlines():
        m = CHRON.match(line.strip())
        if m:
            r["nsteps"], r["sps"] = int(m.group(1)), float(m.group(3))
            continue
        m = TIMING.match(line.strip())
        if m:
            grp, phase, _, _, sps = m.groups()
            if phase != "total_measured":
                r["phases"][f"{grp.split('_')[0]}/{phase}"] = float(sps)
            continue
        m = LEAVES.search(line)
        if m:
            r["leaves"] = int(m.group(1)); continue
        m = SIZES.search(line)
        if m:
            r["sizes"] = tuple(int(x) for x in m.groups()); continue
        m = MB.search(line)
        if m:
            r["mb"] = tuple(float(x) for x in m.groups())
    cfg = read_config(d / "config.ini")
    r["cfg"] = cfg
    r["cells"] = cells_of(cfg, r["leaves"])
    parts = d.name.rsplit("_n", 1)
    r["config"], r["ranks"] = parts[0], int(parts[1]) if len(parts) == 2 else 0
    return r


def fmt(x, w=9, p=4):
    return "-".rjust(w) if x is None else f"{x:{w}.{p}f}"


def main():
    res = Path(sys.argv[1] if len(sys.argv) > 1 else "results")
    runs = [r for r in (parse_run(d) for d in sorted(res.iterdir()) if d.is_dir()) if r]
    runs = [r for r in runs if r["sps"]]
    if not runs:
        print("No completed runs found in", res)
        # A failed matrix must be visible, not silently empty.
        for d in sorted(res.glob("*/run.FAILED.log")):
            print("  FAILED:", d.parent.name)
        return 1

    prov = res / "provenance.txt"
    print("# 2:1 refinement timing + profiling matrix\n")
    if prov.exists():
        print("```"); print(prov.read_text().rstrip()); print("```\n")

    print("## Runs\n")
    print("| run | ranks | leaves | cells (M) | s/step | s per M cells |")
    print("|---|---|---|---|---|---|")
    for r in sorted(runs, key=lambda r: (r["config"], r["ranks"])):
        cM = r["cells"] / 1e6 if r["cells"] else None
        per = r["sps"] / cM if cM else None
        print(f"| {r['config']} | {r['ranks']} | {r['leaves'] or '-'} | "
              f"{fmt(cM,7,2)} | {fmt(r['sps'],8,5)} | {fmt(per,9,6)} |")

    # ---- the headline ratios, per rank count ------------------------------
    by = {(r["config"], r["ranks"]): r for r in runs}

    def per_cell(cfg, n):
        r = by.get((cfg, n))
        if not r or not r["cells"]:
            return None
        return r["sps"] / (r["cells"] / 1e6)

    print("\n## Headline ratios (cost per cell, equal rank count)\n")
    print("`refined_big_rect_jacobi` and `rect_jacobi` share level-0 grid, block")
    print("shape and solver, so their ratio is the 2:1 machinery alone.\n")
    print("| ranks | 2:1 overhead (big/rect) | block tax (rect/base) | small refined/rect |")
    print("|---|---|---|---|")
    for n in sorted({r["ranks"] for r in runs}):
        def ratio(a, b):
            x, y = per_cell(a, n), per_cell(b, n)
            return f"{x/y:.3f}" if x and y else "-"
        print(f"| {n} | {ratio('refined_big_rect_jacobi','rect_jacobi')} "
              f"| {ratio('rect_jacobi','base_jacobi')} "
              f"| {ratio('refined_yp82_rect_jacobi','rect_jacobi')} |")

    # ---- strong scaling ---------------------------------------------------
    print("\n## Strong scaling (per config; efficiency vs the smallest rank count)\n")
    for cfg in sorted({r["config"] for r in runs}):
        rs = sorted([r for r in runs if r["config"] == cfg], key=lambda r: r["ranks"])
        if len(rs) < 2:
            continue
        base = rs[0]
        print(f"\n### {cfg}  (reference: {base['ranks']} rank(s), {base['sps']:.5f} s/step)\n")
        print("| ranks | s/step | speedup | efficiency | exchange s/step | exch % |")
        print("|---|---|---|---|---|---|")
        for r in rs:
            speed = base["sps"] / r["sps"]
            eff = speed / (r["ranks"] / base["ranks"])
            exch = sum(v for k, v in r["phases"].items() if k.startswith("exch/"))
            print(f"| {r['ranks']} | {r['sps']:.5f} | {speed:.2f}x | {eff*100:.1f}% "
                  f"| {exch:.5f} | {exch/r['sps']*100:.1f}% |")

    # ---- phases -----------------------------------------------------------
    print("\n## Per-phase seconds per step\n")
    keys = sorted({k for r in runs for k in r["phases"]})
    hdr = sorted(runs, key=lambda r: (r["config"], r["ranks"]))
    print("| phase | " + " | ".join(f"{r['config']}/n{r['ranks']}" for r in hdr) + " |")
    print("|---" * (len(hdr) + 1) + "|")
    for k in keys:
        vals = " | ".join(f"{r['phases'].get(k, 0.0):.5f}" for r in hdr)
        print(f"| {k} | {vals} |")

    # ---- exchange traffic -------------------------------------------------
    print("\n## Exchange traffic and balance\n")
    print("`send pts/rank min..max` is the number to read on any new")
    print("decomposition: MPI_Waitall pays for the slowest peer, so a wide")
    print("spread is a straggler even when every rank owns the same leaf count.\n")
    print("| run | ranks | peers/rank | send pts min | max | spread | total send | local copy | MB/round nv=4 |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in hdr:
        if not r["sizes"]:
            continue
        peers, smin, smax, stot, loc = r["sizes"]
        spread = f"{smax/smin:.2f}x" if smin else "-"
        mb = f"{r['mb'][2]:.3f}" if r["mb"] else "-"
        print(f"| {r['config']} | {r['ranks']} | {peers} | {smin} | {smax} | {spread} "
              f"| {stot} | {loc} | {mb} |")

    failed = sorted(res.glob("*/run.FAILED.log"))
    if failed:
        print("\n## FAILED runs\n")
        for f in failed:
            print(f"- `{f.parent.name}` -- see `{f.name}`")
    return 0


if __name__ == "__main__":
    sys.exit(main())
