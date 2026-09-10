#!/usr/bin/env python3
"""Tables for the op-split and A0 passes of results_exchange/.

    collect_exchange.py <results_dir>

Reads only run.log: the `exchange sizes` / `exchange by op` lines written at
init, the exch_timing buckets (seconds AND calls, so everything is quoted per
CALL = per exchange round) and the chron loop time.

The table that matters is "time per call against points per call": a bucket
whose per-call time barely moves while its points move by an order of magnitude
is paying a fixed per-round cost, and no amount of copying fewer points will
touch it. Nothing here is fitted for the reader -- the two columns are printed
side by side and the fit is stated separately in the report.
"""
import re
import sys
from pathlib import Path

OPS = ("copy", "restrict", "prolong")
BUCKETS = ("pack", "mpi_post", "mpi_wait", "unpack", "local_copy", "copy_cross")


def parse(d):
    """Everything this report needs out of one run directory."""
    r = {"ops": {}, "sec": {}, "calls": {}, "name": d.name}
    m = re.search(r"_[nr](\d+)", d.name)
    r["ranks"] = int(m.group(1)) if m else 1
    for line in (d / "run.log").read_text(errors="replace").splitlines():
        m = re.search(r"exchange sizes: peers/rank\(max\) (\d+).*total send pts (\d+)"
                      r"\s+local copy pts (\d+)", line)
        if m:
            r["peers"], r["send"], r["local"] = (int(m.group(i)) for i in (1, 2, 3))
        m = re.search(r"exchange by op:\s+(\w+)\s+local pts\s+(\d+) entries\s+(\d+)"
                      r"\s+send pts\s+(\d+) entries\s+(\d+)", line)
        if m:
            r["ops"][m.group(1)] = tuple(int(m.group(i)) for i in (2, 3, 4, 5))
        m = re.search(r"exch_timing: (\S+) calls (\d+) nsteps (\d+) seconds\s+(\S+)", line)
        if m:
            r["calls"][m.group(1)] = int(m.group(2))
            r["sec"][m.group(1)] = float(m.group(4))
            r["nsteps"] = int(m.group(3))
        m = re.search(r"timing: nsteps\s+(\d+) loop_seconds\s+\S+ seconds_per_step\s+(\S+)", line)
        if m:
            r["sstep"] = float(m.group(2))
        m = re.search(r"A0 PROBE ACTIVE: (\d+) iterations", line)
        if m:
            r["a0iter"] = int(m.group(1))
        m = re.match(r"\s*\d+\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*$", line)
        if m:
            r["l2div"] = m.group(1)
    return r


def per_call(r, bucket):
    """Seconds per exchange round for one bucket, in microseconds."""
    n = r["calls"].get(bucket, 0)
    return 1e6 * r["sec"].get(bucket, 0.0) / n if n else 0.0


def calls_per_step(r, bucket):
    return r["calls"].get(bucket, 0) / r.get("nsteps", 1)


def main():
    res = Path(sys.argv[1])
    runs = {d.name: parse(d) for d in sorted(res.iterdir())
            if (d / "run.log").is_file()}

    print("# Exchange op split, per-round cost, and the A0 probe at scale\n")
    prov = res / "provenance.txt"
    if prov.is_file():
        print("```\n" + prov.read_text().strip() + "\n```\n")

    print("## Gate -- fields against the pre-change binary\n")
    got_gate = False
    for f in sorted(res.glob("gate_*.txt")):
        got_gate = True
        print(f"**{f.stem}**\n\n```\n{f.read_text().strip()}\n```\n")
    if not got_gate:
        print("_no gate output in this directory_\n")

    ops = [(n, r) for n, r in runs.items() if n.startswith("op_") and r["ops"]]

    print("\n## Pass 1a -- exchange volume by op (points summed over ranks)\n")
    print("| run | peers/rank | local copy | local restrict | local prolong | "
          "send copy | send restrict | send prolong | cross-level % |")
    print("|" + "---|" * 9)
    for name, r in ops:
        lp = [r["ops"].get(o, (0, 0, 0, 0))[0] for o in OPS]
        sp = [r["ops"].get(o, (0, 0, 0, 0))[2] for o in OPS]
        tot = sum(lp) + sum(sp)
        cross = (sum(lp[1:]) + sum(sp[1:])) / tot if tot else 0.0
        print(f"| {name[3:]} | {r.get('peers','?')} | " +
              " | ".join(f"{v:,}" for v in lp + sp) + f" | {100*cross:.2f} % |")

    print("\n## Pass 1b -- entries, and points per entry\n")
    print("| run | copy entries | cross entries | pts/entry copy | pts/entry cross |")
    print("|" + "---|" * 5)
    for name, r in ops:
        ent = [r["ops"].get(o, (0, 0, 0, 0))[1] + r["ops"].get(o, (0, 0, 0, 0))[3]
               for o in OPS]
        pts = [r["ops"].get(o, (0, 0, 0, 0))[0] + r["ops"].get(o, (0, 0, 0, 0))[2]
               for o in OPS]
        ce, cp = sum(ent[1:]), sum(pts[1:])
        print(f"| {name[3:]} | {ent[0]:,} | {ce:,} | "
              f"{pts[0]/ent[0]:.0f} | {cp/ce:.0f} |" if ce and ent[0] else
              f"| {name[3:]} | {ent[0]:,} | {ce:,} | - | - |")

    print("\n## Pass 1c -- microseconds per exchange ROUND, against points per round\n")
    print("Points are PER RANK (the printed totals divided by the rank count);")
    print("times are rank 0's. `local_copy` is same-level only, `copy_cross` is")
    print("the 2:1 restrict/prolong kernel (exactly zero on a single-level grid).\n")
    print("| run | s/step | local pts/rank | send pts/rank | pack us | unpack us | "
          "local_copy us | copy_cross us | mpi_wait us | rounds/step |")
    print("|" + "---|" * 10)
    for name, r in runs.items():
        if not name.startswith("op_") or "sstep" not in r:
            continue
        n = r["ranks"]
        print(f"| {name[3:]} | {r['sstep']:.6f} | {r.get('local',0)//n:,} | "
              f"{r.get('send',0)//n:,} | " +
              " | ".join(f"{per_call(r,b):.1f}" for b in
                         ("pack", "unpack", "local_copy", "copy_cross", "mpi_wait")) +
              f" | {calls_per_step(r,'pack'):.0f} |")

    print("\n### The same, as a share of the step\n")
    print("| run | exchange total s/step | % of step | device-local % | mpi_wait % |")
    print("|" + "---|" * 5)
    for name, r in runs.items():
        if not name.startswith("op_") or "sstep" not in r:
            continue
        tot = sum(r["sec"].get(b, 0.0) for b in BUCKETS) / r["nsteps"]
        dev = sum(r["sec"].get(b, 0.0) for b in
                  ("pack", "unpack", "local_copy", "copy_cross")) / r["nsteps"]
        w = r["sec"].get("mpi_wait", 0.0) / r["nsteps"]
        print(f"| {name[3:]} | {tot:.6f} | {100*tot/r['sstep']:.1f} % | "
              f"{100*dev/r['sstep']:.1f} % | {100*w/r['sstep']:.1f} % |")

    print("\n## Pass 2 -- the A0 overlap probe\n")
    print("If an in-flight transfer progressed while the probe kernel ran, `mpi_wait`")
    print("would collapse toward zero as `a0_probe` grows.\n")
    print("| case | probe us/round | mpi_wait us/round | vs its baseline | s/step |")
    print("|" + "---|" * 5)
    for name, r in sorted(runs.items()):
        if not name.startswith("a0_") or "sstep" not in r:
            continue
        base = runs.get(name.rsplit("_", 1)[0] + "_base", {})
        w = per_call(r, "mpi_wait")
        p = per_call(r, "a0_probe")
        wb = per_call(base, "mpi_wait") if base else 0.0
        rel = f"{w/wb:.2f}x" if wb else "-"
        print(f"| {name[3:]} | {p:.1f} | {w:.1f} | {rel} | {r['sstep']:.6f} |")

    print("\n## L2_div -- the correctness check\n")
    print("| run | last L2_div |")
    print("|---|---|")
    for name, r in sorted(runs.items()):
        if "l2div" in r:
            print(f"| {name} | {r['l2div']} |")


if __name__ == "__main__":
    main()
