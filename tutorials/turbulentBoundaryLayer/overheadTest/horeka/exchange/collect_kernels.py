#!/usr/bin/env python3
"""Per-kernel device time of a step, from the nsys traces of run_timeline.sh.

    collect_kernels.py <results_dir>

For every `nsys_*` run: one table, every kernel, launches per step, device
microseconds per step, share of the device-busy time, mean microseconds per
launch -- measured in the steady window (40-95 % of the trace) and normalised
by the step count read from `jacobi_compute_phi`'s launches (18 per step for
Chebyshev-Jacobi) or, on a red-black run, from `redblack_sweep`'s (36 per
step: 2 colours x 6 iterations x 3 substages). The run.log's `timing:` line
gives the wall step for the share of the step that is device time.
"""
import csv
import re
import sys
from pathlib import Path

KERNEL_RE = re.compile(r"nvkernel_(.+?)__F1L(\d+)_")
PER_STEP = (("jacobi_compute_phi", 18.0), ("redblack_sweep", 36.0))


def short(name):
    m = KERNEL_RE.match(name)
    return (m.group(1) + "@" + m.group(2)) if m else name


def col(row, *names):
    for n in names:
        for k in row:
            if n.lower() in k.lower():
                return k
    return None


def main():
    res = Path(sys.argv[1])
    print("# Device time per kernel, per step\n")
    prov = res / "provenance.txt"
    if prov.is_file():
        print("```\n" + prov.read_text().strip() + "\n```\n")
    for run in sorted(res.glob("nsys_*")):
        f = run / "rep_0_cuda_gpu_trace.csv"
        if not f.is_file():
            print(f"_{run.name}: no trace_\n")
            continue
        with f.open(newline="") as fh:
            rows = list(csv.DictReader(fh))
        if not rows:
            continue
        r0 = rows[0]
        cs, cd, cn = col(r0, "Start"), col(r0, "Duration"), col(r0, "Name")
        ev = []
        for r in rows:
            try:
                ev.append((int(r[cs]), int(r[cd]), r[cn]))
            except (ValueError, KeyError):
                pass
        ev.sort()
        t0, t1 = ev[0][0], ev[-1][0]
        lo, hi = t0 + 0.40*(t1 - t0), t0 + 0.95*(t1 - t0)
        win = [e for e in ev if lo <= e[0] <= hi and not e[2].startswith("[CUDA")]
        steps = 0.0
        for key, per in PER_STEP:
            n = sum(1 for e in win if key in e[2])
            if n:
                steps = n/per
                break
        if steps <= 0:
            print(f"_{run.name}: no step marker kernel in the window_\n")
            continue
        kern = {}
        for _, dur, name in win:
            d = kern.setdefault(short(name), [0, 0])
            d[0] += 1; d[1] += dur
        busy = sum(v[1] for v in kern.values())
        wall = None
        log = run / "run.log"
        if log.is_file():
            m = re.findall(r"^timing:.*seconds_per_step\s+(\S+)", log.read_text(errors="replace"), re.M)
            if m:
                wall = float(m[-1])*1e6
        print(f"## {run.name[5:]}  ({steps:.1f} steps in the window; device busy "
              f"{busy/steps/1e3:.2f} us/step" + (f" = {100*busy/steps/1e3/wall:.1f} % of the {wall/1e3:.2f} ms step (nsys run)" if wall else "") + ")\n")
        print("| kernel@line | launches/step | us/step | share of device | us/launch |")
        print("|---|---|---|---|---|")
        for k, (n, d) in sorted(kern.items(), key=lambda x: -x[1][1]):
            print(f"| `{k}` | {n/steps:.1f} | {d/steps/1e3:,.1f} | {100*d/busy:.1f} % | {d/n/1e3:.1f} |")
        print()


if __name__ == "__main__":
    main()
