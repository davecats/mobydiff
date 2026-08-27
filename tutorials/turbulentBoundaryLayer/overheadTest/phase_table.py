#!/usr/bin/env python3
"""Per-phase table from the profiled runs (runs/*_prof/run.log).

The Phase-0 deliverable: for each timing category, seconds per step in each run
and the nb16/base ratio. The point of the ratio column is the decision rule in
docs/next_session_block_overhead.md --

  the block tax lands on exchange pack/unpack/local_copy  -> traffic-bound
  it is spread over the volume kernels in proportion to
      allocated volume, (nb+2)^3/nb^3 = 1.4238 at nb=16   -> footprint-bound
  it lands on mpi_wait                                    -> re-run at 4 ranks

`ibm_mu` is the control: update_ibm_mu is a pure pointwise pass over
halo-carrying arrays with no stencil and no exchange, so its ratio is what a
purely footprint-bound phase looks like on this machine.

Usage: ./phase_table.py [--markdown]
"""

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
PRED = {"nb16": (18 / 16) ** 3, "nb8": (10 / 8) ** 3}   # allocated-volume predictions

LINE = re.compile(
    r"^(step_timing|proj_timing|exch_timing):\s+(\S+)\s+calls\s+(\d+)\s+nsteps\s+(\d+)"
    r"\s+seconds\s+(\S+)\s+seconds_per_step\s+(\S+)")
CHRON = re.compile(r"timing:\s+nsteps\s+\d+\s+loop_seconds\s+(\S+)\s+seconds_per_step\s+(\S+)")
COVER = re.compile(r"step_timing: coverage measured\s+\S+\s+loop_seconds\s+\S+\s+fraction\s+(\S+)")

# Which base each run is read against, keyed by CONFIG name. The refined case
# pairs with nb16_jacobi (same block tax on both sides) so its column isolates
# the refinement. Run directories carry a tag on top of the config name
# (_prof from PROFILE=1, plus anything from SUFFIX=<tag>); pairing is done on
# the config with the tag carried across, so an A/B run still gets its ratios.
PAIRS = {
    "nb16_redblack": "base_redblack",
    "nb16_jacobi": "base_jacobi",
    "nb8_jacobi": "base_jacobi",
    "rect_jacobi": "base_jacobi",
    "refined_yp100_jacobi": "nb16_jacobi",
}
KNOWN = set(PAIRS) | set(PAIRS.values())
ORDER = ["base_redblack", "nb16_redblack", "base_jacobi",
         "nb8_jacobi", "nb16_jacobi", "rect_jacobi", "refined_yp100_jacobi"]


def stem(name):
    """Config name inside a run-directory name (longest known prefix)."""
    best = None
    for key in KNOWN:
        if name == key or name.startswith(key + "_"):
            if best is None or len(key) > len(best):
                best = key
    return best


def base_of(name):
    s = stem(name)
    if s is None or s not in PAIRS:
        return None
    return PAIRS[s] + name[len(s):]


def read(d):
    log = (d / "run.log")
    if not log.exists():
        return None
    text = log.read_text(errors="replace")
    phases, chron, cover = {}, None, None
    for line in text.splitlines():
        m = LINE.match(line.strip())
        if m:
            group, label, _calls, _n, _sec, per_step = m.groups()
            if label != "total_measured":
                phases[f"{group.split('_')[0]}/{label}"] = float(per_step)
            continue
        m = CHRON.search(line)
        if m:
            chron = float(m.group(2))
        m = COVER.search(line)
        if m:
            cover = float(m.group(1))
    return {"phases": phases, "chron": chron, "coverage": cover} if phases else None


def main():
    md = "--markdown" in sys.argv
    # profiled run dirs, in ORDER first then whatever else is present
    found = [d.name for d in sorted(RUNS.iterdir())
             if d.is_dir() and "_prof" in d.name] if RUNS.is_dir() else []
    found.sort(key=lambda n: (ORDER.index(stem(n)) if stem(n) in ORDER else 99, n))
    runs = {}
    for name in found:
        r = read(RUNS / name)
        if r:
            runs[name] = r
    if not runs:
        sys.exit("no profiled runs -- PROFILE=1 ./run_overhead.sh")

    names = list(runs)
    labels = []
    for r in runs.values():
        for k in r["phases"]:
            if k not in labels:
                labels.append(k)

    short = [n.replace("_prof", "").replace("_yp100", "") for n in names]
    short = [s[-14:] for s in short]
    sep = " | " if md else "  "
    head = f"{'phase':<22}" + sep + sep.join(f"{s:>14}" for s in short)
    ratios = [n for n in names if base_of(n) in runs]
    head += sep + sep.join(
        f"{(n.replace('_prof','').replace('_yp100',''))[-16:]+'/base':>22}"
        for n in ratios)
    print(("| " + head + " |") if md else head)
    if md:
        print("|" + "|".join(["---"] * (1 + len(names) + len(ratios))) + "|")
    else:
        print("-" * len(head))

    for lab in labels + ["TOTAL/chron"]:
        cells = []
        for n in names:
            v = runs[n]["chron"] if lab == "TOTAL/chron" else runs[n]["phases"].get(lab)
            cells.append(f"{v:14.5f}" if v is not None else " " * 14)
        rcells = []
        for n in ratios:
            b = runs[base_of(n)]
            a, bb = (runs[n]["chron"], b["chron"]) if lab == "TOTAL/chron" else \
                    (runs[n]["phases"].get(lab), b["phases"].get(lab))
            # Below ~1 ms/step a phase is noise; a ratio there means nothing.
            rcells.append(f"{a/bb:22.3f}" if a and bb and bb > 1e-3 else " " * 22)
        row = f"{lab:<22}" + sep + sep.join(cells) + (sep + sep.join(rcells) if rcells else "")
        print(("| " + row + " |") if md else row)

    print()
    for k, v in PRED.items():
        print(f"allocated-volume prediction, {k}: {v:.4f}")
    for n in names:
        c = runs[n]["coverage"]
        if c is None:
            print(f"WARNING {n}: no coverage line")
        elif abs(c - 1.0) > 0.03:
            print(f"WARNING {n}: step_timing covers only {c:.3f} of the loop -- "
                  f"the instrumentation is missing time and nothing above is trustworthy")
        else:
            print(f"{n}: coverage {c:.4f}")


if __name__ == "__main__":
    main()
