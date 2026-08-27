#!/usr/bin/env python3
"""Summarise the block-overhead timing runs under runs/.

Reports TWO independent rates per run, as the handout requires:

  chron       the solver's own loop timer ("timing: ... seconds_per_step"),
              which covers the whole main loop including start-up effects of
              the first steps.
  marginal    differenced from runtime.txt. That file's last column is the
              CUMULATIVE average since the run started, so the marginal rate
              between two lines i < j is

                  (sps_j*step_j - sps_i*step_i) / (step_j - step_i)

              which assumes the run started at step 0. That holds for the
              cold-start configs here and is CHECKED (a restart writes absolute
              step numbers while averaging only over steps since the restart, so
              sps_i*step_i would not be elapsed time); read_run warns and drops
              the marginal rate if the first line does not sit at one output
              interval. Anchoring at the first line past WARMUP drops
              allocation, the first device maps and the JIT.

The two should agree to a few percent; a gap means the run was disturbed
(shared GPU, another job) and the numbers should not be used.

Ratios are taken against the matching base_* run, which is the only quantity
that transfers between machines.
"""

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
WARMUP = 100          # steps discarded before the marginal-rate anchor

CHRON_RE = re.compile(
    r"timing:\s+nsteps\s+(\d+)\s+loop_seconds\s+(\S+)\s+seconds_per_step\s+(\S+)")
LEAF_RE = re.compile(r"block refinement:\s+(\d+)\s+leaves,\s+(\d+)\s+refined")

# Which base each run is compared against, and the expected cell count.
# The refined case is read against nb16_jacobi, NOT base_jacobi: both carry the
# same block tax, so the ratio isolates what the refinement itself buys.
PAIRS = {
    "nb16_redblack": "base_redblack",
    "nb16_jacobi": "base_jacobi",
    "nb8_jacobi": "base_jacobi",
    "refined_yp100_jacobi": "nb16_jacobi",
}
CELLS = {
    "base_redblack": 4096 * 176 * 192,
    "nb16_redblack": 4096 * 176 * 192,
    "base_jacobi": 4096 * 176 * 192,
    "nb16_jacobi": 4096 * 176 * 192,
    "nb8_jacobi": 4096 * 176 * 192,
    # 2048*176*96 level-0 cells with 3 of 11 y-tiles refined 2x in x and z.
    "refined_yp100_jacobi": 2048 * 176 * 96 * 8 // 11 + 2048 * 176 * 96 * 3 // 11 * 4,
}


KNOWN = set(PAIRS) | set(PAIRS.values()) | set(CELLS)


def read_run(d):
    out = {"name": d.name, "chron": None, "marginal": None, "leaves": None,
           "l2_div": None, "linf": None, "restarted": False}

    log = d / "run.log"
    if log.exists():
        text = log.read_text(errors="replace")
        m = CHRON_RE.search(text)
        if m:
            out["chron"] = float(m.group(3))
            out["nsteps"] = int(m.group(1))
        m = LEAF_RE.search(text)
        if m:
            out["leaves"] = (int(m.group(1)), int(m.group(2)))

    rt = d / "runtime.txt"
    if rt.exists():
        rows = []
        for line in rt.read_text().splitlines()[1:]:
            f = line.split()
            if len(f) == 8:
                rows.append((int(f[0]), float(f[2]), float(f[4]), float(f[7])))
        if len(rows) >= 2:
            # cold-start check: line k sits at (k+1)*interval, so the first line
            # must equal the spacing. A restart breaks the sps*step identity.
            spacing = rows[1][0] - rows[0][0]
            if spacing <= 0 or rows[0][0] != spacing:
                out["restarted"] = True
                return out
            anchor = next((r for r in rows if r[0] > WARMUP), rows[0])
            last = rows[-1]
            if last[0] > anchor[0]:
                out["marginal"] = ((last[3] * last[0] - anchor[3] * anchor[0])
                                   / (last[0] - anchor[0]))
            # The runtime line is the nb-independence gate: identical physics
            # must come out of every layout of the same grid.
            out["l2_div"], out["linf"] = last[1], last[2]
    return out


def main():
    if not RUNS.is_dir():
        sys.exit(f"no runs/ directory -- run ./run_overhead.sh first ({RUNS})")

    runs = {d.name: read_run(d) for d in sorted(RUNS.iterdir()) if d.is_dir()}
    if not runs:
        sys.exit("runs/ is empty")

    # A profiled matrix lives in <name>_prof/ and is summarised exactly like a
    # plain one: strip the suffix to look a run up, keep it to find its base.
    # A run directory is a known config name plus an arbitrary tag from
    # PROFILE=1 / SUFFIX=<tag> (e.g. nb16_jacobi_prof_opt). Recover the config
    # by longest-prefix match rather than stripping a hardcoded list of tags --
    # an unrecognised SUFFIX must not silently blank the ratio column AND the
    # nb-independence gate, which is exactly what an A/B run is there for.
    def stem(name):
        best = None
        for key in KNOWN:
            if name == key or name.startswith(key + "_"):
                if best is None or len(key) > len(best):
                    best = key
        return best if best is not None else name

    def base_of(name):
        b = PAIRS.get(stem(name))
        if b is None:
            return None
        # pair a run with a base measured the same way (same binary, same
        # profiling setting) -- never across an A/B boundary
        return b + name[len(stem(name)):]

    print(f"{'run':<24} {'chron s/step':>13} {'marginal':>10} {'ratio':>7} "
          f"{'cells':>8} {'leaves':>8}  {'L2_div':>12} {'Linf_vel':>12}")
    print("-" * 104)
    for name, r in runs.items():
        base = runs.get(base_of(name) or "")
        ratio = ""
        if base and base.get("marginal") and r.get("marginal"):
            ratio = f"{r['marginal'] / base['marginal']:.3f}"
        cellrat = ""
        if stem(name) in CELLS and PAIRS.get(stem(name)) in CELLS:
            cellrat = f"{CELLS[stem(name)] / CELLS[PAIRS[stem(name)]]:.3f}"
        leaves = f"{r['leaves'][0]}" if r["leaves"] else "-"
        chron = f"{r['chron']:.4f}" if r["chron"] else "-"
        marg = f"{r['marginal']:.4f}" if r["marginal"] else "-"
        l2 = f"{r['l2_div']:.6e}" if r["l2_div"] is not None else "-"
        li = f"{r['linf']:.6e}" if r["linf"] is not None else "-"
        print(f"{name:<24} {chron:>13} {marg:>10} {ratio:>7} {cellrat:>8} "
              f"{leaves:>8}  {l2:>12} {li:>12}")

    print()
    for name, r in runs.items():
        if r.get("restarted"):
            print(f"WARNING {name}: runtime.txt does not start at one output "
                  f"interval -- not a cold start, marginal rate not computed")
    for name, r in runs.items():
        if r["chron"] and r["marginal"]:
            gap = abs(r["chron"] - r["marginal"]) / r["marginal"]
            if gap > 0.05:
                print(f"WARNING {name}: chron and marginal rates differ by "
                      f"{100*gap:.1f}% -- the run was probably disturbed")

    # nb-independence: the block layout must not change the answer.
    for nbname in [n for n in runs if re.match(r"nb\d+_", stem(n))]:
        basename = base_of(nbname) or ""
        a, b = runs.get(nbname), runs.get(basename)
        if a and b and a["l2_div"] is not None and b["l2_div"] is not None:
            same = (a["l2_div"] == b["l2_div"] and a["linf"] == b["linf"])
            print(f"nb-independence {nbname} vs {basename}: "
                  f"{'OK (runtime lines identical)' if same else 'FAILED -- runtime lines differ'}")


if __name__ == "__main__":
    main()
