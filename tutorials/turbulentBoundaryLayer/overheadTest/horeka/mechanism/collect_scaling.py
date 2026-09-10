#!/usr/bin/env python3
"""Re-emit the 2026-09-07 campaign tables for two binaries side by side.

    collect_scaling.py <ref_results_dir> <new_results_dir> [> scaling.md]

Every table in results_horeka_2026-09-07.md was measured with the old
`device = local_rank mod ndev` mapping, which mismatches GPU affinity classes at
every node boundary. This regenerates them from the same 23-run matrix run twice
in one allocation, so the two columns differ only in the mapping.

Cell counts come from the run itself: `block refinement: N leaves` times the
product of `[blocks] nb`, falling back to nx*ny*nz when nb is unset. That
reproduces the published 60.56 M / 138.41 M / 242.22 M exactly. Stdlib only.
"""
import re, sys
from pathlib import Path

T  = re.compile(r"^timing:\s+nsteps\s+\d+\s+loop_seconds\s+\S+\s+seconds_per_step\s+(\S+)")
W  = re.compile(r"^exch_timing:\s+mpi_wait\s+calls\s+(\d+)\s+nsteps\s+(\d+)\s+seconds\s+(\S+)")
LF = re.compile(r"block refinement:\s+(\d+)\s+leaves")


def cfg_val(text, key):
    m = re.search(r"^\s*%s\s*=\s*(.+?)\s*$" % key, text, re.M)
    return m.group(1) if m else None


def parse(d):
    log = d / "run.log"
    if not log.is_file():
        return None
    t = log.read_text(errors="replace")
    m = re.match(r"(.+)_n(\d+)$", d.name)
    if not m:
        return None
    r = {"cfg": m.group(1), "n": int(m.group(2))}
    for line in t.splitlines():
        s = line.strip()
        mm = T.match(s)
        if mm:
            r["sps"] = float(mm.group(1))
        mm = W.match(s)
        if mm:
            calls, nsteps, sec = int(mm.group(1)), int(mm.group(2)), float(mm.group(3))
            if calls > 0 and nsteps > 0:
                # A 1-rank run has no peers, so mpi_wait is reported with zero
                # calls -- that is a legitimate row, not a broken one.
                r["rounds"] = calls / nsteps
                r["wait"] = sec / calls * 1e6
    if "sps" not in r:
        return None
    ini = (d / "config.ini")
    text = ini.read_text(errors="replace") if ini.is_file() else ""
    nb = cfg_val(text, "nb")
    mm = LF.search(t)
    if nb and mm:
        p = [int(x) for x in nb.split()]
        if len(p) == 1:
            p = p * 3
        r["cells"] = int(mm.group(1)) * p[0] * p[1] * p[2]
    else:
        try:
            r["cells"] = (int(cfg_val(text, "nx")) * int(cfg_val(text, "ny"))
                          * int(cfg_val(text, "nz")))
        except (TypeError, ValueError):
            r["cells"] = None
    return r


def load(p):
    p = Path(p)
    if not p.is_dir():
        return {}
    out = {}
    for d in sorted(p.iterdir()):
        if not d.is_dir():
            continue
        r = parse(d)
        if r:
            out[(r["cfg"], r["n"])] = r
    return out


def main():
    ref = load(sys.argv[1]); new = load(sys.argv[2])
    if not new:
        print("No completed runs in the new-binary matrix."); return 1
    RANKS = [1, 2, 4, 8, 16]
    CFGS = ["base_jacobi", "rect_jacobi", "refined_yp82_rect_jacobi",
            "refined_yp82_rect_redblack", "refined_big_rect_jacobi"]

    def g(src, c, n, k="sps"):
        r = src.get((c, n))
        return r.get(k) if r else None

    print("# The campaign, re-measured at the corrected rank-to-GPU mapping\n")
    print("Same 23-run matrix as `results_horeka_2026-09-07.md` (`run_matrix.sh`,")
    print("`--map-by numa --bind-to core`, 400 steps), run twice in ONE allocation:")
    print("`ref` is the pre-change binary (`device = local_rank mod ndev`), `new`")
    print("the topology-aware mapping. The two columns differ only in the mapping.\n")

    print("## s/step, and what the mapping is worth\n")
    print("| config | ranks | ref | new | gain |")
    print("|---|---|---|---|---|")
    for c in CFGS:
        for n in RANKS:
            a, b = g(ref, c, n), g(new, c, n)
            if b is None:
                continue
            gain = "-" if not a else f"**{100*(a-b)/a:+.1f} %**"
            print(f"| `{c}` | {n} | {a if a is None else f'{a:.5f}'} | {b:.5f} | {gain} |")

    print("\n## Per-round `mpi_wait`\n")
    print("| config | " + " | ".join(f"n={n}" for n in RANKS) + " |")
    print("|---" * (len(RANKS) + 1) + "|")
    for c in CFGS:
        for src, tag in ((ref, "ref"), (new, "new")):
            cells = [g(src, c, n, "wait") for n in RANKS]
            if not any(x is not None for x in cells):
                continue
            row = " | ".join("-" if x is None else f"{x:.1f} us" for x in cells)
            print(f"| `{c}` ({tag}) | {row} |")

    print("\n## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)\n")
    print("| ranks | Mcell/GPU | base ref | base new | rect ref | rect new "
          "| **tax ref** | **tax new** |")
    print("|---|---|---|---|---|---|---|---|")
    for n in RANKS:
        vals = {}
        for src, tag in ((ref, "ref"), (new, "new")):
            for c in ("base_jacobi", "rect_jacobi"):
                r = src.get((c, n))
                vals[(c, tag)] = (r["sps"] * n / r["cells"] * 1e9
                                  if r and r.get("cells") else None)
        rb = new.get(("rect_jacobi", n))
        mc = f"{rb['cells']/n/1e6:.2f}" if rb and rb.get("cells") else "-"
        f = lambda k: "-" if vals.get(k) is None else f"{vals[k]:.3f}"
        tax = lambda tg: ("-" if not (vals.get(("rect_jacobi", tg)) and vals.get(("base_jacobi", tg)))
                          else f"**{vals[('rect_jacobi',tg)]/vals[('base_jacobi',tg)]:.3f}**")
        print(f"| {n} | {mc} | {f(('base_jacobi','ref'))} | {f(('base_jacobi','new'))} "
              f"| {f(('rect_jacobi','ref'))} | {f(('rect_jacobi','new'))} "
              f"| {tax('ref')} | {tax('new')} |")

    print("\n## Strong-scaling efficiency (vs each config's smallest rank count)\n")
    print("| config | " + " | ".join(f"n={n}" for n in RANKS) + " |")
    print("|---" * (len(RANKS) + 1) + "|")
    for c in CFGS:
        for src, tag in ((ref, "ref"), (new, "new")):
            have = [n for n in RANKS if g(src, c, n)]
            if not have:
                continue
            n0 = have[0]; t0 = g(src, c, n0)
            row = []
            for n in RANKS:
                t = g(src, c, n)
                row.append("-" if t is None else f"{100*t0*n0/(t*n):.0f} %")
            print(f"| `{c}` ({tag}) | " + " | ".join(row) + " |")

    print("\n## Headline 1 — the 2:1 machinery, like for like\n")
    print("`refined_big_rect_jacobi` shares grid, block shape and solver with")
    print("`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells")
    print("against the coarse cells beside them:\n")
    print("| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |")
    print("|---|---|---|---|---|---|")
    for n in RANKS:
        for src, tag in ((ref, "ref"), (new, "new")):
            a, b = src.get(("rect_jacobi", n)), src.get(("refined_big_rect_jacobi", n))
            if not (a and b and a.get("cells") and b.get("cells")):
                continue
            extra = b["cells"] - a["cells"]
            ns = (b["sps"] - a["sps"]) / extra * 1e9
            base = a["sps"] / a["cells"] * 1e9
            print(f"| {n} | {tag} | {a['sps']:.5f} | {b['sps']:.5f} | {ns:.4f} | "
                  f"**{ns/base:.3f}x** |")

    print("\n## Headline 3 — red-black against Jacobi\n")
    print("| binary | " + " | ".join(f"n={n}" for n in RANKS) + " |")
    print("|---" * (len(RANKS) + 1) + "|")
    for src, tag in ((ref, "ref"), (new, "new")):
        row = []
        for n in RANKS:
            a = g(src, "refined_yp82_rect_redblack", n)
            b = g(src, "refined_yp82_rect_jacobi", n)
            row.append("-" if not (a and b) else f"{a/b:.3f}")
        print(f"| {tag} | " + " | ".join(row) + " |")

    fails = sorted(Path(sys.argv[2]).glob("*/run.FAILED.log"))
    if fails:
        print("\n## FAILED (new)\n")
        for f in fails:
            print(f"- `{f.parent.name}`")
    return 0


if __name__ == "__main__":
    sys.exit(main())
