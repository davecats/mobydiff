#!/usr/bin/env python3
"""Read the placement matrix and decide between the mechanism hypotheses.

    collect_mechanism.py <results_dir> [> mechanism.md]

The predictions are PRE-REGISTERED (see README.md) so the verdict is read off
the data rather than argued after it. Standard library only.
"""
import re, sys
from pathlib import Path

T = re.compile(r"^(exch_timing):\s+(\S+)\s+calls\s+(\d+)\s+nsteps\s+(\d+)"
               r"\s+seconds\s+\S+\s+seconds_per_step\s+(\S+)")
C = re.compile(r"^timing:\s+nsteps\s+\d+\s+loop_seconds\s+\S+\s+seconds_per_step\s+(\S+)")
S = re.compile(r"exchange sizes: peers/rank\(max\)\s+(\d+)\s+send pts/rank min\s+(\d+)"
               r"\s+max\s+(\d+)\s+total send pts\s+(\d+)\s+local copy pts\s+(\d+)")

def parse(d):
    log = d/"run.log"
    if not log.is_file(): return None
    t = log.read_text(errors="replace")
    r = {"dir": d.name}
    m = re.match(r"(.+)_r(\d+)_N(\d+)$", d.name)
    if not m: return None
    r["cfg"], r["ranks"], r["nodes"] = m.group(1), int(m.group(2)), int(m.group(3))
    r["per_node"] = r["ranks"]//r["nodes"]
    # cross-node links in a 1-D chain of `ranks` ranks laid out `per_node` per node
    r["xlinks"] = r["nodes"]-1
    for l in t.splitlines():
        mm = C.match(l.strip())
        if mm: r["sps"] = float(mm.group(1))
        mm = T.match(l.strip())
        if mm and mm.group(2) == "mpi_wait":
            r["wait"] = float(mm.group(5)); r["rounds"] = int(mm.group(3))/int(mm.group(4))
        if mm and mm.group(2) == "local_copy":
            r["copy_t"] = float(mm.group(5))
    mm = S.search(t)
    if mm:
        r["peers"], r["smin"], r["smax"], r["stot"], r["loc"] = [int(x) for x in mm.groups()]
    hosts = (d/"hosts.txt")
    r["hosts"] = hosts.read_text().split() if hosts.is_file() else \
                 sorted(set(re.findall(r"Data for node:\s+(\S+)", t)))
    # a run OpenMPI did not lay out as requested must not be read as if it had
    r["void"] = (d/"VOID").is_file() or len(r["hosts"]) != r["nodes"]
    r["us"] = r.get("wait",0)/r["rounds"]*1e6 if r.get("rounds") else None
    return r if "sps" in r else None

def main():
    res = Path(sys.argv[1] if len(sys.argv)>1 else "results_mechanism")
    runs = [x for x in (parse(d) for d in sorted(res.iterdir()) if d.is_dir()) if x]
    if not runs:
        print("No completed runs in", res)
        for f in sorted(res.glob("*/run.FAILED.log")): print("  FAILED:", f.parent.name)
        return 1
    by = {(r["cfg"], r["ranks"], r["nodes"]): r for r in runs}

    print("# Mechanism probe: is the blocked exchange's stall caused by node crossing?\n")
    print("## All runs\n")
    print("| run | ranks | nodes | ranks/node | hosts | s/step | wait/round | peers | send max | local copy |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(runs, key=lambda r:(r["cfg"], r["ranks"], r["nodes"])):
        flag = "  **VOID**" if r["void"] else ""
        print(f"| {r['cfg']}{flag} | {r['ranks']} | {r['nodes']} | {r['per_node']} | {len(r['hosts'])} "
              f"| {r['sps']:.5f} | {r['us']:.1f} us | {r.get('peers','-')} "
              f"| {r.get('smax','-')} | {r.get('loc','-')} |")

    def line(cfg, ranks, nodes):
        r = by.get((cfg, ranks, nodes))
        return None if not r or r["void"] else r["us"]   # never quote a void run

    nvoid = sum(1 for r in runs if r["void"])
    if nvoid:
        print(f"\n> **WARNING: {nvoid} run(s) were not laid out as requested and are")
        print("> excluded from every comparison below.** Fix the launcher before")
        print("> drawing any conclusion -- see README 'Reading it', item 1.\n")

    # ---- TIER 1: the decisive pair -----------------------------------------
    print("\n## TIER 1 — one link, no chain: 2 ranks on 1 node vs 2 nodes\n")
    print("At 2 ranks each rank has exactly ONE peer and the pair is symmetric,")
    print("so no stall can propagate. Only the link's locality changes.\n")
    print("| config | 2x1 (intra) | 2x2 (cross) | ratio |")
    print("|---|---|---|---|")
    verdicts = {}
    for cfg in ("rect_jacobi", "base_jacobi", "refined_yp82_rect_jacobi"):
        a, b = line(cfg,2,1), line(cfg,2,2)
        if a and b:
            verdicts[cfg] = b/a
            print(f"| {cfg} | {a:.1f} us | {b:.1f} us | **{b/a:.2f}x** |")
    print()
    # The discriminator is the BLOCKED-vs-UNBLOCKED gap on the SAME single
    # cross-node link, not each config's own jump: the unblocked control is
    # itself expected to rise when it starts crossing (38.9 -> 121 us in the
    # main matrix), so its jump carries no information on its own.
    rc, bc = line("rect_jacobi",2,2), line("base_jacobi",2,2)
    ri, bi = line("rect_jacobi",2,1), line("base_jacobi",2,1)
    if rc and bc:
        gap = rc/bc
        print(f"\nBlocked / unblocked on ONE cross-node link: **{gap:.2f}x** "
              f"(rect {rc:.1f} us vs base {bc:.1f} us)")
        if ri and bi:
            print(f"Same ratio intra-node: {ri/bi:.2f}x (rect {ri:.1f} vs base {bi:.1f})")
        print()
        if gap > 3:
            print("**VERDICT: H-transport.** A SINGLE cross-node link is intrinsically")
            print("slow in the blocked path -- two symmetric ranks, one peer each, no")
            print("chain to propagate a stall and no third rank to be imbalanced")
            print("against, while the unblocked control over the same link is fine.")
            print("The cause is in the blocked path's buffers or message structure.")
            print("Chain propagation is NOT required, so re-partitioning would NOT")
            print("fix it. Next: the blk8 probe below says whether the device-local")
            print("copy is the co-factor.")
        elif gap < 1.6:
            print("**VERDICT: H-chain.** One cross-node link is as cheap for the blocked")
            print("path as for the unblocked one, so the 10x at 8+ ranks is COLLECTIVE:")
            print("it needs several ranks and links. That is chain propagation or load")
            print("skew, and the target is the PARTITIONING (the 2-peer Morton chain).")
            print("Confirm with the link-count series: saturating after the first")
            print("crossing is the chain signature.")
        else:
            print(f"**VERDICT: intermediate ({gap:.2f}x).** Both effects plausibly present.")
            print("Weight the link-count series (proportional = per-link transport,")
            print("saturating = chain) and the blk8 probe before choosing a fix.")
    else:
        print("\n**Tier 1 incomplete** -- the decisive pair did not both produce")
        print("valid runs. Nothing below can settle the mechanism without it.")

    # ---- link-count series -------------------------------------------------
    print("\n## Link count at FIXED rank count\n")
    print("Rank count, decomposition, peer count and volume held fixed; only the")
    print("number of node boundaries the chain crosses changes.\n")
    print("| config | ranks | 1 node | 2 nodes | 4 nodes |")
    print("|---|---|---|---|---|")
    for cfg in ("rect_jacobi","base_jacobi","refined_yp82_rect_jacobi","nb16_jacobi","blk8_jacobi"):
        for ranks in (2,4,8):
            cells = [line(cfg,ranks,n) for n in (1,2,4)]
            if not any(cells): continue
            f = lambda v: "-" if v is None else f"{v:.1f} us"
            print(f"| {cfg} | {ranks} | {f(cells[0])} | {f(cells[1])} | {f(cells[2])} |")
    print("\nProportional in the number of crossings => per-link transport cost.")
    print("Saturating after the FIRST crossing => a single slow link stalls the chain.")

    # ---- blk8: node crossing without the local-copy load --------------------
    print("\n## Discriminator — node crossing WITHOUT device-local copy (`blk8_jacobi`)\n")
    print("One block per rank at 8 ranks: same 2-peer chain, same cells, but the")
    print("same-rank copy volume collapses.\n")
    print("| run | wait/round | local copy pts |")
    print("|---|---|---|")
    for key in (("rect_jacobi",8,2),("blk8_jacobi",8,2),("base_jacobi",8,2),
                ("rect_jacobi",4,1),("blk8_jacobi",4,1)):
        r = by.get(key)
        if r: print(f"| {key[0]} {key[1]}x{key[2]} | {r['us']:.1f} us | {r.get('loc','-')} |")
    a, b = by.get(("blk8_jacobi",8,2)), by.get(("rect_jacobi",8,2))
    if a and b:
        if a["us"] < 0.4*b["us"]:
            print("\n**Local-copy / transfer CONTENTION confirmed**: removing the same-rank")
            print("copy load removes most of the stall at the same placement and topology.")
        else:
            print("\n**Local copy is INNOCENT**: the stall survives without it, so the cause")
            print("is the 2-peer chain topology or the blocked path's message structure.")

    # ---- volume ------------------------------------------------------------
    print("\n## Latency or bandwidth — exchange volume at a fixed placement\n")
    print("| run | send pts (max/rank) | wait/round | us per MB |")
    print("|---|---|---|---|")
    for key in (("rect_jacobi",4,2),("nb16_jacobi",4,2),("blk8_jacobi",8,2),("rect_jacobi",8,2)):
        r = by.get(key)
        if r and r.get("smax"):
            mb = r["smax"]*32/1e6
            print(f"| {key[0]} {key[1]}x{key[2]} | {r['smax']} | {r['us']:.1f} us | {r['us']/mb:.1f} |")
    print("\nWait flat against a large change in bytes => latency/serialization.")
    print("Wait tracking the bytes => bandwidth.")

    failed = sorted(res.glob("*/run.FAILED.log"))
    if failed:
        print("\n## FAILED\n")
        for f in failed: print(f"- `{f.parent.name}`")
    return 0

if __name__ == "__main__":
    sys.exit(main())
