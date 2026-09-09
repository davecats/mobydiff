#!/usr/bin/env python3
"""What is a rank waiting for? Read the be48d44 exchange diagnostics.

    collect_balance.py <normal_results_dir> [<barrier_results_dir>]

Two independent questions, one per diagnostic:

  `exchange balance:`   mpi_wait min/mean/max ACROSS RANKS plus the argmax.
                        max/min ~ 1 -> every rank waits equally (fabric).
                        max/min >> 1 with a stable argmax -> everyone waits for
                        one late rank (skew / partitioning).

  exchange_barrier      an MPI_Barrier before every Waitall. The barrier absorbs
                        the ARRIVAL SKEW into bucket `skew_barrier`; what stays
                        in mpi_wait is TRANSFER. Serialising: its step times are
                        meaningless and are never printed as such here.

The decision table is the one pre-registered in HANDOUT_cluster_session.md, so
the verdict is read off the data rather than argued after it. Stdlib only.
"""
import re, sys
from pathlib import Path

T = re.compile(r"^exch_timing:\s+(\S+)\s+calls\s+(\d+)\s+nsteps\s+(\d+)\s+seconds\s+(\S+)")
C = re.compile(r"^timing:\s+nsteps\s+\d+\s+loop_seconds\s+\S+\s+seconds_per_step\s+(\S+)")
B = re.compile(r"exchange balance: mpi_wait seconds min\s+(\S+)\s+mean\s+(\S+)\s+max\s+(\S+)"
               r"\s+max/min\s+(\S+)\s+slowest rank\s+(\d+)")
S = re.compile(r"exchange sizes: peers/rank\(max\)\s+(\d+)\s+send pts/rank min\s+(\d+)"
               r"\s+max\s+(\d+)\s+total send pts\s+(\d+)\s+local copy pts\s+(\d+)")


def parse(d):
    log = d / "run.log"
    if not log.is_file():
        return None
    t = log.read_text(errors="replace")
    m = re.match(r"(.+)_r(\d+)_N(\d+)$", d.name)
    if not m:
        return None
    r = {"dir": d.name, "cfg": m.group(1), "ranks": int(m.group(2)), "nodes": int(m.group(3))}
    r["per_node"] = r["ranks"] // r["nodes"]
    r["buckets"] = {}
    for line in t.splitlines():
        s = line.strip()
        mm = C.match(s)
        if mm:
            r["sps"] = float(mm.group(1))
        mm = T.match(s)
        if mm:
            r["buckets"][mm.group(1)] = (int(mm.group(2)), int(mm.group(3)), float(mm.group(4)))
    mm = B.search(t)
    if mm:
        r["bmin"], r["bmean"], r["bmax"] = (float(mm.group(i)) for i in (1, 2, 3))
        r["bratio"], r["bargmax"] = float(mm.group(4)), int(mm.group(5))
    mm = S.search(t)
    if mm:
        r["peers"], r["smin"], r["smax"], r["stot"], r["loc"] = [int(x) for x in mm.groups()]
    # last printed L2_div: the rank-independence check
    l2 = [LL.split()[2] for LL in t.splitlines()
          if re.match(r"^\s*\d+\s+\S+\s+\S+\s+\S+", LL) and LL.split()[0].isdigit()]
    r["l2"] = l2[-1] if l2 else None
    hosts = d / "hosts.txt"
    r["hosts"] = hosts.read_text().split() if hosts.is_file() else []
    r["void"] = (d / "VOID").is_file() or len(r["hosts"]) != r["nodes"]
    if "mpi_wait" in r["buckets"]:
        calls, nsteps, sec = r["buckets"]["mpi_wait"]
        r["rounds"] = calls / nsteps
        r["calls"] = calls
        r["wait_us"] = sec / calls * 1e6          # rank-0 wait per round
        for k in ("bmin", "bmean", "bmax"):
            if k in r:
                r[k + "_us"] = r[k] / calls * 1e6
    if "skew_barrier" in r["buckets"]:
        calls, _, sec = r["buckets"]["skew_barrier"]
        r["skew_us"] = sec / calls * 1e6 if calls else 0.0
    return r if "sps" in r else None


def load(p):
    if p is None or not Path(p).is_dir():
        return {}
    runs = [x for x in (parse(d) for d in sorted(Path(p).iterdir()) if d.is_dir()) if x]
    return {(r["cfg"], r["ranks"], r["nodes"]): r for r in runs}


def main():
    A = load(sys.argv[1] if len(sys.argv) > 1 else "results_balance")
    Bp = load(sys.argv[2] if len(sys.argv) > 2 else None)
    if not A:
        print("No completed runs in the normal pass.")
        return 1

    print("# What a rank is waiting for — the be48d44 exchange diagnostics\n")

    print("## Pass A — timing and per-rank wait balance\n")
    print("`wait r0` is the rank-0 bucket every earlier report quoted; min/mean/max")
    print("are that same accumulated `mpi_wait` reduced across ALL ranks, per round.\n")
    print("| run | r x N | /node | s/step | wait r0 | min | mean | max | max/min | argmax | peers | copy pts |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for k in sorted(A, key=lambda k: (k[0], k[1], k[2])):
        r = A[k]
        f = lambda x: f"{r[x]:.1f}" if x in r else "-"
        flag = " **VOID**" if r["void"] else ""
        print(f"| {r['cfg']}{flag} | {r['ranks']}x{r['nodes']} | {r['per_node']} | {r['sps']:.5f} "
              f"| {f('wait_us')} | {f('bmin_us')} | {f('bmean_us')} | {f('bmax_us')} "
              f"| {r.get('bratio', float('nan')):.2f} | {r.get('bargmax', '-')} "
              f"| {r.get('peers', '-')} | {r.get('loc', '-')} |")

    print("\n### Rank independence (`L2_div` must match across placements)\n")
    bycfg = {}
    for (cfg, ranks, _), r in A.items():
        bycfg.setdefault((cfg, ranks), set()).add(r["l2"])
    bad = {k: v for k, v in bycfg.items() if len(v) > 1}
    if bad:
        for k, v in bad.items():
            print(f"- **MISMATCH** `{k[0]}` {k[1]} ranks: {sorted(v)}")
    else:
        print("All configs: identical `L2_div` across every placement at each rank count.")

    if Bp:
        print("\n## Pass B — barrier: arrival SKEW vs TRANSFER\n")
        print("`exchange_barrier = true` puts an MPI_Barrier before every Waitall, so")
        print("the barrier absorbs the arrival skew (`skew_barrier`) and what is left")
        print("in `mpi_wait` is transfer. **Serialising — its step times mean nothing.**\n")
        print("| run | r x N | skew/round | transfer/round | skew share | pass-A wait r0 |")
        print("|---|---|---|---|---|---|")
        for k in sorted(Bp, key=lambda k: (k[0], k[1], k[2])):
            r = Bp[k]
            sk, tr = r.get("skew_us"), r.get("wait_us")
            share = f"{100*sk/(sk+tr):.0f} %" if sk is not None and tr else "-"
            a = A.get(k)
            ref = f"{a['wait_us']:.1f} us" if a and "wait_us" in a else "-"
            print(f"| {r['cfg']} | {r['ranks']}x{r['nodes']} | {sk:.1f} us "
                  f"| {tr:.1f} us | {share} | {ref} |")

    print("\n## Verdict, against the pre-registered table\n")
    print("| run | max/min | skew vs transfer | reading |")
    print("|---|---|---|---|")
    for k in sorted(A, key=lambda k: (k[0], k[1], k[2])):
        r = A[k]
        if "bratio" not in r:
            continue
        b = Bp.get(k)
        sk = b.get("skew_us") if b else None
        tr = b.get("wait_us") if b else None
        if sk is None or tr is None:
            sv, verdict = "-", "no barrier pass"
        else:
            sv = f"{sk:.0f} / {tr:.0f} us"
            skewy = sk > 2 * tr
            spread = r["bratio"] > 2
            if not spread and not skewy:
                verdict = "**transport** — every rank waits equally, and the wait is transfer"
            elif not spread and skewy:
                verdict = "**collective skew** — all ranks late together (upstream imbalance in lockstep)"
            elif spread and skewy:
                verdict = "**skew** — everyone waits for one late rank; fix upstream (load balance)"
            else:
                verdict = "**per-rank transfer** — ranks differ in volume, not arrival"
        print(f"| {r['cfg']} {r['ranks']}x{r['nodes']} | {r['bratio']:.2f} (rank {r['bargmax']}) | {sv} | {verdict} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
