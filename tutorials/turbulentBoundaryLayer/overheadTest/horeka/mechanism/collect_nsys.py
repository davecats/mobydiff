#!/usr/bin/env python3
"""Reduce the Nsight Systems GPU traces to what the exchange timers cannot say.

    collect_nsys.py <results_nsys_dir> [> nsys.md]

NO MPI EVENTS. nsys `--trace=mpi` captures nothing from this solver: OpenMPI's
Fortran mpi_f08 bindings reach the C layer as PMPI_*, so the interception of the
MPI_* symbols never fires ("does not contain MPI event data"). An LD_PRELOAD
shim misses it for the same reason, and solver-side NVTX is a gated change. The
GPU timeline is enough for the two questions that matter, because the exchange
has a kernel signature:

    pack -> copy_local -> [ gap = MPI post + Waitall ] -> unpack

  1. IS THE WAIT DEVICE WORK?  GPU busy/idle over the steady window. A rank
     blocked on an unfinished CUDA stream keeps the GPU BUSY; a rank blocked on
     a peer leaves it IDLE. This is the one hypothesis no timer in the project
     can test, since every such stall is booked as `mpi_wait` either way.

  2. ALL ROUNDS, OR A FEW?  Per-round wait is the pack/copy -> unpack gap. The
     share held by the slowest 1 % says whether a mean over 7800 calls is a fair
     description or is hiding a handful of catastrophic rounds.

READ SHAPE, NOT MAGNITUDE. Tracing inflates the wait unevenly (measured: base
120 -> ~690 us, rect 738 -> ~1255 us), so absolute values here must never be
quoted against an untraced run. Ratios WITHIN one traced run, and the
distribution shape, are what this file is for. Stdlib only.
"""
import csv, re, sys, statistics
from pathlib import Path


def load(path):
    rows = list(csv.reader(open(path, newline="", errors="replace")))
    if not rows:
        return []
    idx = {c.strip(): i for i, c in enumerate(rows[0])}
    try:
        ni, di, si = idx["Name"], idx["Duration (ns)"], idx["Start (ns)"]
    except KeyError:
        return []
    ev = []
    for row in rows[1:]:
        if len(row) <= max(ni, di, si):
            continue
        try:
            ev.append((float(row[si]), float(row[di]), row[ni]))
        except ValueError:
            pass
    ev.sort()
    return ev


def steady(ev, frac=0.4):
    """Drop the first `frac` of the trace: init, and the one ~600 MB map that
    would otherwise dominate every memcpy statistic."""
    t0 = ev[0][0]
    t1 = max(s + d for s, d, _ in ev)
    cut = t0 + frac * (t1 - t0)
    return [e for e in ev if e[0] >= cut], (t1 - cut)


def merge(ev):
    out = []
    for s, d, _ in ev:
        e = s + d
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def analyse(path):
    ev = load(path)
    if not ev:
        return None
    ss, span = steady(ev)
    if not ss:
        return None
    busy = merge(ss)
    bt = sum(e - s for s, e in busy)
    gaps = sorted((busy[i+1][0] - busy[i][1]) for i in range(len(busy)-1))
    big = [g for g in gaps if g > 50e3]
    r = {"span_ms": span/1e6, "busy_pct": 100*bt/span,
         "idle_ms": sum(gaps)/1e6, "nbig": len(big),
         "big_ms": sum(big)/1e6,
         "big_med_us": statistics.median(big)/1e3 if big else 0.0,
         "big_p90_us": big[int(0.9*len(big))]/1e3 if big else 0.0}
    r["big_share"] = 100*sum(big)/sum(gaps) if gaps else 0.0
    # per-round wait: last kernel before the Waitall -> the unpack that closes it
    W, pre = [], None
    for s, d, n in ss:
        if "pack_entries" in n and "unpack" not in n:
            pre = s + d
        elif "copy_local" in n:
            pre = s + d
        elif "unpack" in n and pre is not None:
            if s >= pre:
                W.append(s - pre)
            pre = None
    if W:
        W.sort()
        tot = sum(W)
        r.update(n=len(W), mean_us=tot/len(W)/1e3, p50_us=W[len(W)//2]/1e3,
                 p90_us=W[int(0.9*len(W))]/1e3, max_us=W[-1]/1e3,
                 top1=100*sum(W[-max(1, len(W)//100):])/tot)
    for key, tag in (("h2d", "Host-to-Device"), ("p2p", "Peer-to-Peer"),
                     ("copy", "copy_local")):
        v = [d for _, d, n in ss if tag in n]
        r[key + "_n"], r[key + "_ms"] = len(v), sum(v)/1e6
    return r


def main():
    res = Path(sys.argv[1] if len(sys.argv) > 1 else "results_nsys")
    print("# Timeline probe — is the wait device work, and is it every round?\n")
    print("No MPI events: OpenMPI's Fortran mpi_f08 bindings reach the C layer as")
    print("PMPI_*, so nsys's MPI_* interception never fires. The exchange's kernel")
    print("signature carries the same information: `pack -> copy_local -> [gap =")
    print("post + Waitall] -> unpack`.\n")
    print("**READ SHAPE, NOT MAGNITUDE.** Tracing inflates the wait unevenly")
    print("(measured: base 120 -> ~690 us, rect 738 -> ~1255 us), so no absolute")
    print("value here may be quoted against an untraced run. Ratios WITHIN one")
    print("traced run, and the distribution shape, are what this file is for.\n")
    for run in sorted(d for d in res.iterdir() if d.is_dir()):
        csvs = sorted(run.glob("rep_*_cuda_gpu_trace.csv"))
        if not csvs:
            continue
        hosts = (run/"hosts.txt").read_text().split() if (run/"hosts.txt").is_file() else []
        print(f"\n## `{run.name}` — {len(hosts)} node(s): {' '.join(hosts)}\n")
        print("| rank | span ms | **GPU busy** | idle ms | gaps>50us | median | p90 | "
              "share of idle | rounds | mean wait | p50 | top 1 % | H2D n | P2P ms | copy ms |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        rows = []
        for c in csvs:
            m = re.search(r"rep_(\d+)", c.name)
            a = analyse(c)
            if not a:
                continue
            a["rank"] = int(m.group(1)) if m else -1
            rows.append(a)
        for a in sorted(rows, key=lambda a: a["rank"]):
            g = lambda k, f="{:.1f}": f.format(a[k]) if k in a else "-"
            print(f"| {a['rank']} | {a['span_ms']:.0f} | **{a['busy_pct']:.1f} %** | "
                  f"{a['idle_ms']:.0f} | {a['nbig']} | {a['big_med_us']:.0f} us | "
                  f"{a['big_p90_us']:.0f} | {a['big_share']:.0f} % | {g('n','{:.0f}')} | "
                  f"{g('mean_us')} us | {g('p50_us')} | {g('top1','{:.0f}')} % | "
                  f"{a['h2d_n']} | {a['p2p_ms']:.1f} | {a['copy_ms']:.1f} |")
        if rows:
            bp = [a["busy_pct"] for a in rows]
            print(f"\n- GPU busy {min(bp):.1f}–{max(bp):.1f} % across ranks.")
            if max(bp) < 85:
                print("  **The GPU is idle through the wait — the stall is NOT unfinished")
                print("  device work.** A rank blocked on a CUDA stream would show the GPU busy.")
            t1 = [a["top1"] for a in rows if "top1" in a]
            if t1 and max(t1) < 12:
                print(f"- The slowest 1 % of rounds hold only {min(t1):.0f}–{max(t1):.0f} % of the")
                print("  wait: it is EVERY round, so a per-round mean is a fair description.")
            elif t1:
                print(f"- The slowest 1 % of rounds hold {min(t1):.0f}–{max(t1):.0f} % of the wait.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
