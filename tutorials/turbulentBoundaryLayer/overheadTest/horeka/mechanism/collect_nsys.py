#!/usr/bin/env python3
"""Reduce the Nsight Systems traces to the question: what is a rank waiting for?

    collect_nsys.py <results_nsys_dir> [> nsys.md]

Three things the aggregate `exch_timing` buckets cannot answer, one section each:

  1. IS THE WAIT EVEN MPI?  Every MPI_Waitall interval is intersected with the
     merged CUDA kernel intervals from the SAME report (same process, same
     clock). A Waitall that overlaps device activity is waiting for the GPU, not
     the wire, and no partitioning or transport change would touch it.

  2. IS IT ALL ROUNDS OR A FEW?  752 us is a mean over 7800 calls; thirty-eight
     cheap rounds and one catastrophic one give the same mean as uniform
     slowness and need a different fix. Reported as the share of total wait held
     by the slowest 1 % and 10 % of calls.

  3. WHO WAITS, AND AFTER WHAT?  Per-rank wait totals (the balance line,
     decomposed) and the gap between a rank's last posted request and its entry
     into Waitall.

DELIBERATELY NOT DONE: comparing absolute timestamps ACROSS ranks. nsys aligns
clocks within a report; across nodes that alignment is not good enough to call
one rank "late" by microseconds. Everything here is a per-rank duration or a
same-report overlap. Stdlib only.
"""
import csv, re, sys
from pathlib import Path


def _num(s):
    try:
        return float(str(s).replace(",", ""))
    except (TypeError, ValueError):
        return None


def read_csv(path):
    """Return (header, rows) with the header sniffed -- nsys column names move
    between versions, so nothing here hard-codes them."""
    with open(path, newline="", errors="replace") as fh:
        r = list(csv.reader(fh))
    if not r:
        return [], []
    return r[0], r[1:]


def col(header, *wants):
    """First column whose name contains all of `wants` (case-insensitive)."""
    for i, h in enumerate(header):
        hl = h.lower()
        if all(w in hl for w in wants):
            return i
    return None


def intervals(path, start_key=("start",), dur_key=("dur",), name_filter=None, name_key=None):
    header, rows = read_csv(path)
    si, di = col(header, *start_key), col(header, *dur_key)
    ni = col(header, *name_key) if name_key else None
    if si is None or di is None:
        return []
    out = []
    for row in rows:
        if len(row) <= max(si, di, ni if ni is not None else 0):
            continue
        if name_filter is not None and ni is not None:
            if name_filter not in row[ni]:
                continue
        s, d = _num(row[si]), _num(row[di])
        if s is None or d is None:
            continue
        out.append((s, s + d, row[ni] if ni is not None else ""))
    out.sort()
    return out


def merge(iv):
    out = []
    for s, e, _ in iv:
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return out


def overlap_total(waits, merged):
    """Total intersection of `waits` with the merged busy set. Both sorted."""
    tot, j = 0.0, 0
    for s, e, _ in waits:
        while j < len(merged) and merged[j][1] <= s:
            j += 1
        k = j
        while k < len(merged) and merged[k][0] < e:
            tot += min(e, merged[k][1]) - max(s, merged[k][0])
            k += 1
    return tot


def analyse_rank(mpi_csv, gpu_csv):
    header, _ = read_csv(mpi_csv)
    ev = col(header, "event") or col(header, "name")
    if ev is None:
        return None
    waits = intervals(mpi_csv, name_filter="Waitall", name_key=("event",) if col(header, "event") is not None else ("name",))
    if not waits:
        return None
    posts = intervals(mpi_csv, name_filter="Isend", name_key=("event",) if col(header, "event") is not None else ("name",))
    posts += intervals(mpi_csv, name_filter="Irecv", name_key=("event",) if col(header, "event") is not None else ("name",))
    posts.sort()
    d = sorted(e - s for s, e, _ in waits)
    tot = sum(d)
    r = {"n": len(d), "total_s": tot / 1e9, "mean_us": tot / len(d) / 1e3,
         "p50_us": d[len(d)//2] / 1e3, "p90_us": d[int(0.9*len(d))] / 1e3,
         "max_us": d[-1] / 1e3}
    top1 = max(1, len(d)//100)
    r["top1pct"] = 100 * sum(d[-top1:]) / tot if tot else 0.0
    r["top10pct"] = 100 * sum(d[-max(1, len(d)//10):]) / tot if tot else 0.0
    if gpu_csv and Path(gpu_csv).is_file():
        k = intervals(gpu_csv, name_key=("name",))
        if k:
            ov = overlap_total(waits, merge(k))
            r["gpu_overlap_pct"] = 100 * ov / tot if tot else 0.0
            r["gpu_busy_s"] = sum(e - s for s, e in merge(k)) / 1e9
    # gap from the last posted request to entry into the Waitall that follows it
    gaps, i = [], 0
    for s, e, _ in waits:
        while i + 1 < len(posts) and posts[i+1][0] < s:
            i += 1
        if posts and posts[i][1] <= s:
            gaps.append(s - posts[i][1])
    if gaps:
        gaps.sort()
        r["post_gap_us"] = gaps[len(gaps)//2] / 1e3
    return r


def main():
    res = Path(sys.argv[1] if len(sys.argv) > 1 else "results_nsys")
    runs = sorted(d for d in res.iterdir() if d.is_dir())
    print("# Timeline probe — what a rank is waiting for inside MPI_Waitall\n")
    print("Nsight Systems, `--trace=mpi,cuda`, one report per rank. Durations and")
    print("overlaps are per-report (same process, same clock); no cross-rank")
    print("timestamp is compared, because nsys clock alignment across nodes is not")
    print("good enough to call one rank late by microseconds.\n")

    for run in runs:
        mpis = sorted(run.glob("*mpi_event_trace*.csv"))
        if not mpis:
            continue
        hosts = (run / "hosts.txt").read_text().split() if (run / "hosts.txt").is_file() else []
        void = " **VOID**" if (run / "VOID").is_file() else ""
        print(f"\n## `{run.name}`{void} — {len(hosts)} node(s): {' '.join(hosts)}\n")
        print("| rank | Waitall calls | total s | mean | p50 | p90 | max | top 1 % share | top 10 % | **GPU overlap** | post->wait |")
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        rows = []
        for m in mpis:
            mm = re.search(r"rep_(\d+)", m.name)
            rank = int(mm.group(1)) if mm else -1
            g = sorted(run.glob(f"rep_{rank}_*cuda_gpu_trace*.csv"))
            a = analyse_rank(m, g[0] if g else None)
            if not a:
                continue
            a["rank"] = rank
            rows.append(a)
        for a in sorted(rows, key=lambda a: a["rank"]):
            ov = f"{a['gpu_overlap_pct']:.1f} %" if "gpu_overlap_pct" in a else "-"
            pg = f"{a['post_gap_us']:.1f} us" if "post_gap_us" in a else "-"
            print(f"| {a['rank']} | {a['n']} | {a['total_s']:.3f} | {a['mean_us']:.1f} us "
                  f"| {a['p50_us']:.1f} | {a['p90_us']:.1f} | {a['max_us']:.1f} "
                  f"| {a['top1pct']:.0f} % | {a['top10pct']:.0f} % | **{ov}** | {pg} |")
        if rows:
            tots = [a["total_s"] for a in rows]
            lo, hi = min(tots), max(tots)
            arg = max(rows, key=lambda a: a["total_s"])["rank"]
            print(f"\n- wait spread across ranks: min {lo:.3f} s, max {hi:.3f} s, "
                  f"**max/min {hi/lo if lo else float('nan'):.2f}**, slowest rank {arg}")
            ovs = [a["gpu_overlap_pct"] for a in rows if "gpu_overlap_pct" in a]
            if ovs:
                print(f"- GPU overlap of the wait: {min(ovs):.1f} – {max(ovs):.1f} % across ranks")
                if min(ovs) > 50:
                    print("  - **the wait is device work, not the wire**: the Waitall is")
                    print("    running concurrently with CUDA kernels on its own GPU.")
                elif max(ovs) < 10:
                    print("  - the GPU is idle through the wait: this is genuinely MPI.")
            t1 = [a["top1pct"] for a in rows]
            if min(t1) > 40:
                print(f"- **concentrated**: the slowest 1 % of calls hold "
                      f"{min(t1):.0f}–{max(t1):.0f} % of all wait — the mean is not the shape.")
            elif max(t1) < 10:
                print("- uniform across rounds: no single round dominates.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
