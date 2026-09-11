#!/usr/bin/env python3
"""What every kernel launch copies to the device, from an nsys sqlite export.

    analyse_launch_traffic.py <rep.sqlite> [<rep.sqlite> ...]

`nsys stats` summaries say how much host-to-device traffic there is; they do not
say WHICH launch pays for it. This walks the HOST-side CUDA runtime timeline --
the sequence is `cuMemcpyHtoDAsync`* -> `cuLaunchKernel` -> `cuStreamSynchronize`
for every OpenMP target region -- and attributes each copy to the launch that
follows it, then reports the SIZE HISTOGRAM per kernel. The size is what names
the mechanism: 3-8 B is a scalar, 152/200 B is an nvfortran rank-1/rank-2 array
descriptor, and one large block is a whole derived-type object.

Host-side attribution, not GPU-side: kernel execution is asynchronous and
reordered against issue order, so attributing by GPU timestamp mixes regions.
Counts and byte sizes are exact and unaffected by tracing overhead; the elapsed
host time per launch is NOT, and is deliberately not reported here.

Everything is measured in a steady window (40-95 % of the trace) and normalised
by the steps in it, counted from jacobi_compute_phi's 18 launches per step.
"""
import collections
import sqlite3
import sys

PHI_PER_STEP = 18.0


def analyse(path):
    c = sqlite3.connect(path)
    names = dict(c.execute("select id, value from StringIds"))
    kern = {cid: (names.get(sn) or names.get(dn)) for cid, sn, dn in c.execute(
        "select correlationId, shortName, demangledName from CUPTI_ACTIVITY_KIND_KERNEL")}
    nbytes = dict(c.execute(
        "select correlationId, bytes from CUPTI_ACTIVITY_KIND_MEMCPY"))
    ev = list(c.execute("select start, nameId, correlationId "
                        "from CUPTI_ACTIVITY_KIND_RUNTIME order by start"))
    if not ev:
        return None
    t0, t1 = ev[0][0], ev[-1][0]
    lo, hi = t0 + 0.40 * (t1 - t0), t0 + 0.95 * (t1 - t0)

    hist = collections.defaultdict(collections.Counter)
    launches = collections.Counter()
    pend = []
    for start, nid, cid in ev:
        if start > hi:
            break
        if start < lo:
            continue
        nm = names[nid]
        if nm.startswith("cuMemcpyHtoDAsync"):
            pend.append(nbytes.get(cid, 0))
        elif nm.startswith("cuLaunchKernel"):
            k = (kern.get(cid) or "?").replace("nvkernel_", "").split("__F1L")[0]
            hist[k].update(pend)
            launches[k] += 1
            pend = []
    steps = sum(n for k, n in launches.items()
                if "jacobi_compute_phi" in k) / PHI_PER_STEP
    return hist, launches, steps


def main():
    for path in sys.argv[1:]:
        got = analyse(path)
        print(f"\n## {path}\n")
        if not got:
            print("_no runtime events_")
            continue
        hist, launches, steps = got
        if steps <= 0:
            print("_no jacobi_compute_phi launches; cannot normalise_")
            continue
        print(f"Steady window = {steps:.1f} steps\n")
        print("| kernel | launches/step | H2D copies/launch | bytes/launch | size histogram per launch |")
        print("|" + "---|" * 5)
        for k, n in sorted(launches.items(), key=lambda x: -x[1]):
            if n < 5:
                continue
            h = hist[k]
            tot = sum(h.values())
            b = sum(s * m for s, m in h.items())
            sizes = ", ".join(f"{s:,} B x{m/n:.0f}" for s, m in sorted(h.items()))
            print(f"| `{k}` | {n/steps:.1f} | {tot/n:.1f} | {b/n:,.0f} | {sizes} |")


if __name__ == "__main__":
    main()
