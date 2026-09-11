#!/usr/bin/env python3
"""What the per-launch fixed cost of the exchange kernels consists of.

    collect_timeline.py <results_dir>

Pairs pass A (untraced `exch_timing`/`proj_timing` brackets, i.e. the host-side
cost per launch) with pass B (nsys: device execution time, the CUDA API mix, and
the host-to-device copies that precede each launch).

The load-bearing quantity is
        untraced host bracket  -  traced DEVICE duration,
because nsys inflates host timings but not GPU execution. Anything in that
difference is host-side work that is not the kernel. API-call DURATIONS from the
trace are quoted only to say which calls are involved; their magnitudes are
inflated and are marked as such. Call COUNTS are exact.
"""
import csv
import re
import sys
from pathlib import Path

# exch_timing / proj_timing bucket -> the kernels that run inside it, with how
# many times each runs per step. 21 velocity rounds and 18 scalar rounds make up
# the 39; the cross-level kernel runs in 6 full velocity exchanges + 18 scalar.
BUCKETS = {
    "pack":       [("comm_pack_entries", 21), ("comm_pack_scalar_entries", 18)],
    "unpack":     [("comm_unpack_entries", 21), ("comm_unpack_scalar_entries", 18)],
    "local_copy": [("comm_copy_local_same_level", 21),
                   ("comm_copy_local_scalar_same_level", 18)],
    "copy_cross": [("comm_copy_local_cross_level", 6),
                   ("comm_copy_local_scalar_entries", 18)],
    "sweep":      [("pressure_solver_jacobi_compute_phi", 18)],
    "apply":      [("pressure_solver_jacobi_apply", 36),
                   ("pressure_solver_interface_correct", 18)],
}
KERNEL_RE = re.compile(r"nvkernel_(.+?)__F1L(\d+)_")


def short(name):
    m = KERNEL_RE.match(name)
    return m.group(1) if m else name


def read_csv(path):
    if not path.is_file():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def col(row, *names):
    """First column whose header contains one of `names` (case-insensitive)."""
    for n in names:
        for k in row:
            if n.lower() in k.lower():
                return k
    return None


def host_brackets(run):
    """Per-call microseconds from an untraced run's profiler lines."""
    out = {}
    log = run / "run.log"
    if not log.is_file():
        return out
    for line in log.read_text(errors="replace").splitlines():
        m = re.search(r"(?:exch|proj)_timing: (\S+) calls (\d+) nsteps (\d+) "
                      r"seconds\s+\S+ seconds_per_step\s+\S+ seconds_per_call\s+(\S+)", line)
        if m and int(m.group(2)):
            out[m.group(1)] = (float(m.group(4)) * 1e6, int(m.group(2)) // int(m.group(3)))
    return out


def analyse_trace(run):
    """Per-kernel device time, and the memcpys that precede each launch.

    Everything is measured inside a STEADY WINDOW (40-95 % of the trace) and
    normalised by the number of steps in it, counted from jacobi_compute_phi's
    18 launches per step -- so init-time kernels cannot leak into a per-step
    figure.
    """
    rows = read_csv(run / "rep_0_cuda_gpu_trace.csv")
    if not rows:
        return None
    r0 = rows[0]
    cs, cd, cn = col(r0, "Start"), col(r0, "Duration"), col(r0, "Name")
    cb = col(r0, "Bytes")
    ev = []
    for r in rows:
        try:
            ev.append((int(r[cs]), int(r[cd]), r[cn],
                       float(r[cb]) if cb and r[cb] else 0.0))
        except (ValueError, KeyError):
            continue
    ev.sort()
    t0, t1 = ev[0][0], ev[-1][0]
    lo, hi = t0 + 0.40 * (t1 - t0), t0 + 0.95 * (t1 - t0)
    win = [e for e in ev if lo <= e[0] <= hi]

    nphi = sum(1 for e in win if "jacobi_compute_phi" in e[2])
    steps = nphi / 18.0
    if steps <= 0:
        return None

    kern, pend_n, pend_b, pending = {}, {}, {}, [0, 0.0]
    for start, dur, name, mb in win:
        if name.startswith("[CUDA"):
            if "Host-to-Device" in name:
                pending[0] += 1
                pending[1] += mb
            continue
        k = short(name)
        d = kern.setdefault(k, [0, 0])
        d[0] += 1
        d[1] += dur
        pend_n[k] = pend_n.get(k, 0) + pending[0]
        pend_b[k] = pend_b.get(k, 0.0) + pending[1]
        pending[0], pending[1] = 0, 0.0

    h2d = [e for e in win if "Host-to-Device" in e[2]]
    return dict(steps=steps, kern=kern, pend_n=pend_n, pend_b=pend_b,
                h2d_count=len(h2d),
                h2d_mb=sum(e[3] for e in h2d),
                h2d_ns=sum(e[1] for e in h2d),
                gpu_busy=sum(e[1] for e in win if not e[2].startswith("[CUDA")),
                span=hi - lo)


def main():
    res = Path(sys.argv[1])
    print("# What the per-launch fixed cost is made of\n")
    prov = res / "provenance.txt"
    if prov.is_file():
        print("```\n" + prov.read_text().strip() + "\n```\n")

    for nrun in sorted(res.glob("nsys_*")):
        cfg = nrun.name[len("nsys_"):]
        hrun = res / ("host_" + cfg)
        tr = analyse_trace(nrun)
        hb = host_brackets(hrun)
        print(f"\n## {cfg}\n")
        if tr is None:
            print("_no usable GPU trace_\n")
            continue
        print(f"Steady window {tr['span']/1e6:.0f} ms = {tr['steps']:.1f} steps; "
              f"GPU busy {100*tr['gpu_busy']/tr['span']:.1f} % of it.\n")

        print("### Per kernel: what the host pays, and what the device does\n")
        print("| kernel | launches/step | device us/launch | H2D copies before each | "
              "bucket | host us/launch (untraced) | host − device |")
        print("|" + "---|" * 7)
        for bucket, parts in BUCKETS.items():
            for kname, per_step in parts:
                d = next((v for k, v in tr["kern"].items() if kname in k), None)
                if d is None:
                    continue
                n, tot = d
                dev = tot / n / 1e3
                pn = tr["pend_n"].get(next(k for k in tr["kern"] if kname in k), 0) / n
                hostcall = hb.get(bucket, (0.0, 0))[0]
                print(f"| `{kname}` | {n/tr['steps']:.1f} | **{dev:.1f}** | {pn:.1f} | "
                      f"`{bucket}` | {hostcall:.1f} | "
                      f"{hostcall - dev:.1f} |")
        print("\n(The `bucket` column repeats for the velocity and scalar variants that")
        print("share one profiler bucket; its host figure is their average, so compare")
        print("it against the two device figures together, not against either alone.)\n")

        print("### Host-to-device copies\n")
        print(f"- **{tr['h2d_count']/tr['steps']:.0f} per step**, "
              f"{1e3*tr['h2d_mb']/tr['h2d_count']:.2f} kB each on average "
              f"({tr['h2d_mb']/tr['steps']:.2f} MB/step in total)")
        print(f"- {tr['h2d_ns']/tr['h2d_count']/1e3:.2f} us each on the device, "
              f"{tr['h2d_ns']/tr['steps']/1e6:.2f} ms/step of device time\n")

        api = read_csv(nrun / "rep_0_cuda_api_sum.csv")
        if api:
            cn2 = col(api[0], "Name")
            cc = col(api[0], "Num Calls", "Calls", "Instances")
            ct = col(api[0], "Total Time", "Total")
            cavg = col(api[0], "Avg")
            print("### CUDA API calls per step (counts exact; DURATIONS INFLATED BY TRACING)\n")
            print("| API | calls/step | avg us (inflated) | share of API time |")
            print("|" + "---|" * 4)
            tot = sum(float(r[ct]) for r in api if r.get(ct))
            for r in api[:12]:
                try:
                    c = float(r[cc])
                except (TypeError, ValueError):
                    continue
                print(f"| `{r[cn2]}` | {c/tr['steps']:.0f} | "
                      f"{float(r[cavg])/1e3:.2f} | {100*float(r[ct])/tot:.1f} % |")
            print()


if __name__ == "__main__":
    main()
