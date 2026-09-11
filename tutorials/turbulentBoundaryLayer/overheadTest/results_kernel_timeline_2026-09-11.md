# The per-launch cost is an 8 kB derived-type descriptor blob, re-sent every launch

The nsys per-kernel timeline that `results_horeka_exchange_2026-09-10.md` §8 made
the precondition for any implementation. Job 5141872, HoreKa `dev_accelerated`,
hkn0402, 1 node, 4 ranks, 30 steps, commit `f7a597a`. Two passes over the same
two configs in one allocation: **A** untraced (the `exch_timing`/`proj_timing`
brackets) and **B** under `nsys --trace=cuda`.

**Answer: every launch of an exchange kernel copies the entire `comm_type` object
— 7 952 bytes, all 45 of its array descriptors — from host to device, plus 14–23
separate descriptor and scalar copies, whatever the kernel actually names.** The
projection kernels, which touch no such object, pay 8–148 bytes and are 3–17x
cheaper per launch. That is the 56–114 us.

## 1 — Method, and why tracing overhead does not reach the conclusion

nsys inflates host-side timings unevenly (`results_horeka_2026-09-10.md` §5,
"READ SHAPE, NOT MAGNITUDE"). It does not inflate a kernel's GPU execution time.
So the load-bearing quantity is

> **(untraced host bracket) − (traced device duration)**,

whose second term is tracing-independent: whatever the difference is, it is
host-side work that is not the kernel. Byte sizes and call counts are exact.
Host-side *durations* from the trace are used nowhere in this report's
conclusions.

The method's own noise floor is measurable from the control: `jacobi_compute_phi`
gives −15.1 us on `rect_jacobi` (a 1 553 us bracket) and +9.7 us on
`refined_yp82` (699 us) — so about **1 % of the device time, or ~15 us, whichever
is larger**. For the exchange kernels, whose device time is 6–22 us, that is
under 4 us against differences of 62–88.

## 2 — The per-launch host-side cost, measured a second way

Device times are launch-weighted over each bucket's velocity and scalar variants
(21 : 18 rounds per step; 6 : 18 for the cross-level kernel):

| bucket | host us/launch (untraced) | device us/launch (traced) | **host − device** | 2026-09-10 fit intercept |
|---|---|---|---|---|
| `pack` (refined / rect) | 100.3 / 102.9 | 12.6 / 14.6 | **87.7 / 88.3** | 92–114 |
| `unpack` | 87.0 / 92.8 | 9.3 / 15.1 | **77.7 / 77.7** | 83–100 |
| `local_copy` (same level) | 199.0 / 426.3 | 136.7 / 363.6 | **62.3 / 62.7** | 56–67 |
| `copy_cross` (refined) | 116.7 | 35.3 | **81.4** | 84.6–87.3 |
| `sweep` | 698.5 / 1553.4 | 688.8 / 1568.5 | 9.7 / −15.1 | 24–26 |

**Two configurations whose device times differ by 3x give the same host-side cost
to within 1 %** (62.3/62.7, 77.7/77.7, 87.7/88.3), and those costs reproduce the
affine-fit intercepts of 2026-09-10 — derived from an entirely different
measurement, the slope of per-round time against points across rank counts. The
"fixed cost per launch" reading is now established by two independent methods.

The third pre-registered outcome — *device duration ≈ host bracket, the kernels
are simply inefficient at small sizes, and sections 5 and 8 are wrong* — **did not
occur.**

The `apply` bucket needs its own arithmetic because it contains five launches per
call on a refined grid (2 × `jacobi_apply` + **3** × `interface_correct`, not the
one that report assumed): device 1 986.8 us against a 2 095.8 us bracket = 109.0
us per call, against 9.6 us per call for the single-level case's two launches.
**99 us for three extra launches, i.e. ~33 us each.**

## 3 — What is actually copied, per launch

Host-side attribution (`horeka/exchange/analyse_launch_traffic.py`): the CUDA
runtime sequence for an OpenMP target region is `cuMemcpyHtoDAsync`* →
`cuLaunchKernel` → `cuStreamSynchronize`, so each copy belongs to the launch that
follows it. `refined_yp82`, steady window, 30 steps:

| kernel | launches/step | H2D copies/launch | bytes/launch | size histogram per launch |
|---|---|---|---|---|
| `comm_pack_entries` | 21.0 | 14.0 | **10 264** | 152 B ×6, 200 B ×7, **7 952 B ×1** |
| `comm_pack_scalar_entries` | 18.0 | 14.0 | **10 264** | 152 B ×6, 200 B ×7, **7 952 B ×1** |
| `comm_unpack_entries` | 21.0 | 23.0 | **10 192** | 8 B ×9, 152 B ×9, 200 B ×4, **7 952 B ×1** |
| `comm_unpack_scalar_entries` | 18.0 | 20.0 | **9 400** | 8 B ×11, 152 B ×5, 200 B ×3, **7 952 B ×1** |
| `comm_copy_local_same_level` | 21.0 | 21.0 | **9 264** | 8 B ×13, 152 B ×4, 200 B ×3, **7 952 B ×1** |
| `comm_copy_local_scalar_same_level` | 18.0 | 21.0 | **9 264** | 8 B ×13, 152 B ×4, 200 B ×3, **7 952 B ×1** |
| `comm_copy_local_cross_level` | 6.0 | 22.0 | **10 472** | 8 B ×7, 152 B ×7, 200 B ×7, **7 952 B ×1** |
| `comm_copy_local_scalar_entries` | 18.0 | 19.0 | **9 968** | 8 B ×7, 152 B ×5, 200 B ×6, **7 952 B ×1** |
| `pressure_solver_jacobi_compute_phi` | 18.0 | 11.8 | 95 | 8 B ×12 |
| `pressure_solver_jacobi_apply` | 36.0 | 2.0 | 8 | 3 B ×2 |
| `pressure_solver_interface_correct` | 54.0 | 3.0 | 9 | 3 B ×3 |
| `boundary_apply_bc` | 21.0 | 1.0 | 12 | 12 B ×1 |
| `step_momentum` | 5.9 | 2.3 | 27 | 8 B ×2 |

**Every `comm_*` kernel pays exactly one 7 952-byte copy. No other kernel in the
solver pays one at all.** The single-level `rect_jacobi` trace is identical in
this respect (7 952 B ×1 on `pack`/`unpack`/`copy_local_same_level`, none on the
projection kernels), so it is not a refinement effect.

### The 7 952 bytes are `comm_type` itself

`comm_type` declares **45 allocatable array components** — 26 of rank 1 and 19 of
rank 2. The trace's other two sizes identify nvfortran's descriptors directly:
the rank-1 arrays it copies individually are **152 B** and the rank-2 ones
**200 B**. So

> 26 × 152 B + 19 × 200 B = **7 752 B** of array descriptors,

and the remaining ~200 B is the type's scalars, logicals, MPI handles and the
`dims`/`activeVars`/op-counter arrays. **The 7 952-byte block is the whole
`comm_type` object, descriptor table included, re-sent on every launch — even
though `copy_local_same_level` names 8 of those 45 components and
`pack_entries` 15.**

The three sizes also explain the 14–23 separate small copies: they are the
individually-mapped components (152/200 B each) and the mapped scalars (8 B).
A kernel that maps more `c%` components pays more of them — which is why the
first addendum of `horeka/exchange/PREREGISTERED.md` saw a map-count correlation
at all, and why counting map items alone could never reproduce it: the dominant
term is a per-launch constant that does not depend on how many components the
kernel names.

Across the whole step: **2 736 (single level) to 3 484 (refined) host-to-device
copies per step**, ~1 MB in total, against 298–407 `cuLaunchKernel` calls — about
**12 copies per launch**. The ~1 300/step that `results_horeka_2026-09-10.md` §5
saw and left unchased were this, undercounted because that run had fewer
exchange kernels on a coarser configuration.

## 4 — The pre-registered scorecard (third addendum)

| prediction | outcome |
|---|---|
| control: `jacobi_compute_phi` shows ~25 us of host-side cost per launch | **PARTIAL** — 9.7 and −15.1 us, i.e. 0–25 us. Consistent with the 24–26 us intercept but not resolved; the method's floor is ~15 us here. |
| the outcome that would sink §5 and §8 (device ≈ host for the exchange kernels) | **did NOT occur** — 62–88 us, 20x the floor |
| H2D per launch tracks the count of components/module arrays: 4–6 sweep, 9–11 copy, 15–17 pack/unpack, 16–18 cross | **REFUTED in detail** — measured 11.8, 21, 14/23, 22. The counts are higher and do not track the component count |
| if the counts track → the fixed cost is map-clause marshalling | **fires in substance** — it *is* marshalling, but the variable I failed to predict is that nvfortran ships the whole parent object, not one copy per named component |

So the mechanism named in the first addendum was right in kind and wrong in every
quantitative detail, twice. The trace was needed.

## 5 — What follows, and what it is worth

The `c%` arrays are **already resident** — `init_block_exchange` maps all of them
in a `target enter data` region that lives for the run. Nothing in the 10 kB per
launch is data the device needs; it is the descriptor table being re-materialised
because the kernel body writes `c%lOff`, `c%sExt`, … and nvfortran therefore
marshals the parent object.

Recoverable, estimated from the floor the projection kernels already achieve
(5–22 us per launch against the exchange kernels' 62–88):

| | launches/step | us/launch now | recoverable | ms/step |
|---|---|---|---|---|
| `pack` | 39 | 88 | ~66–78 | 2.6–3.0 |
| `unpack` | 39 | 78 | ~56–68 | 2.2–2.7 |
| `local_copy` | 39 | 62 | ~40–52 | 1.6–2.0 |
| `copy_cross` | 24 | 81 | ~59–71 | 1.4–1.7 |
| **total** | 141 | | | **7.8–9.4 ms/step** |

Against the 16-rank steps measured on 2026-09-10 that is **19–23 % of the refined
step** and 12–14 % of the single-level one — i.e. essentially the whole
launch-fixed bill that report attributed but could not explain.

**The fix direction is to stop referencing `c%component` inside the target
regions** — bind what each kernel needs to local arrays/pointers outside the
region, or pass them as dummy arguments, so the map list contains plain arrays
that are already present rather than a derived type. Which exact form nvfortran
rewards is an empirical question, and **the next step is one kernel, not six**:
take `copy_local_same_level` (11 map items, the simplest), hoist its `c%`
references, and re-measure (a) its H2D histogram — the 7 952 B block must
disappear — and (b) its `local_copy` bracket, which should fall from 62 us
toward the ~20 us floor. It is a scheduling change, so it must be **bit-exact**:
the 7-case suite plus the production-case Pass G, CPU and GPU.

## 6 — What this does NOT support

- **It does not show that the descriptor traffic is the whole 62–88 us.** It
  shows 10 kB in ~20 transfers per launch and a cost consistent with that; the
  proof is the A/B above, which has not been run.
- **It does not measure the host time of those copies.** Their *device* time is
  ~2 us each (5.5–6.8 ms/step), and the API durations in the trace are inflated.
  The 62–88 us figure comes from the untraced brackets, not from summing API
  calls.
- **It says nothing about `mpi_wait`**, which is a separate 2–10 % of the step and
  was closed as a partitioning/overlap target on 2026-09-10.
- **4 ranks, one node, one machine.** The fixed cost is a per-launch constant
  measured at 2, 4, 8 and 16 ranks, so the dissection should carry — but it was
  dissected at 4.
