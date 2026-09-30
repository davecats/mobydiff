# Pre-registered: the predictor register guard (written before job submission, 2026-09-30)

**Change.** `step.f90 momentum`: the three skew-symmetric convection
corrections (u, v, w) are back inside `if (skew)`, where `skew` is a local
logical set from runtime data (`dns%re > 0`, always true since the S3
lockdown hardwired skew convection) and mapped `to:` the kernel. One
never-taken branch per component; no expression moves.

**Why.** `results_horeka_2026-09-26.md` finding 1 + ncu job 5163917: removing
the guards let nvfortran cut the fused predictor from 128 to 100 registers,
occupancy did NOT move (23.76 -> 23.69 %, the kernel is capped by grid/block
shape under `collapse(4)`), and both utilisation axes fell (DRAM 34.51 ->
31.02 % of peak, SM 40.33 -> 34.09 %): the kernel time rose 15,226 -> 16,938 us
(+11.2 %) with traffic byte-identical, i.e. the shorter live ranges lengthened
the dependency chains that 15 of 48 warps cannot hide. The step lost
1.3-2.0 % on every config and rank count.

**Predictions (ref = the head without the guard, new = with it).**

1. ncu on `refined_yp82_rect_jacobi`, predictor kernel (`step_momentum` F1L<n>_6):
   registers **100 -> 128** (the 8fa0fc2 schedule), occupancy unchanged within
   0.1 point (~23.7 %), doubles/cell and sectors/request identical, DRAM %peak
   **31 -> ~34.5**, SM %peak **34 -> ~40**, kernel time **-9 to -11 %**. The
   second momentum kernel (F1L<n>_8, 66 registers) is the control: unmoved.
   If registers stay at 100 the compiler folded the guard and the change is
   inert -- then `skew` needs a source the compiler sees even less of, not a
   different guard.
2. Step at 4 and 8 ranks, 200 steps: `momentum` bucket **-6 to -7 %**
   (the 09-26 delta was +2.644 ms on 41.4 at `rect` 4 ranks), step **-1.3 to
   -2.0 %** on base / rect / refined_yp82 alike, flat in rank count (it is a
   per-cell effect on a kernel whose share of the step is rank-independent).
   Every other bucket within noise (±0.3 %).
3. Bit-exact under nofma on the 7-case and 9-case suites, CPU and GPU, against
   `~/step8_ref_binaries` (max_abs 0 on every dataset): the guard adds no
   arithmetic. Under PRODUCTION flags the comparison may read ulps (FMA
   contraction can follow the new schedule) -- not a gate, not a defect.

**What would falsify it.** Registers back at 128 but the time unmoved: the
09-26 mechanism (register-driven scheduling) is wrong and the +11 % has
another cause -- then look at the SASS, not at registers. Registers unmoved:
the guard was folded (see 1).
